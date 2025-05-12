#![feature(
    slice_from_ptr_range, 
    new_range_api, 
    hash_set_entry, 
    assert_matches, 
    likely_unlikely,
    cold_path,
)]

mod alignment;
mod position;
mod task;
mod utils;

use position::{fill_positions, Position};
use rmp_serde::from_read;
use task::{scan_alignment_segments, TaskResult};
use utils::asc2dnacomp;

use std::{hint::cold_path, path::Path, sync::{mpsc, LazyLock}};
use anyhow::Result;
use memmap2::{Advice, Mmap};
use rayon::{iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator}, ThreadPoolBuilder};
use clap::Parser;

use std::{fs::File, path::PathBuf};
use std::collections::HashMap;
use ascii::{AsciiString, ToAsciiChar};
use std::io::{Write, BufReader};
use ahash::AHashMap;
use crate::task::{Task2, TaskIter2};

#[derive(clap::Parser, Debug)]
#[command(version, about)]
struct Arguments {
    #[arg(
        long = "alignments",
        value_name = "alignmentFile",
        help = "SORTED SAM filename. Please enter '-' for standard input."
    )]
    alignment_file: PathBuf,
    #[arg(
        long = "refIndex",
        value_name = "refFileIndex",
        help = "reference file (should be dna_index's output for an FASTA format reference file)."
    )]
    reference_file_index: PathBuf,
    #[arg(
        long,
        value_name = "outputFile",
        help = "file name to save the 3n table (tsv format). By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console)."
    )]
    output_name: PathBuf,
    #[arg(
        long, 
        value_parser = |s: &str| -> Result<((u8, u8), (u8, u8)), String> {
            let s = Vec::from_iter(s.trim().split(','));
            if s.len() != 2 || !s.iter().all(|b| b.len() == 1) {
                return Err("format error".to_owned())                
            }
            let bases = utils::BASE_CHARS;
            let from = s[0].chars().next().unwrap().to_ascii_char().unwrap().to_ascii_uppercase();
            let from = u8::try_from(from).unwrap();
            let from_comp = asc2dnacomp(from);
            let to = s[1].chars().next().unwrap().to_ascii_char().unwrap().to_ascii_uppercase();
            let to = u8::try_from(to).unwrap();
            let to_comp = asc2dnacomp(to);
            if !bases.contains(&from) || !bases.contains(&to) {
                return Err("no such base (or use uppercase)".to_owned());
            }
            Ok(((from, from_comp), (to, to_comp)))
        },
        help = "the char1 is the nucleotide converted from, the char2 is the nucleotide converted to."
    )]
    /// ((convert_from, complement), (convert_to, convert_to_complement))
    base_change: ((u8, u8), (u8, u8)),
    #[arg(
        short,
        long,
        default_value_t = false,
        help = "only count the base which is in unique mapped reads."
    )]
    unique_only: bool,
    #[arg(
        short,
        long,
        default_value_t = false,
        help = "only count the base which is in multiple mapped reads."
    )]
    multiple_only: bool,
    #[arg(
        short,
        long,
        default_value_t = false,
        help = "only count CG and ignore CH in reference."
    )]
    cg_only: bool,
    #[arg(
        short,
        long,
        default_value_t = false,
        help = "please add this option if you use --add-chrname during HISAT-3N alignment."
    )]
    added_chrname: bool,
    #[arg(
        short,
        long,
        default_value_t = false,
        help = "please add this option if you use --remove-chrname during HISAT-3N alignment."
    )]
    removed_chrname: bool,
    #[arg(
        short = 'p',
        long,
        default_value_t = 1,
        help = "number of threads to launch (1)."
    )]
    threads: usize,
    #[arg(
        long,
        default_value_t = 20000000,
        help = "max number of Alignment record lines in a Task (20000)",
    )]
    align_block_size: usize,
    #[arg(
        long,
        default_value_t = 20000000,
        help = "max number of chromosome Position s in a Task (20000)",
    )]
    ref_block_size: usize,
}

static ARGS: LazyLock<Arguments> = LazyLock::new(|| { Arguments::parse() });

fn static_mmap_str(p: &Path) -> &'static [u8] {
    let alignment_file = Box::new(File::open(p).unwrap());
    let alignment_file: &'static File = Box::leak(alignment_file);
    let alignment_map = Box::new(unsafe {
        let mmap = Mmap::map(alignment_file).unwrap();
        mmap.advise(Advice::Sequential).unwrap();
        mmap
    });
    Box::leak(alignment_map)
}

// a comprehensive survey shows that LazyLock has no sync overhead after init
// deref ops after init is just like normal deref ops
static ALIGN_FILE: LazyLock<&'static [u8]> = LazyLock::new(|| static_mmap_str(&ARGS.alignment_file));
static DNAS: LazyLock<AHashMap<&'static [u8], &'static [u8]>> = LazyLock::new(|| {
    let ref_index_file = BufReader::new(File::open(&ARGS.reference_file_index).unwrap());
    let by_ascii: HashMap::<AsciiString, AsciiString> = from_read(ref_index_file).unwrap();
    let dnas: AHashMap<_, _> = by_ascii
      .into_iter()
      .map(|(k, v)| { (Box::leak(k.into_boxed_ascii_str()).as_bytes(), Box::leak(v.into_boxed_ascii_str()).as_bytes()) })
      .collect();
    dnas
});

#[inline(never)]
fn worker2(task: Task2<'static>) -> Vec<Position<'static>> {
    let mut positions = Vec::new();
    Vec::reserve(&mut positions, task.position_range.len());
    let dna_name = task.dna_name;
    // let ulen = DNAS.get(dna_name).unwrap().len();
    // eprintln!("{}, {}", str::from_utf8(dna_name).unwrap(), ulen);
    fill_positions(&mut positions, DNAS.get(dna_name).unwrap(), dna_name, task.position_range.start, task.position_range.end);

    for alignment in task.alignments {
        debug_assert_eq!(alignment.dna, task.dna_name);
        if !alignment.mapped || alignment.bases.is_empty() {
            continue;
        }
        // int firstPos = refPositions[0]->location;
        //         return targetPos - firstPos;
        for base in &alignment.bases {
            if base.remove {
                continue;
            }

            let index = (alignment.location as usize) - task.position_range.start + (TryInto::<usize>::try_into(base.ref_pos).unwrap());
            if index >= positions.len() {
                continue;
            }
            // eprintln!("index: {}, positions: {}", index, positions.len());
            let position = &mut positions[index as usize];
            assert_eq!(position.location, alignment.location + base.ref_pos);

            if position.strand.is_none() {
                continue;
            }

            position.append_base(base, &alignment);
        }
    }

    positions
}

fn main() -> Result<()> {
    ThreadPoolBuilder::new().num_threads(ARGS.threads).build_global()?;

    let (tx, rx) = mpsc::channel();
    let dna_align_segments: Vec<(&[u8], std::ops::Range<usize>)> = scan_alignment_segments(&ALIGN_FILE);

    std::thread::spawn(move || {
        let tx = tx.clone();
        dna_align_segments
            .par_iter()
            .flat_map(|(_, r)| TaskIter2::new(&ALIGN_FILE[r.start..r.end]).par_bridge())
            .map(worker2)
            .for_each(|positions| { tx.send(Some(positions)).unwrap(); });
        
        tx.send(None).unwrap();
    });

    let mut output = std::io::BufWriter::with_capacity(1 * 1024 * 1024, File::create(&ARGS.output_name)?);

    writeln!(output, "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount")?;

    loop {
        let res = rx.recv()?;
        match res {
            TaskResult::Some(positions) => {
                for p in positions {
                    if p.converted_qualities.is_empty() && p.unconverted_qualities.is_empty() {
                        continue;
                    }
                    let len1 = p.converted_qualities.len();
                    let len2 = p.unconverted_qualities.len();
                    writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}", str::from_utf8(p.dna).unwrap(), p.location, char::from(p.strand.unwrap_or(b'?')), String::from_utf8(p.converted_qualities).unwrap(), len1, String::from_utf8(p.unconverted_qualities).unwrap(), len2)?;
                }
            }
            TaskResult::None => {
                cold_path();
                break;
            }
        }
    }

    Ok(())
}

#[test]
fn test_arguments() {
    use clap::CommandFactory;
    Arguments::command().debug_assert();
}
