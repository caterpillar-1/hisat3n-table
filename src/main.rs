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

use position::Position;
use task::{Task, TaskIter, TaskResult};

use std::{path::Path, sync::{mpsc, LazyLock}};
use memmap2::{Advice, Mmap};
use rayon::{iter::{ParallelBridge, ParallelIterator}, ThreadPoolBuilder};
use clap::Parser;

use std::{fs::File, path::PathBuf};
use ascii::{AsAsciiStr, AsciiChar, AsciiStr};
use std::io::Write;

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
        long = "ref",
        value_name = "refFile",
        help = "reference file (should be FASTA format)."
    )]
    reference_file: PathBuf,
    #[arg(
        long,
        value_name = "outputFile",
        help = "file name to save the 3n table (tsv format). By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console)."
    )]
    output_name: PathBuf,
    #[arg(
        long, 
        value_parser = |s: &str| -> Result<(char, char), String> {
            let s = Vec::from_iter(s.trim().split(','));
            if s.len() != 2 || !s.iter().all(|b| b.len() == 1) {
                return Err("format error".to_owned())                
            }
            let bases = utils::BASE_CHARS;
            let from = s[0].chars().next().unwrap();
            let to = s[1].chars().next().unwrap();
            if !bases.contains(&from) || !bases.contains(&to) {
                return Err("no such base (or use uppercase)".to_owned());
            }
            Ok((from, to))
        },
        help = "the char1 is the nucleotide converted from, the char2 is the nucleotide converted to."
    )]
    /// (convert_from, convert_to)
    base_change: (char, char),
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
        default_value_t = 1000000,
        help = "max number of Alignment record lines in a Task (1000000)",
    )]
    align_block_size: usize,
    #[arg(
        long,
        default_value_t = 2000000,
        help = "max number of chromosome Position s in a Task (2000000)",
    )]
    ref_block_size: usize,
}

static ARGS: LazyLock<Arguments> = LazyLock::new(|| { Arguments::parse() });

fn static_mmap_asciistr(p: &Path) -> &'static AsciiStr {
    let alignment_file = Box::new(File::open(p).unwrap());
    let alignment_file: &'static File = Box::leak(alignment_file);
    let alignment_map = Box::new(unsafe {
        let mmap = Mmap::map(alignment_file).unwrap();
        mmap.advise(Advice::Sequential).unwrap();
        mmap
    });
    let alignment_map: &'static Mmap = Box::leak(alignment_map);
    &alignment_map.as_ascii_str().unwrap()
}

static ALIGN_FILE: LazyLock<&'static AsciiStr> = LazyLock::new(|| static_mmap_asciistr(&ARGS.alignment_file));
static REF_FILE: LazyLock<&'static AsciiStr> = LazyLock::new(|| static_mmap_asciistr(&ARGS.reference_file));

fn worker(task: Task) -> Vec<Position> {
    let mut positions = Vec::from_iter(task.positions);
    assert!(positions.len() > 0);

    let dna_location = positions[0].location;
    for alignment in task.alignments {
        if !alignment.mapped || alignment.bases.is_empty() {
            continue;
        }

        let align_location = alignment.location;

        for base in &alignment.bases {
            if base.remove {
                continue;
            }

            let index = align_location - dna_location + base.ref_pos;
            assert!(0 <= index && index < positions.len() as isize);
            let position = &mut positions[index as usize];
            
            if position.strand.is_none() {
                continue;
            }

            position.append_base(base, &alignment);
        }
    }

    positions
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    ThreadPoolBuilder::new().num_threads(ARGS.threads).build_global().unwrap();

    let (tx, rx) = mpsc::channel();
    let tasks = TaskIter::new(&ALIGN_FILE, &REF_FILE);

    std::thread::spawn(move || {
        let tx = tx.clone();
        tasks.into_iter().par_bridge().map(worker).for_each(|positions| {
            tx.send(TaskResult::Some(positions)).unwrap();
        });
        
        tx.send(TaskResult::None).unwrap();
    });

    let mut output = std::io::BufWriter::new(File::create(&ARGS.output_name).unwrap());
    writeln!(output, "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount")?;

    loop {
        let res = rx.recv().unwrap();
        match res {
            TaskResult::Some(positions) => {
                for p in positions {
                    if p.empty() || p.strand.is_none() {
                        continue;
                    }
                    writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}", p.dna, p.location, p.strand.unwrap(), p.converted_qualities, p.converted_qualities.len(), p.unconverted_qualities, p.unconverted_qualities.len())?;
                }
            }
            TaskResult::None => {
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
