use std::{collections::HashMap, fs::File, path::PathBuf};

use anyhow::Result;
use ascii::{AsAsciiStr, AsciiChar, AsciiStr, AsciiString, IntoAsciiString};
use memmap2::Mmap;
use serde::{Serialize, Deserialize};
use rmp_serde::{Serializer, Deserializer};

use clap::Parser;

#[derive(Parser, Debug)]
struct Arguments {
    #[arg(
        short = 'r',
    )]
    reference_file: PathBuf,
    #[arg(
        short = 'i',
    )]
    index_file: PathBuf,
}

fn get_dna_name(info_line: &AsciiStr) -> AsciiString {
    assert_eq!(info_line.first().unwrap(), AsciiChar::GreaterThan);
    let info_line = &info_line[1..];
    info_line
        .as_str()
        .split_ascii_whitespace()
        .next()
        .unwrap()
        .into_ascii_string()
        .unwrap()
}

fn main() -> Result<()> {
    let args = Arguments::parse();
    let dna_file = {
        let ref_file = Box::new(File::open(args.reference_file)?);
        let ref_file: &'static File = Box::leak(ref_file);
        let ref_map = Box::new(unsafe {
            let mmap = Mmap::map(ref_file)?;
            mmap.advise(memmap2::Advice::Sequential)?;
            mmap
        });
        Box::leak(ref_map).as_ascii_str()?
    };

    let mut dnas = HashMap::new();
    let mut lines = dna_file.lines().peekable();
    loop {
        if let Some(line) = lines.next() {
            match line.first() {
                Some(AsciiChar::GreaterThan) => {
                    let dna = get_dna_name(line);
                    let mut text = AsciiString::new();

                    while let Some(line) = lines.peek() {
                        match line.first() {
                            Some(AsciiChar::GreaterThan) => break,
                            Some(_) => {
                                text.push_str(line);
                                lines.next();
                            }
                            None => {
                                break;
                            }
                        }
                    }

                    dnas.insert(dna, text);
                }
                _ => (),
            }
        } else {
            break;
        }
    }

    let mut index_file = File::create(args.index_file)?;
    dnas.serialize(&mut Serializer::new(&mut index_file))?;

    Ok(())
}
