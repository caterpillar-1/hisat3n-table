use ascii::{AsciiChar, AsciiStr, AsciiString};
use std::hint::cold_path;
use std::{cell::Cell, hash::Hash};

use crate::{
    ARGS,
    alignment::{Alignment, PosQuality},
};

#[derive(Default, Debug, Clone)]
pub struct UniqueID {
    pub read_name_id: usize,
    pub converted: bool,
    pub quality: AsciiChar,
    pub removed: Cell<bool>,
}

impl UniqueID {
    fn new(read_name_id: usize, converted: bool, quality: AsciiChar) -> Self {
        Self {
            read_name_id,
            converted,
            quality,
            removed: Cell::new(false),
        }
    }
}

pub struct Position {
    pub dna: &'static AsciiStr,
    pub location: isize,
    pub strand: Option<AsciiChar>,
    pub converted_qualities: AsciiString,
    pub unconverted_qualities: AsciiString,
    pub unique_ids: Vec<UniqueID>,
}

impl Position {
    pub fn new() -> Self {
        Self {
            dna: Default::default(),
            location: -1,
            strand: None,
            converted_qualities: AsciiString::new(),
            unconverted_qualities: AsciiString::new(),
            unique_ids: Vec::default(),
        }
    }

    pub fn empty(&self) -> bool {
        self.converted_qualities.is_empty() && self.unconverted_qualities.is_empty()
    }

    fn search_read_name_id(&mut self, read_name_id: usize, start: usize, end: usize) -> usize {
        if self.unique_ids.is_empty() {
            return 0;
        }
        if start <= end {
            let middle = (start + end) / 2;
            if self.unique_ids[middle].read_name_id == read_name_id {
                return middle;
            }
            if self.unique_ids[middle].read_name_id > read_name_id {
                return self.search_read_name_id(read_name_id, start, middle - 1);
            }
            return self.search_read_name_id(read_name_id, middle + 1, end);
        }
        return start; // return the bigger one
    }

    fn append_read_name_id(&mut self, base: &PosQuality, a: &Alignment) -> bool {
        let id_count = self.unique_ids.len();
        if id_count == 0 || a.read_name_id > self.unique_ids.last().unwrap().read_name_id {
            self.unique_ids
                .push(UniqueID::new(a.read_name_id, base.converted, base.qual));
            return true;
        }
        let index = self.search_read_name_id(a.read_name_id, 0, id_count);
        if self.unique_ids[index].read_name_id == a.read_name_id {
            // if the new base is consistent with exist base's conversion status, ignore
            // otherwise, delete the exist conversion status
            if self.unique_ids[index].removed.get() {
                return false;
            }
            if self.unique_ids[index].converted != base.converted {
                self.unique_ids[index].removed.set(true);
                if self.unique_ids[index].converted {
                    for i in 0..self.converted_qualities.len() {
                        if self.converted_qualities[i] == base.qual {
                            self.converted_qualities.remove(i);
                            return false;
                        }
                    }
                } else {
                    for i in 0..self.converted_qualities.len() {
                        if self.unconverted_qualities[i] == base.qual {
                            self.unconverted_qualities.remove(i);
                            return false;
                        }
                    }
                }
            }
            return false;
        } else {
            self.unique_ids.insert(
                index,
                UniqueID::new(a.read_name_id, base.converted, base.qual),
            );
            return true;
        }

        // let id = UniqueID::new(a.read_name_id, base.converted, base.qual);
        // let e = self.unique_ids.get_or_insert(id.clone());
        // if e.converted != id.converted {
        //     if e.removed.get() {
        //         return false;
        //     }
        //     if e.converted != id.converted {
        //         e.removed.set(true);
        //         if e.converted {
        //             let i = match self
        //                 .converted_qualities
        //                 .as_slice()
        //                 .iter()
        //                 .position(|ch| *ch == base.qual)
        //             {
        //                 Some(i) => i,
        //                 None => return false,
        //             };
        //             let _ = self.converted_qualities.remove(i);
        //             return false;
        //         } else {
        //             let i = match self
        //                 .unconverted_qualities
        //                 .as_slice()
        //                 .iter()
        //                 .position(|ch| *ch == base.qual)
        //             {
        //                 Some(i) => i,
        //                 None => return false,
        //             };
        //             let _ = self.unconverted_qualities.remove(i);
        //             return false;
        //         }
        //     }
        //     false
        // } else {
        //     true
        // }
    }

    pub fn append_base(&mut self, input: &PosQuality, a: &Alignment) {
        if self.append_read_name_id(input, a) {
            if input.converted {
                self.converted_qualities.push(input.qual);
            } else {
                self.unconverted_qualities.push(input.qual);
            }
        }
    }
}

fn search_read_name_id(read_name_id: usize, arg: i32, id_count: usize) -> usize {
    todo!()
}

pub struct PositionIter {
    dna: &'static AsciiStr,
    location: isize,
    seq: &'static AsciiStr, // file text, can contain blank spaces
    last_base: Option<AsciiChar>,
    last_pos: Option<Position>,
}

impl PositionIter {
    // location: pass in 0-based location
    //           internally use 1-based location
    pub fn new(dna: &'static AsciiStr, location: isize, seq: &'static AsciiStr) -> Self {
        Self {
            dna,
            location: location + 1,
            seq,
            last_base: None,
            last_pos: None,
        }
    }
}

impl Iterator for PositionIter {
    type Item = Position;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.seq.first() {
                Some(AsciiChar::CarriageReturn | AsciiChar::LineFeed | AsciiChar::Space) => {
                    self.seq = &self.seq[1..];
                    continue;
                }
                Some(AsciiChar::At) => {
                    cold_path();
                    panic!("wrong ref parser impl")
                }
                Some(ch) => {
                    let ch = ch.to_ascii_uppercase();
                    let mut p = Position::new();
                    p.dna = self.dna;
                    p.location = self.location;
                    self.location += 1;
                    if ARGS.cg_only {
                        if matches!(self.last_base, Some(AsciiChar::C))
                            && matches!(ch, AsciiChar::G)
                        {
                            self.last_pos.as_mut().unwrap().strand = Some(AsciiChar::Plus);
                            p.strand = Some(AsciiChar::Minus);
                        }
                    } else {
                        if ch == ARGS.base_change.0.0 {
                            p.strand = Some(AsciiChar::Plus);
                        } else if ch == ARGS.base_change.1.0 {
                            p.strand = Some(AsciiChar::Minus);
                        }
                    }
                    self.last_base = Some(ch);
                    self.seq = &self.seq[1..];
                    let ret = std::mem::replace(&mut self.last_pos, Some(p));
                    if ret.is_none() {
                        continue;
                    } else {
                        return ret;
                    }
                }
                None => {
                    cold_path();
                    match &self.last_pos {
                        Some(_) => {
                            let _ = std::mem::replace(&mut self.last_base, None);
                            return std::mem::replace(&mut self.last_pos, None);
                        }
                        None => return None,
                    }
                }
            }
        }
    }
}
