use std::ops::Index;
use std::{cell::Cell, collections::HashSet, hash::Hash};
use std::hint::{likely, unlikely, cold_path};

use ascii::{AsAsciiStr, AsciiChar, AsciiStr, AsciiString};

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

impl Hash for UniqueID {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.read_name_id.hash(state);
    }
}

impl PartialEq for UniqueID {
    fn eq(&self, other: &Self) -> bool {
        self.read_name_id == other.read_name_id
    }
}

impl Eq for UniqueID {}

pub struct Position {
    pub dna: &'static AsciiStr,
    pub location: isize,
    pub strand: Option<AsciiChar>,
    pub converted_qualities: AsciiString,
    pub unconverted_qualities: AsciiString,
    pub unique_ids: HashSet<UniqueID>,
}

impl Position {
    pub fn new() -> Self {
        Self {
            dna: Default::default(),
            location: 0,
            strand: None,
            converted_qualities: AsciiString::new(),
            unconverted_qualities: AsciiString::new(),
            unique_ids: HashSet::default(),
        }
    }

    pub fn empty(&self) -> bool {
        self.converted_qualities.is_empty() && self.unconverted_qualities.is_empty()
    }

    fn append_read_name_id(&mut self, base: &PosQuality, a: &Alignment) -> bool {
        let id = UniqueID::new(a.read_name_id, base.converted, base.qual);
        let e = self.unique_ids.get_or_insert(id.clone());
        if e.converted != id.converted {
            if e.removed.get() {
                return false;
            }
            if e.converted != id.converted {
                e.removed.set(true);
                if e.converted {
                    let i = match self.converted_qualities.as_slice().iter().position(|ch| *ch == base.qual) {
                        Some(i) => i,
                        None => return false,
                    };
                    let _ = self.converted_qualities.remove(i);
                    return false;
                } else {
                    let i = match self.unconverted_qualities.as_slice().iter().position(|ch| *ch == base.qual) {
                        Some(i) => i,
                        None => return false,
                    };
                    let _ = self.unconverted_qualities.remove(i);
                    return false;
                }
            }
            false
        } else {
            true
        }
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
        Self { dna, location: location + 1, seq, last_base: None, last_pos: None }
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
                Some(AsciiChar::At) => { cold_path(); panic!("wrong ref parser impl") },
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
