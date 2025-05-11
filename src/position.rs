use std::{collections::BTreeMap, hint::cold_path};

use crate::{
    ARGS,
    alignment::{Alignment, PosQuality},
};

#[derive(Default, Debug, Clone)]
pub struct UniqueID {
    pub read_name_id: u64,
    pub converted: bool,
    pub quality: u8,
    pub removed: bool,
}

impl UniqueID {
    fn new(read_name_id: u64, converted: bool, quality: u8) -> Self {
        Self {
            read_name_id,
            converted,
            quality,
            removed: false,
        }
    }
}

pub struct Position<'a> {
    pub dna: &'a [u8],
    pub location: isize,
    pub strand: Option<u8>,
    pub converted_qualities: Vec<u8>,
    pub unconverted_qualities: Vec<u8>,
    pub unique_ids: BTreeMap<u64, UniqueID>,
}

impl<'a> Position<'a> {
    pub fn new(dna: &'a [u8], location: isize) -> Self {
        Self {
            dna,
            location,
            strand: None,
            converted_qualities: Vec::new(),
            unconverted_qualities: Vec::new(),
            unique_ids: BTreeMap::new(),
        }
    }

    fn append_read_name_id(&mut self, in_base: &PosQuality, in_align: &Alignment) -> bool {
        match self.unique_ids.entry(in_align.read_name_id) {
            std::collections::btree_map::Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(UniqueID::new(in_align.read_name_id, in_base.converted, in_base.qual));
                return true;
            },
            std::collections::btree_map::Entry::Occupied(mut occupied_entry) => {
                let ent = occupied_entry.get_mut();
                // if the new base is consistent with exist base's conversion status, ignore
                // otherwise, delete the exist conversion status
                if ent.removed {
                    return false;
                }
                if ent.converted != in_base.converted {
                    ent.removed = true;
                    if ent.converted {
                        for i in 0..self.converted_qualities.len() {
                            if self.converted_qualities[i] == in_base.qual {
                                self.converted_qualities.remove(i);
                                return false;
                            }
                        }
                    } else {
                        for i in 0..self.converted_qualities.len() {
                            if self.unconverted_qualities[i] == in_base.qual {
                                self.unconverted_qualities.remove(i);
                                return false;
                            }
                        }
                    }
                }
                return false;
            },
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

pub fn fill_positions<'a>(positions: &mut Vec<Position<'a>>, text: &'a [u8], dna: &'a [u8],
                          start_pos: usize, end_pos: usize) {
    positions.reserve(end_pos - start_pos);
    let mut last_base = 0u8;
    for i in start_pos..end_pos {
        if i >= text.len() {
            break;
        }
        let ch = text[i - 1];
        assert!(ch.is_ascii_alphabetic());
        let mut p = Position::new(dna, i as isize);
        if ARGS.cg_only {
            cold_path();
            if last_base == b'C' && ch == b'G' {
                positions.last_mut().unwrap().strand = Some(b'+');
                p.strand = Some(b'-');
            }
        } else {
            if ch == ARGS.base_change.0.0 {
                p.strand = Some(b'+');
            } else if ch == ARGS.base_change.0.1 {
                p.strand = Some(b'-');
            }
        }
        positions.push(p);
        last_base = ch;
    }
}
