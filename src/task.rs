// we use term dna instead of chromosome in this module

use core::slice;
use std::{assert_matches::assert_matches, collections::HashMap, hint::{cold_path, likely, unlikely}, ops::Range};

use ascii::{AsciiChar, AsciiStr, AsciiString, IntoAsciiString};

use crate::{alignment::AlignmentIter, position::{Position, PositionIter}, ARGS};


pub struct Task {
    pub dna: &'static AsciiStr,
    pub alignments: AlignmentIter,
    pub positions: PositionIter,
}

pub type TaskResult = Option<Vec<Position>>;

pub struct TaskIter {
    /* files */
    align_lines: Box<dyn Iterator<Item = &'static AsciiStr>>,
    dnas: HashMap<AsciiString, &'static AsciiStr>,
    /* states */
    // it acts as a cache to the dna entry in HashMap dnas,
    // we don't need to scan the same dna text from start for every alignment record 
    // when the alignment lines are sorted based on dna name and location in the dna.
    // (name, remaining seq text (lines), peeked dna line, location of first bp in remaining seq text)
    current_dna: Option<(&'static AsciiStr, Box<dyn Iterator<Item = &'static AsciiStr>>, Option<&'static AsciiStr>, isize)>,
    // the next align_line is stored only when next is returning the current Task and need
    // to place the peeked align_line
    peeked_align_line: Option<&'static AsciiStr>,
}

unsafe impl Send for TaskIter {}

impl TaskIter {
    fn get_dna_location(align_line: &AsciiStr) -> Option<(&AsciiStr, isize)> {
        let mut res: (&AsciiStr, isize) = Default::default();
        let mut parts = align_line.split(AsciiChar::Tab);
        match parts.nth(2) {
            Some(chr) => {
                res.0 = chr;
            },
            None => return None,
        }
        match parts.next() {
            Some(chr) => {
                if chr.as_str() == "*" {
                    return None
                } else {
                    res.1 = chr.as_str().parse().unwrap();
                }
            },
            None => return None,
        }
        Some(res)
    }

    fn get_dna_name(info_line: &AsciiStr) -> AsciiString {
        assert_eq!(info_line.first().unwrap(), AsciiChar::GreaterThan);
        let info_line = &info_line[1..];
        info_line.as_str().split_ascii_whitespace().next().unwrap().into_ascii_string().unwrap()
    }

    pub fn new(align_file: &'static AsciiStr, dna_file: &'static AsciiStr) -> Self {
        let mut dnas = HashMap::new();
        let mut dna_seq_start: *const AsciiChar;
        let mut dna_seq_end: *const AsciiChar;
        let mut lines = dna_file.lines().peekable();
        loop {
            if let Some(line) = lines.next() {
                match line.first() {
                    Some(AsciiChar::GreaterThan) => {
                        let info_line = line;
                        let dna = Self::get_dna_name(info_line);

                        // set seq to empty slice
                        dna_seq_start = info_line.as_slice().as_ptr_range().end;
                        dna_seq_end = dna_seq_start;

                        while let Some(line) = lines.peek() {
                            match line.first() {
                                Some(AsciiChar::GreaterThan) => break,
                                Some(_) => { 
                                    dna_seq_end = line.as_slice().as_ptr_range().end;
                                    lines.next();
                                }
                                None => {
                                    break;
                                },
                            }
                        }

                        let text_lines: &'static AsciiStr = (unsafe { std::slice::from_ptr_range(dna_seq_start..dna_seq_end) }).into();

                        dnas.insert(dna, text_lines);
                    },
                    _ => (),
                }
            } else {
                break;
            }
        }
        
        Self { align_lines: Box::new(align_file.lines()), dnas, current_dna: None, peeked_align_line: None }
    }
}

impl Iterator for TaskIter {
    type Item = Task;

    /// we perform task partitioning here, scan through two files quickly
    fn next(&mut self) -> Option<Self::Item> {
        // we construct a Task from self.dna and these two components
        // 3 components
        let mut align_t: Option<Range<*const AsciiChar>> = None;
        let mut dna_t: Option<Range<*const AsciiChar>> = None; // dna text slice
        let mut dna_l: Option<isize> = None; // dna seq location

        // Alignment (line) s in this chunk
        let mut align_count = 0usize;
        // Position (non-blank char) s in this chunk
        let mut dna_count = 0usize;

        let advance_dna_location_to_align_location = |align_location: isize, dna_lines: &mut Box<dyn Iterator<Item = &AsciiStr>>, dna_location: &mut isize, dna_count: &mut usize| -> (Option<Range<*const AsciiChar>>, Option<isize>) {
            // seek to align_location
            while *dna_location < align_location {
                let dna_line = dna_lines.next().expect("alignment refer to a non-exist location to chromosome");
                *dna_count += 1;
                let next_dna_location = *dna_location + dna_line.len() as isize;

                if *dna_location <= align_location && align_location < next_dna_location {
                    // the start dna line, set to dna_t
                    let ret = (Some(dna_line.as_slice().as_ptr_range()), Some(*dna_location));
                    *dna_location = next_dna_location;
                    return ret;
                }
                *dna_location = next_dna_location;
            }
            (None, None)
        };

        // the loop basically tries to partition the alignment file into chunks by different reference dna name (hard split)
        // or when the chunk is too lengthy (soft split)
        loop {
            // put in the loop to process no such chromosome condition
            if unlikely(self.peeked_align_line.is_some()) {
                // process the first line in the new chunk (peeked by processing last chunk)
                // this line is always in self.peeked_align_line

                // 3 components must be None
                assert_matches!(align_t, None);

                let (align_dna, align_location) = match Self::get_dna_location(self.peeked_align_line.unwrap()) {
                    Some(x) => x,
                    None => {
                        // strange record
                        self.peeked_align_line = None;
                        continue;
                    },
                };

                match &mut self.current_dna {
                    Some((current_dna, dna_lines, peeked_dna_line, dna_location)) => {
                        // cache valid

                        if align_dna == *current_dna {
                            // cache hit

                            match peeked_dna_line {
                                Some(dna_line) => {
                                    let prev_dna_location = *dna_location - dna_line.len() as isize;

                                    if align_location < prev_dna_location {
                                        cold_path();
                                        panic!("")
                                    } else if prev_dna_location <= align_location && align_location < *dna_location {
                                        // this chunk starts on the peeked dna line
                                        dna_l = Some(prev_dna_location);
                                        dna_t = Some(dna_line.as_slice().as_ptr_range());

                                        // no peeked dna line now
                                        *peeked_dna_line = None;

                                        align_t = Some(self.peeked_align_line.unwrap().as_slice().as_ptr_range());

                                        self.peeked_align_line = None;
                                        continue;
                                    } else {
                                        // we can safely discard the peeked dna line and process as normal
                                        *peeked_dna_line = None;
                                        continue;
                                    }

                                }
                                None => {
                                    // normally
                                    (dna_t, dna_l) = advance_dna_location_to_align_location(align_location, dna_lines, dna_location, &mut dna_count);

                                    assert_matches!(dna_t, Some(_));

                                    align_t = Some(self.peeked_align_line.unwrap().as_slice().as_ptr_range());

                                    self.peeked_align_line = None;

                                    continue;
                                }
                            }
                        } else {
                            // cache invalid

                            // we can simply set cache to empty and enter next loop
                            self.current_dna = None;
                            continue;
                        }
                    }
                    None => {
                        // cache miss
                        // temporary cache of self.current_dna
                        let (current_dna, mut dna_lines, mut dna_location) = {
                            let dna_text = match self.dnas.get(align_dna) {
                                Some(t) => *t,
                                None => {
                                    eprintln!("no such chromosome: {}", align_dna);
                                    self.peeked_align_line = None;
                                    continue;
                                },
                            };

                            let b: Box<dyn Iterator<Item = &AsciiStr>> = Box::new(dna_text.lines());

                            // we just set current_dna to its complete sequence
                            (align_dna, b, 0)
                        };

                        (dna_t, dna_l) = advance_dna_location_to_align_location(align_location, &mut dna_lines, &mut dna_location, &mut dna_count);

                        // now dna_t must be set
                        assert_matches!(dna_t, Some(_));

                        // cache write through
                        self.current_dna = Some((current_dna, dna_lines, None, dna_location));

                        // set alignment range to first line and set no peeked lines
                        align_t = Some(
                            std::mem::replace(
                                &mut self.peeked_align_line, 
                                None
                            ).unwrap().as_slice().as_ptr_range()
                        );
                        continue;
                    }
                }
            }
            
            align_count += 1;
            match self.align_lines.next() {
                Some(align_line) => {
                    // a new line in alignment file
                    if unlikely(align_line.first().is_none_or(|ch| ch == AsciiChar::At)) {
                        // this line is empty or is header
                        continue;
                    }

                    let align_dna = Self::get_dna_location(align_line);
                    if align_dna.is_none() {
                        // strange record
                        continue;
                    }
                    let (align_dna, align_location) = align_dna.unwrap();

                    // cache lookup
                    match &mut self.current_dna {
                        Some((current_dna, dna_lines, peeked_dna_line, dna_location)) => {
                            // cache hit

                            if likely(align_dna == current_dna.as_slice()) {
                                // cache valid (in the same chunk)

                                // process peeked
                                if let Some(dna_line) = peeked_dna_line {
                                    let prev_dna_location = *dna_location - dna_line.len() as isize;

                                    assert!(align_location < prev_dna_location, "the input file is not sorted.");

                                    if prev_dna_location <= align_location && align_location < *dna_location {
                                        match &dna_t {
                                            Some(_) => {
                                                dna_t.as_mut().unwrap().end = dna_line.as_slice().as_ptr_range().end;
                                            }
                                            None => {
                                                dna_t = Some(dna_line.as_slice().as_ptr_range());
                                                dna_l = Some(prev_dna_location);
                                            }
                                        }
                                    } else {
                                        *peeked_dna_line = None;
                                    }

                                    continue;
                                }

                                // no peeked line

                                assert!(align_location < *dna_location, "the input file is not sorted.");

                                // check whether the chunk size is too large (a soft split)
                                if align_count >= ARGS.align_block_size || dna_count >= ARGS.ref_block_size {
                                    // reasoning is the same as hard split
                                    self.peeked_align_line = Some(align_line);
                                    break;
                                } else {
                                    // update ranges!

                                    let (last_dna_t, last_dna_l) = advance_dna_location_to_align_location(align_location, dna_lines, dna_location, &mut dna_count);

                                    match &dna_t {
                                        Some(_) => {
                                            dna_t.as_mut().unwrap().end = last_dna_t.unwrap().end;
                                        }
                                        None => {
                                            // this should not happen
                                            eprintln!("unexpected branch");
                                            cold_path();
                                            (dna_t, dna_l) = (last_dna_t, last_dna_l);
                                        }
                                    }
                                }
                            } else {
                                // cache invalid (a hard split by different dna name)

                                // save the peeked line
                                self.peeked_align_line = Some(align_line);

                                // we postpone setting self.current_dna to next call in checking self.peeked_align_line

                                // now that former loops have maintained align_r and dna_r,
                                // we break the loop to return current Task
                                break;
                            }
                        }
                        None => {
                            cold_path();
                            // cache miss on the first record line of the alignment file,
                            // or after the first record line referring to an nonexist dna

                            // we process this case like a hard split

                            // 3 components must be None
                            assert_matches!(align_t, None);
                            self.peeked_align_line = Some(align_line);
                            continue;
                        }
                    }
                }
                None => {
                    cold_path();
                    if matches!(align_t, Some(_)) {
                        // remaining chunk, return it
                        break;
                    } else {
                        return None
                    }
                }
            }
        }

        // when the loop exits, we found a new chunk
        assert_matches!(align_t, Some(_));
        assert_matches!(dna_t, Some(_));
        assert_matches!(dna_l, Some(_));

        return Some(Task { 
            dna: self.current_dna.as_mut().unwrap().0,
            alignments: AlignmentIter::new(unsafe {
                slice::from_ptr_range(dna_t.unwrap())
            }.into()),
            positions: PositionIter::new(self.current_dna.as_mut().unwrap().0, dna_l.unwrap(), unsafe { slice::from_ptr_range(align_t.unwrap()) }.into()),  
        })
    }
}
