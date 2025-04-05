// we use term dna instead of chromosome in this module

use core::slice;
use std::{
    assert_matches::assert_matches, collections::HashMap, fmt::Debug, hint::{cold_path, likely, unlikely}, ops::Range, ptr::dangling, task
};

use ascii::{AsciiChar, AsciiStr, AsciiString, IntoAsciiString};
use named_tuple::named_tuple;

use crate::{
    ARGS,
    alignment::AlignmentIter,
    position::{Position, PositionIter},
};

pub struct Task {
    pub dna: &'static AsciiStr,
    pub alignments: AlignmentIter,
    pub alignment_count: usize,
    pub positions: PositionIter,
    pub position_line_count: usize,
}

impl Debug for Task {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Task").field("dna", &self.dna).field("alignment_count", &self.alignment_count).field("position_line_count", &self.position_line_count).finish()
    }
}

pub type TaskResult = Option<Vec<Position>>;

named_tuple!(
    struct DnaCacheResp {
        hit: bool,
        // indicate the start position of line in the dna
        location: isize,
        scan_line_count: usize,
        line: &'static AsciiStr,
    }
);

named_tuple!(
    struct DnaCacheLine {
        current_dna: &'static AsciiStr,
        dna_lines: Box<dyn Iterator<Item = &'static AsciiStr>>,
        peeked_dna_line: Option<&'static AsciiStr>,
        location: isize,
    }
);

pub struct TaskIter {
    /* files */
    align_lines: Box<dyn Iterator<Item = &'static AsciiStr>>,
    // the next align_line is stored only when next is returning the current Task and need
    // to place the peeked align_line
    peeked_align_line: Option<&'static AsciiStr>,
    dnas: HashMap<AsciiString, &'static AsciiStr>,
    /* states */
    // it acts as a cache to the dna entry in HashMap dnas,
    // we don't need to scan the same dna text from start for every alignment record
    // when the alignment lines are sorted based on dna name and location in the dna.
    // (name, remaining seq text (lines), peeked dna line, location of first bp in remaining seq text)
    current_dna: Option<DnaCacheLine>,
}

unsafe impl Send for TaskIter {}

impl TaskIter {
    fn get_dna_location(align_line: &AsciiStr) -> Option<(&AsciiStr, isize)> {
        if align_line.is_empty() || align_line.first().unwrap() == '@' {
            return None;
        }
        let mut res: (&AsciiStr, isize) = Default::default();
        let mut parts = align_line.split(AsciiChar::Tab);
        match parts.nth(2) {
            Some(chr) => {
                res.0 = chr;
            }
            None => return None,
        }
        match parts.next() {
            Some(chr) => {
                if chr.as_str() == "*" {
                    return None;
                } else {
                    res.1 = chr.as_str().parse().unwrap();
                }
            }
            None => return None,
        }
        Some(res)
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
                                }
                            }
                        }

                        let text_lines: &'static AsciiStr =
                            (unsafe { std::slice::from_ptr_range(dna_seq_start..dna_seq_end) })
                                .into();

                        dnas.insert(dna, text_lines);
                    }
                    _ => (),
                }
            } else {
                break;
            }
        }

        // assert!(
        //     dnas.values().all(|t| t
        //         .chars()
        //         .all(|ch| ch.is_ascii_alphabetic() || ch.is_ascii_whitespace())),
        //     "wrong dnas parser impl."
        // );

        Self {
            align_lines: Box::new(align_file.lines()),
            dnas,
            current_dna: None,
            peeked_align_line: None,
        }
    }

    fn current_dna_update(
        &mut self,
        align_dna: &'static AsciiStr,
        align_location: isize,
    ) -> Option<DnaCacheResp> {
        let mut cache_hit = false;
        let mut scan_line_count = 0;
        // cache lookup
        if self.current_dna.is_some()
            && self.current_dna.as_ref().unwrap().field_values().0 == align_dna
        {
            // cache hit
            cache_hit = true;
        } else {
            cold_path();
            // cache invalid or miss
            // refill
            match self.dnas.get(align_dna) {
                Some(dna_text) => {
                    // we should navigate to align_location first to reuse the code below
                    let mut dna_lines = Box::new(dna_text.lines());
                    let mut dna_location = 0;
                    let mut peeked_dna_line = None;

                    loop {
                        let dna_line = match dna_lines.next() {
                            Some(l) => l,
                            None => {
                                eprintln!(
                                    "alignment refer to chromosome {} at position {}, which is out of range",
                                    align_dna, align_location
                                );
                                self.current_dna = None;
                                return None
                            }
                        };
                        let next_dna_location = dna_location + dna_line.len() as isize;
                        if (dna_location..next_dna_location).contains(&align_location) {
                            peeked_dna_line = Some(dna_line);
                            let ret = DnaCacheResp::new(cache_hit, dna_location, scan_line_count, dna_line);
                            dna_location = next_dna_location;
                            self.current_dna = Some(DnaCacheLine::new(
                                align_dna,
                                dna_lines,
                                peeked_dna_line,
                                dna_location,
                            ));
                            return Some(ret);
                        }
                        dna_location = next_dna_location;
                    }
                }
                None => {
                    cold_path();
                    return None;
                }
            }
        }

        let DnaCacheLine((dna_dna, dna_lines, peeked_dna_line, dna_location)) =
            self.current_dna.as_mut().unwrap();
        assert_eq!(*dna_dna, align_dna);

        match peeked_dna_line {
            Some(dna_line) => {
                let prev_dna_location = *dna_location - dna_line.len() as isize;
                assert!(
                    align_location >= prev_dna_location,
                    "alignment file is not sorted"
                );

                if (prev_dna_location..*dna_location).contains(&align_location) {
                    return Some(DnaCacheResp::new(true, prev_dna_location, scan_line_count, dna_line));
                } else {
                    *peeked_dna_line = None;
                }
            }
            None => (),
        }

        assert!(align_location >= *dna_location);
        assert!(peeked_dna_line.is_none());

        loop {
            let dna_line = match dna_lines.next() {
                Some(l) => l,
                None => {
                    eprintln!(
                        "alignment refer to chromosome {} at position {}, which is out of range",
                        align_dna, align_location
                    );
                    return None;
                }
            };
            scan_line_count += 1;
            let next_dna_location = *dna_location + dna_line.len() as isize;
            if (*dna_location..next_dna_location).contains(&align_location) {
                *peeked_dna_line = Some(dna_line);
                let ret = DnaCacheResp::new(cache_hit, *dna_location, scan_line_count, dna_line);
                *dna_location = next_dna_location;
                return Some(ret);
            }
            *dna_location = next_dna_location;
        }
    }
}

impl Iterator for TaskIter {
    type Item = Task;

    /// we perform task partitioning here, scan through two files quickly
    fn next(&mut self) -> Option<Self::Item> {
        // we construct a Task from self.dna and these two components
        // 3 components
        let mut align_text: Option<Range<*const AsciiChar>> = None;
        let mut dna_text: Option<Range<*const AsciiChar>> = None; // dna text slice
        let mut dna_location: Option<isize> = None; // dna seq location

        // Alignment (line) s in this chunk
        let mut align_count = 0usize;
        // Position (non-blank char) s in this chunk
        let mut dna_count = 0usize;

        fn task_from_components(
            current_dna: &'static AsciiStr,
            align_text: Option<Range<*const AsciiChar>>,
            dna_text: Option<Range<*const AsciiChar>>,
            dna_location: Option<isize>,
            align_count: usize,
            dna_count: usize
        ) -> Task {
            let t = Task {
                dna: current_dna,
                alignments: AlignmentIter::new(
                    unsafe { slice::from_ptr_range(align_text.unwrap()) }.into(),
                ),
                positions: PositionIter::new(
                    current_dna,
                    dna_location.unwrap(),
                    unsafe { slice::from_ptr_range(dna_text.unwrap()) }.into(),
                ),
                alignment_count: align_count,
                position_line_count: dna_count,
            };
            eprintln!("{:?}", t);
            t
        }

        loop {
            // this condition is satisfied when next lines of code fills self.current_dna with some content
            if align_count >= ARGS.align_block_size || dna_count >= ARGS.ref_block_size {
                cold_path();
                return Some(task_from_components(
                    self.current_dna.as_ref().unwrap().current_dna(),
                    align_text,
                    dna_text,
                    dna_location,
                    align_count,
                    dna_count
                ));
            }

            let align_line = match self.peeked_align_line {
                Some(l) => {
                    // process first line of alignment chunks rather than the first Alignment of the file
                    let t = l;                    
                    self.peeked_align_line = None;
                    t
                }
                None => {
                    match self.align_lines.next() {
                        Some(l) => l,
                        None => {
                            cold_path();
                            match align_text {
                                Some(_) => {
                                    return Some(task_from_components(
                                        self.current_dna.as_ref().unwrap().current_dna(),
                                        align_text,
                                        dna_text,
                                        dna_location,
                                        align_count,
                                        dna_count
                                    ));
                                }
                                None => {
                                    eprintln!("Assignment Done!");
                                    return None;
                                }
                            }
                        }
                    }
                }
            };

            let (align_dna, align_location) = match Self::get_dna_location(align_line) {
                Some(t) => t,
                None => continue,
            };

            let DnaCacheResp((dna_hit, new_dna_location, dna_scan_line_count, dna_line)) =
                match self.current_dna_update(align_dna, align_location) {
                    Some(t) => t,
                    None => {
                        cold_path();
                        continue;
                    }
                };

            if dna_hit {
                // cache hit, add align_line and ..dna_line.end to this chunk

                if align_text.is_none() {
                    // we are processing the first line of a new chunk rather than the first in the file
                    cold_path();
                    align_text = Some(align_line.as_slice().as_ptr_range());
                    dna_text = Some(dna_line.as_slice().as_ptr_range());
                    dna_location = Some(new_dna_location);

                    align_count += 1;
                    dna_count += 1; // hack, dna_scan_line_count is 0 in this case
                } else {
                    // we are in the middle of the text chunk

                    align_text.as_mut().unwrap().end = align_line.as_slice().as_ptr_range().end;
                    dna_text.as_mut().unwrap().end = dna_line.as_slice().as_ptr_range().end;

                    align_count += 1;
                    dna_count += dna_scan_line_count;
                }
            } else {
                // cache miss, hard split
                cold_path();
                if align_text.is_none() {
                    // the most rare case: its the first alignment line
                    cold_path();
                    align_text = Some(align_line.as_slice().as_ptr_range());
                    dna_text = Some(dna_line.as_slice().as_ptr_range());
                    dna_location = Some(new_dna_location);

                    align_count += 1;
                    dna_count += dna_scan_line_count;
                } else {
                    // on the boundry of areas split by different dna
                    self.peeked_align_line = Some(align_line);
                    return Some(task_from_components(
                        self.current_dna.as_ref().unwrap().current_dna(),
                        align_text,
                        dna_text,
                        dna_location,
                        align_count,
                        dna_count
                    ));
                }
            }
        }
    }
}
