// we use term dna instead of chromosome in this module

use core::slice;
use std::{
    assert_matches::assert_matches, collections::HashMap, fmt::Debug, fs::File, hint::{cold_path, likely, unlikely}, ops::Range, ptr::dangling, task
};

use ascii::{AsciiChar, AsciiStr, AsciiString, IntoAsciiString};
use named_tuple::named_tuple;

use crate::{
    ARGS,
    alignment::AlignmentIter,
    position::{Position, PositionIter},
};

use crate::DNAS;

pub struct Task {
    pub dna_name: &'static AsciiStr,
    pub alignments: AlignmentIter,
    pub alignment_count: usize,
    pub position_range: Range<usize>,
}

impl Debug for Task {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Task").field("dna", &self.dna_name).field("alignment_count", &self.alignment_count).field("position_range", &self.position_range).finish()
    }
}

pub type TaskResult = Option<Vec<Position>>;

pub type DnaIndex = HashMap<AsciiString, AsciiString>;

pub struct TaskIter {
    /* files */
    align_lines: Box<dyn Iterator<Item = &'static AsciiStr> + Send>,
    // the next align_line is stored only when next is returning the current Task and need
    // to place the peeked align_line
    peeked_align_line: Option<&'static AsciiStr>,
    current_dna: Option<&'static AsciiStr>,
}


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

    pub fn new(align_file: &'static AsciiStr) -> Self {
        Self {
            align_lines: Box::new(align_file.lines()),
            peeked_align_line: None,
            current_dna: None,
        }
    }
}

impl Iterator for TaskIter {
    type Item = Task;

    fn next(&mut self) -> Option<Self::Item> {
        // we construct a Task from self.dna and these two components
        // 3 components
        let mut task_dna: Option<&'static AsciiStr> = None;
        let mut align_text: Option<Range<*const AsciiChar>> = None;

        // Alignment (line) s in this chunk
        let mut align_count = 0usize;
        // Position (non-blank char) s in this chunk
        let mut position_range: Range<isize> = 0..0;

        fn task_from_components(task_dna: Option<&'static AsciiStr>, align_text: Option<Range<*const AsciiChar>>, align_count: usize, position_range: Range<isize>) -> Task {
            let align_iter = AlignmentIter::new(unsafe {slice::from_ptr_range(align_text.unwrap()).into() });
            let task = Task {
                dna_name: task_dna.unwrap(),
                alignments: align_iter,
                alignment_count: align_count,
                position_range: (position_range.start as usize .. position_range.end as usize),
            };
            // eprintln!("{:?}", task);
            return task;
        }

        loop {
            // this condition is satisfied when next lines of code fills self.current_dna with some content
            if align_count >= ARGS.align_block_size || position_range.len() >= ARGS.ref_block_size {
                // soft split
                cold_path();
                return Some(task_from_components(task_dna, align_text, align_count, position_range));
            }

            let align_line = match self.peeked_align_line {
                Some(l) => {
                    // process first line of alignment chunks rather than the first Alignment of the file
                    self.peeked_align_line = None;
                    l
                }
                None => {
                    match self.align_lines.next() {
                        Some(l) => l,
                        None => {
                            cold_path();
                            match align_text {
                                Some(_) => {
                                    return Some(task_from_components(task_dna, align_text, align_count, position_range));
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

            if task_dna.is_none() {
                cold_path();

                task_dna = Some(align_dna);
                align_text = Some(align_line.as_slice().as_ptr_range());
                align_count = 1;
                position_range = align_location..align_location + 1;
            } else {
                if task_dna.unwrap() != align_dna {
                    // hard split                    
                    self.peeked_align_line = Some(align_line);

                    return Some(task_from_components(task_dna, align_text, align_count, position_range));
                } else {
                    align_text.as_mut().unwrap().end = align_line.as_slice().as_ptr_range().end;
                    align_count += 1;
                    position_range.end = align_location + 1;
                }
            }
        }
    }
}
