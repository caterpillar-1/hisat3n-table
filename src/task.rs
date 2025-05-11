// we use term dna instead of chromosome in this module

use std::hint::cold_path;
use std::ops::Range;

use crate::alignment::Alignment;
use crate::{
    position::Position,
    ARGS,
};

pub struct Task2<'a> {
    pub dna_name: &'a [u8],
    pub alignments: Vec<Alignment<'a>>,
    pub position_range: Range<usize>,
}

pub type TaskResult<'a> = Option<Vec<Position<'a>>>;

pub struct TaskIter2<'a> {
    src: &'a [u8],
    current_position: usize,
}

impl<'a> TaskIter2<'a> {
    pub fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            current_position: 0,
        }
    }
}

impl<'a> Iterator for TaskIter2<'a> {
    type Item = Task2<'a>;

    #[inline(never)]
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_position >= self.src.len() {
            return None;
        }
        let chunk_start = self.current_position;
        let lines = memchr::memchr_iter(b'\n', &self.src[chunk_start..]);
        let mut line_start = chunk_start;
        let mut current_dna_name = &self.src[0..0];
        let mut current_chunk_beginning_pos = usize::MAX;
        let mut current_chunk_end_pos = usize::MAX;
        let mut n = 0;
        let mut chunk_end: usize = 0;
        let mut alignments: Vec<Alignment<'a>> = Vec::new();
        for line_feed_pos in lines {
            let actual_feed_pos = chunk_start + line_feed_pos;
            let line = &self.src[line_start..actual_feed_pos];
            line_start = actual_feed_pos + 1;
            let alignment = match Alignment::from_file(line) {
                Ok(alignment) => alignment,
                Err(_) => continue,
            };
            let seq_len: usize = alignment.bases.iter().map(|it| it.ref_pos).max().unwrap_or(alignment.sequence.len() as isize).try_into().unwrap();
            let pos = alignment.location as usize;

            if current_dna_name.len() == 0 {
                cold_path();
                current_dna_name = alignment.dna;
            } else if current_dna_name != alignment.dna {
                break; // 更换 ref 文件，放回当前行
            }
            if current_chunk_beginning_pos == usize::MAX {
                cold_path();
                current_chunk_beginning_pos = pos;
                current_chunk_end_pos = pos + seq_len + 1;
            } else if pos - current_chunk_beginning_pos > ARGS.ref_block_size && pos > current_chunk_end_pos {
                break; // 当前 chunk 过大，放回当前行
            }
            if n >= ARGS.align_block_size && pos > current_chunk_end_pos {
                break; // // 当前 chunk 过大，放回当前行
                // 注意必须保证各个段之间即使算上 location ~bases~ 延申之后还没有任何重叠！
                // 并且还不能紧密连接，因此这里是大于不是大于等于，因为下一个碱基可能影响上一个的 strand
            }
            current_chunk_end_pos = std::cmp::max(current_chunk_end_pos, pos + seq_len + 1); // 因此如果还在重叠区间内就不能分割

            n += 1;
            chunk_end = line_start;
            alignments.push(alignment);
            if line_start >= self.src.len() {
                break;
            }
        }
        self.current_position = chunk_end;
        if current_dna_name.len() == 0 || current_chunk_beginning_pos == usize::MAX {
            None
        } else {
            // eprintln!("fn {} position range {} - {}, size {}", str::from_utf8(&current_dna_name).unwrap(), current_chunk_beginning_pos, current_chunk_end_pos, current_chunk_end_pos - current_chunk_beginning_pos);
            Some(Task2 {
                dna_name: current_dna_name,
                alignments,
                position_range: current_chunk_beginning_pos .. current_chunk_end_pos,
            })
        }
    }
}
