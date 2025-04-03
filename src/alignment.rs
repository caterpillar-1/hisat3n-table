use ascii::{AsciiChar, AsciiStr, AsciiString, IntoAsciiString, ToAsciiChar};

use crate::ARGS;

#[derive(Debug, Default)]
pub struct PosQuality {
    pub read_pos: isize,
    pub ref_pos: isize,
    pub qual: AsciiChar,
    pub converted: bool,
    pub remove: bool,
}

impl PosQuality {
    pub fn new(pos: isize) -> Self {
        PosQuality {
            read_pos: pos,
            ref_pos: pos,
            remove: true,
            ..Default::default()
        }
    }

    pub fn set_qual(&mut self, qual: AsciiChar, converted: bool) {
        self.qual = qual;
        self.converted = converted;
        self.remove = false;
    }
}

pub struct Alignment {
    pub dna: &'static AsciiStr,
    pub location: isize,
    pub mate_location: isize,
    pub flag: i32,
    pub mapped: bool,
    pub strand: AsciiChar,
    pub sequence: &'static AsciiStr,
    pub quality: &'static AsciiStr,
    pub unique: bool,
    pub map_q: &'static AsciiStr,
    pub nh: i32,
    pub bases: Vec<PosQuality>,
    pub cigar: &'static AsciiStr,
    pub md: &'static AsciiStr,
    pub read_name_id: usize,
    pub sequence_covered_length: usize,
    pub overlap: bool,
    pub paired: bool,
}

impl TryFrom<&'static AsciiStr> for Alignment {
    type Error = ();

    fn try_from(s: &'static AsciiStr) -> Result<Self, Self::Error> {
        let mut a = Self::new();

        let mut s = s.split('\t'.to_ascii_char().map_err(|_| ())?);
        // 0
        a.read_name_id = Self::hash_name(s.next().ok_or(())?);
        // 1
        a.flag = s.next().ok_or(())?.as_str().parse().map_err(|_| ())?;
        // 2
        a.dna = s.next().ok_or(())?;
        // 3
        a.location = s.next().ok_or(())?.as_str().parse().map_err(|_| ())?;
        // 4
        a.map_q = s.next().ok_or(())?;
        a.unique = a.map_q != "1".into_ascii_string().map_err(|_| ())?;
        // 5
        a.cigar = s.next().ok_or(())?;
        // 6
        s.next().ok_or(())?;
        // 7
        a.mate_location = s.next().ok_or(())?.as_str().parse().map_err(|_| ())?;
        // 8
        s.next().ok_or(())?;
        // 9
        a.sequence = s.next().ok_or(())?;
        // 10
        a.quality = s.next().ok_or(())?;
        // > 10
        while let Some(s) = s.next() {
            let ss = s.as_str();
            if ss.starts_with("MD") {
                a.md = &s[5..];
            } else if ss.starts_with("NM") {
                a.nh = (&ss[5..]).parse().map_err(|_| ())?;
            } else if ss.starts_with("YZ") {
                a.strand = s.last().ok_or(())?;
            }
        }

        if (ARGS.unique_only && !a.unique) || (ARGS.multiple_only && a.unique) {
            return Ok(a);
        }

        a.append_base();

        Ok(a)
    }
}

impl Alignment {
    fn new() -> Self {
        Self {
            dna: Default::default(),
            location: -1,
            mate_location: -1,
            flag: -1,
            mapped: false,
            md: Default::default(),
            sequence: Default::default(),
            quality: Default::default(),
            unique: false,
            map_q: Default::default(),
            nh: -1,
            bases: Vec::new(),
            read_name_id: 0,
            sequence_covered_length: 0,
            overlap: false,
            paired: false,
            strand: Default::default(),
            cigar: Default::default(),
        }
    }

    fn hash_name(name: &AsciiStr) -> usize {
        let mut r = 0;
        let a = 63689;
        for i in 0..name.len() {
            r = (r * a) + name[i] as usize;
        }
        r
    }

    fn cigar_get_next_segment(
        cigar_string: &AsciiStr,
        start: &mut usize,
        str_len: usize,
        len: &mut usize,
        symbol: &mut AsciiChar,
    ) -> bool {
        if *start == str_len {
            return false;
        }
        *len = 0;
        let mut current_index = *start;
        loop {
            if cigar_string[current_index].is_ascii_alphabetic() {
                *len = (&cigar_string[*start..current_index - *start])
                    .as_str()
                    .parse()
                    .unwrap();
                *symbol = cigar_string[current_index];
                *start = current_index + 1;
                return true;
            }
            current_index += 1;
        }
    }

    fn md_get_next_segment(
        md_string: &AsciiStr,
        start: &mut usize,
        str_len: usize,
        seg: &mut AsciiString,
    ) -> bool {
        if *start >= str_len {
            return false;
        }
        seg.clear();
        let mut current_index = *start;
        let mut deletion = false;

        loop {
            if current_index >= str_len {
                *start = current_index + 1;
                return !seg.is_empty();
            }
            if seg.is_empty() && md_string[current_index] == '0' {
                current_index += 1;
                continue;
            }
            if md_string[current_index].is_alphabetic() {
                if seg.is_empty() {
                    *seg = md_string[current_index].into();
                    *start = current_index + 1;
                    return true;
                } else {
                    if deletion {
                        (*seg).push(md_string[current_index]);
                    } else {
                        *start = current_index;
                        return true;
                    }
                }
            } else if md_string[current_index] == '^' {
                if seg.is_empty() {
                    *seg = md_string[current_index].into();
                    deletion = true;
                } else {
                    *start = current_index;
                    return true;
                }
            } else {
                // number
                if seg.is_empty() {
                    *seg = md_string[current_index].into();
                } else {
                    if deletion || seg.last().unwrap().is_alphabetic() {
                        *start = current_index;
                        return true;
                    } else {
                        seg.push(md_string[current_index]);
                    }
                }
            }
            current_index += 1;
        }
    }

    fn adjust_pos(&mut self) -> isize {
        let mut read_pos = 0;
        let mut return_pos = 0;
        let mut seq_length = self.sequence.len();

        let mut cigar_symbol: AsciiChar = Default::default();
        let mut cigar_len = 0;
        self.sequence_covered_length = 0;

        let mut cigar_start = 0;
        let cigar_str_len = self.cigar.len();

        while Self::cigar_get_next_segment(
            self.cigar,
            &mut cigar_start,
            cigar_str_len,
            &mut cigar_len,
            &mut cigar_symbol,
        ) {
            self.sequence_covered_length += cigar_len;
            match cigar_symbol {
                AsciiChar::S => {
                    if read_pos == 0 {
                        return_pos = cigar_len;
                        for i in cigar_len..seq_length {
                            self.bases[i].ref_pos -= cigar_len as isize;
                        }
                    } else {
                    }
                    read_pos += cigar_len;
                }
                AsciiChar::N => {
                    for i in read_pos..seq_length {
                        self.bases[i].ref_pos += cigar_len as isize;
                    }
                }
                AsciiChar::M => {
                    for i in read_pos..read_pos + cigar_len {
                        self.bases[i].remove = false;
                    }
                    read_pos += cigar_len;
                }
                AsciiChar::I => {
                    for i in read_pos + cigar_len..seq_length {
                        self.bases[i].ref_pos -= cigar_len as isize;
                    }
                    read_pos += cigar_len;
                }
                AsciiChar::D => {
                    for i in read_pos..seq_length {
                        self.bases[i].ref_pos += cigar_len as isize;
                    }
                }
                _ => (),
            }
        }
        return return_pos as isize;
    }

    fn append_base(&mut self) {
        // TODO: check understanding
        // original impl checks sequence_covered_length, which should
        // always be 0 at this time
        if !self.mapped || self.sequence_covered_length > 500000 {
            return;
        }

        self.bases.reserve_exact(self.sequence.len());
        for i in 0..self.sequence.len() {
            self.bases.push(PosQuality::new(i as isize));
        }

        let mut pos = self.adjust_pos() as usize;

        let mut md_start = 0;
        let md_str_len = self.md.len();
        let mut seg = AsciiString::new();
        while Self::md_get_next_segment(self.md, &mut md_start, md_str_len, &mut seg) {
            let ch = seg.first().unwrap();
            if ch.is_ascii_digit() {
                let len: usize = seg.as_str().parse().unwrap();
                for i in 0..len {
                    while self.bases[pos].remove {
                        pos += 1;
                    }
                    if self.strand == '+' && self.sequence[pos] == ARGS.base_change.0.0
                        || self.strand == '-' && self.sequence[pos] == ARGS.base_change.0.1
                    {
                        self.bases[pos].set_qual(self.quality[pos], false);
                    } else {
                        self.bases[pos].remove = true;
                    }
                    pos += 1;
                }
            } else if ch.is_alphabetic() {
                while self.bases[pos].remove {
                    pos += 1;
                }

                if (self.strand == '+'
                    && ch == ARGS.base_change.0.0
                    && self.sequence[pos] == ARGS.base_change.1.0)
                    || (self.strand == '-'
                        && ch == ARGS.base_change.0.1
                        && self.sequence[pos] == ARGS.base_change.1.1)
                {
                    self.bases[pos].set_qual(self.quality[pos], true);
                } else {
                    self.bases[pos].remove = true;
                }
                pos += 1;
            } else {

            }
        }
    }
}

pub struct AlignmentIter {
    lines: Box<dyn DoubleEndedIterator<Item = &'static AsciiStr>>,
}

impl AlignmentIter {
    pub fn new(text: &'static AsciiStr) -> Self {
        Self {
            lines: Box::new(text.lines()),
        }
    }
}

impl Iterator for AlignmentIter {
    type Item = Alignment;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.lines.next() {
                Some(line) => {
                    match line.first() {
                        None | Some(AsciiChar::At) => continue,
                        _ => (),
                    }
                    let a = match Alignment::try_from(line) {
                        Ok(a) => a,
                        Err(e) => {
                            eprintln!("'{}' is not a valid Alignment.", line);
                            continue;
                        }
                    };
                    if !a.mapped || a.bases.is_empty() {
                        continue;
                    }
                    return Some(a);
                }
                None => return None,
            }
        }
    }
}

unsafe impl Send for AlignmentIter {}
