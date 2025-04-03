use ascii::{AsciiChar, AsciiStr, IntoAsciiString, ToAsciiChar};

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

        let mut s = s.split('\t'.to_ascii_char().unwrap());
        // 0
        a.read_name_id = Self::hash_name(s.next().unwrap());
        // 1
        a.flag = s.next().unwrap().as_str().parse().unwrap();
        // 2
        a.dna = s.next().unwrap();
        // 3
        a.location = s.next().unwrap().as_str().parse().unwrap();
        // 4
        a.map_q = s.next().unwrap();
        a.unique = a.map_q != "1".into_ascii_string().unwrap();
        // 5
        a.cigar = s.next().unwrap();
        // 6
        s.next().unwrap();
        // 7
        a.mate_location = s.next().unwrap().as_str().parse().unwrap();
        // 8
        s.next().unwrap();
        // 9
        a.sequence = s.next().unwrap();
        // 10
        a.quality = s.next().unwrap();
        // > 10
        while let Some(s) = s.next() {
            let ss = s.as_str();
            if ss.starts_with("MD") {
                a.md = &s[5..];
            } else if ss.starts_with("NM") {
                a.nh = (&ss[5..]).parse().unwrap();
            } else if ss.starts_with("YZ") {
                a.strand = s.last().unwrap();
            }
        }

        Ok(a)
    }
}

impl Alignment {
    pub fn new() -> Self {
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
                    let a = Alignment::try_from(line).unwrap();
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
