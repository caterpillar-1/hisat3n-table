use crate::utils::{md_get_next_segment, ChunkIterator, CigarIterator, StringSearchState};
use crate::ARGS;

#[derive(Debug, Default)]
pub struct PosQuality {
    pub ref_pos: isize,
    pub qual: u8,
    pub converted: bool,
    pub remove: bool,
}

impl PosQuality {
    pub fn new(pos: isize) -> Self {
        PosQuality {
            ref_pos: pos,
            remove: true,
            ..Default::default()
        }
    }

    pub fn set_qual(&mut self, qual: u8, converted: bool) {
        self.qual = qual;
        self.converted = converted;
        self.remove = false;
    }
}

pub struct Alignment<'a> {
    pub dna: &'a [u8],
    pub location: isize,
    pub mate_location: isize,
    pub flag: i32,
    pub mapped: bool,
    pub strand: u8,
    pub sequence: &'a [u8],
    pub quality: &'a [u8],
    pub unique: bool,
    pub map_q: &'a [u8],
    pub nh: i32,
    pub bases: Vec<PosQuality>,
    pub cigar: &'a [u8],
    pub md: &'a [u8],
    pub read_name_id: u64,
    pub sequence_covered_length: usize,
    pub overlap: bool,
    pub paired: bool,
}

// static debugfile: std::sync::LazyLock<std::sync::Mutex<File>> = std::sync::LazyLock::new(|| std::sync::Mutex::new(File::create("test2.check").unwrap()));

impl<'a> Alignment<'a> {
    pub(crate) fn from_file(data: &'a [u8])  -> Result<Self, ()> {
        if data.len() < 1 || data[0] == b'@' {
            return Err(());
        }
        let mut a = Self::new();

        let iter = memchr::memchr_iter(b'\t', data);
        let mut s = ChunkIterator::new(data, iter);
        // 0
        a.read_name_id = Self::name_hash_str(s.next().ok_or(())?);
        // 1
        a.flag = atoi_simd::parse(s.next().ok_or(())?).map_err(|_| ())?;
        a.mapped = (a.flag & 4) == 0;
        a.paired = (a.flag & 1) != 0;

        // 2
        a.dna = s.next().ok_or(())?;
        // 3
        a.location = atoi_simd::parse(s.next().ok_or(())?).map_err(|_| ())?;
        // 4
        a.map_q = s.next().ok_or(())?;
        a.unique = a.map_q != b"1";
        // 5
        a.cigar = s.next().ok_or(())?;
        // 6
        s.next().ok_or(())?;
        // 7
        a.mate_location = atoi_simd::parse(s.next().ok_or(())?).map_err(|_| ())?;
        // 8
        s.next().ok_or(())?;
        // 9
        a.sequence = s.next().ok_or(())?;
        // 10
        a.quality = s.next().ok_or(())?;
        // > 10
        while let Some(s) = s.next() {
            if s.starts_with(b"MD") {
                a.md = &s[5..];
            } else if s.starts_with(b"NM") {
                a.nh = atoi_simd::parse(&s[5..]).map_err(|_| ())?;
            } else if s.starts_with(b"YZ") {
                a.strand = *s.last().ok_or(())?;
            }
        }

        if (ARGS.unique_only && !a.unique) || (ARGS.multiple_only && a.unique) {
            return Ok(a);
        }
        a.append_base();
        Ok(a)
    }

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

    pub fn name_hash_str(name: &[u8]) -> u64 {
        let mut hash: u64 = 0;
        let a: u64 = 63689;
        for byte in name {
            let byte_val = *byte as u64;
            hash = hash.wrapping_mul(a).wrapping_add(byte_val);
        }
        hash
    }

    fn adjust_pos(&mut self) -> usize {
        let mut read_pos = 0;
        let mut return_pos = 0;
        let seq_length = self.sequence.len();
        self.sequence_covered_length = 0;
        for (cigar_len, symbol) in CigarIterator::new(self.cigar) {
            self.sequence_covered_length += cigar_len;
            match symbol {
                b'S' => {
                    if read_pos == 0 {
                        return_pos = cigar_len;
                        for i in cigar_len..seq_length {
                            self.bases[i].ref_pos -= cigar_len as isize;
                        }
                    }
                    read_pos += cigar_len;
                }
                b'N' => {
                    for i in read_pos..seq_length {
                        self.bases[i].ref_pos += cigar_len as isize;
                    }
                }
                b'M' => {
                    for i in read_pos..(read_pos + cigar_len) {
                        self.bases[i].remove = false;
                    }
                    read_pos += cigar_len;
                }
                b'I' => {
                    for i in (read_pos + cigar_len)..seq_length {
                        self.bases[i].ref_pos -= cigar_len as isize;
                    }
                    read_pos += cigar_len;
                }
                b'D' => {
                    for i in read_pos..seq_length {
                        self.bases[i].ref_pos += cigar_len as isize;
                    }
                },
                _ => {}
            }
        }
        return_pos
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

        let mut pos = self.adjust_pos();
        let mut search = StringSearchState::new(self.md);
        let mut seg = Vec::<u8>::new();
        while md_get_next_segment(&mut search, &mut seg) {
            let ref_base = seg.first().unwrap();
            if ref_base.is_ascii_digit() {
                let len: usize = atoi_simd::parse(seg.as_slice()).unwrap();
                for _ in 0..len {
                    while self.bases[pos].remove {
                        pos += 1;
                    }
                    if self.strand == b'+' && self.sequence[pos] == ARGS.base_change.0.0
                      || self.strand == b'-' && self.sequence[pos] == ARGS.base_change.0.1
                    {
                        self.bases[pos].set_qual(self.quality[pos], false);
                    } else {
                        self.bases[pos].remove = true;
                    }
                    pos += 1;
                }
            } else if ref_base.is_ascii_alphabetic() {
                while self.bases[pos].remove {
                    pos += 1;
                }
                if self.strand == b'+' && *ref_base == ARGS.base_change.0.0 && self.sequence[pos] == ARGS.base_change.1.0
                  || self.strand == b'-' && *ref_base == ARGS.base_change.0.1 && self.sequence[pos] == ARGS.base_change.1.1
                {
                    self.bases[pos].set_qual(self.quality[pos], true);
                } else {
                    self.bases[pos].remove = true;
                }
                pos += 1;
            }
        }
    }
}
