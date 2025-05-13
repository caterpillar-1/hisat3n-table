pub static BASE_CHARS: [u8; 4] = [b'A', b'T', b'C', b'G'];

#[inline]
pub fn asc2dnacomp(ch: u8) -> u8 {
    match ch {
        b'-' => b'-', // Gap
        b'A' => b'T', // Adenine -> Thymine
        b'B' => b'V', // B (C,G,T) -> V (A,C,G)
        b'C' => b'G', // Cytosine -> Guanine
        b'D' => b'H', // D (A,G,T) -> H (A,C,T)
        b'G' => b'C', // Guanine -> Cytosine
        b'H' => b'D', // H (A,C,T) -> D (A,G,T)
        b'K' => b'M', // K (G,T) -> M (A,C)
        b'M' => b'K', // M (A,C) -> K (G,T)
        b'N' => b'N', // N (any) -> N (any)
        b'R' => b'Y', // R (A,G) -> Y (C,T) (Purine -> Pyrimidine)
        b'S' => b'S', // S (C,G) -> S (C,G) (Strong -> Strong)
        b'T' => b'A', // Thymine -> Adenine
        b'U' => b'A', // Uracil -> Adenine (RNA context, often handled similarly to T)
        b'V' => b'B', // V (A,C,G) -> B (C,G,T)
        b'W' => b'W', // W (A,T) -> W (A,T) (Weak -> Weak)
        b'Y' => b'R',
        _ => b'\0',
    }
}

pub struct ChunkIterator<'a, I: DoubleEndedIterator<Item = usize>> {
    buffer: &'a [u8],
    separator_indices: I,
    current_start: usize,
    finished: bool,
}

impl<'a, I: DoubleEndedIterator<Item = usize>> ChunkIterator<'a, I> {
    pub fn new(buffer: &'a [u8], separator_indices: I) -> Self {
        ChunkIterator {
            buffer,
            separator_indices,
            current_start: 0,
            finished: false,
        }
    }
}

impl<'a, I: DoubleEndedIterator<Item = usize>> Iterator for ChunkIterator<'a, I> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        match self.separator_indices.next() {
            Some(sep_index) => {
                let start = self.current_start;
                self.current_start = sep_index + 1;
                Some(&self.buffer[start..sep_index])
            }
            None => {
                let slice = &self.buffer[self.current_start..];
                self.finished = true;
                Some(slice)
            }
        }
    }
}

pub struct StringSearchState<'a> {
    s: &'a [u8],
    start: usize,
}

pub struct CigarIterator<'a> {
    s: &'a [u8],
    start: usize,
}

impl<'a> StringSearchState<'a> {
    pub fn new(s: &'a [u8]) -> Self {
        Self { s, start: 0 }
    }
}

impl<'a> CigarIterator<'a> {
    pub fn new(input_string: &'a [u8]) -> Self {
        Self { s: input_string, start: 0 }
    }
}

impl<'a> Iterator for CigarIterator<'a> {
    type Item = (usize, u8);

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.s.len() {
            return None;
        }
        let mut current_index = self.start;
        while current_index < self.s.len() {
            let current_char_byte = self.s[current_index];
            if current_char_byte.is_ascii_alphabetic() {
                let num_slice = &self.s[self.start..current_index];
                return match atoi_simd::parse::<usize>(num_slice) {
                    Ok(len) => {
                        self.start = current_index + 1;
                        Some((len, current_char_byte))
                    }
                    Err(_) => None
                }
            } else if !current_char_byte.is_ascii_digit() {
                // Malformed CIGAR: Character in the length part is not a digit.
                return None;
            }
            current_index += 1;
        }
        None
    }
}

/// Parses the next MD tag segment from the state.
/// An MD tag segment can be a number of matching bases (as Vec<u8> representing the digits),
/// a mismatched base (Vec<u8> with one char), or a deletion (Vec<u8> starting with '^').
pub fn md_get_next_segment(state: &mut StringSearchState, seg: &mut Vec<u8>) -> bool {
    if state.start >= state.s.len() {
        return false;
    }
    seg.clear();
    let mut current_index = state.start;
    let mut deletion = false; // True if the current segment under construction is a deletion
    loop {
        if current_index >= state.s.len() {
            state.start = current_index + 1;
            return if seg.is_empty() { false } else { true };
        }
        let current_char_byte = state.s[current_index];
        if seg.is_empty() && current_char_byte == b'0' { // skip zero-prefixes?
            current_index += 1;
            continue;
        }
        if current_char_byte.is_ascii_alphabetic() {
            if seg.is_empty() { // Segment is a single mismatched base (e.g., 'A').
                seg.push(current_char_byte);
                state.start = current_index + 1;
                return true;
            } else { // Segment not empty. Alphabetic char may end current segment or extend deletion.
                if deletion { // Current segment is a deletion, like "^AC". Append current char.
                    seg.push(current_char_byte);
                } else { // Current segment is numeric (e.g. "123"). Alpha char starts the next segment.
                    state.start = current_index;
                    return true;
                }
            }
        } else if current_char_byte == b'^' {
            if seg.is_empty() {
                // Start of a new deletion segment.
                seg.push(current_char_byte);
                deletion = true;
            } else {
                // Segment not empty (e.g., numeric "123").
                // '^' signifies the start of the next segment, so end this one.
                state.start = current_index;
                return true;
            }
        } else if current_char_byte.is_ascii_digit() {
            if seg.is_empty() { // Start of a new numeric segment.
                seg.push(current_char_byte);
            } else { // Segment not empty. `seg.last()` is safe.
                let last_char_in_seg = *seg.last().unwrap();
                if deletion || last_char_in_seg.is_ascii_alphabetic() {
                    // Current segment is deletion (e.g., "^A") or mismatch (e.g., "A").
                    // A digit signifies the end of that segment and start of new numeric one.
                    state.start = current_index;
                    return true;
                } else { // Current segment is numeric (e.g., "12"), and current char is another digit. Append.
                    seg.push(current_char_byte);
                }
            }
        } else {
            // Invalid character in MD string (not alphabet, not '^', not digit).
            // SAM spec: MD tags only contain [0-9A-Z^].
            // C++ code's `else { // number }` might misinterpret.
            // This Rust version: end current segment (if any) and stop.
            state.start = current_index;
            return if seg.is_empty() { false } else { true };
        }

        current_index += 1;
    }
}

