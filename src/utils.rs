use ascii::AsciiChar;

pub static BASE_CHARS: [char; 4] = ['A', 'T', 'C', 'G'];

#[inline]
pub fn asc2dnacomp(ch: AsciiChar) -> AsciiChar {
    use AsciiChar::*;
    match ch {
        Minus => Minus,
        A => T,
        B => V,
        C => G,
        D => H,
        G => C,
        H => D,
        K => M,
        M => K,
        N => N,
        R => Y,
        S => S,
        T => A,
        V => B,
        W => W,
        Y => R,
        _ => Null,
    }
}
