//! Deterministic exact-output Experimental regression across independent builds.
//!
//! The pinned main and candidate must compile this *same* source with the same
//! flags on the same host. Every generated input's bits and its result bits
//! are emitted; the CI compares complete output files (not aggregate sums).
//! This is finite evidence and does not replace the continuous IEEE certificate.
use implied_vol::solver::{BlackSolver, Experimental};
use implied_vol::{DefaultSpecialFn, PriceBlackScholesNormalised};
use std::io::{self, BufWriter, Write};

#[inline]
fn random_word(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

#[inline]
fn unit(state: &mut u64) -> f64 {
    ((random_word(state) >> 11) as f64) * (1.0 / 9_007_199_254_740_992.0)
}

#[inline]
fn exp_between(state: &mut u64, lower_log: f64, upper_log: f64) -> f64 {
    (lower_log + (upper_log - lower_log) * unit(state)).exp()
}

fn generate(case: usize, state: &mut u64) -> (f64, f64) {
    let a = match case {
        0 => 0.0,
        1 => 0.45 * unit(state),
        2 => 0.5 + 1.1 * unit(state),
        3 => 1.6 + 8.4 * unit(state),
        4 => 0.0625 + 9.9375 * unit(state),
        5 => 0.5 + 9.5 * unit(state),
        6 => 10.0 + 1_405.0 * unit(state),
        7 => exp_between(state, -800.0, -10.0),
        8 => 0.25 * unit(state),
        9 => {
            let guard = [0.0625, 0.25, 0.45, 0.5, 1.6, 4.0, 10.0];
            let k = (random_word(state) as usize) % guard.len();
            let bits = guard[k].to_bits();
            f64::from_bits(bits + (random_word(state) % 3) - 1)
        }
        10 => 10.0 * unit(state),
        11 => exp_between(state, -700.0, 2.3),
        12 => 0.01 + 9.99 * unit(state),
        13 => 1e-5 + 1.0 * unit(state),
        14 => 2.0 + 8.0 * unit(state),
        _ => exp_between(state, -100.0, 2.3),
    };
    let cap = (-0.5 * a).exp();
    let b = match case {
        4 => cap * exp_between(state, -680.0, -5.0),
        5 => {
            let k = 1 + (random_word(state) % 56) as i32;
            cap * (1.0 - 2.0_f64.powi(-k))
        }
        6 => (-0.5 * a - 0.01 - 300.0 * unit(state)).exp(),
        7 => a * exp_between(state, -700.0, 7.0),
        8 => cap * (1.0 - 0.99 * unit(state)),
        9 => {
            let b0 = if (random_word(state) & 1) == 0 {
                cap * (0.001 + 0.998 * unit(state))
            } else {
                cap * 2.0_f64.powi(-((random_word(state) % 55) as i32 + 1))
            };
            f64::from_bits(b0.to_bits().saturating_add(random_word(state) % 3).saturating_sub(1))
        }
        11 => a * exp_between(state, -700.0, 4.0),
        12 => cap * unit(state),
        13 => cap * (0.5 + 0.49999999 * unit(state)),
        14 => cap * exp_between(state, -250.0, 0.0),
        _ => {
            let s = exp_between(state, -7.0, 1.8);
            PriceBlackScholesNormalised::builder()
                .log_moneyness(-a)
                .total_volatility(s)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>()
        }
    };
    (a, b)
}

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let mut state = 0x8d37_a5e4_2b9c_1f06_u64;
    let mut count = 0_usize;
    for i in 0..(16 * 8192) {
        let case = i % 16;
        let (a, b) = generate(case, &mut state);
        let got = Experimental::implied_total_volatility(-a, b);
        match got {
            Some(s) => writeln!(
                out,
                "{case:02x}|{:016x}|{:016x}|S|{:016x}",
                a.to_bits(),
                b.to_bits(),
                s.to_bits()
            )?,
            None => writeln!(
                out,
                "{case:02x}|{:016x}|{:016x}|N",
                a.to_bits(),
                b.to_bits()
            )?,
        }
        count += 1;
    }
    assert_eq!(count, 131_072);
    out.flush()
}
