//! Scalar/whole-book A/B benchmark: identical prepared quotes and one process.
//! Use with benches/compare_volfi.py, not the historical cross-host README values.
use implied_vol::solver::{
    BlackSolver, Experimental, VolfiBook, VolfiBookCompensated, VolfiBookRsqrt, VolfiBookSinh,
    VolfiBookSelective,
};
use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised,
    PriceBlackScholesNormalised, volfi_book,
};
use rand::{RngExt, SeedableRng};
use std::{hint::black_box, time::Instant};

const CASES: &[&str] = &[
    "book_a",
    "book_b",
    "book_seam",
    "fallback",
    "mixed_book",
    "mixed_market",
    "mixed_full",
];

fn from_intrinsic(h: f64, a: f64) -> (f64, f64) {
    let c = h.exp_m1() / a.exp_m1();
    (-h, c * (-0.5 * h).exp())
}

fn price(x: f64, v: f64) -> f64 {
    PriceBlackScholesNormalised::builder()
        .log_moneyness(x)
        .total_volatility(v)
        .build()
        .expect("valid price parameters")
        .calculate::<DefaultSpecialFn>()
}

fn make_cases(name: &str) -> Vec<(f64, f64)> {
    match name {
        "book_a" => vec![from_intrinsic(0.2, 1.4)],
        "book_b" => vec![from_intrinsic(1.0, 9.0)],
        "book_seam" => vec![from_intrinsic(0.8, std::f64::consts::TAU)],
        "fallback" => vec![(-5.0, 0.05 * (-2.5_f64).exp())],
        "mixed_book" => {
            let mut rng = rand::rngs::StdRng::from_seed([124; 32]);
            let mut inputs = Vec::with_capacity(4096);
            while inputs.len() < 4096 {
                let h: f64 = rng.random_range(0.03..1.4);
                let a: f64 = if inputs.len() % 2 == 0 {
                    rng.random_range((h / 0.3 + 0.1)..6.15)
                } else {
                    rng.random_range(8.0..40.0)
                };
                let (x, b) = from_intrinsic(h, a);
                if b.is_normal() && volfi_book::chart_total_volatility(-x, b).is_some() {
                    inputs.push((x, b));
                }
            }
            inputs
        }
        "mixed_market" | "mixed_full" => {
            let mut rng = rand::rngs::StdRng::from_seed([132; 32]);
            let mut inputs = Vec::with_capacity(4096);
            while inputs.len() < 4096 {
                let h: f64 = rng.random_range(0.005..4.0);
                let v: f64 = rng.random_range(0.07..3.0);
                let x = -h;
                let b = price(x, v);
                if b >= f64::MIN_POSITIVE && b < (-0.5 * h).exp() {
                    inputs.push((x, b));
                }
            }
            inputs
        }
        _ => panic!("unknown case: {name}"),
    }
}

#[inline(never)]
fn time_core<S: BlackSolver>(cases: &[ImpliedBlackVolatilityNormalised], n: usize) -> (f64, f64) {
    let now = Instant::now();
    let mut checksum = 0.0;
    for i in 0..n {
        checksum += black_box(
            black_box(&cases[i % cases.len()])
                .calculate_with::<S>()
                .expect("finite benchmark output"),
        );
    }
    #[allow(clippy::cast_precision_loss)]
    let elapsed = now.elapsed().as_secs_f64() * 1e9 / n as f64;
    (elapsed, black_box(checksum))
}

#[inline(never)]
fn time_full<S: BlackSolver>(cases: &[ImpliedBlackVolatility], n: usize) -> (f64, f64) {
    let now = Instant::now();
    let mut checksum = 0.0;
    for i in 0..n {
        checksum += black_box(
            black_box(&cases[i % cases.len()])
                .calculate_with::<S>()
                .expect("finite benchmark output"),
        );
    }
    #[allow(clippy::cast_precision_loss)]
    let elapsed = now.elapsed().as_secs_f64() * 1e9 / n as f64;
    (elapsed, black_box(checksum))
}

fn run<S: BlackSolver>(name: &str, n: usize, rounds: usize) {
    let inputs = make_cases(name);
    let covered = inputs
        .iter()
        .filter(|&&(x, b)| volfi_book::chart_total_volatility(x.abs(), b).is_some())
        .count();
    let norm = inputs
        .iter()
        .map(|&(x, b)| {
            ImpliedBlackVolatilityNormalised::builder()
                .log_moneyness(x)
                .normalised_price(b)
                .build()
                .expect("valid normalized input")
        })
        .collect::<Vec<_>>();
    let full = inputs
        .iter()
        .map(|&(x, b)| {
            let f: f64 = 100.0;
            let k = f * (-x).exp();
            let p = b * (f * k).sqrt();
            ImpliedBlackVolatility::builder()
                .forward(f)
                .strike(k)
                .expiry(1.0)
                .is_call(true)
                .option_price(p)
                .build()
                .expect("valid full input")
        })
        .collect::<Vec<_>>();
    for &(x, b) in &inputs {
        let baseline = Experimental::implied_total_volatility(x, b);
        let candidate = VolfiBook::implied_total_volatility(x, b);
        assert!(matches!(baseline, Some(v) if v.is_finite() && v > 0.0));
        assert!(matches!(candidate, Some(v) if v.is_finite() && v > 0.0));
    }
    let full_api = name == "mixed_full";
    if full_api {
        black_box(time_full::<S>(&full, 10000));
    } else {
        black_box(time_core::<S>(&norm, 10000));
    }
    println!("case,round,ns_per_call,checksum,covered,total");
    for round in 0..rounds {
        let (ns, checksum) = if full_api {
            time_full::<S>(&full, n)
        } else {
            time_core::<S>(&norm, n)
        };
        println!(
            "{name},{round},{ns:.6},{checksum:.17e},{covered},{}",
            inputs.len()
        );
    }
}

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|x| x == "--test") {
        return;
    }
    let mut solver = "experimental";
    let mut case = "mixed_book";
    let mut numbers = Vec::new();
    let mut it = args.iter().filter(|x| *x != "--bench");
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--solver" => solver = it.next().expect("--solver value"),
            "--case" => case = it.next().expect("--case value"),
            _ => numbers.push(arg),
        }
    }
    assert!(CASES.contains(&case), "unknown case: {case}");
    let n = numbers
        .first()
        .map_or(150000, |v| v.parse().expect("iterations"));
    let rounds = numbers.get(1).map_or(2, |v| v.parse().expect("rounds"));
    assert!(n > 0 && rounds > 0 && numbers.len() <= 2);
    match solver {
        "experimental" => run::<Experimental>(case, n, rounds),
        "volfi-book" => run::<VolfiBook>(case, n, rounds),
        "volfi-sinh" => run::<VolfiBookSinh>(case, n, rounds),
        "volfi-rsqrt" => run::<VolfiBookRsqrt>(case, n, rounds),
        "volfi-compensated" => run::<VolfiBookCompensated>(case, n, rounds),
        "volfi-selective" => run::<VolfiBookSelective>(case, n, rounds),
        _ => panic!("unknown solver: {solver}"),
    }
}
