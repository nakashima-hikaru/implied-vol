//! Stable benchmark of fixed solver paths and a reproducible mixed workload.
//! Run with `cargo bench --bench implied_black_vol_paths -- 200000 5 --solver hybrid`.

use implied_vol::solver::{BlackSolver, Hybrid, Jaeckel};
use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, PriceBlackScholes,
    PriceBlackScholesNormalised,
};
use rand::{RngExt, SeedableRng};
use std::{hint::black_box, time::Instant};

fn time_normalised<S: BlackSolver>(
    cases: &[ImpliedBlackVolatilityNormalised],
    n: usize,
) -> (f64, f64) {
    let start = Instant::now();
    let mut sum = 0.0;
    for i in 0..n {
        sum += black_box(
            black_box(&cases[i % cases.len()])
                .calculate_with::<S>()
                .unwrap(),
        );
    }
    #[allow(clippy::cast_precision_loss)]
    let ns = start.elapsed().as_secs_f64() * 1e9 / n as f64;
    (ns, black_box(sum))
}

fn time_full<S: BlackSolver>(cases: &[ImpliedBlackVolatility], n: usize) -> (f64, f64) {
    let start = Instant::now();
    let mut sum = 0.0;
    for i in 0..n {
        sum += black_box(
            black_box(&cases[i % cases.len()])
                .calculate_with::<S>()
                .unwrap(),
        );
    }
    #[allow(clippy::cast_precision_loss)]
    let ns = start.elapsed().as_secs_f64() * 1e9 / n as f64;
    (ns, black_box(sum))
}

fn normalised_price(x: f64, s: f64) -> f64 {
    PriceBlackScholesNormalised::builder()
        .log_moneyness(x)
        .total_volatility(s)
        .build()
        .unwrap()
        .calculate::<DefaultSpecialFn>()
}

fn main() {
    // Cargo also passes --bench or --test to harness-free benchmarks.
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|arg| arg == "--test") {
        return;
    }
    let mut args = args.iter().filter(|arg| *arg != "--bench");
    let mut solver = "hybrid";
    let mut numbers = Vec::new();
    while let Some(arg) = args.next() {
        if arg == "--solver" {
            solver = args.next().expect("--solver requires a solver name");
        } else {
            numbers.push(arg);
        }
    }
    let n: usize = numbers.first().map_or(200_000, |arg| {
        arg.parse()
            .expect("first argument must be an iteration count")
    });
    let rounds: usize = numbers.get(1).map_or(5, |arg| {
        arg.parse().expect("second argument must be a round count")
    });
    assert!(n > 0 && rounds > 0 && numbers.len() <= 2);

    // Select one monomorphized runner before constructing or timing inputs.
    match solver {
        "hybrid" => run::<Hybrid>(n, rounds),
        "jaeckel" => run::<Jaeckel>(n, rounds),
        "flashiv" => {
            #[cfg(feature = "flashiv")]
            run::<implied_vol::solver::FlashIv>(n, rounds);
            #[cfg(not(feature = "flashiv"))]
            usage_error("solver flashiv requires the flashiv Cargo feature");
        }
        "experimental" => {
            #[cfg(feature = "experimental")]
            run::<implied_vol::solver::Experimental>(n, rounds);
            #[cfg(not(feature = "experimental"))]
            usage_error("solver experimental requires the experimental Cargo feature");
        }
        _ => usage_error(&format!(
            "unknown solver {solver}; expected hybrid, jaeckel, flashiv, or experimental"
        )),
    }
}

fn usage_error(message: &str) -> ! {
    eprintln!("{message}");
    std::process::exit(2);
}

fn run<S: BlackSolver>(n: usize, rounds: usize) {
    let mut fixed = Vec::new();
    for (name, x, s) in [
        ("atm", 0.0, 0.2),
        ("lowest", -0.5, 0.125),
        ("lower_middle", -0.5, 0.7),
        ("upper_middle", -0.5, 1.6),
        ("highest", -0.5, 4.0),
        ("near_atm", -0.01, 0.5),
        ("near_atm_wider", -0.1, 0.6),
    ] {
        fixed.push((
            name,
            vec![
                ImpliedBlackVolatilityNormalised::builder()
                    .log_moneyness(x)
                    .normalised_price(normalised_price(x, s))
                    .build()
                    .unwrap(),
            ],
        ));
    }

    // Prepare all prices and builders before timing. All solvers use the
    // same seeded log-moneyness and total-volatility distribution.
    let mut rng = rand::rngs::StdRng::from_seed([42; 32]);
    let mut normalised = Vec::with_capacity(4096);
    let mut full = Vec::with_capacity(4096);
    while normalised.len() < 4096 {
        let x: f64 = -rng.random_range(1e-4..5.0);
        let s: f64 = rng.random_range(1e-3..5.0);
        let beta = normalised_price(x, s);
        if beta < f64::MIN_POSITIVE || beta >= (0.5 * x).exp() {
            continue;
        }
        normalised.push(
            ImpliedBlackVolatilityNormalised::builder()
                .log_moneyness(x)
                .normalised_price(beta)
                .build()
                .unwrap(),
        );
        let f = 100.0;
        let k = f * (-x).exp();
        let expiry: f64 = rng.random_range(0.01..2.0);
        let price = PriceBlackScholes::builder()
            .forward(f)
            .strike(k)
            .volatility(s / expiry.sqrt())
            .expiry(expiry)
            .is_call(true)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        full.push(
            ImpliedBlackVolatility::builder()
                .forward(f)
                .strike(k)
                .expiry(expiry)
                .is_call(true)
                .option_price(price)
                .build()
                .unwrap(),
        );
    }

    for (_, cases) in &fixed {
        black_box(time_normalised::<S>(cases, 10_000));
    }
    black_box(time_normalised::<S>(&normalised, 10_000));
    black_box(time_full::<S>(&full, 10_000));

    println!("case,round,ns_per_call,checksum");
    for round in 0..rounds {
        for (name, cases) in &fixed {
            let (ns, sum) = time_normalised::<S>(cases, n);
            println!("{name},{round},{ns:.6},{sum:.16e}");
        }
        let (ns, sum) = time_normalised::<S>(&normalised, n);
        println!("mixed_normalised,{round},{ns:.6},{sum:.16e}");
        let (ns, sum) = time_full::<S>(&full, n);
        println!("mixed_full,{round},{ns:.6},{sum:.16e}");
    }
}
