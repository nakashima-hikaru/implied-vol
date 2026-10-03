//! Stable benchmark of fixed solver paths and a reproducible mixed workload.
//! Run with `cargo bench --bench implied_black_vol_paths -- 200000 5 --solver hybrid`.

use implied_vol::solver::{BlackSolver, Hybrid, Jaeckel};
use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, PriceBlackScholes,
    PriceBlackScholesNormalised,
};
use rand::{RngExt, SeedableRng};
use std::{hint::black_box, marker::PhantomData, time::Instant};

#[path = "support/mod.rs"]
mod support;

// Keep the comparison adapter local to the benchmark. Rust retains its original
// prebuilt case layout, while C++ receives the same inputs as raw arguments.
trait BenchmarkSolver {
    type NormalisedCase;
    type FullCase;

    fn prepare_normalised(x: f64, beta: f64) -> Self::NormalisedCase;
    fn prepare_full(price: f64, f: f64, k: f64, expiry: f64) -> Self::FullCase;
    fn normalised_volatility(case: &Self::NormalisedCase) -> Option<f64>;
    fn full_volatility(case: &Self::FullCase) -> Option<f64>;
}

struct RustSolver<S>(PhantomData<S>);

impl<S: BlackSolver> BenchmarkSolver for RustSolver<S> {
    type NormalisedCase = ImpliedBlackVolatilityNormalised;
    type FullCase = ImpliedBlackVolatility;

    fn prepare_normalised(x: f64, beta: f64) -> Self::NormalisedCase {
        ImpliedBlackVolatilityNormalised::builder()
            .log_moneyness(x)
            .normalised_price(beta)
            .build()
            .unwrap()
    }

    fn prepare_full(price: f64, f: f64, k: f64, expiry: f64) -> Self::FullCase {
        ImpliedBlackVolatility::builder()
            .forward(f)
            .strike(k)
            .expiry(expiry)
            .is_call(true)
            .option_price(price)
            .build()
            .unwrap()
    }

    #[inline(always)]
    fn normalised_volatility(case: &Self::NormalisedCase) -> Option<f64> {
        case.calculate_with::<S>()
    }

    #[inline(always)]
    fn full_volatility(case: &Self::FullCase) -> Option<f64> {
        case.calculate_with::<S>()
    }
}

#[cfg(feature = "cxx_bench")]
struct CppSolver;

#[cfg(feature = "cxx_bench")]
impl BenchmarkSolver for CppSolver {
    type NormalisedCase = (f64, f64);
    type FullCase = (f64, f64, f64, f64);

    fn prepare_normalised(x: f64, beta: f64) -> Self::NormalisedCase {
        (x, beta)
    }

    fn prepare_full(price: f64, f: f64, k: f64, expiry: f64) -> Self::FullCase {
        (price, f, k, expiry)
    }

    #[inline(always)]
    fn normalised_volatility(&(x, beta): &Self::NormalisedCase) -> Option<f64> {
        // Every generated case is an OTM call (x <= 0), so no intrinsic-value
        // conversion is needed in the normalized C++ entry point.
        Some(implied_vol::cxx::ffi::NormalisedImpliedBlackVolatility(
            beta, x, 1.0,
        ))
    }

    #[inline(always)]
    fn full_volatility(&(price, f, k, expiry): &Self::FullCase) -> Option<f64> {
        Some(implied_vol::cxx::ffi::ImpliedBlackVolatility(
            price, f, k, expiry, 1.0,
        ))
    }
}

fn time_normalised<S: BenchmarkSolver>(cases: &[S::NormalisedCase], n: usize) -> (f64, f64) {
    let start = Instant::now();
    let mut sum = 0.0;
    for i in 0..n {
        sum += black_box(S::normalised_volatility(black_box(&cases[i % cases.len()])).unwrap());
    }
    #[allow(clippy::cast_precision_loss)]
    let ns = start.elapsed().as_secs_f64() * 1e9 / n as f64;
    (ns, black_box(sum))
}

fn time_full<S: BenchmarkSolver>(cases: &[S::FullCase], n: usize) -> (f64, f64) {
    let start = Instant::now();
    let mut sum = 0.0;
    for i in 0..n {
        sum += black_box(S::full_volatility(black_box(&cases[i % cases.len()])).unwrap());
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
    let mut experimental_paths = false;
    let mut numbers = Vec::new();
    while let Some(arg) = args.next() {
        if arg == "--solver" {
            solver = args.next().expect("--solver requires a solver name");
        } else if arg == "--experimental-paths" {
            experimental_paths = true;
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
        "hybrid" => run::<RustSolver<Hybrid>>(n, rounds, experimental_paths),
        "jaeckel" => run::<RustSolver<Jaeckel>>(n, rounds, experimental_paths),
        "flashiv" => {
            #[cfg(feature = "flashiv")]
            run::<RustSolver<implied_vol::solver::FlashIv>>(n, rounds, experimental_paths);
            #[cfg(not(feature = "flashiv"))]
            usage_error("solver flashiv requires the flashiv Cargo feature");
        }
        "experimental" => {
            #[cfg(feature = "experimental")]
            run::<RustSolver<implied_vol::solver::Experimental>>(n, rounds, experimental_paths);
            #[cfg(not(feature = "experimental"))]
            usage_error("solver experimental requires the experimental Cargo feature");
        }
        "cpp" => {
            #[cfg(feature = "cxx_bench")]
            run::<CppSolver>(n, rounds, experimental_paths);
            #[cfg(not(feature = "cxx_bench"))]
            usage_error("solver cpp requires the cxx_bench Cargo feature");
        }
        _ => usage_error(&format!(
            "unknown solver {solver}; expected hybrid, jaeckel, flashiv, experimental, or cpp"
        )),
    }
}

fn usage_error(message: &str) -> ! {
    eprintln!("{message}");
    std::process::exit(2);
}

fn run<S: BenchmarkSolver>(n: usize, rounds: usize, experimental_paths: bool) {
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
        fixed.push((name, vec![S::prepare_normalised(x, normalised_price(x, s))]));
    }

    if experimental_paths {
        // Exact represented (a,b) inputs exercise paths absent from the
        // ordinary a<5 mixed workload. No price evaluation is timed.
        for (name, a, b) in [
            ("central_deferred", 0.1, 0.05),
            ("wing_seed", 1.0, 1e-10),
            ("finite_seed", 5.0, 0.005),
            ("large_lower", 20.0, (-10.0_f64).exp() * 0.1),
            ("large_upper", 20.0, (-10.0_f64).exp() * 0.9),
            ("large_asymptotic", 100.0, f64::MIN_POSITIVE),
            ("large_near_cap", 300.0, (-150.0_f64).exp() * (1.0 - 1e-12)),
            ("microscopic", 1e-200, 1e-201),
        ] {
            fixed.push((name, vec![S::prepare_normalised(-a, b)]));
        }
        let mut rng = rand::rngs::StdRng::from_seed([73; 32]);
        let mut large = Vec::with_capacity(4096);
        while large.len() < 4096 {
            let a: f64 = rng.random_range(10.0..1416.0);
            let z: f64 = rng.random_range(0.01..300.0);
            #[allow(clippy::suboptimal_flops)] // Keep workload construction arithmetic explicit.
            let b = (-0.5 * a - z).exp();
            if b >= f64::MIN_POSITIVE {
                large.push(S::prepare_normalised(-a, b));
            }
        }
        fixed.push(("mixed_large", large));
        let mut near_atm = Vec::with_capacity(4096);
        let mut rng = rand::rngs::StdRng::from_seed([91; 32]);
        for _ in 0..4096 {
            let x: f64 = -rng.random_range(1e-6..0.25);
            let s: f64 = rng.random_range(0.1..0.5);
            near_atm.push(S::prepare_normalised(x, normalised_price(x, s)));
        }
        fixed.push(("mixed_near_atm", near_atm));
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
        normalised.push(S::prepare_normalised(x, beta));
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
        full.push(S::prepare_full(price, f, k, expiry));
    }

    // Match the concrete inputs in the legacy and split Actions benchmarks.
    // Prices and both adapters' case containers are prepared before timing.
    let mut fixed_full = Vec::new();
    for (name, source) in [
        ("deep_otm_full", support::DEEP_OTM_CALL_LONG),
        ("near_atm_short_full", support::NEAR_ATM_CALL_SHORT),
    ] {
        let case = support::black_implied_case(source);
        fixed_full.push((
            name,
            vec![S::prepare_full(
                case.option_price,
                case.forward,
                case.strike,
                case.expiry,
            )],
        ));
    }
    let mut legacy_rng = rand::rngs::StdRng::from_seed([13; 32]);
    let (r, r2, r3): (f64, f64, f64) = legacy_rng.random();
    fixed_full.push((
        "legacy_otm_full",
        vec![S::prepare_full(r * r2, r, 1.0, 1e5 * r3)],
    ));

    // Validate every prepared input before warming up or reporting timings.
    // Numerical failure must not be mistaken for a fast solver result.
    for (name, cases) in &fixed {
        preflight_normalised::<S>(name, cases);
    }
    preflight_normalised::<S>("mixed_normalised", &normalised);
    for (index, case) in full.iter().enumerate() {
        preflight_output("mixed_full", index, S::full_volatility(case));
    }
    for (name, cases) in &fixed_full {
        for (index, case) in cases.iter().enumerate() {
            preflight_output(name, index, S::full_volatility(case));
        }
    }

    for (_, cases) in &fixed {
        black_box(time_normalised::<S>(cases, 10_000));
    }
    black_box(time_normalised::<S>(&normalised, 10_000));
    black_box(time_full::<S>(&full, 10_000));
    for (_, cases) in &fixed_full {
        black_box(time_full::<S>(cases, 10_000));
    }

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
        for (name, cases) in &fixed_full {
            let (ns, sum) = time_full::<S>(cases, n);
            println!("{name},{round},{ns:.6},{sum:.16e}");
        }
    }
}

fn preflight_normalised<S: BenchmarkSolver>(name: &str, cases: &[S::NormalisedCase]) {
    for (index, case) in cases.iter().enumerate() {
        preflight_output(name, index, S::normalised_volatility(case));
    }
}

fn preflight_output(name: &str, index: usize, output: Option<f64>) {
    assert!(
        output.is_some_and(|value| value.is_finite() && value > 0.0 && value < f64::MAX),
        "solver failed preflight for {name}[{index}]: {output:?}"
    );
}
