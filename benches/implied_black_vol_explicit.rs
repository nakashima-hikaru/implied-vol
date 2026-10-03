#![feature(test)]

extern crate test;

use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, PriceBlackScholes,
    PriceBlackScholesNormalised, SpecialFn,
};
use test::Bencher;

const PAPER_DELTAS: [f64; 8] = [0.05, 0.20, 0.30, 0.45, 0.55, 0.70, 0.80, 0.95];

fn paper_total_vols() -> Vec<f64> {
    let mut vols = Vec::with_capacity(41);
    vols.push(0.01);
    let mut current = 0.05;
    while current <= 2.0 + f64::EPSILON {
        vols.push(current);
        current += 0.05;
    }
    vols
}

fn paper_grid_normalised_cases() -> Vec<ImpliedBlackVolatilityNormalised> {
    let mut cases = Vec::with_capacity(PAPER_DELTAS.len() * 41);
    for total_volatility in paper_total_vols() {
        for delta in PAPER_DELTAS {
            let log_moneyness = total_volatility
                * (DefaultSpecialFn::inverse_norm_cdf(delta) - 0.5 * total_volatility);
            let normalised_price = PriceBlackScholesNormalised::builder()
                .log_moneyness(log_moneyness)
                .total_volatility(total_volatility)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>();
            cases.push(
                ImpliedBlackVolatilityNormalised::builder()
                    .log_moneyness(log_moneyness)
                    .normalised_price(normalised_price)
                    .build_unchecked(),
            );
        }
    }
    cases
}

fn paper_grid_full_cases() -> Vec<ImpliedBlackVolatility> {
    let mut cases = Vec::with_capacity(PAPER_DELTAS.len() * 41);
    for total_volatility in paper_total_vols() {
        for delta in PAPER_DELTAS {
            let log_moneyness = total_volatility
                * (DefaultSpecialFn::inverse_norm_cdf(delta) - 0.5 * total_volatility);
            let forward = 1.0;
            let strike = (-log_moneyness).exp();
            let option_price = PriceBlackScholes::builder()
                .forward(forward)
                .strike(strike)
                .volatility(total_volatility)
                .expiry(1.0)
                .is_call(true)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>();
            cases.push(
                ImpliedBlackVolatility::builder()
                    .option_price(option_price)
                    .forward(forward)
                    .strike(strike)
                    .expiry(1.0)
                    .is_call(true)
                    .build_unchecked(),
            );
        }
    }
    cases
}

fn bench_normalised_current(b: &mut Bencher) {
    let cases = paper_grid_normalised_cases();
    b.iter(|| {
        let mut acc = 0.0;
        for case in test::black_box(&cases) {
            acc += test::black_box(case.calculate::<DefaultSpecialFn>().unwrap());
        }
        test::black_box(acc)
    });
}

fn bench_normalised_explicit(b: &mut Bencher) {
    let cases = paper_grid_normalised_cases();
    b.iter(|| {
        let mut acc = 0.0;
        for case in test::black_box(&cases) {
            acc += test::black_box(case.calculate_explicit::<DefaultSpecialFn>().unwrap());
        }
        test::black_box(acc)
    });
}

fn bench_full_current(b: &mut Bencher) {
    let cases = paper_grid_full_cases();
    b.iter(|| {
        let mut acc = 0.0;
        for case in test::black_box(&cases) {
            acc += test::black_box(case.calculate::<DefaultSpecialFn>().unwrap());
        }
        test::black_box(acc)
    });
}

fn bench_full_explicit(b: &mut Bencher) {
    let cases = paper_grid_full_cases();
    b.iter(|| {
        let mut acc = 0.0;
        for case in test::black_box(&cases) {
            acc += test::black_box(case.calculate_explicit::<DefaultSpecialFn>().unwrap());
        }
        test::black_box(acc)
    });
}

#[bench]
fn paper_grid_normalised_current(b: &mut Bencher) {
    bench_normalised_current(b);
}

#[bench]
fn paper_grid_normalised_explicit(b: &mut Bencher) {
    bench_normalised_explicit(b);
}

#[bench]
fn paper_grid_full_current(b: &mut Bencher) {
    bench_full_current(b);
}

#[bench]
fn paper_grid_full_explicit(b: &mut Bencher) {
    bench_full_explicit(b);
}
