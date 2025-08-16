#![feature(test)]

use rand::Rng;

extern crate test;

use implied_vol::{DefaultSpecialFn, ImpliedNormalVolatility};
use test::Bencher;

#[bench]
fn call_atm(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = true;
    let iv_builder = ImpliedNormalVolatility::builder()
        .option_price(price)
        .forward(f)
        .strike(k)
        .expiry(t)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>().unwrap());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_atm_cpp(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = 1.0;
    b.iter(|| implied_vol::cxx::ffi::ImpliedNormalVolatility(f, k, t, q, price));
}

#[bench]
fn call_itm(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
    let f = 1.0;
    let k = 1.0 * r;
    let t = 1e5 * r3;
    let q = true;
    let iv_builder = ImpliedNormalVolatility::builder()
        .option_price(price)
        .forward(f)
        .strike(k)
        .expiry(t)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>().unwrap());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_itm_cpp(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
    let f = 1.0;
    let k = 1.0 * r;
    let t = 1e5 * r3;
    let q = 1.0;
    b.iter(|| implied_vol::cxx::ffi::ImpliedNormalVolatility(f, k, t, q, price));
}

#[bench]
fn call_otm(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * r * r2;
    let f = 1.0 * r;
    let k = 1.0;
    let t = 1e5 * r3;
    let q = true;
    let iv_builder = ImpliedNormalVolatility::builder()
        .option_price(price)
        .forward(f)
        .strike(k)
        .expiry(t)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>().unwrap());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_otm_cpp(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * r * r2;
    let f = 1.0 * r;
    let k = 1.0;
    let t = 1e5 * r3;
    let q = 1.0;
    b.iter(|| implied_vol::cxx::ffi::ImpliedNormalVolatility(f, k, t, q, price));
}

#[bench]
fn put_atm(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = false;
    let iv_builder = ImpliedNormalVolatility::builder()
        .option_price(price)
        .forward(f)
        .strike(k)
        .expiry(t)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>().unwrap());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn put_atm_cpp(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = -1.0;
    b.iter(|| implied_vol::cxx::ffi::ImpliedNormalVolatility(f, k, t, q, price));
}

#[bench]
fn put_itm(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * r * r2;
    let f = 1.0;
    let k = 1.0 * r;
    let t = 1e5 * r3;
    let q = false;
    let iv_builder = ImpliedNormalVolatility::builder()
        .option_price(price)
        .forward(f)
        .strike(k)
        .expiry(t)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>().unwrap());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn put_itm_cpp(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * r * r2;
    let f = 1.0;
    let k = 1.0 * r;
    let t = 1e5 * r3;
    let q = -1.0;
    b.iter(|| implied_vol::cxx::ffi::ImpliedNormalVolatility(f, k, t, q, price));
}

#[bench]
fn put_otm(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
    let f = 1.0 * r;
    let k = 1.0;
    let t = 1e5 * r3;
    let q = false;
    let iv_builder = ImpliedNormalVolatility::builder()
        .option_price(price)
        .forward(f)
        .strike(k)
        .expiry(t)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>().unwrap());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn put_otm_cpp(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
    let f = 1.0 * r;
    let k = 1.0;
    let t = 1e5 * r3;
    let q = -1.0;
    b.iter(|| implied_vol::cxx::ffi::ImpliedNormalVolatility(f, k, t, q, price));
}
