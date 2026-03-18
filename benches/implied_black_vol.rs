#![feature(test)]

use rand::RngExt;

extern crate test;

use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};
use test::Bencher;

#[bench]
fn call_atm(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = true;

    b.iter(|| {
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(test::black_box(price))
            .forward(test::black_box(f))
            .strike(test::black_box(k))
            .expiry(test::black_box(t))
            .is_call(test::black_box(q))
            .build()
            .unwrap();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>().unwrap())
    });
}

#[bench]
fn call_itm(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0f64.mul_add(1.0 - r, 1.0 * r * r2);
    let f = 1.0;
    let k = 1.0 * r;
    let t = 1e5 * r3;
    let q = true;
    b.iter(|| {
        let price = test::black_box(price);
        let f = test::black_box(f);
        let k = test::black_box(k);
        let t = test::black_box(t);
        let q = test::black_box(q);
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build_unchecked();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>().unwrap())
    });
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
    b.iter(|| {
        let price = test::black_box(price);
        let f = test::black_box(f);
        let k = test::black_box(k);
        let t = test::black_box(t);
        let q = test::black_box(q);
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build_unchecked();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>().unwrap())
    });
}

#[bench]
fn put_atm(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = false;
    b.iter(|| {
        let price = test::black_box(price);
        let f = test::black_box(f);
        let k = test::black_box(k);
        let t = test::black_box(t);
        let q = test::black_box(q);
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build_unchecked();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>().unwrap())
    });
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
    b.iter(|| {
        let price = test::black_box(price);
        let f = test::black_box(f);
        let k = test::black_box(k);
        let t = test::black_box(t);
        let q = test::black_box(q);
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build_unchecked();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>().unwrap())
    });
}

#[bench]
fn put_otm(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0f64.mul_add(1.0 - r, 1.0 * r * r2);
    let f = 1.0 * r;
    let k = 1.0;
    let t = 1e5 * r3;
    let q = false;
    b.iter(|| {
        let price = test::black_box(price);
        let f = test::black_box(f);
        let k = test::black_box(k);
        let t = test::black_box(t);
        let q = test::black_box(q);
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build_unchecked();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>().unwrap())
    });
}
#[cfg(feature = "cxx_bench")]
use implied_vol::cxx::ffi::ImpliedBlackVolatility;

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_atm_cpp(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = 1.0;
    b.iter(|| test::black_box(ImpliedBlackVolatility(price, f, k, t, q)));
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_itm_cpp(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0f64.mul_add(1.0 - r, 1.0 * r * r2);
    let f = 1.0;
    let k = 1.0 * r;
    let t = 1e5 * r3;
    let q = 1.0;
    b.iter(|| test::black_box(ImpliedBlackVolatility(price, f, k, t, q)));
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
    b.iter(|| test::black_box(ImpliedBlackVolatility(price, f, k, t, q)));
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn put_atm_cpp(b: &mut Bencher) {
    let price = 0.01;
    let f = 100.0;
    let k = f;
    let t = 1.0;
    let q = -1.0;
    b.iter(|| test::black_box(ImpliedBlackVolatility(price, f, k, t, q)));
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
    b.iter(|| test::black_box(ImpliedBlackVolatility(price, f, k, t, q)));
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn put_otm_cpp(b: &mut Bencher) {
    let seed: [u8; 32] = [13; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
    let (r, r2, r3): (f64, f64, f64) = rng.random();
    let price = 1.0f64.mul_add(1.0 - r, 1.0 * r * r2);
    let f = 1.0 * r;
    let k = 1.0;
    let t = 1e5 * r3;
    let q = -1.0;
    b.iter(|| test::black_box(ImpliedBlackVolatility(price, f, k, t, q)));
}
