#![feature(test)]

extern crate test;

use implied_vol::{DefaultSpecialFn, PriceBlackScholes};
use test::Bencher;

#[bench]
fn call_atm(b: &mut Bencher) {
    let f = 100.0;
    let k = f;
    let vol = 0.2;
    let t = 1.0;
    let q = true;
    let iv_builder = PriceBlackScholes::builder()
        .forward(f)
        .strike(k)
        .expiry(t)
        .volatility(vol)
        .is_call(q)
        .build()
        .unwrap();
    b.iter(|| iv_builder.calculate::<DefaultSpecialFn>());
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_atm_cpp(b: &mut Bencher) {
    let f = 100.0;
    let k = f;
    let vol = 0.2;
    let t = 1.0;
    let q = 1.0;
    b.iter(|| implied_vol::cxx::ffi::Black(f, k, t, vol, q));
}
