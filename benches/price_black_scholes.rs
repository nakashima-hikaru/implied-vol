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
    b.iter(|| {
        let f = test::black_box(f);
        let k = test::black_box(k);
        let vol = test::black_box(vol);
        let t = test::black_box(t);
        let q = test::black_box(q);
        let iv_builder = PriceBlackScholes::builder()
            .forward(f)
            .strike(k)
            .expiry(t)
            .volatility(vol)
            .is_call(q)
            .build_unchecked();
        test::black_box(iv_builder.calculate::<DefaultSpecialFn>())
    });
}

#[cfg(feature = "cxx_bench")]
#[bench]
fn call_atm_cpp(b: &mut Bencher) {
    let f = 100.0;
    let k = f;
    let vol = 0.2;
    let t = 1.0;
    let q = 1.0;
    b.iter(|| {
        test::black_box(implied_vol::cxx::ffi::Black(
            test::black_box(f),
            test::black_box(k),
            test::black_box(vol),
            test::black_box(t),
            test::black_box(q),
        ))
    });
}
