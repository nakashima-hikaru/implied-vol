#![feature(test)]

extern crate test;

use implied_vol::{DefaultSpecialFn, PriceBlackScholes};
use test::Bencher;

#[path = "support/mod.rs"]
mod support;

use support::PriceCase;

fn bench_build_validated(b: &mut Bencher, case: PriceCase) {
    b.iter(|| {
        test::black_box(
            PriceBlackScholes::builder()
                .forward(test::black_box(case.forward))
                .strike(test::black_box(case.strike))
                .volatility(test::black_box(case.volatility))
                .expiry(test::black_box(case.expiry))
                .is_call(test::black_box(case.is_call))
                .build(),
        )
    });
}

fn bench_build_unchecked(b: &mut Bencher, case: PriceCase) {
    b.iter(|| {
        test::black_box(
            PriceBlackScholes::builder()
                .forward(test::black_box(case.forward))
                .strike(test::black_box(case.strike))
                .volatility(test::black_box(case.volatility))
                .expiry(test::black_box(case.expiry))
                .is_call(test::black_box(case.is_call))
                .build_unchecked(),
        )
    });
}

fn bench_calculate_only(b: &mut Bencher, case: PriceCase) {
    let price = PriceBlackScholes::builder()
        .forward(case.forward)
        .strike(case.strike)
        .volatility(case.volatility)
        .expiry(case.expiry)
        .is_call(case.is_call)
        .build_unchecked();

    b.iter(|| test::black_box(test::black_box(&price).calculate::<DefaultSpecialFn>()));
}

fn bench_end_to_end_validated(b: &mut Bencher, case: PriceCase) {
    b.iter(|| {
        let price = PriceBlackScholes::builder()
            .forward(test::black_box(case.forward))
            .strike(test::black_box(case.strike))
            .volatility(test::black_box(case.volatility))
            .expiry(test::black_box(case.expiry))
            .is_call(test::black_box(case.is_call))
            .build()
            .unwrap();
        let price = test::black_box(price);
        test::black_box(price.calculate::<DefaultSpecialFn>())
    });
}

fn bench_end_to_end_unchecked(b: &mut Bencher, case: PriceCase) {
    b.iter(|| {
        let price = PriceBlackScholes::builder()
            .forward(test::black_box(case.forward))
            .strike(test::black_box(case.strike))
            .volatility(test::black_box(case.volatility))
            .expiry(test::black_box(case.expiry))
            .is_call(test::black_box(case.is_call))
            .build_unchecked();
        let price = test::black_box(price);
        test::black_box(price.calculate::<DefaultSpecialFn>())
    });
}

macro_rules! define_case_benches {
    (
        $validated:ident,
        $unchecked:ident,
        $calculate:ident,
        $e2e_validated:ident,
        $e2e_unchecked:ident,
        $cpp:ident,
        $case:expr
    ) => {
        #[bench]
        fn $validated(b: &mut Bencher) {
            bench_build_validated(b, $case);
        }

        #[bench]
        fn $unchecked(b: &mut Bencher) {
            bench_build_unchecked(b, $case);
        }

        #[bench]
        fn $calculate(b: &mut Bencher) {
            bench_calculate_only(b, $case);
        }

        #[bench]
        fn $e2e_validated(b: &mut Bencher) {
            bench_end_to_end_validated(b, $case);
        }

        #[bench]
        fn $e2e_unchecked(b: &mut Bencher) {
            bench_end_to_end_unchecked(b, $case);
        }

        #[cfg(feature = "cxx_bench")]
        #[bench]
        fn $cpp(b: &mut Bencher) {
            let case = $case;
            b.iter(|| {
                test::black_box(implied_vol::cxx::ffi::Black(
                    test::black_box(case.forward),
                    test::black_box(case.strike),
                    test::black_box(case.volatility),
                    test::black_box(case.expiry),
                    test::black_box(support::cxx_q(case.is_call)),
                ))
            });
        }
    };
}

define_case_benches!(
    atm_call_build_validated,
    atm_call_build_unchecked,
    atm_call_calculate_only,
    atm_call_end_to_end_validated,
    atm_call_end_to_end_unchecked,
    atm_call_cpp_direct,
    support::ATM_CALL
);

define_case_benches!(
    deep_itm_put_long_build_validated,
    deep_itm_put_long_build_unchecked,
    deep_itm_put_long_calculate_only,
    deep_itm_put_long_end_to_end_validated,
    deep_itm_put_long_end_to_end_unchecked,
    deep_itm_put_long_cpp_direct,
    support::DEEP_ITM_PUT_LONG
);

define_case_benches!(
    deep_otm_call_long_build_validated,
    deep_otm_call_long_build_unchecked,
    deep_otm_call_long_calculate_only,
    deep_otm_call_long_end_to_end_validated,
    deep_otm_call_long_end_to_end_unchecked,
    deep_otm_call_long_cpp_direct,
    support::DEEP_OTM_CALL_LONG
);

define_case_benches!(
    near_atm_call_short_build_validated,
    near_atm_call_short_build_unchecked,
    near_atm_call_short_calculate_only,
    near_atm_call_short_end_to_end_validated,
    near_atm_call_short_end_to_end_unchecked,
    near_atm_call_short_cpp_direct,
    support::NEAR_ATM_CALL_SHORT
);

define_case_benches!(
    edge_zero_expiry_call_build_validated,
    edge_zero_expiry_call_build_unchecked,
    edge_zero_expiry_call_calculate_only,
    edge_zero_expiry_call_end_to_end_validated,
    edge_zero_expiry_call_end_to_end_unchecked,
    edge_zero_expiry_call_cpp_direct,
    support::EDGE_ZERO_EXPIRY_CALL
);
