use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, PriceBlackScholes,
    PriceBlackScholesNormalised,
};

#[test]
fn adjacent_forward_and_strike_remain_distinguishable() {
    for forward in [1e-300_f64, 1e-100, 1e100, 1e300] {
        let strike = forward.next_up();
        for is_call in [true, false] {
            let price = PriceBlackScholes::builder()
                .forward(forward)
                .strike(strike)
                .volatility(0.2)
                .expiry(1.0)
                .is_call(is_call)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>();
            assert!(price.is_finite() && price > 0.0);
            let recovered = ImpliedBlackVolatility::builder()
                .option_price(price)
                .forward(forward)
                .strike(strike)
                .expiry(1.0)
                .is_call(is_call)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>()
                .unwrap();
            assert!((recovered - 0.2).abs() <= 8.0 * f64::EPSILON);
        }
    }
}

#[test]
fn zero_expiry_only_has_intrinsic_price() {
    for (forward, strike) in [(100.0, 100.0), (120.0, 100.0), (100.0, 120.0)] {
        for is_call in [true, false] {
            let intrinsic = f64::max(
                if is_call {
                    forward - strike
                } else {
                    strike - forward
                },
                0.0,
            );
            for price in [intrinsic, intrinsic + 1.0, forward, strike] {
                let iv = ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(forward)
                    .strike(strike)
                    .expiry(0.0)
                    .is_call(is_call)
                    .build()
                    .unwrap();
                let expected = (price == intrinsic).then_some(0.0);
                assert_eq!(iv.calculate::<DefaultSpecialFn>(), expected);
                assert_eq!(iv.calculate_explicit::<DefaultSpecialFn>(), expected);
            }
        }
    }
}

#[test]
fn zero_and_overflowed_total_volatility_have_finite_black_prices() {
    for (forward, strike) in [(100.0, 100.0), (120.0, 100.0), (100.0, 120.0)] {
        for is_call in [true, false] {
            for (volatility, expiry) in
                [(0.0, f64::INFINITY), (f64::INFINITY, 0.0), (f64::MAX, 4.0)]
            {
                let price = PriceBlackScholes::builder()
                    .forward(forward)
                    .strike(strike)
                    .volatility(volatility)
                    .expiry(expiry)
                    .is_call(is_call)
                    .build()
                    .unwrap()
                    .calculate::<DefaultSpecialFn>();
                let expected = if volatility == 0.0 || expiry == 0.0 {
                    f64::max(
                        if is_call {
                            forward - strike
                        } else {
                            strike - forward
                        },
                        0.0,
                    )
                } else if is_call {
                    forward
                } else {
                    strike
                };
                assert_eq!(price, expected);
            }
        }
    }
}

#[test]
fn tiny_volatility_does_not_produce_nan() {
    let volatility = f64::from_bits(1);
    for (forward, strike, expected) in [(100.0, 110.0, 0.0), (110.0, 100.0, 10.0)] {
        let price = PriceBlackScholes::builder()
            .forward(forward)
            .strike(strike)
            .volatility(volatility)
            .expiry(1.0)
            .is_call(true)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        assert_eq!(price, expected);
    }
    let price = PriceBlackScholesNormalised::builder()
        .log_moneyness(-1.0)
        .total_volatility(volatility)
        .build()
        .unwrap()
        .calculate::<DefaultSpecialFn>();
    assert_eq!(price, 0.0);
}

#[test]
fn subnormal_log_moneyness_can_roundtrip() {
    for x in [f64::from_bits(1), f64::MIN_POSITIVE] {
        for v in [0.2, 3.0] {
            let price = PriceBlackScholesNormalised::builder()
                .log_moneyness(x)
                .total_volatility(v)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>();
            let recovered = ImpliedBlackVolatilityNormalised::builder()
                .log_moneyness(x)
                .normalised_price(price)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>()
                .unwrap();
            assert!((recovered - v).abs() <= 16.0 * f64::EPSILON * v);
        }
    }
}
