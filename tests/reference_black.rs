use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatilityNormalised};

fn implied_volatility(x: f64, beta: f64) -> f64 {
    ImpliedBlackVolatilityNormalised::builder()
        .log_moneyness(x)
        .normalised_price(beta)
        .build()
        .unwrap()
        .calculate::<DefaultSpecialFn>()
        .unwrap()
}

fn conditioned_accuracy(x: f64, beta: f64, s: f64) -> f64 {
    let h = x / s;
    let t = 0.5 * s;
    let inverse_vega = (2.0 * std::f64::consts::PI).sqrt() * (0.5 * h.mul_add(h, t * t)).exp();
    f64::EPSILON * (1.0 + beta * inverse_vega / s)
}

#[test]
fn normalised_solver_matches_high_precision_roots() {
    // Roots of exp(x/2) Phi(x/s+s/2) - exp(-x/2) Phi(x/s-s/2) = beta,
    // evaluated with 100 decimal digits, using the exact binary64 x and beta.
    // These include cases where both the original C++ solver and the Rust port
    // exceed a one-epsilon conditioning estimate without final refinement.
    let cases = [
        (
            -1.576_447_743_634_572,
            6.276_826_885_269_499e-3,
            7.878_278_148_479_239e-1,
        ),
        (
            -8.513_962_679_747_584e-1,
            2.328_973_812_973_88e-2,
            6.198_729_138_022_668e-1,
        ),
        (
            -2.324_962_302_403_283e-1,
            1.125_940_199_629_584_8e-1,
            5.281_269_791_297_323e-1,
        ),
        (-1e-8, 1.470_539_346_550_702e-10, 6.266_099_313_387_908e-9),
        (-1e-12, 1.471_060_431_232_883e-14, 6.266_565_974_016_432e-13),
    ];
    for (x, beta, reference_s) in cases {
        let actual_s = implied_volatility(x, beta);
        let conditioning = conditioned_accuracy(x, beta, reference_s);
        assert!(
            (actual_s - reference_s).abs() / reference_s <= 4.0 * conditioning,
            "x={x:.17e}, beta={beta:.17e}, actual={actual_s:.17e}, root={reference_s:.17e}"
        );
    }
}

#[test]
fn historical_regressions_match_high_precision_roots() {
    // Equation: exp(x/2) Phi(x/s+s/2) - exp(-x/2) Phi(x/s-s/2) = beta,
    // the normalized Black formula in source/lets_be_rational.cpp. Phi(z) was
    // evaluated independently as erfc(-z/sqrt(2))/2 with mpmath 1.4.1. Inputs
    // were converted through integer ratios, preserving every binary64 bit.
    // Bisection at both 100 and 200 decimal digits gave identical rounded roots
    // and correction bits; the latter retain each root's sub-ulp remainder.
    // This avoids using a floating-price plateau as the inverse's oracle.
    let cases = [
        (
            0xbfdd_33e9_3b98_5c83,
            0x1414_c40c_5d29_351e,
            0x3f8e_5705_cb46_0f03,
            0xbc27_7d48_2038_8541,
        ),
        (
            0xbfe4_9fdc_5054_7152,
            0x3faf_1129_b89a_7174,
            0x3fe5_cf46_c34e_3fb5,
            0x3c81_7421_af65_a18b,
        ),
        (
            0xbfe4_9fdc_5054_7152,
            0x3faf_1129_b89a_716c,
            0x3fe5_cf46_c34e_3fb3,
            0x3c79_08c6_9e9e_6fc2,
        ),
    ];
    for (x_bits, beta_bits, root_bits, correction_bits) in cases {
        let x = f64::from_bits(x_bits);
        let beta = f64::from_bits(beta_bits);
        let reference_s = f64::from_bits(root_bits);
        let root_correction = f64::from_bits(correction_bits);
        let actual_s = implied_volatility(x, beta);
        let root_error = ((actual_s - reference_s) - root_correction).abs() / reference_s;
        let conditioning = conditioned_accuracy(x, beta, reference_s);
        // Preserve the historical one-estimate limit for these regressions.
        assert!(
            root_error <= conditioning,
            "x={x:.17e}, beta={beta:.17e}, actual={actual_s:.17e}, root={reference_s:.17e}, root_error={root_error:.3e}, conditioning={conditioning:.3e}"
        );
    }
}

#[cfg(feature = "cxx_bench")]
mod cpp_reference {
    use super::*;
    use implied_vol::{PriceBlackScholesNormalised, SpecialFn, cxx::ffi};
    use rand::RngExt;

    fn normalised_price(x: f64, s: f64) -> f64 {
        PriceBlackScholesNormalised::builder()
            .log_moneyness(x)
            .total_volatility(s)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>()
    }

    fn assert_matches_reference(x: f64, beta: f64) {
        let actual_s = implied_volatility(x, beta);
        let reference_s = ffi::NormalisedImpliedBlackVolatility(beta, x, 1.0);
        assert!(actual_s.is_finite() && actual_s > 0.0);
        assert!(reference_s.is_finite() && reference_s > 0.0);

        // The attainable-accuracy formula is a first-order conditioning estimate,
        // not a bound on all floating-point evaluations. Independent 100-digit
        // roots above show errors above one estimate in both implementations.
        // Allow each solver the four-estimate envelope used by roundtrip tests.
        // Evaluate the same formula from beta so boundary rounding at DBL_MIN
        // cannot invoke C++'s subnormal fallback of a 100% relative tolerance.
        let conditioning = conditioned_accuracy(x, beta, reference_s);
        let reference_conditioning = ffi::ImpliedVolatilityAttainableAccuracy(x, reference_s, 1.0);
        assert!(conditioning.is_finite() && conditioning > 0.0);
        assert!(reference_conditioning.is_finite() && reference_conditioning > 0.0);
        let tolerance = 8.0 * conditioning.min(reference_conditioning);
        assert!(
            (actual_s - reference_s).abs() / reference_s <= tolerance,
            "x={x:.17e}, beta={beta:.17e}, actual={actual_s:.17e}, reference={reference_s:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn normalised_solver_matches_reference_branches_and_boundaries() {
        for x in [
            -1e-12, -1e-8, -1e-4, -1e-2, -0.5, -3.0, -8.0, -190.0, -191.0, -580.0, -581.0, -1000.0,
        ] {
            let s_c = (-2.0_f64 * x).sqrt();
            let ome = DefaultSpecialFn::one_minus_erfcx((-x).sqrt());
            let half_sqrt_pi = (0.5 * std::f64::consts::PI).sqrt();
            let s_l = s_c - half_sqrt_pi * ome;
            let s_u = s_c + half_sqrt_pi * (2.0 - ome);
            let b_max = (0.5 * x).exp();

            for s in [
                0.5 * s_l,
                s_l.next_down(),
                s_l,
                s_l.next_up(),
                s_l.midpoint(s_c),
                s_c.next_down(),
                s_c,
                s_c.next_up(),
                s_c.midpoint(s_u),
                s_u.next_down(),
                s_u,
                s_u.next_up(),
                s_u + 0.25,
            ] {
                // Use prices from each implementation as shared solver inputs.
                for beta in [normalised_price(x, s), ffi::NormalisedBlack(x, s, 1.0)] {
                    if beta >= f64::MIN_POSITIVE && beta < b_max {
                        assert_matches_reference(x, beta);
                    }
                }
            }
            assert_matches_reference(x, f64::MIN_POSITIVE);
            assert_matches_reference(x, b_max.next_down());
        }
    }

    #[test]
    fn normalised_solver_matches_reference_randomised() {
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed([42; 32]);
        let mut checked = 0;
        for _ in 0..20_000 {
            let x = -rng.random_range(1e-4..5.0);
            let s = rng.random_range(1e-3..5.0);
            let beta = normalised_price(x, s);
            if beta >= f64::MIN_POSITIVE && beta < (0.5 * x).exp() {
                assert_matches_reference(x, beta);
                checked += 1;
            }
        }
        assert!(checked > 19_000);
    }
}
