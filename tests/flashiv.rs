//! The default hybrid retains its conditioned precision qualification. Pure
//! `FlashIV` has separate Algorithm 1 evaluation-count and objective-roundoff
//! regressions: direct erfcx subtraction need not achieve Jaeckel's bound.

#[cfg(feature = "flashiv")]
#[path = "support/flashiv.rs"]
mod paper_error;
#[cfg(feature = "flashiv")]
use paper_error::precision_limit as paper_precision_limit;

#[cfg(feature = "flashiv")]
use implied_vol::solver::FlashIv;
use implied_vol::solver::{BlackSolver, Hybrid};
use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, PriceBlackScholes,
    PriceBlackScholesNormalised, SpecialFn,
};
use std::marker::PhantomData;

fn price(x: f64, s: f64) -> f64 {
    PriceBlackScholesNormalised::builder()
        .log_moneyness(x)
        .total_volatility(s)
        .build()
        .unwrap()
        .calculate::<DefaultSpecialFn>()
}

const fn inverse(x: f64, beta: f64) -> ImpliedBlackVolatilityNormalised {
    ImpliedBlackVolatilityNormalised::builder()
        .log_moneyness(x)
        .normalised_price(beta)
        .build()
        .unwrap()
}

fn conditioning(x: f64, beta: f64, s: f64) -> f64 {
    f64::EPSILON * (1.0 + price_error_as_volatility_error(x, beta, s))
}

// Match the common solver-mode bounds without changing their rounding order.
#[allow(clippy::suboptimal_flops)]
fn price_error_as_volatility_error(x: f64, error: f64, s: f64) -> f64 {
    if error == 0.0 {
        return 0.0;
    }
    let h = x / s;
    let t = 0.5 * s;
    // Evaluate error/(s*vega) together: inverse vega alone can overflow for
    // normal prices near the cap even when the condition number is finite.
    (error.ln() + 0.5 * (2.0 * std::f64::consts::PI).ln() + 0.5 * h.mul_add(h, t * t) - s.ln())
        .exp()
}

fn check_solver_price_input<S: BlackSolver>(x: f64, beta: f64) -> f64 {
    check_price_with_limit::<S>(x, beta, hybrid_precision_limit)
}

fn hybrid_precision_limit(x: f64, beta: f64, s: f64) -> f64 {
    4.0 * conditioning(x, beta, s)
}

fn check_price_with_limit<S: BlackSolver>(
    x: f64,
    beta: f64,
    precision_limit: fn(f64, f64, f64) -> f64,
) -> f64 {
    assert!(beta >= f64::MIN_POSITIVE && beta < (-0.5 * x.abs()).exp());
    let s = inverse(x, beta).calculate_with::<S>().unwrap();
    assert!(s.is_finite() && s > 0.0, "x={x}, beta={beta}, s={s}");
    let reprice = price(x, s);
    let residual_error = price_error_as_volatility_error(x, (reprice - beta).abs(), s);
    let limit = precision_limit(x, beta, s);
    assert!(residual_error.is_finite() && limit.is_finite());
    assert!(
        residual_error <= limit,
        "x={x:.17e}, beta={beta:.17e}, s={s:.17e}, reprice={reprice:.17e}, \
         residual_error={residual_error:.3e}, limit={limit:.3e}"
    );
    assert_eq!(
        s.to_bits(),
        inverse(-x, beta).calculate_with::<S>().unwrap().to_bits(),
        "normalized time-value inversion must preserve reciprocal-strike symmetry"
    );
    s
}

fn check_solver_roundtrip<S: BlackSolver>(x: f64, expected_s: f64) {
    let beta = price(x, expected_s);
    if beta < f64::MIN_POSITIVE {
        return;
    }
    let s = check_solver_price_input::<S>(x, beta);
    let input_error = (s - expected_s).abs() / expected_s;
    let limit = 4.0 * conditioning(x, beta, expected_s);
    assert!(
        input_error <= limit,
        "x={x:.17e}, beta={beta:.17e}, expected_s={expected_s:.17e}, \
         s={s:.17e}, input_error={input_error:.3e}, limit={limit:.3e}"
    );
}

#[cfg(feature = "flashiv")]
fn check_paper_roundtrip(x: f64, expected_s: f64) {
    let beta = price(x, expected_s);
    if beta < f64::MIN_POSITIVE || beta >= (-0.5 * x.abs()).exp() {
        return;
    }
    let s = check_price_with_limit::<FlashIv>(x, beta, paper_precision_limit);
    let error = (s - expected_s).abs() / expected_s;
    let limit = paper_precision_limit(x, beta, expected_s);
    assert!(
        error <= limit,
        "FlashIV x={x}, beta={beta}, expected_s={expected_s}, s={s}, error={error:.3e}, limit={limit:.3e}"
    );
}

fn check_price_input(x: f64, beta: f64) {
    check_solver_price_input::<Hybrid>(x, beta);
    #[cfg(feature = "flashiv")]
    check_price_with_limit::<FlashIv>(x, beta, paper_precision_limit);
}

fn check_roundtrip(x: f64, expected_s: f64) {
    check_solver_roundtrip::<Hybrid>(x, expected_s);
    #[cfg(feature = "flashiv")]
    check_paper_roundtrip(x, expected_s);
}

#[test]
fn near_atm_low_price_roots_match_each_solver_accuracy_contract() {
    // These are actual Black prices, not arbitrary iteration states. The price
    // can be exponentially smaller than |x|, where the Bachelier seed is close
    // to |x| even though the true total volatility can be as small as |x|/37.
    check_near_atm_roots(&[-1e-8, -1e-7, -1e-6, -1e-5, -1e-4, -0.001, -0.009, -0.01]);
}

#[cfg(feature = "flashiv")]
#[test]
fn flashiv_microscopic_terminal_branch_handles_tiny_normal_prices() {
    for x in [-1e-300_f64, -1e-200, -1e-100, -1e-32, -1e-16, -1e-12] {
        for ratio in [
            0.8, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 37.0,
        ] {
            check_paper_roundtrip(x, x.abs() / ratio);
        }
    }
}

fn check_near_atm_roots(xs: &[f64]) {
    for &x in xs {
        for ratio in [
            0.8, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 37.0,
        ] {
            check_roundtrip(x, x.abs() / ratio);
        }
    }
}

#[test]
fn low_price_dispatch_boundary_respects_solver_accuracy_contracts() {
    // Obtain the mathematical lower-branch boundary from the tangent at the
    // Black inflection point. Independently price neighboring volatilities;
    // do not duplicate the implementation's fitted b_l/b_max polynomial.
    for x in [
        -1e-8_f64, -1e-6, -1e-4, -0.009, -0.01, -0.5, -3.0, -10.0, -191.0,
    ] {
        let s_c = (-2.0 * x).sqrt();
        let s_l = (0.5 * std::f64::consts::PI)
            .sqrt()
            .mul_add(-DefaultSpecialFn::one_minus_erfcx((-x).sqrt()), s_c);
        assert!(s_l > 0.0);
        for s in [s_l.next_down(), s_l, s_l.next_up()] {
            check_roundtrip(x, s);
        }
        let boundary_price = price(x, s_l);
        for beta in [
            boundary_price.next_down(),
            boundary_price,
            boundary_price.next_up(),
        ] {
            check_price_input(x, beta);
        }
    }
}

#[test]
fn price_inputs_across_seed_and_domain_seams_remain_accurate() {
    // c is the actual out-of-the-money Black call price divided by the forward.
    // Neighboring representable inputs exercise seed/domain seams without
    // depending on which implementation branch is taken.
    for abs_x in [
        0.01_f64.next_down(),
        0.01,
        0.01_f64.next_up(),
        3.0_f64.next_down(),
        3.0,
        3.0_f64.next_up(),
    ] {
        let b_max = (-0.5 * abs_x).exp();
        for c in [0.0005_f64.next_down(), 0.0005, 0.0005_f64.next_up()] {
            check_price_input(-abs_x, b_max * c);
        }
    }
    for abs_x in [1e-8_f64.next_down(), 1e-8, 1e-8_f64.next_up()] {
        let b_max = (-0.5 * abs_x).exp();
        for c in [1e-6_f64.next_down(), 1e-6, 1e-6_f64.next_up()] {
            check_price_input(-abs_x, b_max * c);
        }
    }
    // Complementary-price iteration must also span the former upper guard.
    for abs_x in [0.01_f64, 0.5, 3.0] {
        let b_max = (-0.5 * abs_x).exp();
        for c in [0.99_f64.next_down(), 0.99, 0.99_f64.next_up()] {
            check_price_input(-abs_x, b_max * c);
        }
    }
}

#[test]
fn full_call_and_put_inputs_recover_low_price_volatility() {
    for (abs_x, total_volatility) in [(0.5_f64, 0.125_f64), (0.009, 0.009 / 30.0)] {
        for scale in [1.0_f64, 1e100] {
            for is_call in [true, false] {
                let smaller = scale * (-abs_x).exp();
                let (forward, strike) = if is_call {
                    (smaller, scale)
                } else {
                    (scale, smaller)
                };
                let expiry: f64 = 0.25;
                let sigma = total_volatility / expiry.sqrt();
                let observed = PriceBlackScholes::builder()
                    .forward(forward)
                    .strike(strike)
                    .expiry(expiry)
                    .is_call(is_call)
                    .volatility(sigma)
                    .build()
                    .unwrap()
                    .calculate::<DefaultSpecialFn>();
                assert!(observed.is_finite() && observed > 0.0);
                let input = ImpliedBlackVolatility::builder()
                    .forward(forward)
                    .strike(strike)
                    .expiry(expiry)
                    .is_call(is_call)
                    .option_price(observed)
                    .build()
                    .unwrap();
                let actuals = [
                    input.calculate_with::<Hybrid>().unwrap(),
                    input
                        .calculate_with::<implied_vol::solver::Jaeckel>()
                        .unwrap(),
                ];
                #[cfg(feature = "flashiv")]
                {
                    let actual = input.calculate_with::<FlashIv>().unwrap();
                    let beta = observed / (forward.sqrt() * strike.sqrt());
                    let error = (actual - sigma).abs() / sigma;
                    let limit = paper_precision_limit(-abs_x, beta, total_volatility);
                    assert!(
                        error <= limit,
                        "FlashIV sigma={sigma}, actual={actual}, error={error:.3e}, limit={limit:.3e}"
                    );
                }
                for actual in actuals {
                    let beta = observed / (forward.sqrt() * strike.sqrt());
                    let limit = 4.0 * conditioning(-abs_x, beta, total_volatility);
                    let error = (actual - sigma).abs() / sigma;
                    assert!(
                        error <= limit,
                        "F={forward}, K={strike}, call={is_call}, sigma={sigma}, \
                     actual={actual}, error={error:.3e}, limit={limit:.3e}"
                    );
                }
            }
        }
    }
}

// A borrowed type intentionally has no 'static constraint: calculate's generic
// SpecialFn contract permits providers containing a caller's lifetime.
struct BorrowedSpecialFn<'a>(PhantomData<&'a ()>);

macro_rules! delegate {
    ($name:ident) => {
        #[inline]
        fn $name(x: f64) -> f64 {
            DefaultSpecialFn::$name(x)
        }
    };
}

impl SpecialFn for BorrowedSpecialFn<'_> {
    delegate!(erf);
    delegate!(erfc);
    delegate!(erfcx);
    delegate!(erfinv);
    delegate!(inverse_norm_cdf);
    delegate!(norm_cdf);
    delegate!(one_minus_erfcx);
}

// The reference intentionally binds the provider's lifetime to a local scope.
#[allow(clippy::trivially_copy_pass_by_ref)]
fn calculate_with_borrowed_provider<'a>(
    _scope: &'a (),
    iv: &ImpliedBlackVolatilityNormalised,
) -> f64 {
    iv.calculate_with::<Hybrid<BorrowedSpecialFn<'a>>>()
        .unwrap()
}

#[test]
fn feature_preserves_generic_special_function_providers() {
    let scope = ();
    for (x, s) in [(-0.5, 0.125), (-0.009, 0.009 / 30.0), (0.0, 0.2)] {
        let iv = inverse(x, price(x, s));
        let default = iv.calculate::<DefaultSpecialFn>().unwrap();
        let borrowed = calculate_with_borrowed_provider(&scope, &iv);
        assert_eq!(borrowed.to_bits(), default.to_bits());
    }
}

// Count full-precision provider calls independently of the solver's internal
// helpers. The FAST pre-step uses its own A&S polynomial and contributes none.
#[cfg(feature = "flashiv")]
thread_local! {
    static SPECIAL_CALLS: std::cell::Cell<[usize; 3]> = const { std::cell::Cell::new([0; 3]) };
}

#[cfg(feature = "flashiv")]
struct CountingSpecialFn;

#[cfg(feature = "flashiv")]
impl SpecialFn for CountingSpecialFn {
    fn erf(x: f64) -> f64 {
        SPECIAL_CALLS.with(|calls| {
            let mut n = calls.get();
            n[0] += 1;
            calls.set(n);
        });
        DefaultSpecialFn::erf(x)
    }
    fn erfc(x: f64) -> f64 {
        SPECIAL_CALLS.with(|calls| {
            let mut n = calls.get();
            n[1] += 1;
            calls.set(n);
        });
        DefaultSpecialFn::erfc(x)
    }
    fn erfcx(x: f64) -> f64 {
        SPECIAL_CALLS.with(|calls| {
            let mut n = calls.get();
            n[2] += 1;
            calls.set(n);
        });
        DefaultSpecialFn::erfcx(x)
    }
    delegate!(erfinv);
    delegate!(inverse_norm_cdf);
    delegate!(norm_cdf);
    delegate!(one_minus_erfcx);
}

#[cfg(feature = "flashiv")]
fn counted_flashiv(x: f64, beta: f64) -> (f64, [usize; 3]) {
    SPECIAL_CALLS.with(|calls| calls.set([0; 3]));
    let s = inverse(x, beta)
        .calculate_with::<FlashIv<CountingSpecialFn>>()
        .unwrap();
    let calls = SPECIAL_CALLS.with(std::cell::Cell::get);
    (s, calls)
}

#[cfg(feature = "flashiv")]
#[test]
fn flashiv_ordinary_chain_uses_two_exact_steps_and_a_conditional_third() {
    let (s, ordinary_calls) = counted_flashiv(-0.5, price(-0.5, 0.7));
    assert_eq!(
        ordinary_calls,
        [0, 0, 4],
        "two exact steps, two erfcx values each"
    );
    assert!((s - 0.7).abs() <= 4.0 * f64::EPSILON);

    // The small-price seed is deliberately rough for this deep near-ATM tail;
    // the residual entering the second exact step triggers the published third.
    let (_, safety_calls) = counted_flashiv(-0.009, price(-0.009, 0.0009));
    assert_eq!(safety_calls, [0, 0, 6]);
}

#[cfg(feature = "flashiv")]
#[test]
fn flashiv_upper_chain_uses_three_complementary_halley_steps() {
    for x in [-0.01_f64, -0.5, -3.0] {
        let b_max = (0.5 * x).exp();
        let (_, calls) = counted_flashiv(x, 0.995 * b_max);
        assert_eq!(
            calls,
            [0, 0, 6],
            "upper steps evaluate the erfcx sum directly"
        );
    }
}

#[cfg(feature = "flashiv")]
#[test]
fn flashiv_microscopic_guard_is_terminal_before_the_general_volatility_floor() {
    // A central Bachelier seed receives at most two expansion Newton steps,
    // which call erfc and never the ordinary erfcx objective.
    let (central, calls) = counted_flashiv(-1e-200, price(-1e-200, 1e-200));
    assert!(central > 0.0 && central < 1e-10);
    assert!(calls[1] >= 1 && calls[1] <= 2);
    assert_eq!(calls[0], 0);
    assert_eq!(calls[2], 0);

    // The Mills branch solves the normal tail directly in at most four steps.
    let (tail, calls) = counted_flashiv(-1e-100, price(-1e-100, 1e-101));
    assert!(tail > 0.0 && tail < 1e-10);
    assert_eq!(calls[0], 0);
    assert_eq!(calls[1], 0);
    assert!(calls[2] > 0 && calls[2] <= 4);
}

#[cfg(feature = "flashiv")]
#[test]
fn flashiv_direct_objective_retains_its_published_near_atm_rounding_limit() {
    // This input is just outside the microscopic guard. The ordinary paper
    // objective subtracts nearly equal erfcx values, losing digits that the
    // default hybrid's stable expansion preserves. Keep the accuracy policy
    // explicit rather than claiming an attainable-precision guarantee here.
    let x = -1e-7;
    let expected_s = 1.25e-7;
    let beta = price(x, expected_s);
    let s = check_price_with_limit::<FlashIv>(x, beta, paper_precision_limit);
    assert!((s - expected_s).abs() / expected_s < 1e-8);
    check_solver_roundtrip::<Hybrid>(x, expected_s);
}
