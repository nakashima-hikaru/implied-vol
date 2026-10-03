//! Solver selection, provider ownership, and additive feature availability.

#[cfg(feature = "experimental")]
use implied_vol::solver::Experimental;
#[cfg(feature = "flashiv")]
use implied_vol::solver::FlashIv;
use implied_vol::solver::{BlackSolver, Hybrid, Jaeckel};
use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, SpecialFn,
};
use std::marker::PhantomData;

const fn normalized(x: f64, beta: f64) -> ImpliedBlackVolatilityNormalised {
    ImpliedBlackVolatilityNormalised::builder()
        .log_moneyness(x)
        .normalised_price(beta)
        .build()
        .unwrap()
}

const fn full(f: f64, k: f64, t: f64, p: f64, call: bool) -> ImpliedBlackVolatility {
    ImpliedBlackVolatility::builder()
        .forward(f)
        .strike(k)
        .expiry(t)
        .option_price(p)
        .is_call(call)
        .build()
        .unwrap()
}

fn assert_same_result(actual: Option<f64>, expected: Option<f64>) {
    assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
}

#[test]
fn default_calculation_is_always_hybrid_in_both_public_apis() {
    for (x, beta) in [
        (0.0, 0.125),
        (-0.5, 1e-6),
        (-0.5, 0.25),
        (-3.0, 0.2),
        (-1416.0, f64::from_bits(0x0017_c8ab_2288_c9ab)),
        (-1.0, f64::from_bits(1)),
        (-1.0, 0.0),
        (-1.0, 0.9),
    ] {
        let input = normalized(x, beta);
        assert_same_result(
            input.calculate::<DefaultSpecialFn>(),
            input.calculate_with::<Hybrid>(),
        );
    }
    for call in [false, true] {
        for (f, k, t, p) in [
            (1.0, 1.0, 0.25, 0.125),
            (1.0, 2.0, 4.0, 0.125),
            (1.0, 1.0, 1.0, f64::from_bits(1)),
            (1.0, 1.0, 0.0, 0.0),
        ] {
            let input = full(f, k, t, p, call);
            assert_same_result(
                input.calculate::<DefaultSpecialFn>(),
                input.calculate_with::<Hybrid>(),
            );
        }
    }
}

fn check_selected_solver<S: BlackSolver>() {
    // The same containers can be evaluated repeatedly with different solver types.
    for (x, beta) in [(0.0, 0.125), (-0.5, 1e-6), (-0.5, 0.25), (-3.0, 0.2)] {
        let expected = S::implied_total_volatility(x, beta);
        assert!(expected.is_some_and(|s| s.is_finite() && s > 0.0));
        assert_same_result(normalized(x, beta).calculate_with::<S>(), expected);
    }
    for call in [false, true] {
        // F=K=1 retains beta exactly, and sqrt(T)=1/2 makes annualization exact.
        let expected = S::implied_total_volatility(0.0, 0.125).map(|s| 2.0 * s);
        assert_same_result(
            full(1.0, 1.0, 0.25, 0.125, call).calculate_with::<S>(),
            expected,
        );
    }
}

#[test]
fn all_available_solvers_coexist_in_one_build() {
    check_selected_solver::<Hybrid>();
    check_selected_solver::<Jaeckel>();
    #[cfg(feature = "flashiv")]
    check_selected_solver::<FlashIv>();
    #[cfg(feature = "experimental")]
    check_selected_solver::<Experimental>();
}

fn check_shared_boundaries<S: BlackSolver>() {
    for (x, beta) in [
        (f64::NAN, 0.1),
        (f64::INFINITY, 0.1),
        (-0.5, f64::NAN),
        (-0.5, f64::INFINITY),
        (-0.5, -0.1),
    ] {
        assert_eq!(S::implied_total_volatility(x, beta), None);
    }
    assert_eq!(normalized(-0.5, 0.0).calculate_with::<S>(), Some(0.0));
    assert_eq!(
        normalized(0.0, 1.0).calculate_with::<S>(),
        Some(f64::INFINITY)
    );
    assert_eq!(
        normalized(0.0, 1.0_f64.next_up()).calculate_with::<S>(),
        None
    );
    assert_eq!(normalized(-0.5, 0.9).calculate_with::<S>(), None);
    for call in [false, true] {
        let (f, k) = if call { (2.0, 1.0) } else { (1.0, 2.0) };
        assert_eq!(full(f, k, 1.0, 0.5, call).calculate_with::<S>(), None);
        assert_eq!(full(f, k, 1.0, 1.0, call).calculate_with::<S>(), Some(0.0));
        assert_eq!(
            full(f, k, 1.0, 2.0, call).calculate_with::<S>(),
            Some(f64::INFINITY)
        );
        assert_eq!(
            full(f, k, 1.0, 2.0_f64.next_up(), call).calculate_with::<S>(),
            None
        );
        assert_eq!(full(f, k, 0.0, 1.0, call).calculate_with::<S>(), Some(0.0));
        assert_eq!(full(f, k, 0.0, 1.25, call).calculate_with::<S>(), None);
    }
}

#[test]
fn invalid_prices_and_terminal_contracts_are_checked_for_each_solver() {
    check_shared_boundaries::<Hybrid>();
    check_shared_boundaries::<Jaeckel>();
    #[cfg(feature = "flashiv")]
    check_shared_boundaries::<FlashIv>();
    #[cfg(feature = "experimental")]
    check_shared_boundaries::<Experimental>();
}

struct BorrowedSpecialFn<'a>(PhantomData<&'a ()>);
macro_rules! delegate {
    ($name:ident) => {
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

#[allow(clippy::trivially_copy_pass_by_ref)]
fn check_borrowed_provider<'a>(_scope: &'a ()) {
    let normalized = normalized(-0.5, 0.125);
    let full = full(1.0, 1.0, 0.25, 0.125, true);
    assert_same_result(
        normalized.calculate::<BorrowedSpecialFn<'a>>(),
        normalized.calculate_with::<Hybrid>(),
    );
    assert_same_result(
        normalized.calculate_with::<Hybrid<BorrowedSpecialFn<'a>>>(),
        normalized.calculate_with::<Hybrid>(),
    );
    assert_same_result(
        normalized.calculate_with::<Jaeckel<BorrowedSpecialFn<'a>>>(),
        normalized.calculate_with::<Jaeckel>(),
    );
    assert_same_result(
        full.calculate::<BorrowedSpecialFn<'a>>(),
        full.calculate_with::<Hybrid>(),
    );
    assert_same_result(
        full.calculate_with::<Hybrid<BorrowedSpecialFn<'a>>>(),
        full.calculate_with::<Hybrid>(),
    );
    assert_same_result(
        full.calculate_with::<Jaeckel<BorrowedSpecialFn<'a>>>(),
        full.calculate_with::<Jaeckel>(),
    );
    #[cfg(feature = "flashiv")]
    {
        assert_same_result(
            normalized.calculate_with::<FlashIv<BorrowedSpecialFn<'a>>>(),
            normalized.calculate_with::<FlashIv>(),
        );
        assert_same_result(
            full.calculate_with::<FlashIv<BorrowedSpecialFn<'a>>>(),
            full.calculate_with::<FlashIv>(),
        );
    }
}

#[test]
fn solver_provider_types_accept_a_callers_lifetime() {
    let scope = ();
    check_borrowed_provider(&scope);
}

#[cfg(feature = "experimental")]
#[test]
fn experimental_exact_cap_does_not_replace_other_solver_boundaries() {
    let below_exact = normalized(-1416.0, f64::from_bits(0x0017_c8ab_2288_c9ab));
    assert_eq!(
        below_exact.calculate::<DefaultSpecialFn>(),
        Some(f64::INFINITY)
    );
    assert_eq!(below_exact.calculate_with::<Jaeckel>(), Some(f64::INFINITY));
    #[cfg(feature = "flashiv")]
    assert_eq!(below_exact.calculate_with::<FlashIv>(), Some(f64::INFINITY));
    assert!(
        below_exact
            .calculate_with::<Experimental>()
            .is_some_and(f64::is_finite)
    );

    let above_exact = normalized(-1.0, f64::from_bits(0x3fe3_68b2_fc6f_960a));
    assert_eq!(
        above_exact.calculate::<DefaultSpecialFn>(),
        Some(f64::INFINITY)
    );
    assert_eq!(above_exact.calculate_with::<Jaeckel>(), Some(f64::INFINITY));
    #[cfg(feature = "flashiv")]
    assert_eq!(above_exact.calculate_with::<FlashIv>(), Some(f64::INFINITY));
    assert_eq!(above_exact.calculate_with::<Experimental>(), None);
}

#[test]
fn full_interior_price_does_not_become_a_rounded_normalized_cap() {
    let strike = 3.0_f64;
    let p = 1.0_f64.next_down();
    let x = (1.0 / strike).ln();
    let beta = p / strike.sqrt();
    assert_eq!(beta, (0.5 * x).exp());
    let input = full(1.0, strike, 1.0, p, true);
    assert_eq!(input.calculate::<DefaultSpecialFn>(), None);
    assert_eq!(input.calculate_with::<Jaeckel>(), None);
    #[cfg(feature = "flashiv")]
    assert_eq!(input.calculate_with::<FlashIv>(), None);
    #[cfg(feature = "experimental")]
    assert!(
        input
            .calculate_with::<Experimental>()
            .is_some_and(f64::is_finite)
    );
}
