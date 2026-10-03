//! Public contracts for independently selected Black solvers.

use implied_vol::SpecialFn;
#[cfg(feature = "flashiv")]
use implied_vol::solver::FlashIv;
use implied_vol::solver::{BlackSolver, Hybrid, Jaeckel};
use implied_vol::{
    DefaultSpecialFn, ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised, PriceBlackScholes,
    PriceBlackScholesNormalised,
};

#[cfg(feature = "flashiv")]
#[path = "support/flashiv.rs"]
mod paper_error;

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

// Keep the logarithmic conditioning bound's evaluation order identical across
// solver modes rather than contracting different terms into fused operations.
#[allow(clippy::suboptimal_flops)]
fn price_error_as_volatility_error(x: f64, error: f64, s: f64) -> f64 {
    if error == 0.0 {
        return 0.0;
    }
    let h = x / s;
    let t = 0.5 * s;
    (error.ln() + 0.5 * (2.0 * std::f64::consts::PI).ln() + 0.5 * h.mul_add(h, t * t) - s.ln())
        .exp()
}

fn conditioning(x: f64, beta: f64, s: f64) -> f64 {
    f64::EPSILON * (1.0 + price_error_as_volatility_error(x, beta, s))
}

fn attainable_precision_limit(x: f64, beta: f64, s: f64) -> f64 {
    4.0 * conditioning(x, beta, s)
}

fn check_price_input<S: BlackSolver>(x: f64, beta: f64) {
    check_price_with_limit::<S>(x, beta, attainable_precision_limit);
}

fn check_price_with_limit<S: BlackSolver>(
    x: f64,
    beta: f64,
    precision_limit: fn(f64, f64, f64) -> f64,
) {
    let s = inverse(x, beta).calculate_with::<S>().unwrap();
    assert!(s.is_finite() && s > 0.0, "x={x}, beta={beta}, s={s}");
    let repriced = price(x, s);
    let error = price_error_as_volatility_error(x, (repriced - beta).abs(), s);
    let limit = precision_limit(x, beta, s);
    assert!(error.is_finite() && limit.is_finite());
    assert!(
        error <= limit,
        "x={x:.17e}, beta={beta:.17e}, s={s:.17e}, repriced={repriced:.17e}, error={error:.3e}, limit={limit:.3e}"
    );
    assert_eq!(
        s.to_bits(),
        inverse(-x, beta).calculate_with::<S>().unwrap().to_bits()
    );
}

fn check_family_price_input(x: f64, beta: f64) {
    check_price_input::<Hybrid>(x, beta);
    check_price_input::<Jaeckel>(x, beta);
    #[cfg(feature = "flashiv")]
    check_price_with_limit::<FlashIv>(x, beta, paper_error::precision_limit);
}

#[test]
fn normalized_boundaries_are_shared_by_solver_modes() {
    for x in [0.0_f64, -1e-16, -0.5, -3.0, -1000.0] {
        let b_max = (-0.5 * x.abs()).exp();
        assert_eq!(inverse(x, 0.0).calculate::<DefaultSpecialFn>(), Some(0.0));
        assert_eq!(
            inverse(x, b_max).calculate::<DefaultSpecialFn>(),
            Some(f64::INFINITY)
        );
        assert_eq!(
            inverse(x, b_max.next_up()).calculate::<DefaultSpecialFn>(),
            None
        );
        check_family_price_input(x, b_max.next_down());
    }
}

#[test]
fn complementary_prices_respect_each_solver_accuracy_contract() {
    for x in [-1e-8_f64, -0.01, -0.5, -3.0, -190.0, -580.0, -1000.0] {
        let b_max = (0.5 * x).exp();
        for fraction in [0.5_f64.next_down(), 0.5, 0.5_f64.next_up(), 0.75, 0.99] {
            check_family_price_input(x, b_max * fraction);
        }
    }
    // Inverse vega overflows here; beta/(s*vega), and thus the bound, do not.
    for x in [-1380.0_f64, -1400.0, -1415.0] {
        check_family_price_input(x, (0.5 * x).exp().next_down());
    }
}

#[cfg(feature = "flashiv")]
#[test]
fn flashiv_inverts_normal_prices_near_the_exponent_limit() {
    // These ordinary interior prices require finite inverses even when a
    // separately evaluated reciprocal vega or derivative would overflow.
    // At this exponent the public repricer itself can exceed four estimates;
    // independently qualify the inverse against exact binary64 input roots.
    // mpmath 1.4.1 bisection at 100/200 digits agrees on every root/correction.
    // Fractions .1/.3/.7/.98 at x=-1380/-1400, and .7/.98 at x=-1415, are normal.
    for fixture in [
        (
            0xc095_9000_0000_0000,
            0x0182_9dc5_2fef_0bd3,
            0x4049_a4ef_d2d0_d02e,
            0xbce7_4829_f898_1216,
        ),
        (
            0xc095_9000_0000_0000,
            0x019b_eca7_c7e6_91bb,
            0x404a_0434_7ad8_342c,
            0xbcd4_9c9c_c7bb_27e0,
        ),
        (
            0xc095_9000_0000_0000,
            0x01b0_4a0c_89f1_2a58,
            0x404a_8a79_cfd3_f09c,
            0x3cde_b717_6d3e_d9f7,
        ),
        (
            0xc095_9000_0000_0000,
            0x01b6_ce11_8deb_3b48,
            0x404b_5311_1bb6_4c93,
            0xbcd3_5e6a_5310_4369,
        ),
        (
            0xc095_e000_0000_0000,
            0x009b_b1de_7f8e_b3ff,
            0x4049_d575_7dbb_94d3,
            0xbcca_5317_61f9_1077,
        ),
        (
            0xc095_e000_0000_0000,
            0x00b4_c566_dfab_06ff,
            0x404a_34bd_226c_95f8,
            0x3ccf_8aa8_ff95_d7e7,
        ),
        (
            0xc095_e000_0000_0000,
            0x00c8_3ba2_af9c_dd7e,
            0x404a_bb02_60ae_81ae,
            0x3ce8_debe_9daa_f11e,
        ),
        (
            0xc095_e000_0000_0000,
            0x00d0_f68b_7aed_ce3f,
            0x404b_8390_c0ea_a094,
            0x3c8a_b798_742d_de56,
        ),
        (
            0xc096_1c00_0000_0000,
            0x001b_72f6_4af7_0544,
            0x404a_df2e_b630_ffcb,
            0x3ced_c03e_488d_dfe5,
        ),
        (
            0xc096_1c00_0000_0000,
            0x0023_36df_9ae0_1d4a,
            0x404b_a7b6_85a7_194c,
            0x3ceb_6b98_54d5_1b20,
        ),
    ] {
        check_exact_root_with_limit::<FlashIv>(
            fixture.0,
            fixture.1,
            fixture.2,
            fixture.3,
            paper_error::precision_limit,
        );
    }
}

fn check_exact_root<S: BlackSolver>(
    x_bits: u64,
    beta_bits: u64,
    root_bits: u64,
    correction_bits: u64,
) {
    check_exact_root_with_limit::<S>(
        x_bits,
        beta_bits,
        root_bits,
        correction_bits,
        attainable_precision_limit,
    );
}

fn check_exact_root_with_limit<S: BlackSolver>(
    x_bits: u64,
    beta_bits: u64,
    root_bits: u64,
    correction_bits: u64,
    precision_limit: fn(f64, f64, f64) -> f64,
) {
    let x = f64::from_bits(x_bits);
    let beta = f64::from_bits(beta_bits);
    let root = f64::from_bits(root_bits);
    let correction = f64::from_bits(correction_bits);
    let actual = inverse(x, beta).calculate_with::<S>().unwrap();
    assert!(actual.is_finite() && actual > 0.0);
    assert_eq!(
        actual.to_bits(),
        inverse(-x, beta).calculate_with::<S>().unwrap().to_bits()
    );
    let error = ((actual - root) - correction).abs() / root;
    let limit = precision_limit(x, beta, root);
    assert!(limit.is_finite());
    assert!(
        error <= limit,
        "x={x:.17e}, beta={beta:.17e}, actual={actual:.17e}, root={root:.17e}, error={error:.3e}, limit={limit:.3e}"
    );
}

#[test]
fn upper_tail_matches_independent_exact_input_roots() {
    // The normalized Black formula from source/lets_be_rational.cpp, evaluated
    // independently using mpmath 1.4.1 and exact integer-ratio binary64 inputs.
    // Bisection at 100/200 digits agrees on root and sub-ulp correction bits.
    for fixture in [
        (
            0xbfe0_0000_0000_0000,
            0x3fe8_ac22_fbce_9f76,
            0x4015_477a_ca58_0b9a,
            0x3cb5_1276_903c_0a0f,
        ),
        (
            0xc008_0000_0000_0000,
            0x3fcc_4669_ee97_67fe,
            0x4018_3587_d488_7bbb,
            0x3cb0_d6f5_e37f_e21a,
        ),
        (
            0xc095_9000_0000_0000,
            0x01b7_4536_7bea_cec6,
            0x404e_b820_2842_f2d7,
            0x3cd9_8466_c8ed_3156,
        ),
        (
            0xc095_e000_0000_0000,
            0x00d1_4f2b_0fb9_307e,
            0x404e_d7de_7d42_e981,
            0xbce5_6f8b_2616_1efa,
        ),
    ] {
        check_exact_root::<Hybrid>(fixture.0, fixture.1, fixture.2, fixture.3);
        check_exact_root::<Jaeckel>(fixture.0, fixture.1, fixture.2, fixture.3);
        #[cfg(feature = "flashiv")]
        check_exact_root_with_limit::<FlashIv>(
            fixture.0,
            fixture.1,
            fixture.2,
            fixture.3,
            paper_error::precision_limit,
        );
    }
}

#[test]
fn solvers_recover_normal_roots_at_microscopic_moneyness() {
    // All beta/root values here are normal. Exact Black roots agree at two
    // independent precisions. The dimensionless Jaeckel correction must avoid
    // overflow even when its unscaled second derivative is unrepresentable;
    // x=1e-300 needs 400/500 digits to resolve the formula's cancellation.
    for fixture in [
        (
            0xbe7a_d7f2_9abc_af48,
            0x3e50_3a53_0d6b_3fb8,
            0x3e80_d154_947f_a38b,
            0xbb24_4b4a_60d7_df6f,
        ),
        (
            0xb949_f623_d5a8_a733,
            0x38cc_3717_4909_39c0,
            0x3939_f623_d5a8_a733,
            0x35c7_bce1_78fe_d4c0,
        ),
        (
            0x9668_7e92_154e_f7ac,
            0x15ea_9eeb_20a4_3d7a,
            0x1658_7e92_154e_f7ac,
            0x12e5_8fde_d2c0_f043,
        ),
        (
            0x81a5_6e1f_c2f8_f359,
            0x0127_4a5f_c7cb_2062,
            0x0195_6e1f_c2f8_f359,
            0x0000_0000_002d_a5b9,
        ),
    ] {
        check_exact_root::<Hybrid>(fixture.0, fixture.1, fixture.2, fixture.3);
        check_exact_root::<Jaeckel>(fixture.0, fixture.1, fixture.2, fixture.3);
        #[cfg(feature = "flashiv")]
        check_exact_root_with_limit::<FlashIv>(
            fixture.0,
            fixture.1,
            fixture.2,
            fixture.3,
            paper_error::precision_limit,
        );
    }
}

#[test]
fn full_otm_call_and_put_use_the_selected_solver_across_price_regions() {
    for (abs_x, s) in [
        (1e-7_f64, 0.001_f64),
        (0.5, 0.125),
        (0.5, 0.7),
        (0.5, 1.6),
        (3.0, 6.0),
    ] {
        for scale in [1e-100_f64, 1.0, 1e100] {
            for is_call in [true, false] {
                let smaller = scale * (-abs_x).exp();
                let (f, k) = if is_call {
                    (smaller, scale)
                } else {
                    (scale, smaller)
                };
                let expiry: f64 = 0.25;
                let sigma = s / expiry.sqrt();
                let observed = PriceBlackScholes::builder()
                    .forward(f)
                    .strike(k)
                    .expiry(expiry)
                    .is_call(is_call)
                    .volatility(sigma)
                    .build()
                    .unwrap()
                    .calculate::<DefaultSpecialFn>();
                let input = ImpliedBlackVolatility::builder()
                    .forward(f)
                    .strike(k)
                    .expiry(expiry)
                    .is_call(is_call)
                    .option_price(observed)
                    .build()
                    .unwrap();
                let actuals = [
                    (
                        input.calculate_with::<Hybrid>().unwrap(),
                        attainable_precision_limit as fn(f64, f64, f64) -> f64,
                    ),
                    (
                        input.calculate_with::<Jaeckel>().unwrap(),
                        attainable_precision_limit,
                    ),
                    #[cfg(feature = "flashiv")]
                    (
                        input.calculate_with::<FlashIv>().unwrap(),
                        paper_error::precision_limit,
                    ),
                ];
                for (actual, precision_limit) in actuals {
                    let beta = observed / (f.sqrt() * k.sqrt());
                    let error = (actual - sigma).abs() / sigma;
                    let limit = precision_limit(-abs_x, beta, s);
                    assert!(error.is_finite() && limit.is_finite());
                    assert!(
                        error <= limit,
                        "F={f}, K={k}, call={is_call}, sigma={sigma}, actual={actual}, error={error:.3e}, limit={limit:.3e}"
                    );
                }
            }
        }
    }
}

#[cfg(feature = "flashiv")]
struct NoJaeckelInterpolation;

#[cfg(feature = "flashiv")]
impl SpecialFn for NoJaeckelInterpolation {
    fn one_minus_erfcx(_: f64) -> f64 {
        panic!("pure FlashIV must not enter the Jäckel interpolation seed")
    }
}

#[cfg(feature = "flashiv")]
#[test]
fn pure_flashiv_handles_all_price_regions_without_jaeckel_interpolation() {
    for (x, s) in [
        (-1e-16, 0.001),
        (-0.5, 0.125),
        (-0.5, 0.7),
        (-0.5, 1.6),
        (-0.5, 5.0),
    ] {
        let beta = price(x, s);
        let actual = inverse(x, beta)
            .calculate_with::<FlashIv<NoJaeckelInterpolation>>()
            .unwrap();
        let limit = paper_error::precision_limit(x, beta, s);
        assert!((actual - s).abs() / s <= limit);
    }
}

struct ObserveJaeckelInterpolation;

static JAECKEL_SEEDS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

impl SpecialFn for ObserveJaeckelInterpolation {
    fn one_minus_erfcx(x: f64) -> f64 {
        JAECKEL_SEEDS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        DefaultSpecialFn::one_minus_erfcx(x)
    }
}

#[test]
fn pure_jaeckel_retains_its_seed_in_every_non_atm_price_region() {
    for (x, s) in [(-0.5, 0.125), (-0.5, 0.7), (-0.5, 1.6), (-0.5, 5.0)] {
        let beta = price(x, s);
        JAECKEL_SEEDS.store(0, std::sync::atomic::Ordering::Relaxed);
        let actual = inverse(x, beta)
            .calculate_with::<Jaeckel<ObserveJaeckelInterpolation>>()
            .unwrap();
        assert!(JAECKEL_SEEDS.load(std::sync::atomic::Ordering::Relaxed) > 0);
        let limit = 4.0 * conditioning(x, beta, s);
        assert!((actual - s).abs() / s <= limit);
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum HybridProviderCall {
    Erfcx(u64),
    JaeckelSeed(u64),
}

#[derive(Default)]
struct HybridProviderTrace {
    fail_erfcx_at: Option<u64>,
    failures: usize,
    calls: Vec<HybridProviderCall>,
}

thread_local! {
    static HYBRID_PROVIDER_TRACE: std::cell::RefCell<HybridProviderTrace> =
        const { std::cell::RefCell::new(HybridProviderTrace {
            fail_erfcx_at: None,
            failures: 0,
            calls: Vec::new(),
        }) };
}

struct TraceHybridProvider;

impl SpecialFn for TraceHybridProvider {
    fn erfcx(x: f64) -> f64 {
        let fail = HYBRID_PROVIDER_TRACE.with(|trace| {
            let mut trace = trace.borrow_mut();
            trace.calls.push(HybridProviderCall::Erfcx(x.to_bits()));
            let fail = trace.fail_erfcx_at == Some(x.to_bits());
            trace.failures += usize::from(fail);
            fail
        });
        if fail {
            f64::NAN
        } else {
            DefaultSpecialFn::erfcx(x)
        }
    }

    fn one_minus_erfcx(x: f64) -> f64 {
        HYBRID_PROVIDER_TRACE.with(|trace| {
            trace
                .borrow_mut()
                .calls
                .push(HybridProviderCall::JaeckelSeed(x.to_bits()));
        });
        DefaultSpecialFn::one_minus_erfcx(x)
    }
}

fn trace_inverse<S: BlackSolver>(
    x: f64,
    beta: f64,
    fail_erfcx_at: Option<u64>,
) -> (f64, HybridProviderTrace) {
    HYBRID_PROVIDER_TRACE.with(|trace| {
        *trace.borrow_mut() = HybridProviderTrace {
            fail_erfcx_at,
            ..HybridProviderTrace::default()
        };
    });
    let actual = inverse(x, beta).calculate_with::<S>().unwrap();
    let trace = HYBRID_PROVIDER_TRACE.with(|trace| std::mem::take(&mut *trace.borrow_mut()));
    (actual, trace)
}

#[test]
fn hybrid_skips_the_jaeckel_seed_only_inside_its_predispatch_guard() {
    for (x, fraction, seed_calls) in [
        (-20.0_f64, 0.0004, 0),
        (-20.0, 0.001, 1),
        ((-0.01_f64).next_up(), 0.0004, 1),
        (-1420.0, 0.0004, 1),
    ] {
        let beta = fraction * (0.5 * x).exp();
        let (hybrid, hybrid_trace) = trace_inverse::<Hybrid<TraceHybridProvider>>(x, beta, None);
        let (jaeckel, jaeckel_trace) = trace_inverse::<Jaeckel<TraceHybridProvider>>(x, beta, None);
        let seeds = |trace: &HybridProviderTrace| {
            trace
                .calls
                .iter()
                .filter(|call| matches!(call, HybridProviderCall::JaeckelSeed(_)))
                .count()
        };
        assert_eq!(seeds(&hybrid_trace), seed_calls, "x={x}, beta={beta}");
        assert_eq!(seeds(&jaeckel_trace), 1, "x={x}, beta={beta}");
        assert_eq!(
            hybrid.to_bits(),
            inverse(x, beta)
                .calculate_with::<Hybrid>()
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            jaeckel.to_bits(),
            inverse(x, beta)
                .calculate_with::<Jaeckel>()
                .unwrap()
                .to_bits()
        );
    }
}

#[test]
fn hybrid_failed_predispatch_restores_the_jaeckel_path_without_retry() {
    let x = -20.0_f64;
    let beta = 0.0004 * (0.5 * x).exp();
    // This deep case uses the provider in the first exact FlashIV step. Record
    // its argument under the current FMA policy, then fail only that argument;
    // all remaining function values are deterministic default Cody values.
    let (_, qualified) = trace_inverse::<Hybrid<TraceHybridProvider>>(x, beta, None);
    assert!(
        qualified
            .calls
            .iter()
            .all(|call| matches!(call, HybridProviderCall::Erfcx(_)))
    );
    let Some(HybridProviderCall::Erfcx(fail_at)) = qualified.calls.first() else {
        panic!("the fixture must reach the exact FlashIV erfcx objective");
    };
    let fail_at = *fail_at;
    let (expected, jaeckel) = trace_inverse::<Jaeckel<TraceHybridProvider>>(x, beta, Some(fail_at));
    assert_eq!(
        jaeckel.failures, 0,
        "the failed argument must be exclusive to FlashIV"
    );
    let (actual, hybrid) = trace_inverse::<Hybrid<TraceHybridProvider>>(x, beta, Some(fail_at));
    assert_eq!(hybrid.failures, 1);
    assert_eq!(actual.to_bits(), expected.to_bits());
    assert!(
        matches!(hybrid.calls.first(), Some(HybridProviderCall::Erfcx(bits)) if *bits == fail_at)
    );
    assert!(matches!(
        hybrid.calls.get(1),
        Some(HybridProviderCall::Erfcx(_))
    ));
    // One failed exact step evaluates its erfcx pair. Every later provider
    // call must match the original Jaeckel path, including its seed call.
    assert_eq!(&hybrid.calls[2..], jaeckel.calls.as_slice());
    assert!(matches!(
        jaeckel.calls.first(),
        Some(HybridProviderCall::JaeckelSeed(_))
    ));
}
