//! Exact-bit regression tests for the power-of-two/scalar DD specializations.
//!
//! The test reference freezes the old common finish expression, so it also
//! guards the original FMA operation order independently of result tolerance.
use super::{
    Pair, add, finish, fma, forward_D, mul, scale as new_scale, scale_negative_half, square_scalar,
};
use crate::{DefaultSpecialFn, PriceBlackScholesNormalised};

#[inline]
fn old_scale(a: Pair, b: f64) -> Pair {
    mul(a, Pair::from(b))
}

#[inline]
fn frozen_finish(a: f64, b: f64, s: f64) -> f64 {
    let hv = a / s;
    let h = Pair::new(hv, fma(-hv, s, a) / s);
    let t = Pair::from(s * 0.5);
    if !(h.hi >= 0.0 && h.hi <= 8.0 && t.hi > 0.0 && t.hi < 5.5 && h.hi + t.hi < 15.5) {
        return f64::NAN;
    }
    let hh = mul(h, h);
    let tt = mul(t, t);
    let d = forward_D(h, t);
    if !(d.hi > 0.0 && d.hi.is_finite()) {
        return f64::NAN;
    }
    let exponent = old_scale(add(hh, tt), -0.5);
    let ev = exponent.hi.exp();
    let v = mul(
        mul(Pair::from(ev), Pair::new(1.0, exponent.lo)),
        Pair::new(
            f64::from_bits(0x3fd9884533d43651),
            -f64::from_bits(0x3c7cbc0d30ebfd15),
        ),
    );
    let sv = old_scale(v, s);
    let residual = (fma(-sv.hi, d.hi, b) - sv.hi * d.lo) - sv.lo * d.hi;
    let n = residual / sv.hi;
    let h2 = hh.hi - tt.hi;
    let q = 3.0 * hh.hi + tt.hi;
    let c2 = -0.5 * h2;
    let c3 = (2.0 * (h2 * h2) + q) / 6.0;
    let ds = n * fma(n, fma(n, c3, c2), 1.0);
    fma(s, ds, s)
}

#[test]
fn scalar_dd_specialization_preserves_bits() {
    let mut state = 0x43d2_c7b9_8293_54a1_u64;
    let mut covered = 0;
    for i in 0..64_000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let frac = f64::from_bits((0x3ff0_0000_0000_0000_u64) | (state & 0x000f_ffff_ffff_ffff));
        let exponent = (i % 980) as i32 - 480;
        let x = frac * 2.0_f64.powi(exponent);
        let low = fma(x, 0.7, -(x * 0.7));
        let pair = add(Pair::from(x), Pair::from(low));
        for k in [-0.5, 2.0, 0.125, 0.6] {
            let a = old_scale(pair, k);
            let b = new_scale(pair, k);
            assert_eq!(a.hi.to_bits(), b.hi.to_bits(), "scale hi {i}");
            assert_eq!(a.lo.to_bits(), b.lo.to_bits(), "scale lo {i}");
        }
        let old_square = mul(Pair::from(x), Pair::from(x));
        let new_square = square_scalar(x);
        assert_eq!(
            old_square.hi.to_bits(),
            new_square.hi.to_bits(),
            "square hi {i}"
        );
        assert_eq!(
            old_square.lo.to_bits(),
            new_square.lo.to_bits(),
            "square lo {i}"
        );
        if pair.hi > f64::from_bits(0x20b0_0000_0000_0000) {
            let old = old_scale(pair, -0.5);
            let new = scale_negative_half(pair);
            assert_eq!(old.hi.to_bits(), new.hi.to_bits(), "half hi {i}");
            assert_eq!(old.lo.to_bits(), new.lo.to_bits(), "half lo {i}");
            covered += 1;
        }
    }
    assert!(covered > 20_000);
}

#[test]
fn common_finish_matches_frozen_expression() {
    let mut compared = 0_usize;
    for i in 0..2048_u32 {
        let s = 0.1 + f64::from((i * 73) % 251) / 60.0;
        let a = s * (0.05 + f64::from((i * 37) % 397) / 110.0);
        let price = PriceBlackScholesNormalised::builder()
            .log_moneyness(-a)
            .total_volatility(s)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        if !(price > 0.0 && price < (-0.5 * a).exp()) {
            continue;
        }
        let perturb = f64::from((i * 17) % 23) - 11.0;
        let seed = s * (1.0 + perturb * 1e-7);
        let old = frozen_finish(a, price, seed);
        let new = finish(a, price, seed);
        if old.is_finite() && new.is_finite() {
            assert_eq!(
                old.to_bits(),
                new.to_bits(),
                "a={a:.17e}, b={price:.17e}, seed={seed:.17e}"
            );
            compared += 1;
        }
    }
    assert!(compared > 1200, "insufficient accepted finish cases");
}

#[test]
fn exact_exponential_pair_matches_generic_multiplication() {
    for i in 0..20_000_u32 {
        let exponent = -50.0 * f64::from((i * 137) % 20_003) / 20_003.0;
        let high = exponent.exp();
        let low = (f64::from((i * 271) % 19_999) / 19_999.0 - 0.5) * 1e-15;
        let original = mul(Pair::from(high), Pair::new(1.0, low));
        let optimized = add(Pair::from(high), Pair::from(high * low));
        assert_eq!(original.hi.to_bits(), optimized.hi.to_bits(), "high {i}");
        assert_eq!(original.lo.to_bits(), optimized.lo.to_bits(), "low {i}");
    }
}

#[test]
fn finish_dispatch_boundaries_preserve_output_bits() {
    // Include both sides of t=0.001, t=0.5, t=1 and h=2, the
    // boundaries at which the existing D evaluator switches formulas.
    let hs = [0.0001, 0.01, 0.04, 0.1, 0.5, 1.999_999, 2.0, 2.000_001, 7.9];
    let ts = [
        0.000_1, 0.000_499, 0.000_5, 0.000_501, 0.249_999, 0.25, 0.250_001, 0.499_999, 0.5,
        0.500_001, 1.0, 2.5, 5.49,
    ];
    let mut checked = 0_usize;
    for h in hs {
        for t in ts {
            let s = 2.0 * t;
            let a = h * s;
            if a > 10.0 {
                continue;
            }
            let price = PriceBlackScholesNormalised::builder()
                .log_moneyness(-a)
                .total_volatility(s)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>();
            if !(price >= f64::MIN_POSITIVE && price < (-0.5 * a).exp()) {
                continue;
            }
            for tweak in [-1e-5_f64, 0.0, 1e-5] {
                let seed = s * (1.0 + tweak);
                let old = frozen_finish(a, price, seed);
                let new = finish(a, price, seed);
                if old.is_finite() {
                    assert_eq!(
                        old.to_bits(),
                        new.to_bits(),
                        "h={h:.17e} t={t:.17e} tweak={tweak}"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(
        checked > 65,
        "insufficient representable boundary fixtures: {checked}"
    );
}

#[test]
fn scalar_square_subnormal_paths_match_original() {
    for exponent in [-530, -520, -512, -511, -510, -500, -1, 0, 1] {
        let x = 2.0_f64.powi(exponent);
        for scale in [1.0, 1.125, 1.5] {
            let y = x * scale;
            let old = mul(Pair::from(y), Pair::from(y));
            let new = square_scalar(y);
            assert_eq!(old.hi.to_bits(), new.hi.to_bits());
            assert_eq!(old.lo.to_bits(), new.lo.to_bits());
        }
    }
}
