//! Exact-bit regression tests for the power-of-two/scalar DD specializations.
//!
//! The test reference freezes the old common finish expression, so it also
//! guards the original FMA operation order independently of result tolerance.
use super::{
    Pair, add, finish, fma, forward_D, mul, scale as new_scale, scale_negative_half, square_scalar,
};
use crate::{DefaultSpecialFn, PriceBlackScholesNormalised};
use rand::{RngExt, SeedableRng};

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
fn profile_existing_finish_residual_distribution() {
    let mut rng = rand::rngs::StdRng::from_seed([42; 32]);
    let mut complete = 0_usize;
    let mut accepted = 0_usize;
    super::FINISH_N_HIST.with(|h| *h.borrow_mut() = [0; 8]);
    while complete < 4096 {
        let a = rng.random_range(1e-4_f64..5.0_f64);
        let total_vol = rng.random_range(1e-3_f64..5.0_f64);
        let b = PriceBlackScholesNormalised::builder()
            .log_moneyness(-a)
            .total_volatility(total_vol)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        if b < f64::MIN_POSITIVE || b >= (-0.5 * a).exp() {
            continue;
        }
        let vol = super::solve(a, b);
        if vol.is_finite() && vol > 0.0 {
            accepted += 1;
        }
        complete += 1;
    }
    let hist = super::FINISH_N_HIST.with(|h| *h.borrow());
    let total: usize = hist.iter().sum();
    assert!(accepted > 3000 && total > 250);
    eprintln!("FINISH_N_PROFILE|input_count={complete}|valid={accepted}|finish_count={total}|hist={hist:?}");
}
