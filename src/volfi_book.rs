//! Rust adaptation of volfi v0.3.1's two whole-book direct charts.
//!
//! Copyright (c) 2026 Wolfgang Schadner, BSD 3-Clause.
//! Ported from wol-fi/volfi commit 7347f68025fea2463b30f67c714a7448a64d1c39.
//! See `src/volfi_book/LICENSE-BSD-3-CLAUSE` and `docs/volfi-book-study.md`.
//!
//! Unlike upstream volfi, uncovered inputs fall back to the pre-existing Rust
//! Experimental solver; the original four-chart fallback is **not** ported.
//! Rust standard-library expm1/log1p replace upstream custom functions, so
//! bit-identical output versus the C++ implementation is *not* asserted.
//!
//! Inputs: h=abs(log(F/K)), b=OTM price/sqrt(F*K).
//! Output: total implied volatility (the upstream book chart returns variance).
#![allow(non_upper_case_globals, clippy::doc_markdown)]
mod tables;
use tables::*;

/// Upstream whole-book guard, 2^-500.
const H_FLOOR: f64 = f64::from_bits(523_u64 << 52);
const TWO_PI: f64 = 6.283_185_307_179_586;
const A_MAX: f64 = 700.0;
const INV_4PI2: f64 = 2.533_029_591_058_444_3e-2;

#[inline]
fn clenshaw(coef: &[f64], x: f64) -> f64 {
    let x2 = 2.0 * x;
    let mut e0 = 0.0;
    let mut e1 = 0.0;
    for &c in coef[1..].iter().rev() {
        let next = x2.mul_add(e0, c) - e1;
        e1 = e0;
        e0 = next;
    }
    x.mul_add(e0, coef[0]) - e1
}

#[inline]
fn poly_a(row: usize, y: f64) -> f64 {
    let offset = NWRowA_OFF[row];
    let degree = NWRowA_DEG[row];
    let mut p = NW_P[offset];
    for j in 1..=degree {
        p = p.mul_add(y, NW_P[offset + j]);
    }
    p
}

#[inline]
fn poly_b(row: usize, w: f64) -> f64 {
    let offset = NWRowB_OFF[row];
    let degree = NWRowB_DEG[row];
    let mut p = NW_D[offset];
    for j in 1..degree {
        p = p.mul_add(w, NW_D[offset + j]);
    }
    p * w
}

#[inline]
fn direct_a(h: f64, a: f64) -> f64 {
    let g = clenshaw(&NW_G0, a.mul_add(0.318_309_886_183_790_67, -1.0));
    let a0 = a * g;
    let t = h / a0;
    let x = t * t;
    let y = a0 * a0;
    let mut sum = 0.0_f64;
    for row in (1..=15).rev() {
        sum = sum.mul_add(x, poly_a(row, y));
    }
    (t * x).mul_add(sum, t)
}

#[inline]
fn direct_b(h: f64, a: f64) -> f64 {
    let u = a.sqrt().recip();
    let f1 = clenshaw(
        &NW_F1,
        u.mul_add(5.537_929_048_854_298, -1.209_314_043_451_270_6),
    );
    let pp = u * f1;
    let t = h * pp;
    let w = pp * pp;
    let z = h * h * INV_4PI2;
    let sp = (1.0 + z).sqrt() + 1.0;
    let q = z / (sp * sp);
    let mut sum = 0.0_f64;
    for row in (1..=16).rev() {
        sum = sum.mul_add(q, poly_b(row, w));
    }
    (t * q).mul_add(sum, t)
}

/// Direct whole-book chart without fallback; returns `None` outside its region.
///
/// The primary accuracy contract for Experimental (`rho_J < 1`) has not been
/// proved for this experimental port. Finite tests are not certification.
#[must_use]
#[inline]
pub fn chart_total_volatility(h: f64, b: f64) -> Option<f64> {
    if !(h >= H_FLOOR && h <= 16.5 && b > 0.0 && b.is_finite()) {
        return None;
    }
    let c = b * (0.5 * h).exp();
    if !(c > 0.0 && c < 1.0) {
        return None;
    }
    // a is volfi's intrinsic coordinate, not log-moneyness h.
    let ratio = h.exp_m1() / c;
    if !ratio.is_finite() {
        return None;
    }
    let a = ratio.ln_1p();
    let v = if a < TWO_PI {
        if h > 0.35 * a {
            return None;
        }
        direct_a(h, a)
    } else {
        if !(a <= A_MAX && h <= 4.0 && h <= a.mul_add(0.16, 0.85)) {
            return None;
        }
        direct_b(h, a)
    };
    (v > 0.0 && v.is_finite()).then_some(v)
}

/// Whole-book direct chart, with Rust Experimental fallback for uncovered cases.
///
/// This covers Experimental's supported input domain while preserving its
/// cap/zero-price handling in the `solver::VolfiBook` adapter.
#[must_use]
#[inline]
pub fn implied_total_volatility(h: f64, b: f64) -> f64 {
    chart_total_volatility(h, b)
        .unwrap_or_else(|| crate::experimental::implied_total_volatility(h, b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_direct_regions_and_fallback_are_reachable() {
        for (h, a) in [(0.2_f64, 1.4_f64), (1.0, 9.0), (2.0, 20.0)] {
            let c = h.exp_m1() / a.exp_m1();
            let b = c * (-h * 0.5).exp();
            let direct = chart_total_volatility(h, b);
            assert!(direct.is_some(), "h={h} a={a}");
            let reference = crate::experimental::implied_total_volatility(h, b);
            let got = direct.unwrap();
            assert!(reference > 0.0 && reference.is_finite());
            assert!(
                (got / reference - 1.0).abs() < 2e-13,
                "h={h} a={a} book={got} experimental={reference}"
            );
        }
        let h = 5.0_f64;
        let b = 0.05 * (-0.5 * h).exp();
        assert!(chart_total_volatility(h, b).is_none());
        assert_eq!(
            implied_total_volatility(h, b).to_bits(),
            crate::experimental::implied_total_volatility(h, b).to_bits()
        );
        assert!(chart_total_volatility(0.0, 0.1).is_none());
    }
}
