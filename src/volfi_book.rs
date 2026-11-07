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
#![allow(
    non_upper_case_globals,
    clippy::doc_markdown,
    clippy::approx_constant,
    clippy::wildcard_imports,
    clippy::manual_range_contains
)]
mod tables;
use tables::*;

/// Upstream whole-book guard, 2^-500.
const H_FLOOR: f64 = f64::from_bits(523_u64 << 52);
const TWO_PI: f64 = 6.283_185_307_179_586;
const A_MAX: f64 = 700.0;
const INV_4PI2: f64 = 2.533_029_591_058_444_3e-2;

// Keep all coefficient indices and row degrees compile-time constants as in
// volfi's C++ templates. The former runtime-indexed loop prevented LLVM from
// seeing the straight-line structure of the original whole-book chart.
#[inline(always)]
fn clenshaw<const N: usize>(coef: &[f64; N], x: f64) -> f64 {
    let x2 = 2.0 * x;
    let mut e0 = 0.0;
    let mut e1 = 0.0;
    for i in (1..N).rev() {
        let next = x2.mul_add(e0, coef[i]) - e1;
        e1 = e0;
        e0 = next;
    }
    x.mul_add(e0, coef[0]) - e1
}

#[inline(always)]
fn poly_a<const ROW: usize>(y: f64) -> f64 {
    let offset = NWRowA_OFF[ROW];
    let degree = NWRowA_DEG[ROW];
    let mut p = NW_P[offset];
    for j in 1..=degree {
        p = p.mul_add(y, NW_P[offset + j]);
    }
    p
}

#[inline(always)]
fn poly_b<const ROW: usize>(w: f64) -> f64 {
    let offset = NWRowB_OFF[ROW];
    let degree = NWRowB_DEG[ROW];
    let mut p = NW_D[offset];
    for j in 1..degree {
        p = p.mul_add(w, NW_D[offset + j]);
    }
    p * w
}

mod compensated;

// Unroll the outer dimension (each row is a separate monomorphization).
// Within each row the Horner FMA sequence and scalar accumulator are identical
// to the original Rust port; only integer indexing/control are removed.
macro_rules! rows_a {
    ($x:expr, $y:expr) => {{
        let mut sum = 0.0_f64;
        sum = sum.mul_add($x, poly_a::<15>($y));
        sum = sum.mul_add($x, poly_a::<14>($y));
        sum = sum.mul_add($x, poly_a::<13>($y));
        sum = sum.mul_add($x, poly_a::<12>($y));
        sum = sum.mul_add($x, poly_a::<11>($y));
        sum = sum.mul_add($x, poly_a::<10>($y));
        sum = sum.mul_add($x, poly_a::<9>($y));
        sum = sum.mul_add($x, poly_a::<8>($y));
        sum = sum.mul_add($x, poly_a::<7>($y));
        sum = sum.mul_add($x, poly_a::<6>($y));
        sum = sum.mul_add($x, poly_a::<5>($y));
        sum = sum.mul_add($x, poly_a::<4>($y));
        sum = sum.mul_add($x, poly_a::<3>($y));
        sum = sum.mul_add($x, poly_a::<2>($y));
        sum = sum.mul_add($x, poly_a::<1>($y));
        sum
    }};
}

macro_rules! rows_b {
    ($q:expr, $w:expr) => {{
        let mut sum = 0.0_f64;
        sum = sum.mul_add($q, poly_b::<16>($w));
        sum = sum.mul_add($q, poly_b::<15>($w));
        sum = sum.mul_add($q, poly_b::<14>($w));
        sum = sum.mul_add($q, poly_b::<13>($w));
        sum = sum.mul_add($q, poly_b::<12>($w));
        sum = sum.mul_add($q, poly_b::<11>($w));
        sum = sum.mul_add($q, poly_b::<10>($w));
        sum = sum.mul_add($q, poly_b::<9>($w));
        sum = sum.mul_add($q, poly_b::<8>($w));
        sum = sum.mul_add($q, poly_b::<7>($w));
        sum = sum.mul_add($q, poly_b::<6>($w));
        sum = sum.mul_add($q, poly_b::<5>($w));
        sum = sum.mul_add($q, poly_b::<4>($w));
        sum = sum.mul_add($q, poly_b::<3>($w));
        sum = sum.mul_add($q, poly_b::<2>($w));
        sum = sum.mul_add($q, poly_b::<1>($w));
        sum
    }};
}

#[inline]
fn direct_a<const COMPENSATED: bool, const SELECTIVE: bool>(h: f64, a: f64) -> f64 {
    let (a0, t) = if SELECTIVE {
        compensated::region_a_scale_selective(h, a)
    } else if COMPENSATED {
        compensated::region_a_scale(h, a)
    } else {
        let g = clenshaw(&NW_G0, a.mul_add(0.318_309_886_183_790_7, -1.0));
        let a0 = a * g;
        (a0, h / a0)
    };
    let x = t * t;
    let y = a0 * a0;
    let sum = rows_a!(x, y);
    (t * x).mul_add(sum, t)
}

#[inline]
#[allow(clippy::suboptimal_flops)] // Retain the analysed unfused correction order.
fn direct_b<const REFINED_RSQRT: bool, const COMPENSATED: bool, const SELECTIVE: bool>(
    h: f64,
    a: f64,
) -> f64 {
    let root = a.sqrt();
    let mut u = root.recip();
    if REFINED_RSQRT {
        // Remove double rounding in sqrt(a) followed by reciprocal:
        //   a - root^2 and 1 - root*u are FMA residuals.
        // One inverse-square-root Newton correction approximates
        //   u*sqrt(a) = 1 in the exact input a, not only rounded root.
        // This is a finite-precision optimization candidate, not a
        // globally certified correctly-rounded reciprocal square root.
        let square_residual = (-root).mul_add(root, a);
        let reciprocal_residual = (-root).mul_add(u, 1.0);
        let correction = reciprocal_residual - (0.5 * square_residual) * (u * u);
        u = u.mul_add(correction, u);
    }
    let pp = if SELECTIVE {
        compensated::region_b_scale_selective(a)
    } else if COMPENSATED {
        compensated::region_b_scale(a)
    } else {
        let f1 = clenshaw(
            &NW_F1,
            u.mul_add(5.537_929_048_854_298, -1.209_314_043_451_270_6),
        );
        u * f1
    };
    let t = h * pp;
    let w = pp * pp;
    let z = h * h * INV_4PI2;
    let sp = (1.0 + z).sqrt() + 1.0;
    let q = z / (sp * sp);
    let sum = rows_b!(q, w);
    (t * q).mul_add(sum, t)
}

/// Direct whole-book chart without fallback; returns `None` outside its region.
///
/// The primary accuracy contract for Experimental (`rho_J < 1`) has not been
/// proved for this experimental port. Finite tests are not certification.
#[must_use]
#[inline]
pub fn chart_total_volatility(h: f64, b: f64) -> Option<f64> {
    chart_impl::<false, false, false, false>(h, b)
}

/// Alternative algebraically equivalent coordinate:
/// a = log1p(2*sinh(h/2)/b). This avoids both exp(h/2) and expm1(h).
///
/// On the stated chart the region inequalities imply a > h, which is exactly
/// the strict Black cap condition b < exp(-h/2) in real arithmetic.
#[must_use]
#[inline]
pub fn chart_total_volatility_sinh(h: f64, b: f64) -> Option<f64> {
    chart_impl::<true, false, false, false>(h, b)
}

/// Research chart variant with a single FMA-based inverse-sqrt correction
/// in region B. This branch may still exceed rho_J = 1.
#[must_use]
#[inline]
pub fn chart_total_volatility_rsqrt(h: f64, b: f64) -> Option<f64> {
    chart_impl::<false, true, false, false>(h, b)
}

/// High-precision chart-coordinate prototype; evaluated only when opted in.
/// Validated on existing finite reference fixtures, not globally certified.
#[must_use]
#[inline]
pub fn chart_total_volatility_compensated(h: f64, b: f64) -> Option<f64> {
    chart_impl::<false, true, true, false>(h, b)
}

/// Localized twofold Clenshaw research variant, retaining the same A/B
/// approximate coefficient tables and the original Experimental fallback.
#[must_use]
#[inline]
pub fn chart_total_volatility_selective(h: f64, b: f64) -> Option<f64> {
    chart_impl::<false, true, false, true>(h, b)
}

#[inline]
fn chart_impl<
    const SINH: bool,
    const REFINED_RSQRT: bool,
    const COMPENSATED: bool,
    const SELECTIVE: bool,
>(
    h: f64,
    b: f64,
) -> Option<f64> {
    // Both shipped direct charts imply h <= 4: region A has h <= 0.35*2pi
    // and region B explicitly imposes h <= 4.
    if !(h >= H_FLOOR && h <= 4.0 && b > 0.0 && b.is_finite()) {
        return None;
    }
    // a is volfi's intrinsic coordinate, not log-moneyness h.
    let ratio = if SINH {
        2.0 * (0.5 * h).sinh() / b
    } else {
        let c = b * (0.5 * h).exp();
        if !(c > 0.0 && c < 1.0) {
            return None;
        }
        h.exp_m1() / c
    };
    if !ratio.is_finite() {
        return None;
    }
    let a = ratio.ln_1p();
    let v = if a < TWO_PI {
        if h > 0.35 * a {
            return None;
        }
        direct_a::<COMPENSATED, SELECTIVE>(h, a)
    } else {
        if !(a <= A_MAX && h <= 4.0 && h <= a.mul_add(0.16, 0.85)) {
            return None;
        }
        direct_b::<REFINED_RSQRT, COMPENSATED, SELECTIVE>(h, a)
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

/// Research variant using the hyperbolic-sine intrinsic coordinate.
#[must_use]
#[inline]
pub fn implied_total_volatility_sinh(h: f64, b: f64) -> f64 {
    chart_total_volatility_sinh(h, b)
        .unwrap_or_else(|| crate::experimental::implied_total_volatility(h, b))
}

/// Same direct region as VolfiBook, with compensated inverse square root in B.
#[must_use]
#[inline]
pub fn implied_total_volatility_rsqrt(h: f64, b: f64) -> f64 {
    chart_total_volatility_rsqrt(h, b)
        .unwrap_or_else(|| crate::experimental::implied_total_volatility(h, b))
}

/// Compensated research route with Rust Experimental fallback.
#[must_use]
#[inline]
pub fn implied_total_volatility_compensated(h: f64, b: f64) -> f64 {
    chart_total_volatility_compensated(h, b)
        .unwrap_or_else(|| crate::experimental::implied_total_volatility(h, b))
}

/// Localized twofold chart with Experimental fallback.
#[must_use]
#[inline]
pub fn implied_total_volatility_selective(h: f64, b: f64) -> f64 {
    chart_total_volatility_selective(h, b)
        .unwrap_or_else(|| crate::experimental::implied_total_volatility(h, b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn specialized_rows_preserve_dynamic_horner_bits() {
        // The original dynamic evaluator is frozen below. Compile-time row
        // specialization must not change a single live FMA rounding point.
        for x in [0.0_f64, f64::EPSILON, 1e-8, 0.005, 0.02, 0.1, 0.4] {
            for y in [0.0_f64, 1e-12, 0.005, 0.1, 0.8, 2.0, 6.0, 14.0] {
                let mut old_a = 0.0_f64;
                for row in (1..=15).rev() {
                    let offset = NWRowA_OFF[row];
                    let degree = NWRowA_DEG[row];
                    let mut p = NW_P[offset];
                    for j in 1..=degree {
                        p = p.mul_add(y, NW_P[offset + j]);
                    }
                    old_a = old_a.mul_add(x, p);
                }
                let new_a = rows_a!(x, y);
                assert_eq!(new_a.to_bits(), old_a.to_bits(), "A x={x}, y={y}");

                let mut old_b = 0.0_f64;
                for row in (1..=16).rev() {
                    let offset = NWRowB_OFF[row];
                    let degree = NWRowB_DEG[row];
                    let mut p = NW_D[offset];
                    for j in 1..degree {
                        p = p.mul_add(y, NW_D[offset + j]);
                    }
                    old_b = old_b.mul_add(x, p * y);
                }
                let new_b = rows_b!(x, y);
                assert_eq!(new_b.to_bits(), old_b.to_bits(), "B x={x}, y={y}");
            }
        }
    }

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
