//! Compensated volfi chart coordinates in a pair of binary64 values.
//!
//! The double-double steps follow the same original polynomial coefficients,
//! with error-free two-sum / fused-product residuals. No correctly-rounded
//! libm or whole-domain price-to-volatility theorem is implied.
use super::tables::{NW_F1, NW_G0};

#[derive(Clone, Copy)]
struct Pair {
    hi: f64,
    lo: f64,
}

#[inline(always)]
fn exact(hi: f64) -> Pair {
    Pair { hi, lo: 0.0 }
}

#[inline(always)]
fn two_sum(a: f64, b: f64) -> Pair {
    let s = a + b;
    let z = s - a;
    Pair {
        hi: s,
        lo: (a - (s - z)) + (b - z),
    }
}

#[inline(always)]
fn add(a: Pair, b: Pair) -> Pair {
    let sum = two_sum(a.hi, b.hi);
    two_sum(sum.hi, (a.lo + b.lo) + sum.lo)
}

#[inline(always)]
fn sub(a: Pair, b: Pair) -> Pair {
    add(a, Pair { hi: -b.hi, lo: -b.lo })
}

#[inline(always)]
fn mul(a: Pair, b: Pair) -> Pair {
    let prod = a.hi * b.hi;
    let mut lo = a.hi.mul_add(b.hi, -prod);
    lo = a.hi.mul_add(b.lo, lo);
    lo = a.lo.mul_add(b.hi, lo);
    lo = a.lo.mul_add(b.lo, lo);
    two_sum(prod, lo)
}

#[inline(always)]
fn clenshaw<const N: usize>(coefficients: &[f64; N], x: f64) -> Pair {
    let x2 = x + x;
    let mut e0 = exact(0.0);
    let mut e1 = exact(0.0);
    for i in (1..N).rev() {
        let next = add(sub(mul(exact(x2), e0), e1), exact(coefficients[i]));
        e1 = e0;
        e0 = next;
    }
    add(sub(mul(exact(x), e0), e1), exact(coefficients[0]))
}

/// Compensated A0 = a*G(a) and t = h/A0 in volfi region A.
///
/// The division is corrected once with the error of the product and the
/// low part of A0. This preserves its original uncorrected chart separately.
#[inline]
pub(super) fn region_a_scale(h: f64, a: f64) -> (f64, f64) {
    let x = a.mul_add(0.318_309_886_183_790_7, -1.0);
    let g = clenshaw(&NW_G0, x);
    let a0 = mul(exact(a), g);
    let initial = h / a0.hi;
    let residual = (-initial).mul_add(a0.hi, h) - initial * a0.lo;
    let t = a0.hi.recip().mul_add(residual, initial);
    (a0.hi, t)
}

/// FMA-refined 1/sqrt(a), twofold Clenshaw F1, then compensated product.
#[inline]
#[allow(clippy::suboptimal_flops)] // Exact operation sequence of the research correction.
pub(super) fn region_b_scale(a: f64) -> f64 {
    let root = a.sqrt();
    let inv = root.recip();
    let square_residual = (-root).mul_add(root, a);
    let reciprocal_residual = (-root).mul_add(inv, 1.0);
    let correction =
        inv * (reciprocal_residual - (0.5 * square_residual) * (inv * inv));
    let u = two_sum(inv, correction);
    let x = u.hi.mul_add(5.537_929_048_854_298, -1.209_314_043_451_270_6);
    let f1 = clenshaw(&NW_F1, x);
    let pp = mul(u, f1);
    pp.hi + pp.lo
}
