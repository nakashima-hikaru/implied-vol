//! Error-free transformations restricted to the last one or two Clenshaw cells.
//!
//! Unlike full twofold Clenshaw, the high-order polynomial is evaluated by the
//! original scalar FMAs. The low terms of the final cells are reconstructed with
//! a fused-product residual and TwoSum. The independent-root tests verify the
//! represented arithmetic; this is not a continuous error certificate.
use super::tables::{NW_F1, NW_G0};

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let z = s - a;
    (s, (a - (s - z)) + (b - z))
}

#[inline(always)]
fn fma_error(a: f64, b: f64, c: f64, rounded: f64) -> f64 {
    // a*b+c-rounded = [a*b-round(a*b)] +
    //                  [round(a*b)+c-rounded].
    // The first bracket uses an exact fused multiplication residual;
    // TwoSum retains the second bracket even under near-cancellation.
    let p = a * b;
    let prod_low = a.mul_add(b, -p);
    let (sum, sum_low) = two_sum(p, c);
    (sum - rounded) + (sum_low + prod_low)
}

#[inline(always)]
fn stage(
    factor: f64,
    value: f64,
    coefficient: f64,
    previous: f64,
    value_low: f64,
    previous_low: f64,
) -> (f64, f64) {
    let fused = factor.mul_add(value, coefficient);
    let (high, subtract_low) = two_sum(fused, -previous);
    let low = fma_error(factor, value, coefficient, fused)
        + subtract_low
        + factor.mul_add(value_low, -previous_low);
    (high, low)
}

#[inline(always)]
fn clenshaw_tail<const N: usize, const LAST: usize>(coefficients: &[f64; N], x: f64) -> (f64, f64) {
    const {
        assert!(LAST < N);
    }
    let x2 = x + x;
    let mut e0 = 0.0_f64;
    let mut e1 = 0.0_f64;
    for i in (LAST + 1..N).rev() {
        let next = x2.mul_add(e0, coefficients[i]) - e1;
        e1 = e0;
        e0 = next;
    }
    let mut e0_low = 0.0_f64;
    let mut e1_low = 0.0_f64;
    for i in (1..=LAST).rev() {
        let (next, next_low) = stage(x2, e0, coefficients[i], e1, e0_low, e1_low);
        e1 = e0;
        e0 = next;
        e1_low = e0_low;
        e0_low = next_low;
    }
    stage(x, e0, coefficients[0], e1, e0_low, e1_low)
}

/// Endpoint-compensated G0 and corrected h/A0 in region A.
#[inline]
#[allow(clippy::suboptimal_flops)] // Preserve the verified FMA/rounding order.
pub(super) fn region_a_scale(h: f64, a: f64) -> (f64, f64) {
    let x = a.mul_add(0.318_309_886_183_790_7, -1.0);
    let (g, g_low) = clenshaw_tail::<23, 0>(&NW_G0, x);
    let a0 = a * g;
    let a0_low = a.mul_add(g, -a0) + a * g_low;
    let initial = h / a0;
    let residual = (-initial).mul_add(a0, h) - initial * a0_low;
    let t = a0.recip().mul_add(residual, initial);
    (a0, t)
}

/// Endpoint-compensated F1 and corrected 1/sqrt(a) in region B.
#[inline]
#[allow(clippy::suboptimal_flops)] // Preserve the verified FMA/rounding order.
pub(super) fn region_b_scale(a: f64) -> f64 {
    let root = a.sqrt();
    let inv = root.recip();
    let square_residual = (-root).mul_add(root, a);
    let reciprocal_residual = (-root).mul_add(inv, 1.0);
    let correction = inv * (reciprocal_residual - (0.5 * square_residual) * (inv * inv));
    let (u, u_low) = two_sum(inv, correction);
    let x = u.mul_add(5.537_929_048_854_298, -1.209_314_043_451_270_6);
    let (f1, f1_low) = clenshaw_tail::<40, 1>(&NW_F1, x);
    let pp = u * f1;
    let pp_low = u.mul_add(f1, -pp) + u * f1_low + u_low * f1;
    pp + pp_low
}
