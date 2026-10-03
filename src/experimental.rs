//! Native Rust implementation of the experimental normalized Black implied volatility solver.
//!
//! Imported from the current Rust port in `implied-black-volatility`, whose
//! numerical modules correspond to the retained research reference archive.
//! The copied modules preserve their coefficients, evaluation order,
//! native vector lanes, and the experimental port's internal Jaeckel translation.
//! `performance/experimental-source-2026-10-03.json` records the source correspondence.
//!
//! The inputs are `a = abs(ln(F/K))` and the undiscounted out-of-the-money
//! option price divided by `sqrt(F*K)`. The result is **total** volatility,
//! `sigma * sqrt(T)`. Market normalization is the caller's responsibility.
//!
//! The numerical contract targets finite `a >= 0` and positive normal prices
//! below the exact cap `exp(-a/2)`. The original specialized tiny-price paths
//! are retained, but a uniform contract for all subnormal prices is not claimed.
//! Uses its own math functions, explicit fused multiply-add, and the platform's
//! standard math library. These fused operations are retained independently of
//! this crate's `fma` feature. Requires round-to-nearest, gradual underflow, and
//! no unsafe reassociation. No C ABI entry points are imported.

// Keep the source's explicit range guards and NaN-rejecting comparisons visible.
#![allow(clippy::manual_range_contains, clippy::neg_cmp_op_on_partial_ord)]
// Preserve the imported literals, guarded index casts, control flow, and exact
// fused/unfused arithmetic. These observed style lints would otherwise suggest
// rewriting the numerical program; the allowances apply only to this module.
#![allow(
    clippy::bool_to_int_with_if,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::collapsible_else_if,
    clippy::doc_markdown,
    clippy::if_not_else,
    clippy::items_after_statements,
    clippy::manual_assert_eq,
    clippy::manual_let_else,
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::redundant_else,
    clippy::redundant_pub_crate,
    clippy::suboptimal_flops,
    clippy::suspicious_operation_groupings,
    clippy::unreadable_literal,
    clippy::useless_let_if_seq,
    clippy::wildcard_imports
)]

mod lanes;
mod large;
mod large_constants;
mod lbr_fallback;
mod math;
mod small;

/// Invert the normalized out-of-the-money Black price to total volatility.
///
/// `a` is absolute log-moneyness, `b` is price divided by `sqrt(F*K)`.
/// Returns `NaN` for unsupported or invalid inputs, following the experimental dispatcher.
/// This function allocates no memory and calls no C++ solver.
#[inline]
pub fn implied_total_volatility(a: f64, b: f64) -> f64 {
    if a <= 10.0 {
        small::solve(a, b)
    } else {
        large::solve(a, b)
    }
}
