//! Empirical error allowance for the fixed-count paper objective.
//!
//! The direct erfcx difference and subtraction of logarithms can lose bits
//! which the stable Black evaluator preserves. This allowance propagates
//! those errors to relative volatility; it is not Jäckel's precision contract
//! or a proof of convergence over the entire input domain.

use implied_vol::{DefaultSpecialFn, SpecialFn};

#[allow(clippy::suboptimal_flops)]
fn relative_volatility_error(x: f64, price_error: f64, s: f64) -> f64 {
    if price_error == 0.0 {
        return 0.0;
    }
    let h = x / s;
    let t = 0.5 * s;
    (price_error.ln() + 0.5 * (2.0 * std::f64::consts::PI).ln() + 0.5 * h.mul_add(h, t * t)
        - s.ln())
    .exp()
}

#[allow(clippy::suboptimal_flops)]
pub fn precision_limit(x: f64, beta: f64, s: f64) -> f64 {
    let m = x.abs();
    let b_max = (-0.5 * m).exp();
    let c = beta / b_max;
    let base = 4.0 * f64::EPSILON * (1.0 + relative_volatility_error(x, beta, s));
    if c <= 1e-6 && m <= 1e-8 {
        // The source rational seed has approximately 1e-14 relative precision.
        return base + 64.0 * f64::EPSILON;
    }
    if c >= 0.99_f64.next_down() {
        let gap = b_max - beta;
        return base
            + 8.0 * f64::EPSILON * (1.0 + gap.ln().abs()) * relative_volatility_error(x, gap, s);
    }
    let h = -m / s;
    let t = 0.5 * s;
    let z1 = -std::f64::consts::FRAC_1_SQRT_2 * (h + t);
    let z2 = -std::f64::consts::FRAC_1_SQRT_2 * (h - t);
    let erfcx_roundoff = (DefaultSpecialFn::erfcx(z1) + DefaultSpecialFn::erfcx(z2))
        / (s * (2.0 / std::f64::consts::PI).sqrt());
    let log_roundoff = beta.ln().abs() * relative_volatility_error(x, beta, s);
    base + 8.0 * f64::EPSILON * (erfcx_roundoff + log_roundoff)
}
