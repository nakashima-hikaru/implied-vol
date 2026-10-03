//! `FlashIV` log-price Householder inversion.
//!
//! Independently implemented from Le Floch and Healy, arXiv:2605.29102v1,
//! equations (logc), (logvega), (d2), (d3), (h3), and Algorithm 1:
//! <https://arxiv.org/abs/2605.29102v1>
//! The Li seed coefficients are numerical values from Li and Lee, MPRA 6867,
//! equation (42), PDF page 24 / printed page 22:
//! <https://mpra.ub.uni-muenchen.de/6867/1/MPRA_paper_6867.pdf>
//! Abramowitz--Stegun 7.1.26 supplies the inexpensive erfcx polynomial.
//!
//! The restricted entry for the default hybrid follows the paper's cheap H3
//! pre-step and two exact H3 steps, with a third when the residual entering the
//! second exact step is at least 1e-4. The pure `flashiv` entry uses its own
//! Li/asymptotic/Mills seeds, safeguarded dimensionless H3 steps, and a
//! complementary objective near the upper price bound. It never invokes the
//! Jaeckel inverse. Exact iterations may continue until floating-point
//! convergence, rather than imposing the paper's fixed hot-path count.
//!
//! Differences from the paper: use sqrt-forward beta space (the x/2 terms
//! cancel); evaluate the exact objective with the existing scaled Black
//! expansions in their asymptotic/small-t regions to avoid cancellation.
//! If `scaled_B=beta/vega`, then `ln(beta)=ln(scaled_B)+ln(vega)` and
//! `d ln(beta)/dv=1/scaled_B`. Elsewhere use the paper erfcx difference.
//!
//! Only the default hybrid retains LBR outside the restricted entry's qualified
//! lowest-price domain, including `|x|<0.01`. The pure entry handles that domain
//! with scaled variables and the stable price evaluator within its own H3
//! iterations. Near a root it evaluates the log residual with `log1p` of the
//! relative price difference to keep nearby represented prices distinguishable.
//! A pure-mode failure returns `None`; it never switches inverse methods.
//! There is no final output correction or `FlashIV+` price alignment.

// Keep the qualified expression ordering and the crate's opt-in FMA policy.
// Clippy's replacements unconditionally use the standard fused operation.
#![allow(clippy::suboptimal_flops)]

use crate::fused_multiply_add::MulAdd;
use crate::lets_be_rational::bs_option_price::scaled_normalised_black_and_ln_vega;
use crate::lets_be_rational::bs_option_price::uses_scaled_expansion;
#[cfg(feature = "flashiv")]
use crate::lets_be_rational::bs_option_price::{normalised_black, normalised_vega};
use crate::lets_be_rational::special_function::SpecialFn;
use std::f64::consts::FRAC_1_SQRT_2;
use std::f64::consts::LN_2;

const SQRT_2_OVER_PI: f64 = 0.797_884_560_802_865_4;
const INV_SQRT_PI: f64 = 0.564_189_583_547_756_3;
const LN_2_PI: f64 = 1.837_877_066_409_345_3;

/// Self-contained `FlashIV` inversion. The boundary dispatcher supplies an
/// interior sqrt-forward price; no Jaeckel inverse or interpolation is used.
///
/// The paper's cheap pre-step and Li/asymptotic seeds are retained. Exact H3
/// steps use the stable Black evaluator throughout the iteration. Unlike the
/// fixed-count paper hot path, this entry safeguards the seed basin and keeps
/// iterating when the residual still requires it. This is an iterative solver,
/// not an additional correction after another solver has returned a result.
#[cfg(feature = "flashiv")]
#[inline]
pub fn normalised_volatility<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> Option<f64> {
    if !(theta_x <= 0.0 && theta_x.is_finite() && beta > 0.0 && beta < b_max && b_max.is_finite()) {
        return None;
    }
    if theta_x == 0.0 {
        return Some((2.0 * std::f64::consts::SQRT_2) * SpFn::erfinv(beta));
    }
    let m = -theta_x;
    let c = beta / b_max;
    let complementary = beta >= (f64::from_bits(0.99_f64.to_bits() - 1) * b_max);
    let target = if complementary { b_max - beta } else { beta };
    let ln_target = target.ln();
    let seed = if complementary {
        // Solve (h^2+t^2)/2=L on its upper-volatility branch. Including m
        // avoids a very poor seed for large log-moneyness near saturation.
        let half_m = 0.5 * m;
        let l = (-ln_target - 0.5 * LN_2_PI).max(half_m + 1.0);
        2.0 * (l + ((l - half_m) * (l + half_m)).sqrt()).sqrt()
    } else if m < 0.01 && c <= 0.0005 {
        near_atm_seed(m, c, beta)
    } else {
        initial_guess(theta_x, c, beta.ln() + 0.5 * m)
    };
    if !(seed > 0.0 && seed.is_finite()) {
        return None;
    }

    // The inexpensive evaluator is only a pre-step. If its two erfcx values
    // are indistinguishable for microscopic t, retain the seed and use the
    // stable exact objective; this is not a cross-method fallback.
    let inflection = (2.0 * m).sqrt();
    let mut v = pure_objective::<SpFn, true>(theta_x, seed, target, ln_target, complementary)
        .and_then(|state| dimensionless_step(seed, state))
        .filter(|next| {
            *next >= 0.25 * seed && *next <= 4.0 * seed && (!complementary || *next >= inflection)
        })
        .unwrap_or(seed);
    let mut state = pure_objective::<SpFn, false>(theta_x, v, target, ln_target, complementary)?;
    if state.residual == 0.0 {
        return Some(v);
    }

    // Bracket with the same FlashIV objective. Relative doubling/halving
    // preserves the microscopic scale and never constructs the cancellation-
    // prone Jaeckel tangent nodes.
    let below = |residual: f64| {
        if complementary {
            residual > 0.0
        } else {
            residual < 0.0
        }
    };
    let (mut lower, mut upper);
    if below(state.residual) {
        lower = v;
        upper = v;
        let mut found = false;
        for _ in 0..64 {
            upper *= 2.0;
            let upper_state =
                pure_objective::<SpFn, false>(theta_x, upper, target, ln_target, complementary)?;
            if !below(upper_state.residual) {
                found = true;
                break;
            }
            lower = upper;
        }
        if !found {
            return None;
        }
    } else {
        upper = v;
        lower = v;
        let mut found = false;
        for _ in 0..64 {
            lower = if complementary {
                (0.5 * lower).max(inflection)
            } else {
                0.5 * lower
            };
            if lower == 0.0 {
                break;
            }
            let lower_state =
                pure_objective::<SpFn, false>(theta_x, lower, target, ln_target, complementary)?;
            if below(lower_state.residual) {
                found = true;
                break;
            }
            upper = lower;
        }
        if !found {
            return None;
        }
    }
    let mut best = v;
    let mut best_residual = state.residual.abs();
    let mut best_derivative = state.derivative.abs();
    let mut previous = f64::NAN;
    for _ in 0..32 {
        if below(state.residual) {
            lower = v;
        } else {
            upper = v;
        }
        let proposed = dimensionless_step(v, state);
        let converged = best_residual <= 2.0 * f64::EPSILON * best_derivative
            || upper - lower <= (4.0 * f64::EPSILON) * best;
        let next = match proposed {
            Some(next) if next > lower && next < upper => next,
            Some(next) if (next == v || next == previous) && converged => return Some(best),
            _ => lower + 0.5 * (upper - lower),
        };
        if next == v || next == lower || next == upper {
            return converged.then_some(best);
        }
        previous = v;
        v = next;
        state = pure_objective::<SpFn, false>(theta_x, v, target, ln_target, complementary)?;
        if state.residual.abs() < best_residual {
            best = v;
            best_residual = state.residual.abs();
            best_derivative = state.derivative.abs();
        }
        if state.residual == 0.0 {
            return Some(v);
        }
    }
    // A stalled floating-point iteration can be accepted only once its
    // remaining Newton displacement is at the representation limit.
    (best_residual <= 2.0 * f64::EPSILON * best_derivative
        || upper - lower <= (4.0 * f64::EPSILON) * best)
        .then_some(best)
}

#[cfg(feature = "flashiv")]
#[inline(always)]
fn near_atm_seed(m: f64, c: f64, beta: f64) -> f64 {
    let ln_ratio = beta.ln() - m.ln();
    if ln_ratio < -20.0 {
        // Bachelier/Mills tail: beta/m ~ phi(mu)/mu^3, mu=m/v.
        // Solve u+3ln(u)=L approximately without squaring beta/m.
        let l = -2.0 * ln_ratio - LN_2_PI;
        let mut u = l - 3.0 * l.ln();
        for _ in 0..2 {
            u -= (u + 3.0 * u.ln() - l) / (1.0 + 3.0 / u);
        }
        m / u.sqrt()
    } else {
        m.hypot((2.0 * std::f64::consts::PI).sqrt() * c)
    }
}

#[cfg(feature = "flashiv")]
#[derive(Clone, Copy)]
struct PureState {
    residual: f64,
    derivative: f64,
    h2: f64,
    h3: f64,
}

#[inline(always)]
#[cfg(feature = "flashiv")]
fn pure_objective<SpFn: SpecialFn, const FAST: bool>(
    x: f64,
    v: f64,
    target: f64,
    ln_target: f64,
    complementary: bool,
) -> Option<PureState> {
    if !(v > 0.0 && v.is_finite()) {
        return None;
    }
    let h = x / v;
    let t = 0.5 * v;
    let h2 = h * h;
    let t2 = t * t;
    let q = 0.5 * t.mul_add2(t, h2);
    let (scaled_price, ln_vega) = if complementary {
        let z1 = FRAC_1_SQRT_2 * (t + h);
        let z2 = FRAC_1_SQRT_2 * (t - h);
        let sum = if FAST {
            fast_erfcx(z1) + fast_erfcx(z2)
        } else {
            SpFn::erfcx(z1) + SpFn::erfcx(z2)
        };
        (
            (0.5 * std::f64::consts::PI).sqrt() * sum,
            -q - 0.5 * LN_2_PI,
        )
    } else if FAST {
        let delta = fast_erfcx(-FRAC_1_SQRT_2 * (h + t)) - fast_erfcx(-FRAC_1_SQRT_2 * (h - t));
        (
            (0.5 * std::f64::consts::PI).sqrt() * delta,
            -q - 0.5 * LN_2_PI,
        )
    } else {
        scaled_normalised_black_and_ln_vega::<SpFn>(0.5 * x, h, t)
    };
    if !complementary && !FAST && scaled_price == f64::INFINITY {
        // At an upper bracketing endpoint the erfcx ratio can overflow even
        // though the represented Black price is finite and saturated. Its
        // monotone sign is sufficient; a zero derivative forces a bracket
        // step instead of applying H3 to this endpoint.
        let model = normalised_black::<SpFn>(0.5 * x, h, t);
        return (model.is_finite() && model > target).then_some(PureState {
            residual: 1.0,
            derivative: 0.0,
            h2: 0.0,
            h3: 0.0,
        });
    }
    if !(scaled_price > 0.0 && scaled_price.is_finite() && ln_vega.is_finite()) {
        return None;
    }
    let model = if complementary || FAST {
        scaled_price * normalised_vega(h, t)
    } else {
        normalised_black::<SpFn>(0.5 * x, h, t)
    };
    let relative_price = (model - target) / target;
    let residual = if relative_price > -0.5 && relative_price < 0.5 {
        // Keep nearby prices distinguishable instead of subtracting rounded
        // logs whose magnitude grows as a microscopic price tends to zero.
        relative_price.ln_1p()
    } else {
        scaled_price.ln() + ln_vega - ln_target
    };
    let derivative = if complementary {
        -v / scaled_price
    } else {
        v / scaled_price
    };
    let a = h2 - t2;
    let second = a - derivative;
    let third = (-3.0 * derivative).mul_add2(
        second,
        a.mul_add2(a, (-3.0).mul_add2(h2, -t2)) - derivative * derivative,
    );
    if !(residual.is_finite()
        && derivative.is_finite()
        && derivative != 0.0
        && second.is_finite()
        && third.is_finite())
    {
        return None;
    }
    Some(PureState {
        residual,
        derivative,
        h2: second,
        h3: third,
    })
}

#[inline(always)]
#[cfg(feature = "flashiv")]
fn dimensionless_step(v: f64, state: PureState) -> Option<f64> {
    if state.derivative == 0.0 {
        return None;
    }
    let z = -state.residual / state.derivative;
    let numerator = z.mul_add2(0.5 * state.h2, 1.0);
    let denominator = z.mul_add2(z.mul_add2(state.h3 / 6.0, state.h2), 1.0);
    let next = v * (z.mul_add2(numerator / denominator, 1.0));
    (next > 0.0 && next.is_finite()).then_some(next)
}

#[inline]
pub fn try_lowest_branch<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> Option<f64> {
    if !(theta_x <= -0.01 && theta_x.is_finite() && beta > 0.0 && b_max > 0.0) {
        return None;
    }
    let c = beta / b_max;
    if !(c > 0.0 && c < 0.5) {
        return None;
    }
    let ln_beta = beta.ln();
    let ln_c = ln_beta - 0.5 * theta_x;
    let seed = initial_guess(theta_x, c, ln_c).max(1e-10);
    let (v1, _) = h3_step::<SpFn, true>(theta_x, seed, ln_beta)?;
    let (v2, _) = h3_step::<SpFn, false>(theta_x, v1, ln_beta)?;
    let (v3, residual_on_second_entry) = h3_step::<SpFn, false>(theta_x, v2, ln_beta)?;
    if residual_on_second_entry.abs() >= 1e-4 {
        h3_step::<SpFn, false>(theta_x, v3, ln_beta).map(|(v, _)| v)
    } else {
        Some(v3)
    }
}

#[inline(always)]
fn initial_guess(x: f64, c: f64, ln_c: f64) -> f64 {
    let abs_x = -x;
    if abs_x < 3.0 && c > 0.0005 {
        li_seed(x, c)
    } else if c <= 0.5 && ln_c < -2.0 {
        let d2 = -2.0 * ln_c - LN_2_PI;
        (-2.0 * x) / (d2.sqrt() + (d2 - 2.0 * x).sqrt())
    } else {
        (2.0 * abs_x).sqrt()
    }
}

#[inline(always)]
fn li_seed(x: f64, c: f64) -> f64 {
    // Collect sum_{i+j<=3} a_ij x^i c^j as Horner polynomials in c,
    // with coefficient polynomials in x, respecting the opt-in FMA choice.
    let m0 = x
        .mul_add2(-0.216_197_632_156_68, 0.089_753_944_048_51)
        .mul_add2(x, -0.406_619_903_654_27)
        .mul_add2(x, -0.000_061_030_981_65);
    let m1 = x
        .mul_add2(3.838_158_853_945_65, -36.194_052_215_990_28)
        .mul_add2(x, 5.339_676_433_576_88);
    let m2 = x.mul_add2(41.217_726_327_328_34, 3.250_234_253_323_60);
    let numerator = c
        .mul_add2(83.845_932_244_177_96, m2)
        .mul_add2(c, m1)
        .mul_add2(c, m0);
    let n0 = x
        .mul_add2(-0.033_269_442_900_44, 0.430_276_195_531_68)
        .mul_add2(x, -0.484_665_363_616_20)
        .mul_add2(x, 1.0);
    let n1 = x
        .mul_add2(-0.047_638_023_588_53, -1.341_022_799_820_50)
        .mul_add2(x, 22.963_021_090_107_94);
    let n2 = x.mul_add2(2.457_825_742_942_44, -0.772_688_245_324_68);
    let denominator = c
        .mul_add2(-5.705_315_006_451_09, n2)
        .mul_add2(c, n1)
        .mul_add2(c, n0);
    numerator / denominator
}

#[inline(always)]
fn h3_step<SpFn: SpecialFn, const FAST: bool>(
    x: f64,
    v: f64,
    ln_beta_target: f64,
) -> Option<(f64, f64)> {
    if !(v > 0.0 && v.is_finite()) {
        return None;
    }
    let inv_v = v.recip();
    let h = x * inv_v;
    let t = 0.5 * v;
    let (log_beta, d1) = if FAST || !uses_scaled_expansion(h, t) {
        let z_plus = -FRAC_1_SQRT_2 * (h + t);
        let z_minus = -FRAC_1_SQRT_2 * (h - t);
        let delta = if FAST {
            fast_erfcx(z_plus) - fast_erfcx(z_minus)
        } else {
            SpFn::erfcx(z_plus) - SpFn::erfcx(z_minus)
        };
        if !(delta > 0.0 && delta.is_finite()) {
            return None;
        }
        (
            (-0.5).mul_add2(t.mul_add2(t, h * h), delta.ln() - LN_2),
            SQRT_2_OVER_PI / delta,
        )
    } else {
        let (scaled_b, ln_vega) = scaled_normalised_black_and_ln_vega::<SpFn>(0.5 * x, h, t);
        if !(scaled_b > 0.0 && scaled_b.is_finite() && ln_vega.is_finite()) {
            return None;
        }
        (scaled_b.ln() + ln_vega, scaled_b.recip())
    };
    let residual = log_beta - ln_beta_target;
    let h2 = h * h;
    let t2 = t * t;
    let h2_minus_t2 = h2 - t2;
    let b2 = h2_minus_t2 * inv_v;
    let d2 = b2 - d1;
    let b3 = h2_minus_t2.mul_add2(h2_minus_t2, (-3.0).mul_add2(h2, -t2)) * (inv_v * inv_v);
    let d3 = (-3.0 * d1).mul_add2(d2, b3 - d1 * d1);
    let nu = -residual / d1;
    let numerator = nu.mul_add2(0.5 * d2, 1.0);
    let denominator = nu.mul_add2(nu.mul_add2(d3 / 6.0, d2), 1.0);
    let next = (nu * (numerator / denominator)) + v;
    if next > 0.0 && next.is_finite() {
        Some((next, residual))
    } else {
        None
    }
}

#[inline(always)]
fn fast_erfcx(z: f64) -> f64 {
    if z < 0.0 {
        return (z * z).exp().mul_add2(2.0, -fast_erfcx(-z));
    }
    if z >= 2.5 {
        let inv_z = z.recip();
        let w = inv_z * inv_z;
        let correction = w.mul_add2(-1.875, 0.75).mul_add2(w, -0.5).mul_add2(w, 1.0);
        return (INV_SQRT_PI * inv_z) * correction;
    }
    let t = z.mul_add2(0.327_591_1, 1.0).recip();
    t.mul_add2(1.061_405_429, -1.453_152_027)
        .mul_add2(t, 1.421_413_741)
        .mul_add2(t, -0.284_496_736)
        .mul_add2(t, 0.254_829_592)
        * t
}
