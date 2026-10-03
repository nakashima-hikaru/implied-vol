use crate::lets_be_rational::special_function::normal_distribution::norm_pdf;
use crate::{SpecialFn, lets_be_rational};
use std::f64::consts::FRAC_1_SQRT_2;

const ACTUAR_MAX_QUANTILE_ITERATIONS: usize = 100;
const ACTUAR_QUANTILE_TOLERANCE: f64 = 1.0e-14;
const LARGE_KAPPA_THRESHOLD: f64 = 1.0e3;
const LEFT_SMALL_PROBABILITY_LOG_THRESHOLD: f64 = -11.51;
const RIGHT_SMALL_PROBABILITY_LOG_THRESHOLD: f64 = -1.0e-5;
const MAX_QUANTILE_ITERATIONS: usize = 64;
const QUANTILE_BRACKET_TOLERANCE: f64 = 32.0 * f64::EPSILON;

#[inline]
pub(crate) fn implied_black_volatility_normalised<SpFn: SpecialFn>(
    log_moneyness: f64,
    normalised_price: f64,
) -> Option<f64> {
    if normalised_price <= 0.0 {
        return (normalised_price == 0.0).then_some(0.0);
    }

    let abs_x = log_moneyness.abs();
    let b_max = (-0.5 * abs_x).exp();
    if normalised_price >= b_max {
        return (normalised_price == b_max).then_some(f64::INFINITY);
    }

    if abs_x == 0.0 {
        return Some(lets_be_rational::implied_normalised_volatility_atm::<SpFn>(
            normalised_price,
        ));
    }

    let otm_call_price = normalised_price / b_max;
    if otm_call_price < 1.0e-5 {
        let x = inverse_gaussian_quantile_from_survival::<SpFn>(otm_call_price, 2.0 / abs_x)?;
        return Some(2.0 / x.sqrt());
    }
    let probability = 1.0 - otm_call_price;
    if probability <= 0.0 {
        return Some(f64::INFINITY);
    }
    if probability >= 1.0 {
        return Some(0.0);
    }

    let x = inverse_gaussian_quantile::<SpFn>(probability, 2.0 / abs_x)?;
    Some(2.0 / x.sqrt())
}

#[inline]
pub(crate) fn implied_black_volatility_input<SpFn: SpecialFn, const IS_CALL: bool>(
    price: f64,
    f: f64,
    k: f64,
    t: f64,
) -> Option<f64> {
    let intrinsic_value = if IS_CALL { f - k } else { k - f };
    if t == 0.0 {
        return (price == intrinsic_value.max(0.0)).then_some(0.0);
    }
    if price >= if IS_CALL { f } else { k } {
        return (price == if IS_CALL { f } else { k }).then_some(f64::INFINITY);
    }

    let normalised_time_value = if intrinsic_value > 0.0 {
        price - intrinsic_value
    } else {
        price
    } / (f.sqrt() * k.sqrt());

    if normalised_time_value <= f64::MIN_POSITIVE {
        return (normalised_time_value >= 0.0).then_some(0.0);
    }

    Some(
        implied_black_volatility_normalised::<SpFn>(
            lets_be_rational::bs_option_price::negative_abs_log_moneyness(f, k),
            normalised_time_value,
        )? / t.sqrt(),
    )
}

#[inline]
fn inverse_gaussian_quantile<SpFn: SpecialFn>(probability: f64, mu: f64) -> Option<f64> {
    if !(probability > 0.0 && probability < 1.0 && mu > 0.0 && mu.is_finite()) {
        return None;
    }

    inverse_gaussian_quantile_actuar::<SpFn>(probability, mu)
        .filter(|x| x.is_finite() && *x > 0.0)
        .or_else(|| inverse_gaussian_quantile_bracketed::<SpFn>(probability, mu))
}

#[inline]
fn inverse_gaussian_quantile_actuar<SpFn: SpecialFn>(probability: f64, mu: f64) -> Option<f64> {
    let log_probability = probability.ln();
    let phi = mu;
    let mode = inverse_gaussian_mode_standardized(phi);
    let mut x = if log_probability < LEFT_SMALL_PROBABILITY_LOG_THRESHOLD {
        let z = SpFn::inverse_norm_cdf(probability);
        1.0 / (phi * z * z)
    } else if log_probability > RIGHT_SMALL_PROBABILITY_LOG_THRESHOLD {
        (inverse_gaussian_quantile_initial_guess::<SpFn>(probability, mu) / mu).max(mode)
    } else {
        mode
    };

    if !x.is_finite() || !(x > 0.0) {
        x = mode.max(f64::MIN_POSITIVE);
    }

    let mut dx =
        inverse_gaussian_nrstep_standardized::<SpFn>(x, probability, log_probability, phi)?;
    let direction = dx.signum();
    x += dx;
    if !(x > 0.0) || !x.is_finite() {
        return None;
    }

    for _ in 1..ACTUAR_MAX_QUANTILE_ITERATIONS {
        if dx.abs() <= ACTUAR_QUANTILE_TOLERANCE {
            return Some(x * mu);
        }

        dx = inverse_gaussian_nrstep_standardized::<SpFn>(x, probability, log_probability, phi)?;
        if dx * direction < 0.0 {
            return Some(x * mu);
        }

        x += dx;
        if !(x > 0.0) || !x.is_finite() {
            return None;
        }
    }

    Some(x * mu)
}

#[inline]
fn inverse_gaussian_mode_standardized(phi: f64) -> f64 {
    let kappa = 1.5 * phi;
    if kappa <= LARGE_KAPPA_THRESHOLD {
        (1.0 + kappa * kappa).sqrt() - kappa
    } else {
        let reciprocal = 0.5 / kappa;
        reciprocal * (1.0 - reciprocal * reciprocal)
    }
}

#[inline]
fn inverse_gaussian_nrstep_standardized<SpFn: SpecialFn>(
    x: f64,
    probability: f64,
    log_probability: f64,
    phi: f64,
) -> Option<f64> {
    let cdf = inverse_gaussian_cdf_standardized::<SpFn>(x, phi);
    let pdf = inverse_gaussian_pdf_standardized(x, phi);
    if !(pdf > 0.0) || !pdf.is_finite() {
        return None;
    }

    let numerator = if cdf > 0.0 {
        let delta_log_probability = log_probability - cdf.ln();
        if delta_log_probability.abs() < 1.0e-5 {
            delta_log_probability * (log_probability + (-0.5 * delta_log_probability).ln_1p()).exp()
        } else {
            probability - cdf
        }
    } else {
        probability
    };

    Some(numerator / pdf)
}

#[inline]
fn inverse_gaussian_cdf_standardized<SpFn: SpecialFn>(x: f64, phi: f64) -> f64 {
    inverse_gaussian_cdf::<SpFn>(x * phi, phi)
}

#[inline]
fn inverse_gaussian_pdf_standardized(x: f64, phi: f64) -> f64 {
    phi * inverse_gaussian_pdf(x * phi, phi)
}

#[inline]
fn inverse_gaussian_quantile_bracketed<SpFn: SpecialFn>(probability: f64, mu: f64) -> Option<f64> {
    let mut x = inverse_gaussian_quantile_initial_guess::<SpFn>(probability, mu);
    if !x.is_finite() || !(x > 0.0) {
        x = mu;
    }

    let cdf_at_guess = inverse_gaussian_cdf::<SpFn>(x, mu);
    let (mut lower, mut upper) = if cdf_at_guess < probability {
        bracket_quantile_from_below::<SpFn>(x, mu, probability)
    } else {
        bracket_quantile_from_above::<SpFn>(x, mu, probability)
    };

    x = x.clamp(
        if lower > 0.0 {
            lower
        } else {
            f64::MIN_POSITIVE
        },
        upper,
    );

    for _ in 0..MAX_QUANTILE_ITERATIONS {
        let cdf = inverse_gaussian_cdf::<SpFn>(x, mu);
        let error = cdf - probability;
        if error == 0.0 {
            return Some(x);
        }
        if error <= 0.0 {
            lower = x;
        } else {
            upper = x;
        }

        if bracket_is_tight(lower, upper) {
            return Some(bracket_midpoint(lower, upper));
        }

        let pdf = inverse_gaussian_pdf(x, mu);
        let newton = if pdf > 0.0 && pdf.is_finite() {
            let step = error / pdf;
            if step.abs() <= QUANTILE_BRACKET_TOLERANCE * x.max(1.0) {
                return Some(x - step);
            }
            x - step
        } else {
            f64::NAN
        };
        let midpoint = bracket_midpoint(lower, upper);
        let next = if newton.is_finite() && newton > lower && newton < upper {
            newton
        } else {
            midpoint
        };

        if next == x {
            return Some(midpoint);
        }
        x = next;
    }

    Some(bracket_midpoint(lower, upper))
}

#[inline]
fn inverse_gaussian_quantile_from_survival<SpFn: SpecialFn>(
    probability: f64,
    mu: f64,
) -> Option<f64> {
    if !(probability > 0.0 && probability < 1.0 && mu > 0.0 && mu.is_finite()) {
        return None;
    }

    // Use the upper-tail probability directly: 1 - probability rounds to one
    // for small positive option prices and cannot specify the quantile.
    let z = -SpFn::inverse_norm_cdf(probability);
    let mut x = inverse_gaussian_quantile_initial_guess_from_z(z, mu);
    if !x.is_finite() || !(x > 0.0) {
        x = mu;
    }

    let mut lower = x;
    let mut upper = x;
    if inverse_gaussian_survival::<SpFn>(x, mu) > probability {
        loop {
            upper = (upper * 2.0).min(f64::MAX);
            if inverse_gaussian_survival::<SpFn>(upper, mu) <= probability {
                break;
            }
            if upper == f64::MAX {
                // The required quantile exceeds the representable range.
                return None;
            }
            lower = upper;
        }
    } else {
        loop {
            let next = lower * 0.5;
            if next <= f64::MIN_POSITIVE {
                lower = 0.0;
                break;
            }
            lower = next;
            if inverse_gaussian_survival::<SpFn>(lower, mu) >= probability {
                break;
            }
            upper = lower;
        }
    }

    for _ in 0..MAX_QUANTILE_ITERATIONS {
        let survival = inverse_gaussian_survival::<SpFn>(x, mu);
        let error = survival - probability;
        if error == 0.0 {
            return Some(x);
        }
        if error > 0.0 {
            lower = x;
        } else {
            upper = x;
        }

        if bracket_is_tight(lower, upper) {
            return Some(bracket_midpoint(lower, upper));
        }

        let pdf = inverse_gaussian_pdf(x, mu);
        // The survival function is decreasing, so its derivative is -pdf.
        let next = if pdf > 0.0 && pdf.is_finite() {
            let step = if survival > 0.0 && (survival / probability - 1.0).abs() > 0.5 {
                (survival.ln() - probability.ln()) * (survival / pdf)
            } else {
                error / pdf
            };
            if step.abs() <= QUANTILE_BRACKET_TOLERANCE * x {
                return Some(x + step);
            }
            x + step
        } else {
            f64::NAN
        };
        x = if next.is_finite() && next > lower && next < upper {
            next
        } else {
            bracket_midpoint(lower, upper)
        };
    }

    Some(bracket_midpoint(lower, upper))
}

#[inline]
fn inverse_gaussian_survival<SpFn: SpecialFn>(x: f64, mu: f64) -> f64 {
    let sqrt_x = x.sqrt();
    // The same stable erfcx difference used for the normalized Black price
    // equals this survival probability times exp(-1 / mu).
    lets_be_rational::bs_option_price::normalised_black::<SpFn>(
        -mu.recip(),
        -sqrt_x / mu,
        sqrt_x.recip(),
    ) / (-mu.recip()).exp()
}

#[inline]
fn inverse_gaussian_quantile_initial_guess<SpFn: SpecialFn>(probability: f64, mu: f64) -> f64 {
    let z = SpFn::inverse_norm_cdf(probability);
    inverse_gaussian_quantile_initial_guess_from_z(z, mu)
}

#[inline]
fn inverse_gaussian_quantile_initial_guess_from_z(z: f64, mu: f64) -> f64 {
    let sqrt_discriminant = (mu.mul_add(mu * z * z, 4.0 * mu)).sqrt();
    let sqrt_x = if z >= 0.0 {
        0.5 * mu.mul_add(z, sqrt_discriminant)
    } else {
        2.0 * mu / (sqrt_discriminant - mu * z)
    };
    sqrt_x * sqrt_x
}

#[inline]
fn bracket_quantile_from_below<SpFn: SpecialFn>(
    initial: f64,
    mu: f64,
    probability: f64,
) -> (f64, f64) {
    let mut lower = initial;
    let mut upper = initial;
    loop {
        let next = upper * 2.0;
        upper = if next.is_finite() { next } else { f64::MAX };
        if inverse_gaussian_cdf::<SpFn>(upper, mu) >= probability || upper == f64::MAX {
            return (lower, upper);
        }
        lower = upper;
    }
}

#[inline]
fn bracket_quantile_from_above<SpFn: SpecialFn>(
    initial: f64,
    mu: f64,
    probability: f64,
) -> (f64, f64) {
    let mut lower = initial;
    let mut upper = initial;
    loop {
        let next = lower * 0.5;
        if next <= f64::MIN_POSITIVE {
            return (0.0, upper);
        }
        if inverse_gaussian_cdf::<SpFn>(next, mu) <= probability {
            return (next, upper);
        }
        upper = next * 2.0;
        lower = next;
    }
}

#[inline]
fn bracket_midpoint(lower: f64, upper: f64) -> f64 {
    debug_assert!(upper.is_finite() && upper > 0.0);
    if lower <= 0.0 {
        0.5 * upper
    } else {
        (0.5 * (lower.ln() + upper.ln())).exp()
    }
}

#[inline]
fn bracket_is_tight(lower: f64, upper: f64) -> bool {
    if lower <= 0.0 {
        upper <= f64::MIN_POSITIVE
    } else {
        upper / lower - 1.0 <= QUANTILE_BRACKET_TOLERANCE
    }
}

#[inline]
fn inverse_gaussian_cdf<SpFn: SpecialFn>(x: f64, mu: f64) -> f64 {
    debug_assert!(x > 0.0 && mu > 0.0);
    let sqrt_x = x.sqrt();
    let inv_sqrt_x = sqrt_x.recip();
    let leading = sqrt_x / mu;
    let a = leading - inv_sqrt_x;
    let reflected = leading + inv_sqrt_x;

    SpFn::norm_cdf(a) + 0.5 * (-0.5 * a * a).exp() * SpFn::erfcx(reflected * FRAC_1_SQRT_2)
}

#[inline]
fn inverse_gaussian_pdf(x: f64, mu: f64) -> f64 {
    debug_assert!(x > 0.0 && mu > 0.0);
    let sqrt_x = x.sqrt();
    let a = sqrt_x / mu - sqrt_x.recip();
    norm_pdf(a * a) / (x * sqrt_x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DefaultSpecialFn, PriceBlackScholesNormalised};

    const PAPER_DELTAS: [f64; 8] = [0.05, 0.20, 0.30, 0.45, 0.55, 0.70, 0.80, 0.95];

    #[test]
    fn inverse_gaussian_density_matches_closed_form() {
        // At x = 4 and mu = 1 the normal argument is 1.5, so the density is
        // exp(-1.5^2 / 2) / (8 * sqrt(2 * pi)).
        let expected = 0.016_189_699_458_236_467;
        assert!((inverse_gaussian_pdf(4.0, 1.0) - expected).abs() <= f64::EPSILON * expected);
    }

    #[test]
    fn explicit_formula_preserves_small_positive_prices() {
        for normalised_price in [1.0e-20, 1.0e-100, 1.0e-300] {
            let recovered =
                implied_black_volatility_normalised::<DefaultSpecialFn>(1.0, normalised_price)
                    .unwrap();
            assert!(recovered.is_finite() && recovered > 0.0);
            let repriced = PriceBlackScholesNormalised::builder()
                .log_moneyness(1.0)
                .total_volatility(recovered)
                .build()
                .unwrap()
                .calculate::<DefaultSpecialFn>();
            assert!(
                (repriced / normalised_price - 1.0).abs() <= 2.0e-12,
                "price={normalised_price}, recovered={recovered}, repriced={repriced}"
            );
        }
    }

    #[test]
    fn explicit_formula_rejects_unrepresentable_quantile() {
        // Although volatility is representable, the intermediate quantile
        // 4 / volatility^2 exceeds f64::MAX. Do not use an invalid bracket.
        assert!(implied_black_volatility_normalised::<DefaultSpecialFn>(1e-300, 1e-300).is_none());
    }

    fn paper_total_vols() -> Vec<f64> {
        let mut vols = Vec::with_capacity(41);
        vols.push(0.01);
        let mut current = 0.05;
        while current <= 2.0 + f64::EPSILON {
            vols.push(current);
            current += 0.05;
        }
        vols
    }

    #[test]
    fn explicit_formula_recovers_paper_grid() {
        for total_volatility in paper_total_vols() {
            for delta in PAPER_DELTAS {
                let log_moneyness = total_volatility
                    * (DefaultSpecialFn::inverse_norm_cdf(delta) - 0.5 * total_volatility);
                let normalised_price = PriceBlackScholesNormalised::builder()
                    .log_moneyness(log_moneyness)
                    .total_volatility(total_volatility)
                    .build()
                    .unwrap()
                    .calculate::<DefaultSpecialFn>();

                let recovered = implied_black_volatility_normalised::<DefaultSpecialFn>(
                    log_moneyness,
                    normalised_price,
                )
                .unwrap();

                assert!(
                    (recovered - total_volatility).abs() <= 1.0e-12 * total_volatility.max(1.0),
                    "delta={delta}, log_moneyness={log_moneyness}, v={total_volatility}, recovered={recovered}"
                );
            }
        }
    }

    #[test]
    fn explicit_formula_matches_boundaries() {
        assert_eq!(
            implied_black_volatility_normalised::<DefaultSpecialFn>(0.25, 0.0).unwrap(),
            0.0
        );

        let b_max = (-0.5_f64).exp();
        assert_eq!(
            implied_black_volatility_normalised::<DefaultSpecialFn>(1.0, b_max).unwrap(),
            f64::INFINITY
        );
        assert!(
            implied_black_volatility_normalised::<DefaultSpecialFn>(1.0, b_max * 1.000_000_1)
                .is_none()
        );
        assert_eq!(
            implied_black_volatility_normalised::<DefaultSpecialFn>(0.0, 1.0),
            Some(f64::INFINITY)
        );
        assert!(implied_black_volatility_normalised::<DefaultSpecialFn>(0.0, 1.01).is_none());
    }
}
