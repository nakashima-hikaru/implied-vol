use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};
use crate::constants::{DENORMALISATION_CUTOFF, FOURTH_ROOT_DBL_EPSILON, FRAC_SQRT_TWO_PI, HALF_OF_LN_TWO_PI, ONE_OVER_SQRT_THREE, SIXTEENTH_ROOT_DBL_EPSILON, SQRT_DBL_MAX, SQRT_MIN_POSITIVE, SQRT_PI_OVER_TWO, SQRT_THREE, SQRT_THREE_OVER_THIRD_ROOT_TWO_PI, SQRT_TWO_OVER_PI, TWO_PI_OVER_SQRT_TWENTY_SEVEN, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC};
use crate::erf_cody::{erfc_cody, erfcx_cody};
use crate::MulAdd;
use crate::normal_distribution::{inverse_norm_cdf, norm_cdf, norm_pdf};
use crate::rational_cubic::{convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side, convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side, rational_cubic_interpolation};


#[inline(always)]
fn householder3_factor(v: f64, h2: f64, h3: f64) -> f64 { (1.0 + 0.5 * h2 * v) / (1.0 + v * (h2 + h3 * v / 6.0)) }

#[inline(always)]
fn householder4_factor(v: f64, h2: f64, h3: f64, h4: f64) -> f64 { (1.0 + v * (h2 + v * h3 / 6.0)) / (1.0 + v * (1.5 * h2 + v * (h2 * h2 / 4.0 + h3 / 3.0 + v * h4 / 24.0))) }

#[inline(always)]
fn normalised_intrinsic(x: f64, q: bool) -> f64 {
    if (q && !x.is_sign_positive()) || (!q && !x.is_sign_negative()) {
        return 0.0;
    }
    let x2 = x.powi(2);
    if x2 < 98.0 * FOURTH_ROOT_DBL_EPSILON {
        let ret = x * (1.0 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))).abs().max(0.0);
        return if q {
            ret
        } else {
            -ret
        };
    }
    let b_max = (0.5 * x).exp();
    let one_over_b_max = b_max.recip();
    let ret = (b_max - one_over_b_max).abs().max(0.0);
    if q {
        ret
    } else {
        -ret
    }
}

#[inline(always)]
fn normalised_intrinsic_call(x: f64) -> f64 {
    normalised_intrinsic(x, true)
}

const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
const SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

#[inline(always)]
fn asymptotic_expansion_of_normalised_black_call_over_vega(h: f64, t: f64) -> f64 {
    assert!((h < -ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD.abs()) && (h + t < -(SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD).abs()));
    let e = (t / h).powi(2);
    let r = (h + t) * (h - t);
    let q = (h / r).powi(2);
    let asymptotic_expansion_sum = 2.0 + q * (-6.0E0 - 2.0 * e + 3.0 * q * (1.0E1 + e * (2.0E1 + 2.0 * e) + 5.0 * q * (-1.4E1 + e * (-7.0E1 + e * (-4.2E1 - 2.0 * e)) + 7.0 * q * (1.8E1 + e * (1.68E2 + e * (2.52E2 + e * (7.2E1 + 2.0 * e))) + 9.0 * q * (-2.2E1 + e * (-3.3E2 + e * (-9.24E2 + e * (-6.6E2 + e * (-1.1E2 - 2.0 * e)))) + 1.1E1 * q * (2.6E1 + e * (5.72E2 + e * (2.574E3 + e * (3.432E3 + e * (1.43E3 + e * (1.56E2 + 2.0 * e))))) + 1.3E1 * q * (-3.0E1 + e * (-9.1E2 + e * (-6.006E3 + e * (-1.287E4 + e * (-1.001E4 + e * (-2.73E3 + e * (-2.1E2 - 2.0 * e)))))) + 1.5E1 * q * (3.4E1 + e * (1.36E3 + e * (1.2376E4 + e * (3.8896E4 + e * (4.862E4 + e * (2.4752E4 + e * (4.76E3 + e * (2.72E2 + 2.0 * e))))))) + 1.7E1 * q * (-3.8E1 + e * (-1.938E3 + e * (-2.3256E4 + e * (-1.00776E5 + e * (-1.84756E5 + e * (-1.51164E5 + e * (-5.4264E4 + e * (-7.752E3 + e * (-3.42E2 - 2.0 * e)))))))) + 1.9E1 * q * (4.2E1 + e * (2.66E3 + e * (4.0698E4 + e * (2.3256E5 + e * (5.8786E5 + e * (7.05432E5 + e * (4.0698E5 + e * (1.08528E5 + e * (1.197E4 + e * (4.2E2 + 2.0 * e))))))))) + 2.1E1 * q * (-4.6E1 + e * (-3.542E3 + e * (-6.7298E4 + e * (-4.90314E5 + e * (-1.63438E6 + e * (-2.704156E6 + e * (-2.288132E6 + e * (-9.80628E5 + e * (-2.01894E5 + e * (-1.771E4 + e * (-5.06E2 - 2.0 * e)))))))))) + 2.3E1 * q * (5.0E1 + e * (4.6E3 + e * (1.0626E5 + e * (9.614E5 + e * (4.08595E6 + e * (8.9148E6 + e * (1.04006E7 + e * (6.53752E6 + e * (2.16315E6 + e * (3.542E5 + e * (2.53E4 + e * (6.0E2 + 2.0 * e))))))))))) + 2.5E1 * q * (-5.4E1 + e * (-5.85E3 + e * (-1.6146E5 + e * (-1.77606E6 + e * (-9.37365E6 + e * (-2.607579E7 + e * (-4.01166E7 + e * (-3.476772E7 + e * (-1.687257E7 + e * (-4.44015E6 + e * (-5.9202E5 + e * (-3.51E4 + e * (-7.02E2 - 2.0 * e)))))))))))) + 2.7E1 * q * (5.8E1 + e * (7.308E3 + e * (2.3751E5 + e * (3.12156E6 + e * (2.003001E7 + e * (6.919458E7 + e * (1.3572783E8 + e * (1.5511752E8 + e * (1.0379187E8 + e * (4.006002E7 + e * (8.58429E6 + e * (9.5004E5 + e * (4.7502E4 + e * (8.12E2 + 2.0 * e))))))))))))) + 2.9E1 * q * (-6.2E1 + e * (-8.99E3 + e * (-3.39822E5 + e * (-5.25915E6 + e * (-4.032015E7 + e * (-1.6934463E8 + e * (-4.1250615E8 + e * (-6.0108039E8 + e * (-5.3036505E8 + e * (-2.8224105E8 + e * (-8.870433E7 + e * (-1.577745E7 + e * (-1.472562E6 + e * (-6.293E4 + e * (-9.3E2 - 2.0 * e)))))))))))))) + 3.1E1 * q * (6.6E1 + e * (1.0912E4 + e * (4.74672E5 + e * (8.544096E6 + e * (7.71342E7 + e * (3.8707344E8 + e * (1.14633288E9 + e * (2.07431664E9 + e * (2.33360622E9 + e * (1.6376184E9 + e * (7.0963464E8 + e * (1.8512208E8 + e * (2.7768312E7 + e * (2.215136E6 + e * (8.184E4 + e * (1.056E3 + 2.0 * e))))))))))))))) + 3.3E1 * (-7.0E1 + e * (-1.309E4 + e * (-6.49264E5 + e * (-1.344904E7 + e * (-1.4121492E8 + e * (-8.344518E8 + e * (-2.9526756E9 + e * (-6.49588632E9 + e * (-9.0751353E9 + e * (-8.1198579E9 + e * (-4.6399188E9 + e * (-1.6689036E9 + e * (-3.67158792E8 + e * (-4.707164E7 + e * (-3.24632E6 + e * (-1.0472E5 + e * (-1.19E3 - 2.0 * e))))))))))))))))) * q))))))))))))))));

    let b_over_vega = (t / r) * asymptotic_expansion_sum;
    b_over_vega.max(0.0).abs()
}

#[inline(always)]
fn small_t_expansion_of_normalised_black_call_over_vega(h: f64, t: f64) -> f64 {
    let w = t.powi(2);
    let h2 = h.powi(2);
    let a = 1_f64 + h * SQRT_PI_OVER_TWO * erfcx_cody(-FRAC_1_SQRT_2 * h);
    let b_over_vega = 2.0 * t * (a + w * ((-1.0 + 3.0 * a + a * h2) / 6.0 + w * ((-7.0 + 15.0 * a + h2 * (-1.0 + 10.0 * a + a * h2)) / 120.0 + w * ((-57.0 + 105.0 * a + h2 * (-18.0 + 105.0 * a + h2 * (-1.0 + 21.0 * a + a * h2))) / 5040.0 + w * ((-561.0 + 945.0 * a + h2 * (-285.0 + 1260.0 * a + h2 * (-33.0 + 378.0 * a + h2 * (-1.0 + 36.0 * a + a * h2)))) / 362880.0 + w * ((-6555.0 + 10395.0 * a + h2 * (-4680.0 + 17325.0 * a + h2 * (-840.0 + 6930.0 * a + h2 * (-52.0 + 990.0 * a + h2 * (-1.0 + 55.0 * a + a * h2))))) / 39916800.0 + ((-89055.0 + 135135.0 * a + h2 * (-82845.0 + 270270.0 * a + h2 * (-20370.0 + 135135.0 * a + h2 * (-1926.0 + 25740.0 * a + h2 * (-75.0 + 2145.0 * a + h2 * (-1.0 + 78.0 * a + a * h2)))))) * w) / 6227020800.0))))));
    b_over_vega.max(0.0).abs()
}

#[inline(always)]
fn normalised_black_call_with_optimal_use_of_codys_functions(x: f64, s: f64) -> f64 {
    const CODYS_THRESHOLD: f64 = 0.46875;
    let h = x / s;
    let t = 0.5 * s;
    let q1 = -FRAC_1_SQRT_2 * (h + t);
    let q2 = -FRAC_1_SQRT_2 * (h - t);
    let two_b: f64 =
        if q1 < CODYS_THRESHOLD {
            if q2 < CODYS_THRESHOLD {
                0.5 * ((0.5 * x).exp() * erfc_cody(q1) - (-0.5 * x).exp() * erfc_cody(q2))
            } else {
                0.5 * ((0.5 * x).exp() * erfc_cody(q1) - (-0.5 * (h.powi(2) + t.powi(2))).exp() * erfcx_cody(q2))
            }
        } else if q2 < CODYS_THRESHOLD {
            0.5 * ((-0.5 * (h.powi(2) + t.powi(2))).exp() * erfcx_cody(q1) - (-0.5 * x).exp() * erfc_cody(q2))
        } else {
            0.5 * ((-0.5 * (h.powi(2) + t.powi(2))).exp() * (erfcx_cody(q1) - erfcx_cody(q2)))
        };
    two_b.abs().max(0.0)
}

#[inline(always)]
fn normalised_vega(x: f64, s: f64) -> f64 {
    let ax = x.abs();
    if ax <= 0.0 {
        FRAC_SQRT_TWO_PI * (-0.125 * s.powi(2)).exp()
    } else if s <= 0.0 || s <= ax * SQRT_MIN_POSITIVE {
        0.0
    } else {
        FRAC_SQRT_TWO_PI * (-0.5 * ((x / s).powi(2) + (0.5 * s).powi(2))).exp()
    }
}

#[inline(always)]
fn ln_normalised_vega(x: f64, s: f64) -> f64 {
    let ax = x.abs();
    if ax <= 0.0 {
        -HALF_OF_LN_TWO_PI - 0.125 * s.powi(2)
    } else if s <= 0.0 {
        f64::MIN
    } else {
        -HALF_OF_LN_TWO_PI - 0.5 * ((x / s).powi(2) + 0.25 * s.powi(2))
    }
}

#[inline(always)]
fn normalised_black_call(x: f64, s: f64) -> f64 {
    if x.is_sign_positive() {
        return normalised_intrinsic_call(x) + normalised_black_call(-x, s);
    }
    if s <= x.abs() * DENORMALISATION_CUTOFF {
        return normalised_intrinsic_call(x);
    }
    if x < s * ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD && (0.5 * s).powi(2) + x < s * (SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD) {
        return asymptotic_expansion_of_normalised_black_call_over_vega(x / s, 0.5 * s) * normalised_vega(x, s);
    }
    if 0.5 * s < SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD {
        return small_t_expansion_of_normalised_black_call_over_vega(x / s, 0.5 * s) * normalised_vega(x, s);
    }
    normalised_black_call_with_optimal_use_of_codys_functions(x, s)
}


#[inline(always)]
fn normalised_black_call_over_vega_and_ln_vega(x: f64, s: f64) -> (f64, f64) {
    if x.is_sign_positive() {
        let (bx, ln_vega) = normalised_black_call_over_vega_and_ln_vega(-x, s);
        return (normalised_intrinsic_call(x) * (-ln_vega).exp() + bx, ln_vega);
    }
    let ln_vega = ln_normalised_vega(x, s);
    if s <= x.abs() * DENORMALISATION_CUTOFF {
        return (normalised_intrinsic_call(x) * (-ln_vega).exp(), ln_vega);
    }
    if x < s * ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD && 0.5 * s.powi(2) + x < s * (SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD) {
        return (asymptotic_expansion_of_normalised_black_call_over_vega(x / s, 0.5 * s), ln_vega);
    }
    if 0.5 * s < SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD {
        return (small_t_expansion_of_normalised_black_call_over_vega(x / s, 0.5 * s), ln_vega);
    }
    (normalised_black_call_with_optimal_use_of_codys_functions(x, s) * (-ln_vega).exp(), ln_vega)
}

#[inline(always)]
fn normalised_black(x: f64, s: f64, theta: bool) -> f64 {
    normalised_black_call(if !theta { -x } else { x }, s)
}

#[inline(always)]
pub(crate) fn black(f: f64, k: f64, sigma: f64, t: f64, q: bool) -> f64 {
    let intrinsic = if !q { k - f } else { f - k }.max(0f64).abs();
    if (q && ((f - k).is_sign_positive())) || (!q && ((f - k).is_sign_negative())) {
        return intrinsic + black(f, k, sigma, t, !q);
    }
    intrinsic.max(((f * k).sqrt()) * normalised_black((f / k).ln(), sigma * t.sqrt(), q))
}

#[inline(always)]
fn compute_f_lower_map_and_first_two_derivatives(x: f64, s: f64) -> (f64, f64, f64) {
    let ax = x.abs();
    let z = ONE_OVER_SQRT_THREE * ax / s;
    let y = z.powi(2);
    let s2 = s.powi(2);
    let phi_m = norm_cdf(-z);
    let phi = norm_pdf(z);
    let fpp = std::f64::consts::FRAC_PI_6 * y / (s2 * s) * phi_m * (8.0 * SQRT_THREE * s * ax + (3.0 * s2 * (s2 - 8.0) - 8.0 * x.powi(2)) * phi_m / phi) * (2.0 * y + 0.25 * s2).exp();
    let (fp, f) = if s.is_subnormal() {
        (1.0, 0.0)
    } else {
        let phi2 = phi_m.powi(2);
        let fp_val = std::f64::consts::TAU * y * phi2 * (y + 0.125 * s.powi(2)).exp();
        let f_val = if x.is_subnormal() {
            0.0
        } else {
            TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (phi2 * phi_m)
        };
        (fp_val, f_val)
    };
    (f, fp, fpp)
}


#[inline(always)]
fn inverse_f_lower_map(x: f64, f: f64) -> f64 {
    if f.is_subnormal() {
        0.0
    } else {
        (x / (SQRT_THREE * inverse_norm_cdf(SQRT_THREE_OVER_THIRD_ROOT_TWO_PI * f.cbrt() / x.abs().cbrt()))).abs()
    }
}

#[inline(always)]
fn compute_f_upper_map_and_first_two_derivatives(x: f64, s: f64) -> (f64, f64, f64) {
    let f = norm_cdf(-0.5 * s);
    let (fp, fpp) = if x.is_subnormal() {
        (-0.5, 0.0)
    } else {
        let w = (x / s).powi(2);
        (-0.5 * (0.5 * w).exp(),
        SQRT_PI_OVER_TWO * ((w + 0.125 * s.powi(2)).exp()) * w / s)
    };

    (f, fp, fpp)
}


#[inline(always)]
fn inverse_f_upper_map(f: f64) -> f64 {
    -2.0 * inverse_norm_cdf(f)
}

#[inline(always)]
fn take_step(x_min: f64, x_max: f64, x: f64, dx: f64) -> (f64, f64) {
    let new_x = (x + dx).clamp(x_min, x_max);
    (new_x, new_x - x)
}

#[inline(always)]
fn unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    mut beta: f64, mut x: f64, q: bool, n: u8,
) -> f64 {
    if (q && (x.is_sign_positive())) || (!q && (x.is_sign_negative())) {
        beta = (beta - normalised_intrinsic(x, q)).max(0.).abs();
    }
    if !q {
        x = -x;
    }
    if beta <= 0. || beta < DENORMALISATION_CUTOFF {
        return 0.0;
    }
    let b_max = (0.5 * x).exp();
    if beta >= b_max {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    let mut iterations = 0;
    let mut f = f64::MIN;
    let mut s;
    let mut ds = f64::MIN;
    let mut s_left = f64::MIN_POSITIVE;
    let mut s_right = f64::MAX;
    let s_c = SQRT_2 * x.abs().sqrt();
    let b_c = normalised_black_call(x, s_c);
    let v_c = normalised_vega(x, s_c);
    if beta < b_c {
        let s1 = s_c - b_c / v_c;
        let b1 = normalised_black_call(x, s1);
        if beta < b1 {
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) = compute_f_lower_map_and_first_two_derivatives(x, s1);
            let r2 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.0, b1, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
            f = rational_cubic_interpolation(beta, 0.0, b1, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, r2);
            if f <= 0.0 {
                let t = beta / b1;
                f = (f_lower_map_l * t + b1 * (1.0 - t)) * t;
            }
            s = inverse_f_lower_map(x, f);
            s_right = s1;
            let ln_beta = beta.ln();

            ds = 1.0_f64;
            while iterations < n && ds.abs() > f64::EPSILON * s {
                let (bx, ln_vega) = normalised_black_call_over_vega_and_ln_vega(x, s);
                let ln_b = bx.ln() + ln_vega;
                let bpob = bx.recip();
                let h = x / s;
                let b_h2 = (h.powi(2) / s) - s / 4.0;
                let nu = (ln_beta - ln_b) * ln_b / ln_beta * bx;
                let lambda = ln_b.recip();
                let otlambda = lambda.mul_add2(2.0, 1.0);
                let h2 = b_h2 - bpob * otlambda;
                let c = 3.0 * (h / s).powi(2);
                let b_h3 = b_h2.powi(2) - c - 0.25;
                let sq_bpob = bpob.powi(2);
                let mu = 6.0 * lambda * (1.0 + lambda);
                let h3 = b_h3 + sq_bpob * (2.0 + mu) - (b_h2 * bpob * 3.0 * otlambda);
                ds = if x < -190.0 {
                    nu * householder4_factor(nu, h2, h3, ((b_h2 * (b_h3 - 0.5)) - ((b_h2 - 2.0 / s) * 2.0 * c)) - (bpob * (sq_bpob * (6.0 + lambda * (22.0 + lambda * (36.0 + lambda * 24.0))) - (b_h2 * bpob * (12.0 + 6.0 * mu))) - (b_h2 * bpob * 3.0 * otlambda) - (b_h3 * bpob * 4.0 * otlambda)))
                } else {
                    nu * householder3_factor(nu, h2, h3)
                };
                (s, ds) = take_step(s_left, s_right, s, ds);
                iterations += 1;
            }
            return s;
        } else {
            let v1_recip = normalised_vega(x, s1).recip();
            let r_im = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b1, b_c, s1, s_c, v1_recip, v_c.recip(), 0.0, false);
            s = rational_cubic_interpolation(beta, b1, b_c, s1, s_c, v1_recip, v_c.recip(), r_im);
            s_left = s1;
            s_right = s_c;
        }
    } else {
        let s_u = if v_c.is_sign_positive() { s_c + (b_max - b_c) / v_c } else { s_c };
        let b_u = normalised_black_call(x, s_u);
        if beta <= b_u {
            let v_u = normalised_vega(x, s_u);
            let r_u_m = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
                b_c, b_u, s_c, s_u, v_c.recip(), v_u.recip(), 0.0, false);
            s = rational_cubic_interpolation(beta, b_c, b_u, s_c, s_u, v_c.recip(), v_u.recip(), r_u_m);
            s_left = s_c;
            s_right = s_u;
        } else {
            let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) = compute_f_upper_map_and_first_two_derivatives(x, s_u);

            if d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX {
                let r_uu = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_u, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, d2_f_upper_map_h_d_beta2, true);
                f = rational_cubic_interpolation(beta, b_u, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, r_uu);
            }
            if !f.is_sign_positive() {
                let h = b_max - b_u;
                let t = (beta - b_u) / h;
                f = (f_upper_map_h * (1.0 - t) + 0.5 * h * t) * (1.0 - t);
            }
            (s, s_left) = (inverse_f_upper_map(f), s_u);
            if beta > 0.5 * b_max {
                let beta_bar = b_max - beta;
                while iterations < n && ds.abs() > f64::EPSILON * s {
                    let h = x / s;
                    let t = s / 2.0;
                    let gp = SQRT_TWO_OVER_PI / (erfcx_cody((t + h) * FRAC_1_SQRT_2) + erfcx_cody((t - h) * FRAC_1_SQRT_2));
                    let b_bar = normalised_vega(x, s) / gp;
                    let g = (beta_bar / b_bar).ln();
                    let x_over_s_square = h.powi(2) / s;
                    let b_h2 = x_over_s_square - s / 4.0;
                    let c = 3.0 * (h / s).powi(2);
                    let b_h3 = b_h2 * b_h2 - c - 0.25;
                    let nu = -g / gp;
                    let h2 = b_h2 + gp;
                    let h3 = b_h3 + gp * (2.0 * gp + 3.0 * b_h2);
                    ds = if x < -580.0 {
                        nu * householder4_factor(nu, h2, h3, (b_h2 * (b_h3 - 0.5) - (b_h2 - 2.0 / s) * 2.0 * c) + gp * (6.0 * gp * b_h2.mul_add2(2.0,  gp) + 3.0 * b_h2 * b_h2 + 4.0 * b_h3))
                    } else {
                        nu * householder3_factor(nu, h2, h3)
                    };
                    (s, ds) = take_step(s_left, s_right, s, ds);
                    iterations += 1;
                }
                return s;
            }
        }
    }
    for _ in 0..n {
        if ds.abs() <= f64::EPSILON * s {
            break;
        }

        let b = normalised_black_call(x, s);
        let bp = normalised_vega(x, s);
        let nu = (beta - b) / bp;
        let h = x / s;
        let h2 = h.powi(2) / s - s / 4.0;
        let h3 = h2 * h2 - 3.0 * (h / s).powi(2) - 0.25;
        ds = nu * householder3_factor(nu, h2, h3);
        // Never leave the branch (or bracket)
        (s, ds) = take_step(s_left, s_right, s, ds);
    }
    s
}

#[inline(always)]
fn implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    mut price: f64,
    f: f64,
    k: f64,
    t: f64,
    mut q: bool,
    n: u8,
) -> f64 {
    let intrinsic = (if !q { k - f } else { f - k }).max(0.0).abs();
    if price < intrinsic {
        return
            VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
    }
    let max_price = if !q { k } else { f };
    if price >= max_price {
        return
            VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    let x = (f / k).ln();
    // Map in-the-money to out-of-the-money
    if (q && (x.is_sign_positive())) || (!q && (x.is_sign_negative())) {
        price = (price - intrinsic).max(0.0).abs();
        q = !q;
    }

    unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
        price / ((f * k).sqrt()),
        x,
        q,
        n,
    ) / t.sqrt()
}

#[inline(always)]
pub(crate) fn implied_black_volatility(price: f64, f: f64, k: f64, t: f64, q: bool) -> f64 {
    implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price, f, k, t, q, 2)
}


#[cfg(test)]
mod tests {
    use rand::Rng;
    use super::*;

    #[test]
    fn reconstruction_call_atm() {
        for i in 1..100 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let q = true;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            assert!((price - reprice).abs() < f64::EPSILON * 100.0);
        }
    }

    #[test]
    fn reconstruction_put_atm() {
        for i in 1..100 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let q = false;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            assert!((price - reprice).abs() < f64::EPSILON * 100.0);
        }
    }

    #[test]
    fn reconstruction_random_call_intrinsic() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1e5 * r2;
            let f = r + 1e5 * r2;
            let k = f - price;
            let t = 1e5 * r3;
            let q = true;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_itm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
            let q = true;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_otm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let q = true;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
        }
    }


    #[test]
    fn reconstruction_random_put_itm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
            let q = false;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            // if (price - reprice).abs() > 1.5 * f64::EPSILON{
            //     println!("{:?}", (price, f, k, t, q, sigma));
            //     println!("{:?}", (price - reprice).abs() / f64::EPSILON);
            // }
            assert!((price - reprice).abs() <= 1.75 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_put_otm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let q = false;
            let sigma = implied_black_volatility(price, f, k, t, q);
            let reprice = black(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
        }
    }
}
