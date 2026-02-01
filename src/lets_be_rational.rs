pub mod bachelier_impl;
pub mod bs_option_price;
mod constants;
mod householder;
mod rational_cubic;
pub mod special_function;

use crate::fused_multiply_add::MulAdd;

use crate::lets_be_rational::constants::{
    FRAC_1_SQRT_3, FRAC_2_PI_SQRT_27, FRAC_SQRT_3_CUBIC_ROOT_2_PI, SQRT_2_OVER_PI, SQRT_2_PI,
    SQRT_3, SQRT_DBL_MAX, SQRT_PI_OVER_2,
};
use crate::lets_be_rational::rational_cubic::{
    convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side,
    convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side,
    rational_cubic_interpolation,
};
use crate::lets_be_rational::special_function::SpecialFn;
use crate::lets_be_rational::special_function::normal_distribution::inv_norm_pdf;
use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};
use std::ops::{Div, Neg};

#[inline(always)]
fn b_u_over_b_max(s_c: f64) -> f64 {
    if s_c >= 2.449_489_742_783_178 {
        let y = s_c.recip();

        let g = y
            .mul_add2(-1.229_189_712_271_654_4, 6.589_280_957_677_407E2)
            .mul_add2(y, 6.169_692_835_129_17E2)
            .mul_add2(y, 2.983_680_162_805_663E2)
            .mul_add2(y, 8.488_089_220_080_239E1)
            .mul_add2(y, 1.455_319_886_249_397_7E1)
            .mul_add2(y, 1.375_163_082_077_259_1)
            .mul_add2(y, -4.605_394_817_212_609E-2)
            / y.mul_add2(5.206_084_752_279_256E2, 8.881_238_333_960_678E2)
                .mul_add2(y, 8.698_830_313_690_185E2)
                .mul_add2(y, 5.079_647_179_123_228E2)
                .mul_add2(y, 2.030_420_459_952_177_3E2)
                .mul_add2(y, 5.436_378_146_588_073E1)
                .mul_add2(y, 9.327_034_903_790_405)
                .mul_add2(y, 1.0);
        y.mul_add2(g, -1.253_314_137_315_500_3)
            .mul_add2(0.113_984_531_941_499_06 * y, 0.894_954_297_278_031_3)
    } else {
        let g = s_c
            .mul_add2(-3.386_756_817_001_176_5E-9, -8.733_991_026_156_887E-4)
            .mul_add2(s_c, -8.143_812_878_548_491E-3)
            .mul_add2(s_c, -3.512_133_741_041_69E-2)
            .mul_add2(s_c, -8.976_383_086_137_545E-2)
            .mul_add2(s_c, -1.416_368_116_424_721E-1)
            .mul_add2(s_c, -1.344_864_378_589_371E-1)
            .mul_add2(s_c, -6.063_099_880_334_851E-2)
            / s_c
                .mul_add2(1.421_206_743_529_177_8E-2, 1.324_801_623_892_073E-1)
                .mul_add2(s_c, 5.959_161_649_351_221E-1)
                .mul_add2(s_c, 1.652_734_794_196_848_7)
                .mul_add2(s_c, 3.018_638_953_766_389_6)
                .mul_add2(s_c, 3.650_335_036_015_884_6)
                .mul_add2(s_c, 2.722_003_340_655_505_5)
                .mul_add2(s_c, 1.0);

        s_c.mul_add2(g, 0.061_461_680_580_514_74)
            .mul_add2(s_c * s_c, 0.789_908_594_556_062_8)
    }
}

#[inline(always)]
fn b_l_over_b_max(s_c: f64) -> f64 {
    if s_c < 0.709_929_573_971_953_9 {
        let g = s_c
            .mul_add2(4.542_510_209_361_606_4E-7, -6.403_639_934_147_98E-6)
            .mul_add2(s_c, 5.971_692_845_958_919E-3)
            .mul_add2(s_c, 3.976_063_144_567_705_5E-2)
            .mul_add2(s_c, 9.807_891_178_635_89E-2)
            .mul_add2(s_c, 8.074_107_237_288_286E-2)
            / s_c
                .mul_add2(6.125_459_704_983_172E-2, 4.613_270_710_865_565E-1)
                .mul_add2(s_c, 1.365_880_147_571_179)
                .mul_add2(s_c, 1.859_497_767_228_766_5)
                .mul_add2(s_c, 1.0);
        (s_c * s_c)
            * s_c.mul_add2(
                s_c.mul_add2(g, -0.096_727_192_813_394_37),
                0.075_609_966_402_963_62,
            )
    } else if s_c < 2.626_785_107_312_739_5 {
        s_c.mul_add2(6.971_140_063_983_471E-4, 6.584_925_270_230_231E-3)
            .mul_add2(s_c, 2.953_705_895_096_301_8E-2)
            .mul_add2(s_c, 6.917_130_174_466_835E-2)
            .mul_add2(s_c, 7.561_014_227_254_904E-2)
            .mul_add2(s_c, -2.708_128_856_468_558_7E-8)
            .mul_add2(s_c, 1.979_573_792_759_858E-9)
            / s_c
                .mul_add2(6.636_197_582_786_12E-3, 7.171_486_244_882_935E-2)
                .mul_add2(s_c, 3.783_162_225_306_046E-1)
                .mul_add2(s_c, 1.157_148_318_717_978_3)
                .mul_add2(s_c, 2.129_710_354_999_518)
                .mul_add2(s_c, 2.194_144_852_558_658)
                .mul_add2(s_c, 1.0)
    } else if s_c < 7.348_469_228_349_534 {
        s_c.mul_add2(1.701_257_940_724_605_5E-3, 1.002_291_337_825_409E-2)
            .mul_add2(s_c, 3.922_517_740_768_760_6E-2)
            .mul_add2(s_c, 7.403_965_818_682_282E-2)
            .mul_add2(s_c, 7.411_485_544_834_501E-2)
            .mul_add2(s_c, 5.311_803_397_279_465E-4)
            .mul_add2(s_c, -9.332_511_535_483_788E-5)
            / s_c
                .mul_add2(1.619_540_589_593_093_7E-2, 1.174_400_591_971_610_1E-1)
                .mul_add2(s_c, 5.323_125_844_350_184E-1)
                .mul_add2(s_c, 1.391_232_364_627_114)
                .mul_add2(s_c, 2.344_181_670_708_740_4)
                .mul_add2(s_c, 2.221_723_813_222_813_4)
                .mul_add2(s_c, 1.0)
    } else {
        s_c.mul_add2(1.693_020_807_842_147_5E-3, 5.183_252_617_163_152E-3)
            .mul_add2(s_c, 2.934_240_565_862_844_5E-2)
            .mul_add2(s_c, 3.921_610_857_820_463_6E-2)
            .mul_add2(s_c, 7.168_217_831_093_633E-2)
            .mul_add2(s_c, -1.511_669_248_501_119_6E-3)
            .mul_add2(s_c, 1.450_007_229_724_060_4E-3)
            / s_c
                .mul_add2(1.611_699_254_678_867_7E-2, 7.126_137_099_644_303E-2)
                .mul_add2(s_c, 3.754_374_213_737_579E-1)
                .mul_add2(s_c, 8.487_830_756_737_222E-1)
                .mul_add2(s_c, 1.682_315_917_528_153_2)
                .mul_add2(s_c, 1.617_631_350_230_541_5)
                .mul_add2(s_c, 1.0)
    }
}

#[inline(always)]
fn compute_f_lower_map_and_first_two_derivatives<SpFn: SpecialFn>(
    theta_x: f64,
    s: f64,
) -> (f64, f64, f64) {
    debug_assert!(theta_x < 0.0);
    let z = -FRAC_1_SQRT_3 * theta_x / s;
    let y = z * z;
    let s2 = s * s;
    let phi_m = 0.5 * SpFn::erfc(FRAC_1_SQRT_2 * z);

    let phi2 = phi_m * phi_m;
    (
        -FRAC_2_PI_SQRT_27 * theta_x * (phi2 * phi_m),
        std::f64::consts::TAU * y * phi2 * s2.mul_add2(0.125, y).exp(),
        std::f64::consts::FRAC_PI_6 * y / (s2 * s)
            * phi_m
            * (-8.0 * SQRT_3 * s).mul_add2(
                theta_x,
                (3.0 * s2).mul_add2(s2 - 8.0, -(8.0 * theta_x * theta_x)) * phi_m * inv_norm_pdf(y),
            )
            * 2.0f64.mul_add2(y, 0.25 * s2).exp(),
    )
}

#[inline(always)]
fn inverse_f_lower_map<SpFn: SpecialFn>(x: f64, f: f64) -> f64 {
    (x * FRAC_1_SQRT_3
        / SpFn::inverse_norm_cdf(FRAC_SQRT_3_CUBIC_ROOT_2_PI * f.cbrt() / x.abs().cbrt()))
    .abs()
}

#[inline(always)]
fn compute_f_upper_map_and_first_two_derivatives<SpFn: SpecialFn>(
    x: f64,
    s: f64,
) -> (f64, f64, f64) {
    let w = (x / s).powi(2);
    (
        0.5 * SpFn::erfc((0.5 * FRAC_1_SQRT_2) * s),
        -0.5 * (0.5 * w).exp(),
        SQRT_PI_OVER_2 * ((0.125 * s).mul_add2(s, w).exp()) * w / s,
    )
}

#[inline(always)]
fn inverse_f_upper_map<SpFn: SpecialFn>(f: f64) -> f64 {
    -2.0 * SpFn::inverse_norm_cdf(f)
}

#[inline(always)]
fn implied_normalised_volatility_atm<SpFn: SpecialFn>(beta: f64) -> f64 {
    2.0 * SQRT_2 * SpFn::erfinv(beta)
}

#[inline(always)]
fn lets_be_rational<SpFn: SpecialFn>(beta: f64, theta_x: f64) -> Option<f64> {
    debug_assert!(theta_x < 0.0);
    debug_assert!(beta > 0.0);
    let b_max = (0.5 * theta_x).exp();
    if beta >= b_max {
        // time value exceeds the supremum of the model
        None
    } else {
        lets_be_rational_unchecked::<SpFn>(beta, theta_x, b_max)
    }
}

#[inline(always)]
fn lets_be_rational_unchecked<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> Option<f64> {
    let mut s;
    let sqrt_ax = theta_x.neg().sqrt();
    let s_c = SQRT_2 * sqrt_ax;
    let ome = SpFn::one_minus_erfcx(sqrt_ax);
    let b_c = 0.5 * b_max * ome;

    // LOWER HALF: s < s_c
    if beta < b_c {
        debug_assert!(theta_x < 0.0);
        let s_l = (-SQRT_PI_OVER_2).mul_add2(ome, s_c);
        debug_assert!(s_l > 0.0);
        let b_l = b_l_over_b_max(s_c) * b_max;

        // LOWEST BRANCH: s < s_l
        if beta < b_l {
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) =
                compute_f_lower_map_and_first_two_derivatives::<SpFn>(theta_x, s_l);
            let r2 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side::<
                true,
            >(
                b_l,
                f_lower_map_l,
                (1.0, d_f_lower_map_l_d_beta),
                d2_f_lower_map_l_d_beta2,
            );
            let mut f = rational_cubic_interpolation(
                beta,
                b_l,
                (0.0, f_lower_map_l),
                (1.0, d_f_lower_map_l_d_beta),
                r2,
            );

            if !(f > 0.0) {
                let t = beta / b_l;
                f = f_lower_map_l.mul_add2(t, b_l.mul_add2(-t, b_l)) * t;
            }

            s = inverse_f_lower_map::<SpFn>(theta_x, f);
            debug_assert!(s > 0.0);

            let ln_beta = beta.ln();

            let mut s_left = f64::MIN_POSITIVE;
            let mut s_right = s_l;
            let mut direction_reversal_count = 0;
            let mut ds_previous: f64 = 0.0;
            let mut ds = 0.0;

            // Iteration limit strictly 2
            for i in 0..2 {
                // Logic to switch to binary nesting if Householder misbehaves
                if i > 0 {
                    if ds * ds_previous < 0.0 {
                        direction_reversal_count += 1;
                    }
                    if direction_reversal_count == 3 || !(s > s_left && s < s_right) {
                        s = 0.5 * (s_left + s_right);
                        if (s_right - s_left) <= f64::EPSILON * s {
                            return Some(s);
                        }
                        direction_reversal_count = 0;
                        ds_previous = 0.0;
                    } else {
                        ds_previous = ds;
                    }
                }
                debug_assert!(s > 0.0);
                debug_assert!(s.is_finite(), "s is not finite: s={}", s);

                let h = theta_x / s;
                let t = 0.5 * s;
                let (bx, ln_vega) = bs_option_price::scaled_normalised_black_and_ln_vega::<SpFn>(
                    0.5 * theta_x,
                    h,
                    t,
                );

                let ln_b = bx.ln() + ln_vega;
                let b = ln_b.exp();
                let bpob = bx.recip();
                let bp = bpob * b;

                // Update brackets
                if b > beta && s < s_right {
                    s_right = s;
                } else if b < beta && s > s_left {
                    s_left = s;
                }

                // Check for numerical underflow/issues leading to valid b, bp
                if !(b > 0.0 && bp > 0.0) {
                    ds = 0.5 * (s_left + s_right) - s;
                } else {
                    let x2_over_s3 = h * h / s;
                    let b_h2 = t.mul_add2(-0.5, x2_over_s3);
                    let v = (ln_beta - ln_b) * ln_b / ln_beta * bx;
                    let lambda = ln_b.recip();

                    #[cfg(feature = "fma")]
                    let ot_lambda = lambda.mul_add2(2.0, 1.0);
                    #[cfg(not(feature = "fma"))]
                    let ot_lambda = lambda + lambda + 1.0;

                    let h2 = ot_lambda.mul_add2(-bpob, b_h2);
                    let c = 3.0 * (x2_over_s3 / s);
                    let b_h3 = b_h2.mul_add2(b_h2, -c) - 0.25;
                    let sq_bpob = bpob * bpob;
                    let bppob_triple = 3.0 * b_h2 * bpob;
                    let mu_plus_2 = (1.0 + lambda).mul_add2(6.0 * lambda, 2.0);
                    let h3 = bppob_triple.mul_add2(-ot_lambda, sq_bpob.mul_add2(mu_plus_2, b_h3));

                    ds = v * if theta_x < -190.0 {
                        householder::householder_4factor(
                            v,
                            h2,
                            h3,
                            b_h2.mul_add2(b_h3 - 0.5, -((b_h2 - 2.0 / s) * 2.0 * c))
                                - (b_h3 * bpob * 4.0).mul_add2(
                                    -ot_lambda,
                                    bpob.mul_add2(
                                        sq_bpob.mul_add2(
                                            lambda
                                                .mul_add2(24.0, 36.0)
                                                .mul_add2(lambda, 22.0)
                                                .mul_add2(lambda, 6.0),
                                            -(2.0 * bppob_triple * mu_plus_2),
                                        ),
                                        -(bppob_triple * ot_lambda),
                                    ),
                                ),
                        )
                    } else {
                        householder::householder_3factor(v, h2, h3)
                    };
                }

                if ds.abs() <= f64::EPSILON * s {
                    return Some(s);
                }

                s += ds;
            }
            return Some(s);
        }

        // Lower middle: s_l <= s < s_c
        let inv_v = (
            SQRT_2_PI / b_max,
            bs_option_price::inv_normalised_vega(theta_x / s_l, 0.5 * s_l),
        );
        let h = b_c - b_l;
        let r_im = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side::<
            false,
        >(h, s_c - s_l, inv_v, 0.0);
        s = rational_cubic_interpolation(beta - b_l, h, (s_l, s_c), inv_v, r_im);
        debug_assert!(s > 0.0);
    } else {
        // UPPER HALF: s_c <= s
        let s_u = SQRT_PI_OVER_2.mul_add2(2.0 - ome, s_c);
        debug_assert!(s_u > 0.0);
        let b_u = b_u_over_b_max(s_c) * b_max;

        if beta <= b_u {
            // UPPER MIDDLE: s_c <= s <= s_u
            let inv_v = (
                SQRT_2_PI / b_max,
                bs_option_price::inv_normalised_vega(theta_x / s_u, 0.5 * s_u),
            );
            let h = b_u - b_c;
            let r_u_m =
                convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side::<
                    false,
                >(h, s_u - s_c, inv_v, 0.0);
            s = rational_cubic_interpolation(beta - b_c, h, (s_c, s_u), inv_v, r_u_m);
            debug_assert!(s > 0.0);
        } else {
            // HIGHEST BRANCH: s_u < s
            let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) =
                compute_f_upper_map_and_first_two_derivatives::<SpFn>(theta_x, s_u);
            let mut f = if (-SQRT_DBL_MAX..SQRT_DBL_MAX).contains(&d2_f_upper_map_h_d_beta2) {
                let h = b_max - b_u;
                let r_uu =
                    convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side::<
                        true,
                    >(
                        h,
                        -f_upper_map_h,
                        (d_f_upper_map_h_d_beta, -0.5),
                        d2_f_upper_map_h_d_beta2,
                    );
                rational_cubic_interpolation(
                    beta - b_u,
                    h,
                    (f_upper_map_h, 0.0),
                    (d_f_upper_map_h_d_beta, -0.5),
                    r_uu,
                )
            } else {
                f64::MIN
            };
            if f <= 0.0 {
                let h = b_max - b_u;
                let t = (beta - b_u) / h;
                f = f_upper_map_h.mul_add2(1.0 - t, 0.5 * h * t) * (1.0 - t);
            }
            s = inverse_f_upper_map::<SpFn>(f);
            if beta > 0.5 * b_max {
                let beta_bar = b_max - beta;
                for _ in 0..2 {
                    let h = theta_x / s;
                    let t = 0.5 * s;
                    let gp = SQRT_2_OVER_PI
                        / (SpFn::erfcx((t + h) * FRAC_1_SQRT_2)
                            + SpFn::erfcx((t - h) * FRAC_1_SQRT_2));
                    debug_assert!(s > 0.0);
                    let g = (beta_bar * gp).ln() + bs_option_price::ln_inv_normalised_vega(h, t);
                    let x2_over_s3 = h * h / s;
                    let b_h2 = t.mul_add2(-0.5, x2_over_s3);
                    let c = 3.0 * x2_over_s3 / s;
                    let b_h3 = b_h2.mul_add2(b_h2, -c - 0.25);
                    let v = -g / gp;
                    let h2 = b_h2 + gp;
                    let h3 = gp.mul_add2(2.0f64.mul_add2(gp, 3.0 * b_h2), b_h3);
                    let ds = v * if theta_x < -580.0 {
                        householder::householder_4factor(
                            v,
                            h2,
                            h3,
                            gp.mul_add2(
                                4.0f64.mul_add2(
                                    b_h3,
                                    (6.0 * gp).mul_add2(b_h2.mul_add2(2.0, gp), 3.0 * b_h2 * b_h2),
                                ),
                                b_h2.mul_add2(b_h3 - 0.5, -((b_h2 - 2.0 / s) * 2.0 * c)),
                            ),
                        )
                    } else {
                        householder::householder_3factor(v, h2, h3)
                    };

                    if ds.abs() <= f64::EPSILON * s {
                        return Some(s);
                    }

                    s += ds;
                    debug_assert!(s > 0.0);
                }
                return Some(s);
            }
        }
    }

    // MIDDLE BRANCHES (ITERATION)
    const MAX_ITERATIONS: usize = 2;
    for _ in 0..MAX_ITERATIONS {
        debug_assert!(s > 0.0);
        debug_assert!(theta_x < 0.0_f64);
        let h = theta_x / s;
        let t = 0.5 * s;
        let b = bs_option_price::normalised_black::<SpFn>(0.5 * theta_x, h, t);
        let bp = bs_option_price::normalised_vega(h, t);
        let nu = (beta - b) / bp;
        let h2 = t.mul_add2(-0.5, h * h / s);
        let h3 = h2.mul_add2(h2, -(3.0 * (h / s).powi(2))) - 0.25_f64;
        let ds = nu * householder::householder_3factor(nu, h2, h3);

        if ds.abs() <= f64::EPSILON * s {
            return Some(s);
        }

        s += ds;
        debug_assert!(s > 0.0);
    }
    Some(s)
}

#[inline(always)]
pub fn implied_black_volatility_input_unchecked<SpFn: SpecialFn, const IS_CALL: bool>(
    price: f64,
    f: f64,
    k: f64,
    t: f64,
) -> Option<f64> {
    if price >= if IS_CALL { f } else { k } {
        return (price == if IS_CALL { f } else { k }).then_some(f64::INFINITY);
    }
    let intrinsic_value = if IS_CALL { f - k } else { k - f };
    let normalised_time_value = if intrinsic_value > 0.0_f64 {
        price - intrinsic_value
    } else {
        price
    } / (f.sqrt() * k.sqrt());
    if normalised_time_value <= f64::MIN_POSITIVE {
        return (normalised_time_value >= 0.0).then_some(0.0);
    }
    Some(if f == k {
        implied_normalised_volatility_atm::<SpFn>(normalised_time_value) / t.sqrt()
    } else {
        lets_be_rational::<SpFn>(normalised_time_value, (f.ln() - k.ln()).abs().neg())?
            .div(t.sqrt())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lets_be_rational::bs_option_price::{
        black_input_unchecked, scaled_normalised_black_and_ln_vega,
    };
    use crate::lets_be_rational::special_function::DefaultSpecialFn;
    use rand::Rng;

    pub(crate) const FOURTH_ROOT_DBL_EPSILON: f64 = f64::from_bits(0x3f20_0000_0000_0000);

    fn normalised_intrinsic(theta_x: f64) -> f64 {
        // if theta_x <= 0.0 {
        //     return 0.0;
        // }
        let x2 = theta_x * theta_x;
        if x2 < 98.0 * FOURTH_ROOT_DBL_EPSILON {
            return x2
                .mul_add2(1.0 / 92897280.0, 1.0 / 322560.0)
                .mul_add2(x2, 1.0 / 1920.0)
                .mul_add2(x2, 1.0 / 120.0)
                .mul_add2(x2, 1.0 / 24.0)
                .mul_add2(x2, 1.0)
                * theta_x;
        }
        (0.5 * theta_x).exp() - (-0.5 * theta_x).exp()
    }
    fn scaled_normalised_black(theta_x: f64, s: f64) -> f64 {
        debug_assert!(s > 0.0 && theta_x != 0.0);
        let h = theta_x / s;
        let t = 0.5 * s;
        (if theta_x > 0.0 {
            normalised_intrinsic(theta_x) * SQRT_2_PI * (0.5 * t.mul_add2(t, h * h)).exp()
        } else {
            0.0
        }) + scaled_normalised_black_and_ln_vega::<DefaultSpecialFn>(0.5 * -theta_x.abs(), h, t).0
    }

    #[allow(unused)]
    fn black_accuracy_factor(x: f64, s: f64, theta: f64 /* θ=±1 */) -> f64 {
        // When x = 0, we have bx(x,s) = b(x,s) / (∂(b(x,s)/∂s)  =  s·(1+s²/12+s⁴/240+O(s⁶)) for small s.
        if x == 0.0 {
            return if s.abs() < f64::EPSILON {
                1.0
            } else {
                s / (DefaultSpecialFn::erf((0.5 * FRAC_1_SQRT_2) * s)
                    * SQRT_2_PI
                    * (0.125 * s * s).exp())
            };
        }
        let theta_x = if theta < 0.0 { -x } else { x };
        if s <= 0.0 {
            return if theta_x > 0.0 { 0.0 } else { f64::MAX };
        }
        s / scaled_normalised_black(theta_x, s)
    }

    #[test]
    fn reconstruction_call_atm() {
        for i in 1..10000 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            debug_assert!(
                (price - reprice).abs() / price < 4.0 * f64::EPSILON,
                "{f},{k},{t},{sigma},{price},{reprice},{}",
                (price - reprice).abs() / price / f64::EPSILON
            );
        }
    }

    #[test]
    fn reconstruction_call_atm2() {
        for i in 1..=10000 {
            let f = 100.0;
            let k = f;
            let t = 1.0;
            const Q: bool = true;
            let sigma = 0.001 * i as f64;
            let price = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            let sigma2 =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            debug_assert!(
                (sigma - sigma2).abs() / sigma
                    <= 1.0
                        + black_accuracy_factor(f.ln() - k.ln(), sigma * t.sqrt(), 1.0).recip()
                            * f64::EPSILON,
                "f: {f}, k: {k}, t: {t}, sigma: {sigma}, sigma2; {sigma2}, price: {price}, {}, {}",
                (sigma - sigma2).abs() / sigma / f64::EPSILON,
                1.0 + black_accuracy_factor(f.ln() - k.ln(), sigma * t.sqrt(), 1.0).recip()
            );
        }
    }

    #[test]
    fn reconstruction_put_atm() {
        for i in 1..100 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            const Q: bool = false;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
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
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
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
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!(
                (price - reprice).abs() <= 1.5 * f64::EPSILON,
                "{f},{k},{t},{sigma},{price},{reprice},{}",
                (price - reprice).abs() / f64::EPSILON
            );
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
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
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
            const Q: bool = false;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            // if (price - reprice).abs() > 1.5 * f64::EPSILON{
            //     println!("{:?}", (price, f, k, t, q, sigma));
            //     println!("{:?}", (price - reprice).abs() / f64::EPSILON);
            // }
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
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
            const Q: bool = false;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
        }
    }

    #[test]
    fn panic_case() {
        {
            let price = 73.425;
            let f = 12173.425;
            let k = 12100.0;
            let t = 0.0077076327759348934;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice);
        }
        {
            let price = 73.425;
            let f = 12173.425;
            let k = 12100.0;
            let t = 0.007705811088032645;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice);
        }
        {
            let price = 73.425;
            let f = 12173.425;
            let k = 12100.0;
            let t = 0.007705808219781035;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice);
        }
        {
            let price = 73.425;
            let f = 12173.425;
            let k = 12100.0;
            let t = 0.007705804818688366;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice);
        }
        {
            let price = 33.55;
            let f = 11633.55;
            let k = 12100.0;
            let t = 0.007705800716005495;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() / price <= 6.0 * f64::EPSILON,);
        }
        {
            let price = 33.55;
            let f = 11633.55;
            let t = 0.0016085064438058978;
            let k = 11600.0;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice, "f: {f}, k: {k}, t: {t}, sigma: {sigma}");
        }
        {
            let price = 33.55;
            let f = 11633.55;
            let t = 0.0016085064438058978;
            let k = 11600.0;
            const Q: bool = true;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice, "f: {f}, k: {k}, t: {t}, sigma: {sigma}");
        }
    }

    #[test]
    fn time_inf() {
        let price = 20.0;
        let f = 100.0;
        let k = 100.0;
        let t = f64::INFINITY;
        const Q: bool = true;
        let sigma = implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
            .unwrap();
        assert_eq!(sigma, 0.0);
    }

    #[test]
    fn check_relative_error_of_implied_vol() {
        let mut max_relative_error = 0.0;

        let mut error_count = 0;
        let mut panic_count = 0;
        let n = 1_000_000;

        let seed: [u8; 32] = [42; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

        for _ in 0..n {
            // theta_x < 0
            let theta_x = -rng.random_range(1e-4..5.0);
            // s > 0
            let s = rng.random_range(1e-3..5.0);

            // Calculate beta (normalised option price)
            let h = theta_x / s;
            let t = 0.5 * s;
            let beta = bs_option_price::normalised_black::<DefaultSpecialFn>(0.5 * theta_x, h, t);

            if beta <= 0.0 {
                continue;
            }

            // Recalculate s
            let result =
                std::panic::catch_unwind(|| lets_be_rational::<DefaultSpecialFn>(beta, theta_x));

            match result {
                Ok(Some(_s_recalc)) => {
                    // Calculate beta (normalised option price)
                    let beta = bs_option_price::normalised_black::<DefaultSpecialFn>(
                        0.5 * theta_x,
                        theta_x / s,
                        0.5 * s,
                    );

                    if beta < f64::MIN_POSITIVE {
                        continue;
                    }

                    // Calculate implied volatility
                    let b_max = (0.5 * theta_x).exp();
                    let implied_s =
                        lets_be_rational_unchecked::<DefaultSpecialFn>(beta, theta_x, b_max);

                    // lets_be_rational_unchecked returns Option<f64>
                    let implied_s = match implied_s {
                        Some(val) => val,
                        None => {
                            // This should technically not happen for valid beta < b_max
                            continue;
                        }
                    };
                    // Calculate relative error
                    let relative_error = (implied_s - s).abs() / s;

                    if relative_error > max_relative_error {
                        max_relative_error = relative_error;
                    }

                    // Allow for attainable accuracy limit
                    // Note: The C++ reference implementation states that the attainable accuracy is
                    // (|b(s)/(s*b'(s))| + 1) * EPSILON.
                    // We allow a small tolerance factor (e.g. 1.05) on top of this theoretical limit.
                    let attainable_accuracy = implied_volatility_attainable_accuracy::<
                        DefaultSpecialFn,
                    >(
                        theta_x, s, -1.0, /* q=-1 because theta_x < 0 means out-of-the-money */
                    );

                    // Using a slightly relaxed tolerance because the theoretical limit is "to lowest order"
                    // We accept 2.0x attainable accuracy to account for minor differences in math library precision (e.g. erfc)
                    // relative to the C++ reference implementation, while strictly rejecting catastrophic errors.
                    #[cfg(feature = "fma")]
                    let tolerance = attainable_accuracy * 2.0;
                    #[cfg(not(feature = "fma"))]
                    let tolerance = attainable_accuracy * 3.0;

                    if relative_error > tolerance {
                        error_count += 1;
                        // Only print the first few errors to avoid flooding output
                        if error_count <= 10 {
                            println!(
                                "High Relative Error: {:.3e}, Attainable: {:.3e} (Ratio: {:.2}), Theta_x: {:.3e}, s: {:.3e}, Beta: {:.3e}",
                                relative_error,
                                attainable_accuracy,
                                relative_error / attainable_accuracy,
                                theta_x,
                                s,
                                beta
                            );
                        }
                    }
                }
                Ok(None) => {}
                Err(_) => {
                    panic_count += 1;
                    if panic_count < 10 {
                        println!("Panic for beta: {}, theta_x: {}", beta, theta_x);
                    }
                }
            }
        }

        println!("Max Relative Error: {:.3e}", max_relative_error);
        println!("Error Count (exceeding attainable): {}", error_count);
        println!("Panic Count: {}", panic_count);

        assert_eq!(
            panic_count, 0,
            "Panics occurred during implied vol calculation"
        );
        assert_eq!(
            error_count, 0,
            "Relative error exceeded attainable accuracy limit"
        );
    }

    // Helper to calculate theoretical attainable accuracy limit
    // Replicates `ImpliedVolatilityAttainableAccuracy` from C++ source
    fn implied_volatility_attainable_accuracy<SpFn: SpecialFn>(x: f64, s: f64, q: f64) -> f64 {
        if x == 0.0 {
            return f64::EPSILON
                * (1.0
                    + (if s.abs() <= f64::EPSILON {
                        1.0
                    } else {
                        (SpFn::erf(0.5 * std::f64::consts::FRAC_1_SQRT_2 * s)
                            * std::f64::consts::SQRT_2.sqrt() * std::f64::consts::PI.sqrt() // SQRT_2_PI approximation or use crate constant
                            * (0.125 * s * s).exp())
                        .abs()
                            / s.abs()
                    })
                    .abs());
        }

        let theta_x = if q < 0.0 { -x } else { x };

        if s <= 0.0 {
            return if theta_x > 0.0 { 1.0 } else { f64::EPSILON };
        }

        let half_theta_x_neg_abs = -0.5 * x.abs();
        let h_neg_abs = -x.abs() / s;
        let t = 0.5 * s;

        let (scaled_b_part, _ln_vega) =
            crate::lets_be_rational::bs_option_price::scaled_normalised_black_and_ln_vega::<SpFn>(
                half_theta_x_neg_abs,
                h_neg_abs,
                t,
            );

        let intrinsic = if theta_x > 0.0 {
            2.0 * (0.5 * theta_x).sinh()
        } else {
            0.0
        };

        let bx_val = if theta_x > 0.0 {
            intrinsic
                * crate::lets_be_rational::constants::SQRT_2_PI
                * (0.5 * ((theta_x / s).powi(2) + 0.25 * s * s)).exp()
        } else {
            0.0
        } + scaled_b_part;

        let vega = crate::lets_be_rational::bs_option_price::normalised_vega(x / s, 0.5 * s);

        if bx_val.abs() * vega >= f64::MIN_POSITIVE {
            f64::EPSILON * (1.0 + (bx_val / s).abs())
        } else {
            1.0
        }
    }
}
