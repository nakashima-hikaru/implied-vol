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
use std::ops::Neg;

const MAX_HOUSEHOLDER_CORRECTIVE_STEPS: usize = 2;

#[inline]
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

#[inline]
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

#[inline]
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

#[inline]
fn inverse_f_lower_map<SpFn: SpecialFn>(x: f64, f: f64) -> f64 {
    (x * FRAC_1_SQRT_3
        / SpFn::inverse_norm_cdf(FRAC_SQRT_3_CUBIC_ROOT_2_PI * f.cbrt() / x.abs().cbrt()))
    .abs()
}

#[inline]
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

#[inline]
fn inverse_f_upper_map<SpFn: SpecialFn>(f: f64) -> f64 {
    -2.0 * SpFn::inverse_norm_cdf(f)
}

#[inline]
pub(crate) fn implied_normalised_volatility_atm<SpFn: SpecialFn>(beta: f64) -> f64 {
    2.0 * SQRT_2 * SpFn::erfinv(beta)
}

#[inline]
#[cfg(test)]
fn lets_be_rational<SpFn: SpecialFn>(beta: f64, theta_x: f64) -> Option<f64> {
    debug_assert!(theta_x <= 0.0);
    debug_assert!(beta > 0.0);
    let b_max = (0.5 * theta_x).exp();
    if beta >= b_max {
        // time value exceeds the supremum of the model
        None
    } else {
        Some(hybrid_unchecked::<SpFn>(beta, theta_x, b_max))
    }
}

#[inline]
pub fn lets_be_rational_unchecked<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> f64 {
    lets_be_rational_with_dispatch::<SpFn, false>(beta, theta_x, b_max)
}

#[inline]
pub fn hybrid_unchecked<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> f64 {
    lets_be_rational_with_dispatch::<SpFn, true>(beta, theta_x, b_max)
}

#[inline]
fn lets_be_rational_with_dispatch<SpFn: SpecialFn, const HYBRID: bool>(
    beta: f64,
    theta_x: f64,
    b_max: f64,
) -> f64 {
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
            if HYBRID
                && let Some(s) = crate::flashiv::try_lowest_branch::<SpFn>(beta, theta_x, b_max)
            {
                return s;
            }
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
            for _ in 0..MAX_HOUSEHOLDER_CORRECTIVE_STEPS {
                debug_assert!(s > 0.0);
                debug_assert!(s.is_finite(), "s is not finite: s={s}");

                let h = theta_x / s;
                let (bx, ln_vega) = bs_option_price::scaled_normalised_black_and_ln_vega::<SpFn>(
                    0.5 * theta_x,
                    h,
                    0.5 * s,
                );

                let ln_b = bx.ln() + ln_vega;
                if ln_b == ln_beta {
                    return s;
                }
                let bpob = bx.recip();
                let x2_over_s3 = h * h / s;
                let b_h2 = s.mul_add2(-0.25, x2_over_s3);
                let nu = (ln_beta - ln_b) * ln_b / ln_beta / bpob;
                let lambda = ln_b.recip();
                let ot_lambda = lambda.mul_add2(2.0, 1.0);
                let h2 = bpob.mul_add2(-ot_lambda, b_h2);
                let c = 3.0 * (x2_over_s3 / s);
                let b_h3 = b_h2.mul_add2(b_h2, -c - 0.25);
                let sq_bpob = bpob * bpob;
                let bppob = b_h2 * bpob;
                let mu = 6.0 * lambda * (lambda + 1.0);
                let h3 = (bppob * 3.0).mul_add2(-ot_lambda, sq_bpob.mul_add2(2.0 + mu, b_h3));

                let ds = nu
                    * if theta_x < -190.0 {
                        householder::householder_4factor(
                            nu,
                            h2,
                            h3,
                            (b_h3 * bpob * 4.0).mul_add2(
                                -ot_lambda,
                                (bppob * b_h2 * 3.0).mul_add2(
                                    -ot_lambda,
                                    bpob.mul_add2(
                                        -sq_bpob.mul_add2(
                                            lambda.mul_add2(
                                                lambda.mul_add2(lambda.mul_add2(24.0, 36.0), 22.0),
                                                6.0,
                                            ),
                                            -(bppob * 6.0f64.mul_add2(mu, 12.0)),
                                        ),
                                        b_h2 * (b_h3 - 0.5) - (b_h2 - 2.0 / s) * 2.0 * c,
                                    ),
                                ),
                            ),
                        )
                    } else {
                        householder::householder_3factor(nu, h2, h3)
                    };

                s += ds;
                debug_assert!(s > 0.0);
            }
            return s;
        }

        // Lower middle: s_l <= s < s_c
        let inv_v = (
            bs_option_price::inv_normalised_vega(theta_x / s_l, 0.5 * s_l),
            SQRT_2_PI / b_max,
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
                for _ in 0..MAX_HOUSEHOLDER_CORRECTIVE_STEPS {
                    let h = theta_x / s;
                    let t = 0.5 * s;
                    let gp = SQRT_2_OVER_PI
                        / (SpFn::erfcx((t + h) * FRAC_1_SQRT_2)
                            + SpFn::erfcx((t - h) * FRAC_1_SQRT_2));
                    debug_assert!(s > 0.0);
                    let g = (beta_bar * gp).ln() + bs_option_price::ln_inv_normalised_vega(h, t);
                    if g == 0.0 {
                        return s;
                    }
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
                    s += ds;
                    debug_assert!(s > 0.0);
                }
                return s;
            }
        }
    }

    // MIDDLE BRANCHES (ITERATION)
    for _ in 0..MAX_HOUSEHOLDER_CORRECTIVE_STEPS {
        debug_assert!(s > 0.0);
        debug_assert!(theta_x < 0.0_f64);
        let h = theta_x / s;
        let t = 0.5 * s;
        let (b, inv_bp) =
            bs_option_price::normalised_black_and_inv_vega::<SpFn>(0.5 * theta_x, h, t);
        if b == beta {
            return s;
        }
        let nu = (beta - b) * inv_bp;
        let x2_over_s3 = h * h / s;
        let h2 = s.mul_add2(-0.25, x2_over_s3);
        let h3 = h2.mul_add2(h2, (-3.0f64).mul_add2(x2_over_s3 / s, -0.25));
        let ds = nu * householder::householder_3factor(nu, h2, h3);
        s += ds;
        debug_assert!(s > 0.0);
    }
    s
}

#[inline(always)]
#[cfg(test)]
pub fn implied_black_volatility_input_unchecked<SpFn: SpecialFn, const IS_CALL: bool>(
    price: f64,
    f: f64,
    k: f64,
    t: f64,
) -> Option<f64> {
    crate::solver::implied_black_volatility::<crate::solver::Hybrid<SpFn>, IS_CALL>(price, f, k, t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lets_be_rational::bs_option_price::{
        black_input_unchecked, scaled_normalised_black_and_ln_vega,
    };
    use crate::lets_be_rational::special_function::DefaultSpecialFn;
    use rand::RngExt;

    const FOURTH_ROOT_DBL_EPSILON: f64 = f64::from_bits(0x3f20_0000_0000_0000);

    fn normalised_intrinsic(theta_x: f64) -> f64 {
        // if theta_x <= 0.0 {
        //     return 0.0;
        // }
        let x2 = theta_x * theta_x;
        if x2 < 98.0 * FOURTH_ROOT_DBL_EPSILON {
            return x2
                .mul_add2(1.0 / 92_897_280.0, 1.0 / 322_560.0)
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
        const Q: bool = true;
        for i in 1..10000 {
            let price = 0.01 * f64::from(i);
            let f = 100.0;
            let k = f;
            let t = 1.0;
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
        const Q: bool = true;
        for i in 1..=10000 {
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let sigma = 0.001 * f64::from(i);
            let price = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            let sigma2 =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            debug_assert!(
                (sigma - sigma2).abs() / sigma
                    <= black_accuracy_factor(f.ln() - k.ln(), sigma * t.sqrt(), 1.0)
                        .recip()
                        .mul_add2(f64::EPSILON, 1.0),
                "f: {f}, k: {k}, t: {t}, sigma: {sigma}, sigma2; {sigma2}, price: {price}, {}, {}",
                (sigma - sigma2).abs() / sigma / f64::EPSILON,
                1.0 + black_accuracy_factor(f.ln() - k.ln(), sigma * t.sqrt(), 1.0).recip()
            );
        }
    }

    #[test]
    fn reconstruction_put_atm() {
        const Q: bool = false;
        for i in 1..100 {
            let price = 0.01 * f64::from(i);
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() < f64::EPSILON * 100.0);
        }
    }

    #[test]
    fn reconstruction_random_call_intrinsic() {
        const Q: bool = true;
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1e5 * r2;
            let f = 1e5f64.mul_add2(r2, r);
            // let f = 1e5f64.mul_add2(r2, r); // I will fix this in Part 2 if needed or here
            let k = f - price;
            let t = 1e5 * r3;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() <= f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_itm() {
        const Q: bool = true;
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0f64.mul_add2(1.0 - r, 1.0 * r * r2);
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
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
        const Q: bool = true;
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_put_itm() {
        const Q: bool = false;
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
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
        const Q: bool = false;
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0f64.mul_add2(1.0 - r, 1.0 * r * r2);
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON);
        }
    }

    #[test]
    fn panic_case() {
        const Q: bool = true;
        {
            let price = 73.425;
            let f = 12173.425;
            let k = 12100.0;
            let t = 0.007_707_632_775_934_893_4;
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
            let t = 0.007_705_811_088_032_645;
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
            let t = 0.007_705_808_219_781_035;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice);
        }
        {
            const Q: bool = true;
            let price = 73.425;
            let f = 12173.425;
            let k = 12100.0;
            let t = 0.007_705_804_818_688_366;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice);
        }
        {
            const Q: bool = true;
            let price = 33.55;
            let f = 11633.55;
            let k = 12100.0;
            let t = 0.007_705_800_716_005_495;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert!((price - reprice).abs() / price <= 6.0 * f64::EPSILON,);
        }
        {
            const Q: bool = true;
            let price = 33.55;
            let f = 11633.55;
            let t = 0.001_608_506_443_805_897_8;
            let k = 11600.0;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice, "f: {f}, k: {k}, t: {t}, sigma: {sigma}");
        }
        {
            const Q: bool = true;
            let price = 33.55;
            let f = 11633.55;
            let t = 0.001_608_506_443_805_897_8;
            let k = 11600.0;
            let sigma =
                implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
                    .unwrap();
            let reprice = black_input_unchecked::<DefaultSpecialFn, Q>(f, k, sigma, t);
            assert_eq!(price, reprice, "f: {f}, k: {k}, t: {t}, sigma: {sigma}");
        }
    }

    #[test]
    fn time_inf() {
        const Q: bool = true;
        let price = 20.0;
        let f = 100.0;
        let k = 100.0;
        let t = f64::INFINITY;
        let sigma = implied_black_volatility_input_unchecked::<DefaultSpecialFn, Q>(price, f, k, t)
            .unwrap();
        assert_eq!(sigma, 0.0);
    }

    #[test]
    fn flashiv_dispatch_boundary_preserves_accuracy() {
        // Use the actual fitted b_l rather than a reprice of its mathematical
        // tangent point, which can differ by multiple representable prices.
        for theta_x in [
            (-0.01_f64).next_up(),
            -0.01,
            (-0.01_f64).next_down(),
            -0.5,
            -3.0,
            -10.0,
            -191.0,
        ] {
            let b_max = (0.5 * theta_x).exp();
            let s_c = SQRT_2 * (-theta_x).sqrt();
            let b_l = b_l_over_b_max(s_c) * b_max;
            for beta in [b_l.next_down(), b_l, b_l.next_up()] {
                let s = hybrid_unchecked::<DefaultSpecialFn>(beta, theta_x, b_max);
                assert!(s.is_finite() && s > 0.0);
                let residual = residual_based_relative_error::<DefaultSpecialFn>(beta, theta_x, s);
                let limit = 4.0
                    * implied_volatility_attainable_accuracy_from_beta_theta_x(beta, theta_x, s);
                assert!(residual.is_finite() && limit.is_finite());
                assert!(
                    residual <= limit,
                    "x={theta_x:.17e}, beta={beta:.17e}, s={s:.17e}, residual={residual:.3e}, limit={limit:.3e}"
                );
            }
        }
    }

    #[test]
    fn normalised_implied_vol_roundtrip_accuracy() {
        // This is a floating-point round trip, not a bound on the exact inverse.
        // Both pricing and inversion round; the bundled two-step C++ reference
        // also exceeds one attainable-accuracy estimate. Four estimates cover
        // the observed round-trip error without adding any output refinement.
        const ROUNDTRIP_TOLERANCE_MULTIPLIER: f64 = 4.0;
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed([42; 32]);
        let mut checked = 0;

        for _ in 0..1_000_000 {
            let theta_x = -rng.random_range(1e-4..5.0);
            let s = rng.random_range(1e-3..5.0);
            let beta = bs_option_price::normalised_black::<DefaultSpecialFn>(
                0.5 * theta_x,
                theta_x / s,
                0.5 * s,
            );
            // Subnormal prices have a different rounding/conditioning limit.
            if beta < f64::MIN_POSITIVE {
                continue;
            }

            let implied_s = lets_be_rational::<DefaultSpecialFn>(beta, theta_x).unwrap();
            assert!(implied_s.is_finite() && implied_s > 0.0);
            let tolerance = ROUNDTRIP_TOLERANCE_MULTIPLIER
                * implied_volatility_attainable_accuracy_from_beta_theta_x(beta, theta_x, s);
            let input_error = (implied_s - s).abs() / s;
            let residual_error =
                residual_based_relative_error::<DefaultSpecialFn>(beta, theta_x, implied_s);
            assert!(
                input_error <= tolerance && residual_error <= tolerance,
                "theta_x={theta_x:.16e}, s={s:.16e}, beta={beta:.16e}, \
                 implied_s={implied_s:.16e}, input_error={input_error:.3e}, \
                 residual_error={residual_error:.3e}, tolerance={tolerance:.3e}"
            );
            checked += 1;
        }
        assert!(checked > 900_000);
    }

    fn implied_volatility_attainable_accuracy_from_beta_theta_x(
        beta: f64,
        theta_x: f64,
        s: f64,
    ) -> f64 {
        debug_assert!(theta_x <= 0.0);
        if s <= 0.0 {
            return f64::EPSILON;
        }
        if beta < f64::MIN_POSITIVE {
            return 1.0;
        }
        let h = theta_x / s;
        let t = 0.5 * s;
        let bx = beta * bs_option_price::inv_normalised_vega(h, t);
        f64::EPSILON * (1.0 + (bx / s).abs())
    }

    fn residual_based_relative_error<SpFn: SpecialFn>(beta: f64, theta_x: f64, s: f64) -> f64 {
        let h = theta_x / s;
        let t = 0.5 * s;
        let residual = bs_option_price::normalised_black::<SpFn>(0.5 * theta_x, h, t) - beta;
        residual.abs() * bs_option_price::inv_normalised_vega(h, t) / s
    }

    #[test]
    fn regression_lowest_branch_tail_case() {
        // Independent 100- and 200-digit roots for these exact binary64 inputs;
        // see tests/reference_black.rs for the pricing equation and provenance.
        let theta_x = f64::from_bits(0xbfdd_33e9_3b98_5c83);
        let beta = f64::from_bits(0x1414_c40c_5d29_351e);
        let reference_s = f64::from_bits(0x3f8e_5705_cb46_0f03);
        let root_correction = f64::from_bits(0xbc27_7d48_2038_8541);
        let implied_s = lets_be_rational::<DefaultSpecialFn>(beta, theta_x).unwrap();
        let attainable =
            implied_volatility_attainable_accuracy_from_beta_theta_x(beta, theta_x, reference_s);
        let root_error = ((implied_s - reference_s) - root_correction).abs() / reference_s;

        assert!(
            root_error <= attainable,
            "theta_x={theta_x:.16e}, beta={beta:.16e}, implied_s={implied_s:.16e}, reference_s={reference_s:.16e}, root_error={root_error:.3e}, attainable={attainable:.16e}"
        );
    }

    #[test]
    fn regression_lower_middle_case() {
        let theta_x = f64::from_bits(0xbfe4_9fdc_5054_7152);
        // Pricing this case produces different beta values with and without FMA.
        // Fix both inputs instead of bisecting floating-point pricing plateaus.
        for (beta_bits, root_bits, correction_bits) in [
            (
                0x3faf_1129_b89a_7174,
                0x3fe5_cf46_c34e_3fb5,
                0x3c81_7421_af65_a18b,
            ),
            (
                0x3faf_1129_b89a_716c,
                0x3fe5_cf46_c34e_3fb3,
                0x3c79_08c6_9e9e_6fc2,
            ),
        ] {
            let beta = f64::from_bits(beta_bits);
            let reference_s = f64::from_bits(root_bits);
            let root_correction = f64::from_bits(correction_bits);
            let implied_s = lets_be_rational::<DefaultSpecialFn>(beta, theta_x).unwrap();
            let attainable = implied_volatility_attainable_accuracy_from_beta_theta_x(
                beta,
                theta_x,
                reference_s,
            );
            let root_error = ((implied_s - reference_s) - root_correction).abs() / reference_s;
            assert!(
                root_error <= attainable,
                "theta_x={theta_x:.16e}, beta={beta:.16e}, implied_s={implied_s:.16e}, reference_s={reference_s:.16e}, root_error={root_error:.3e}, attainable={attainable:.16e}"
            );
        }
    }
}
