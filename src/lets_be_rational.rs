use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};
use std::ops::Neg;
use crate::constants::{DENORMALISATION_CUTOFF, FOURTH_ROOT_DBL_EPSILON, FRAC_SQRT_TWO_PI, HALF_OF_LN_TWO_PI, ONE_OVER_SQRT_THREE, SIXTEENTH_ROOT_DBL_EPSILON, SQRT_DBL_MAX, SQRT_MIN_POSITIVE, SQRT_PI_OVER_TWO, SQRT_THREE, SQRT_THREE_OVER_THIRD_ROOT_TWO_PI, SQRT_TWO_OVER_PI, SQRT_TWO_PI, TWO_PI_OVER_SQRT_TWENTY_SEVEN, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC};
use crate::erf_cody::{erf_cody, erfc_cody, erfcx_cody};
use crate::MulAdd;
use crate::normal_distribution::{inverse_norm_cdf, norm_cdf, norm_pdf};
use crate::rational_cubic::{convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side, convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side, rational_cubic_interpolation};


#[inline(always)]
fn householder3_factor(v: f64, h2: f64, h3: f64) -> f64 { (1.0 + 0.5 * h2 * v) / (1.0 + v * (h2 + h3 * v / 6.0)) }

#[inline(always)]
fn householder4_factor(v: f64, h2: f64, h3: f64, h4: f64) -> f64 { (1.0 + v * (h2 + v * h3 / 6.0)) / (1.0 + v * (1.5 * h2 + v * (h2 * h2 / 4.0 + h3 / 3.0 + v * h4 / 24.0))) }

#[inline(always)]
fn normalised_intrinsic(x: f64) -> f64 {
    if !x.is_sign_positive() {
        return 0.0;
    }
    let x2 = x * x;
    if x2 < 98.0 * FOURTH_ROOT_DBL_EPSILON {
        let ret = x * (1.0 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))).abs().max(0.0);
        return ret
    }
    let b_max = (0.5 * x).exp();
    let one_over_b_max = b_max.recip();
    (b_max - one_over_b_max).abs().max(0.0)
}

const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
const SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

#[inline(always)]
fn asymptotic_expansion_of_normalised_black_call_over_vega(h: f64, t: f64) -> f64 {
    const fn a0(_e: f64) -> f64 {2.0}
    const fn a1(e: f64) -> f64 {-6.0-2.0*e}
    const fn a2(e: f64) -> f64 {30.0+e*(60.0+6.0*e)}
    const fn a3(e: f64) -> f64 { -2.1E2+e*(-1.05E3+e*(-6.3E2-30.0*e))}
    const fn a4(e: f64) -> f64 { 1.89E3+e*(1.764E4+e*(2.646E4+e*(7.56E3+2.1E2*e)))}
    const fn a5(e: f64) -> f64 {-2.079E4+e*(-3.1185E5+e*(-8.7318E5+e*(-6.237E5+e*(-1.0395E5-1.89E3*e))))}
    const fn a6(e: f64) -> f64 {2.7027E5+e*(5.94594E6+e*(2.675673E7+e*(3.567564E7+e*(1.486485E7+e*(1.62162E6+2.079E4*e)))))}
    const fn a7(e: f64) -> f64 {-4.05405E6+e*(-1.2297285E8+e*(-8.1162081E8+e*(-1.73918745E9+e*(-1.35270135E9+e*(-3.6891855E8+e*(-2.837835E7-2.7027E5*e))))))}
    const fn a8(e: f64) -> f64 {6.891885E7+e*(2.756754E9+e*(2.50864614E10+e*(7.88431644E10+e*(9.85539555E10+e*(5.01729228E10+e*(9.648639E9+e*(5.513508E8+4.05405E6*e)))))))}
    const fn a9(e: f64) -> f64 {-1.30945815E9+e*(-6.678236565E10+e*(-8.013883878E11+e*(-3.4726830138E12+e*(-6.3665855253E12+e*(-5.2090245207E12+e*(-1.8699062382E12+e*(-2.671294626E11+e*(-1.178512335E10-6.891885E7*e))))))))}
    const fn a10(e: f64) -> f64 {2.749862115E10+e*(1.7415793395E12+e*(2.664616389435E13+e*(1.52263793682E14+e*(3.848890340295E14+e*(4.618668408354E14+e*(2.664616389435E14+e*(7.10564370516E13+e*(7.83710702775E12+e*(2.749862115E11+1.30945815E9*e)))))))))}
    const fn a11(e: f64) -> f64 {-6.3246828645E11+e*(-4.870005805665E13+e*(-9.2530110307635E14+e*(-6.74147946527055E15+e*(-2.24715982175685E16+e*(-3.71802806872497E16+e*(-3.14602375045959E16+e*(-1.34829589305411E16+e*(-2.77590330922905E15+e*(-2.4350029028325E14+e*(-6.95715115095E12-2.749862115E10*e))))))))))}
    const fn a12(e: f64) -> f64 {1.581170716125E13+e*(1.454677058835E15+e*(3.36030400590885E16+e*(3.04027505296515E17+e*(1.292_116_897_510_188_8E18+e*(2.81916414002223E18+e*(3.289024830025935E18+e*(2.067387036016302E18+e*(6.840_618_869_171_588E17+e*(1.12010133530295E17+e*(8.0007238235925E15+e*(1.89740485935E14+6.3246828645E11*e)))))))))))}
    const fn a13(e: f64) -> f64 {-4.2691609335375E14+e*(-4.624924344665625E16+e*(-1.2764791191277125E18+e*(-1.404_127_031_040_483_7E19+e*(-7.410_670_441_602_553E19+e*(-2.061_513_777_391_255_6E20+e*(-3.171_559_657_525_009E20+e*(-2.748_685_036_521_674_2E20+e*(-1.333_920_679_488_459_6E20+e*(-3.510_317_577_601_209_5E19+e*(-4.680_423_436_801_613E18+e*(-2.774954606799375E17+e*(-5.54990921359875E15-1.581170716125E13*e))))))))))))}
    const fn a14(e: f64) -> f64 {1.238056670725875E16+e*(1.5599514051146025E18+e*(5.069_842_066_622_458E19+e*(6.663_221_001_846_66E20+e*(4.275_566_809_518_273E21+e*(1.477_013_988_742_676E22+e*(2.897_219_747_149_095_4E22+e*(3.311_108_282_456_109E22+e*(2.215_520_983_114_014E22+e*(8.551_133_619_036_546E21+e*(1.832_385_775_507_831_3E21+e*(2.027_936_826_648_983_3E20+e*(1.013_968_413_324_491_6E19+e*(1.733279339016225E17+4.2691609335375E14*e)))))))))))))}
    const fn a15(e: f64) -> f64 {-3.8379756792502125E17+e*(-5.565_064_734_912_808E19+e*(-2.103_594_469_797_041_5E21+e*(-3.255_562_869_923_992_8E22+e*(-2.495_931_533_608_394_5E23+e*(-1.048_291_244_115_525_7E24+e*(-2.553_529_953_614_742E24+e*(-3.720_857_932_410_053E24+e*(-3.283_109_940_361_811E24+e*(-1.747_152_073_525_876E24+e*(-5.491_049_373_938_467_5E23+e*(-9.766_688_609_771_979E22+e*(-9.115_576_035_787_18E21+e*(-3.895_545_314_438_965_5E20+e*(-5.756_963_518_875_318E18-1.238056670725875E16*e))))))))))))))}
    const fn a16(e: f64) -> f64 {1.266_531_974_152_57E19+e*(2.093_999_530_598_916E21+e*(9.108_897_958_105_285E22+e*(1.639_601_632_458_951_2E24+e*(1.480_195_918_192_108_6E25+e*(7.427_892_244_018_582E25+e*(2.199_798_856_882_426E26+e*(3.980_588_407_692_009_4E26+e*(4.478_161_958_653_511E26+e*(3.142_569_795_546_323E26+e*(1.361_780_244_736_74E26+e*(3.552_470_203_661_060_7E25+e*(5.328_705_305_491_592E24+e*(4.250_819_047_115_799E23+e*(1.570_499_647_949_187E22+e*(2.026_451_158_644_112E20+3.8379756792502125E17*e)))))))))))))))}
    const THRESHOLDS: [f64; 12] = [12.347, 12.958, 13.729, 14.718, 16.016, 17.769, 20.221, 23.816, 29.419, 38.93, 57.171, 99.347];

    assert!((h < -ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD.abs()) && (h + t < -(SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD).abs()));
    let e = (t / h).powi(2);
    let r = (h + t) * (h - t);
    let q = (h / r).powi(2);

    let idx = THRESHOLDS.partition_point(|&x| x<=-h - t + SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD);
    let omega = if idx == 12 {
        return a0(e) + q * (a1(e) + q * (a2(e) + q * (a3(e) + q * a4(e))));
    } else {
        let mut omega = 0.0;
        if idx == 0 {
            omega = q * (a16(e) + omega);
        }
        if idx <= 1 {
            omega = q * (a15(e) + omega);
        }
        if idx <= 2 {
            omega = q * (a14(e) + omega);
        }
        if idx <= 3 {
            omega = q * (a13(e) + omega);
        }
        if idx <= 4 {
            omega = q * (a12(e) + omega);
        }
        if idx <= 5 {
            omega = q * (a11(e) + omega);
        }
        if idx <= 6 {
            omega = q * (a10(e) + omega);
        }
        if idx <= 7 {
            omega = q * (a9(e) + omega);
        }
        if idx <= 8 {
            omega = q * (a8(e) + omega);
        }
        if idx <= 9 {
            omega = q * (a7(e) + omega);
        }
        if idx <= 10 {
            omega = q * (a6(e) + omega);
        }
        if idx <= 11 {
            omega = q * (a5(e) + omega);
        }
        omega
    };
    (t / r) * omega
}

#[inline]
fn yprime_tail_expansion_rational_function_part(w: f64) -> f64 {
    w * (-2.9999999999994663866 + w * (-1.7556263323542206288E2 + w * (-3.4735035445495633334E3 + w * (-2.7805745693864308643E4 + w * (-8.3836021460741980839E4 - 6.6818249032616849037E4 * w))))) / (1.0 + w * (6.3520877744831739102E1 + w * (1.4404389037604337538E3 + w * (1.4562545638507033944E4 + w * (6.6886794165651675684E4 + w * (1.2569970380923908488E5 + 6.9286518679803751694E4 * w))))))
}

fn yprime(h: f64) -> f64 {
    // We copied the thresholds of -0.46875 and -4 from Cody.
    if h < -4.0 {
        // Nonlinear-Remez optimized minimax rational function of order (5,6) for g(w) := (Y'(h)/h²-1)/h² with w:=1/h².
        // The relative accuracy of Y'(h) ≈ w·(1+w·g(w)) is better than 9.8E-17 (in perfect arithmetic) on h in [-∞,-4] (i.e., on w in [0,1/16]).
        let w = 1.0 / (h * h);
        w * (1.0 + yprime_tail_expansion_rational_function_part(w))
    }
    else if h <= -0.46875 {
        // Remez-optimized minimax rational function of order (7,7) of relative accuracy better than 1.6E-16 (in perfect arithmetic) on h in [-4,-0.46875].
        (1.0000000000594317229 - h * (6.1911449879694112749E-1 - h * (2.2180844736576013957E-1 - h * (4.5650900351352987865E-2 - h * (5.545521007735379052E-3 - h * (3.0717392274913902347E-4 - h * (4.2766597835908713583E-8 + 8.4592436406580605619E-10 * h))))))) / (1.0 - h * (1.8724286369589162071 - h * (1.5685497236077651429 - h * (7.6576489836589035112E-1 - h * (2.3677701403094640361E-1 - h * (4.6762548903194957675E-2 - h * (5.5290453576936595892E-3 - 3.0822020417927147113E-4 * h)))))))
    }else {
        1.0 + h * SQRT_PI_OVER_TWO * erfcx_cody(-FRAC_1_SQRT_2 * h)
    }
}

fn small_t_expansion_of_scaled_normalised_black(h: f64, t: f64) -> f64 {
    let a = yprime(h);
    let h2 = h * h;
    let t2 = t * t;
    fn b0(a: f64) -> f64 { 2.0 * a}
    fn b1(a: f64, h2: f64) -> f64 {(-1.0+a*(3.0+h2))/3.0}
    fn b2(a: f64, h2: f64) -> f64 {(-7.0+h2*(10.0+h2)+a*(15.0+h2*(10.0+h2)))/60.0}
    fn b3(a: f64, h2: f64) -> f64 {(-57.0+(-18.0-h2)*h2+a*(105.0+h2*(105.0+h2*(21.0+h2))))/2520.0}
    fn b4(a: f64, h2: f64) -> f64 {(-561.0+h2*(-285.0+(-33.0-h2)*h2)+a*(945.0+h2*(1260.0+h2*(378.0+h2*(36.0+h2)))))/181440.0}
    fn b5(a: f64, h2: f64) -> f64 {(-6555.0+h2*(-4680.0+h2*(-840.0+(-52.0-h2)*h2))+a*(10395.0+h2*(17325.0+h2*(6930.0+h2*(990.0+h2*(55.0+h2))))))/19958400.0}
    fn b6(a: f64, h2: f64) -> f64 {(-89055.0+h2*(-82845.0+h2*(-20370.0+h2*(-1926.0+(-75.0-h2)*h2)))+a*(135135.0+h2*(270270.0+h2*(135135.0+h2*(25740.0+h2*(2145.0+h2*(78.0+h2)))))))/3113510400.0}
    t * (b0(a) + t2 * (b1(a, h2) + t2 * (b2(a, h2) + t2 * (b3(a, h2) + t2 * (b4(a, h2) + t2 * (b5(a, h2) + b6(a, h2) * t2))))))
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
    assert!(s > 0.0, "s must be positive, got: {}", s);
    let h = x / s;
    let t = 0.5 * s;
    SQRT_TWO_PI.recip() * (-0.5 * (h *  h + t * t)).exp()
}

#[inline(always)]
fn inv_normalised_vega(x: f64, s: f64) -> f64 {
    assert!(s > 0.0, "s must be positive, got: {}", s);
    let h = x / s;
    let t = 0.5 * s;
    SQRT_TWO_PI * (0.5 * (h *  h + t * t)).exp()
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
fn normalised_black(x: f64, s: f64) -> f64 {
    assert!(x <= 0.0, "x: {}", x);
    assert!(s > 0.0, "s: {}", s);
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
        return (normalised_intrinsic(x) * (-ln_vega).exp() + bx, ln_vega);
    }
    let ln_vega = ln_normalised_vega(x, s);
    if s <= x.abs() * DENORMALISATION_CUTOFF {
        return (normalised_intrinsic(x) * (-ln_vega).exp(), ln_vega);
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
pub(crate) fn black(f: f64, k: f64, sigma: f64, t: f64, q: bool) -> f64 {
    let s = sigma * t.sqrt();
    if k == f{
        f * erf_cody((0.5 / SQRT_2) * s)
    }else{
        (if q { f- k } else {k - f}).max(0.0) + (if s <= 0.0 {0.0} else {
            f.sqrt() * k.sqrt() * normalised_black((f / k).ln().abs().neg(), s)
        })
    }
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

fn one_minus_erfcx(x: f64) -> f64 {
    if x < -1.0 / 5.0 || x > 1.0 / 3.0 {
        1.0 - erfcx_cody(x)
    } else {
        x * (std::f64::consts::FRAC_2_SQRT_PI - x * (1.0000000000000002 + x * (1.1514967181784756 + x * (5.7689001208873741E-1 + x * (1.4069188744609651E-1 + 1.4069285713634565E-2 * x)))) / (1.0 + x * (1.9037494962421563 + x * (1.5089908593742723 + x * (6.2486081658640257E-1 + x * (1.358008134514386E-1 + 1.2463320728346347E-2 * x))))))
    }
}

#[inline(always)]
fn lets_be_rational(
    beta: f64, x: f64, n: u8,
) -> f64 {
    assert!(x  <= 0.0, "x must be non-positive, but got {}", x);
    if beta <= 0. {
        return if beta == 0.0 {0.0} else {VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC};
    }
    let b_max = (0.5 * x).exp();
    if beta >= b_max {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    let mut iterations = 0;
    let mut f = f64::MIN;
    let mut s;
    let mut ds = f64::MIN;

    let sqrt_ax = x.neg().sqrt();
    let s_c = SQRT_2 * sqrt_ax;
    let ome = one_minus_erfcx(sqrt_ax);
    let b_c = 0.5 * b_max * ome;
    if beta < b_c {
        assert!(x < 0.0, "x must be negative, but got {}", x);
        let s_l = s_c - SQRT_PI_OVER_TWO * ome;
        debug_assert!(s_l > 0.0, "s_l must be positive, but got {}", s_l);
        let b_l = normalised_black(x, s_l);
        // let b_l = b_l_over_b_max(s_c) * b_max;
        if beta < b_l {
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) = compute_f_lower_map_and_first_two_derivatives(x, s_l);
            let r2 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
            f = rational_cubic_interpolation(beta, 0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, r2);
            if !(f > 0.0) {
                let t = beta / b_l;
                f = (f_lower_map_l * t + b_l * (1.0 - t)) * t;
            }
            s = inverse_f_lower_map(x, f);
            assert!(s > 0.0, "s must be positive, but got {}", s);
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
                s += ds;
                assert!(s > 0.0, "s must be positive, but got {}", s);
                iterations += 1;
            }
            return s;
        } else {
            let inv_v_c = SQRT_TWO_PI / b_max;
            let inv_v_l = inv_normalised_vega(x, s_l);
            let r_im = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l, b_c, s_l, s_c, inv_v_l, inv_v_c, 0.0, false);
            s = rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, inv_v_l, inv_v_c, r_im);
            assert!(s > 0.0, "s must be positive, but got {}", s);
        }
    } else {
        let s_u = s_c + SQRT_PI_OVER_TWO * (2.0 - ome);
        assert!(s_u > 0.0, "s_u must be positive, but got {}", s_u);
        let b_u = normalised_black(x, s_u);
        if beta <= b_u {
            let inv_v_c = SQRT_TWO_PI / b_max;

            let inv_v_u = inv_normalised_vega(x, s_u);
            let r_u_m = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
                b_c, b_u, s_c, s_u, inv_v_c, inv_v_u, 0.0, false);
            s = rational_cubic_interpolation(beta, b_c, b_u, s_c, s_u, inv_v_c, inv_v_u, r_u_m);
            assert!(s > 0.0, "s must be positive, but got {}", s);
        } else {
            let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) = compute_f_upper_map_and_first_two_derivatives(x, s_u);

            if d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX {
                let r_uu = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_u, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, d2_f_upper_map_h_d_beta2, true);
                f = rational_cubic_interpolation(beta, b_u, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, r_uu);
            }
            if f <= 0.0 {
                let h = b_max - b_u;
                let t = (beta - b_u) / h;
                f = (f_upper_map_h * (1.0 - t) + 0.5 * h * t) * (1.0 - t);
            }
            s = inverse_f_upper_map(f);
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
                    s += ds;
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

        let b = normalised_black(x, s);
        let bp = normalised_vega(x, s);
        let nu = (beta - b) / bp;
        let h = x / s;
        let h2 = h.powi(2) / s - s / 4.0;
        let h3 = h2 * h2 - 3.0 * (h / s).powi(2) - 0.25;
        ds = nu * householder3_factor(nu, h2, h3);
        // Never leave the branch (or bracket)
        s += ds;
    }
    s
}

#[inline(always)]
pub(crate) fn implied_black_volatility(price: f64, f: f64, k: f64, t: f64, q: bool) -> f64 {
    if price >= if q { f } else { k } {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    let mu = if q { f - k } else { k - f };
    lets_be_rational(if mu > 0.0 {price - mu} else {price} / (f.sqrt() * k.sqrt()), (f / k).ln().abs().neg(), 2) / t.sqrt()
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
            assert!((price - reprice).abs() < f64::EPSILON * 100.0, "{}", (price - reprice).abs() / 100.0);
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
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
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
