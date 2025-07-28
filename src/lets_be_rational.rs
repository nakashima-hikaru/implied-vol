use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};
use std::ops::Neg;
use crate::constants::{HALF_OF_LN_TWO_PI, ONE_OVER_SQRT_THREE, SIXTEENTH_ROOT_DBL_EPSILON, SQRT_DBL_MAX, SQRT_PI_OVER_TWO, SQRT_THREE, SQRT_THREE_OVER_THIRD_ROOT_TWO_PI, SQRT_TWO_OVER_PI, SQRT_TWO_PI, TWO_PI_OVER_SQRT_TWENTY_SEVEN, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC};
use crate::erf_cody::{erf_cody, erfc_cody, erfcx_cody};
use crate::MulAdd;
use crate::normal_distribution::{erfinv, inverse_norm_cdf, norm_cdf, norm_pdf};
use crate::rational_cubic::{convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side, convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side, rational_cubic_interpolation};


#[inline(always)]
fn householder3_factor(v: f64, h2: f64, h3: f64) -> f64 { (1.0 + 0.5 * h2 * v) / (1.0 + v * (h2 + h3 * v / 6.0)) }

#[inline(always)]
fn householder4_factor(v: f64, h2: f64, h3: f64, h4: f64) -> f64 { (1.0 + v * (h2 + v * h3 / 6.0)) / (1.0 + v * (1.5 * h2 + v * (h2 * h2 / 4.0 + h3 / 3.0 + v * h4 / 24.0))) }

const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
const TAU: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

#[inline(always)]
fn asymptotic_expansion_of_scaled_normalised_black(h: f64, t: f64) -> f64 {
    #[inline(always)]
    const fn a0(_e: f64) -> f64 {2.0}
    #[inline(always)]
    const fn a1(e: f64) -> f64 {-6.0-2.0*e}
    #[inline(always)]
    const fn a2(e: f64) -> f64 {30.0+e*(60.0+6.0*e)}
    #[inline(always)]
    const fn a3(e: f64) -> f64 { -2.1E2+e*(-1.05E3+e*(-6.3E2-30.0*e))}
    #[inline(always)]
    const fn a4(e: f64) -> f64 { 1.89E3+e*(1.764E4+e*(2.646E4+e*(7.56E3+2.1E2*e)))}
    #[inline(always)]
    const fn a5(e: f64) -> f64 {-2.079E4+e*(-3.1185E5+e*(-8.7318E5+e*(-6.237E5+e*(-1.0395E5-1.89E3*e))))}
    #[inline(always)]
    const fn a6(e: f64) -> f64 {2.7027E5+e*(5.94594E6+e*(2.675673E7+e*(3.567564E7+e*(1.486485E7+e*(1.62162E6+2.079E4*e)))))}
    #[inline(always)]
    const fn a7(e: f64) -> f64 {-4.05405E6+e*(-1.2297285E8+e*(-8.1162081E8+e*(-1.73918745E9+e*(-1.35270135E9+e*(-3.6891855E8+e*(-2.837835E7-2.7027E5*e))))))}
    #[inline(always)]
    const fn a8(e: f64) -> f64 {6.891885E7+e*(2.756754E9+e*(2.50864614E10+e*(7.88431644E10+e*(9.85539555E10+e*(5.01729228E10+e*(9.648639E9+e*(5.513508E8+4.05405E6*e)))))))}
    #[inline(always)]
    const fn a9(e: f64) -> f64 {-1.30945815E9+e*(-6.678236565E10+e*(-8.013883878E11+e*(-3.4726830138E12+e*(-6.3665855253E12+e*(-5.2090245207E12+e*(-1.8699062382E12+e*(-2.671294626E11+e*(-1.178512335E10-6.891885E7*e))))))))}
    #[inline(always)]
    const fn a10(e: f64) -> f64 {2.749862115E10+e*(1.7415793395E12+e*(2.664616389435E13+e*(1.52263793682E14+e*(3.848890340295E14+e*(4.618668408354E14+e*(2.664616389435E14+e*(7.10564370516E13+e*(7.83710702775E12+e*(2.749862115E11+1.30945815E9*e)))))))))}
    #[inline(always)]
    const fn a11(e: f64) -> f64 {-6.3246828645E11+e*(-4.870005805665E13+e*(-9.2530110307635E14+e*(-6.74147946527055E15+e*(-2.24715982175685E16+e*(-3.71802806872497E16+e*(-3.14602375045959E16+e*(-1.34829589305411E16+e*(-2.77590330922905E15+e*(-2.4350029028325E14+e*(-6.95715115095E12-2.749862115E10*e))))))))))}
    #[inline(always)]
    const fn a12(e: f64) -> f64 {1.581170716125E13+e*(1.454677058835E15+e*(3.36030400590885E16+e*(3.04027505296515E17+e*(1.292_116_897_510_188_8E18+e*(2.81916414002223E18+e*(3.289024830025935E18+e*(2.067387036016302E18+e*(6.840_618_869_171_588E17+e*(1.12010133530295E17+e*(8.0007238235925E15+e*(1.89740485935E14+6.3246828645E11*e)))))))))))}
    #[inline(always)]
    const fn a13(e: f64) -> f64 {-4.2691609335375E14+e*(-4.624924344665625E16+e*(-1.2764791191277125E18+e*(-1.404_127_031_040_483_7E19+e*(-7.410_670_441_602_553E19+e*(-2.061_513_777_391_255_6E20+e*(-3.171_559_657_525_009E20+e*(-2.748_685_036_521_674_2E20+e*(-1.333_920_679_488_459_6E20+e*(-3.510_317_577_601_209_5E19+e*(-4.680_423_436_801_613E18+e*(-2.774954606799375E17+e*(-5.54990921359875E15-1.581170716125E13*e))))))))))))}
    #[inline(always)]
    const fn a14(e: f64) -> f64 {1.238056670725875E16+e*(1.5599514051146025E18+e*(5.069_842_066_622_458E19+e*(6.663_221_001_846_66E20+e*(4.275_566_809_518_273E21+e*(1.477_013_988_742_676E22+e*(2.897_219_747_149_095_4E22+e*(3.311_108_282_456_109E22+e*(2.215_520_983_114_014E22+e*(8.551_133_619_036_546E21+e*(1.832_385_775_507_831_3E21+e*(2.027_936_826_648_983_3E20+e*(1.013_968_413_324_491_6E19+e*(1.733279339016225E17+4.2691609335375E14*e)))))))))))))}
    #[inline(always)]
    const fn a15(e: f64) -> f64 {-3.8379756792502125E17+e*(-5.565_064_734_912_808E19+e*(-2.103_594_469_797_041_5E21+e*(-3.255_562_869_923_992_8E22+e*(-2.495_931_533_608_394_5E23+e*(-1.048_291_244_115_525_7E24+e*(-2.553_529_953_614_742E24+e*(-3.720_857_932_410_053E24+e*(-3.283_109_940_361_811E24+e*(-1.747_152_073_525_876E24+e*(-5.491_049_373_938_467_5E23+e*(-9.766_688_609_771_979E22+e*(-9.115_576_035_787_18E21+e*(-3.895_545_314_438_965_5E20+e*(-5.756_963_518_875_318E18-1.238056670725875E16*e))))))))))))))}
    #[inline(always)]
    const fn a16(e: f64) -> f64 {1.266_531_974_152_57E19+e*(2.093_999_530_598_916E21+e*(9.108_897_958_105_285E22+e*(1.639_601_632_458_951_2E24+e*(1.480_195_918_192_108_6E25+e*(7.427_892_244_018_582E25+e*(2.199_798_856_882_426E26+e*(3.980_588_407_692_009_4E26+e*(4.478_161_958_653_511E26+e*(3.142_569_795_546_323E26+e*(1.361_780_244_736_74E26+e*(3.552_470_203_661_060_7E25+e*(5.328_705_305_491_592E24+e*(4.250_819_047_115_799E23+e*(1.570_499_647_949_187E22+e*(2.026_451_158_644_112E20+3.8379756792502125E17*e)))))))))))))))}
    const THRESHOLDS: [f64; 12] = [12.347, 12.958, 13.729, 14.718, 16.016, 17.769, 20.221, 23.816, 29.419, 38.93, 57.171, 99.347];

    assert!((h < -ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD.abs()) && (h + t < -(TAU + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD).abs()));
    let e = (t / h).powi(2);
    let r = (h + t) * (h - t);
    let q = (h / r).powi(2);

    let idx = THRESHOLDS.partition_point(|&x| x<=-h - t + TAU);
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

#[cfg(not(feature = "use_original_taylor_expansion"))]
#[inline(always)]
fn y_prime_tail_expansion_rational_function_part(w: f64) -> f64 {
    w * (-2.999_999_999_999_466 + w * (-1.755_626_332_354_220_6E2 + w * (-3.473_503_544_549_563_2E3 + w * (-2.780_574_569_386_430_8E4 + w * (-8.383_602_146_074_198E4 - 6.681_824_903_261_685E4 * w))))) / (1.0 + w * (6.352_087_774_483_173_6E1 + w * (1.440_438_903_760_433_7E3 + w * (1.456_254_563_850_703_4E4 + w * (6.688_679_416_565_168E4 + w * (1.256_997_038_092_390_9E5 + 6.928_651_867_980_375E4 * w))))))
}

#[cfg(not(feature = "use_original_taylor_expansion"))]
#[inline(always)]
fn y_prime(h: f64) -> f64 {
    // We copied the thresholds of -0.46875 and -4 from Cody.
    if h < -4.0 {
        // Nonlinear-Remez optimized minimax rational function of order (5,6) for g(w) := (Y'(h)/h²-1)/h² with w:=1/h².
        // The relative accuracy of Y'(h) ≈ w·(1+w·g(w)) is better than 9.8E-17 (in perfect arithmetic) on h in [-∞,-4] (i.e., on w in [0,1/16]).
        let w = (h * h).recip();
        w * (1.0 + y_prime_tail_expansion_rational_function_part(w))
    } else if h <= -0.46875 {
        // Remez-optimized minimax rational function of order (7,7) of relative accuracy better than 1.6E-16 (in perfect arithmetic) on h in [-4,-0.46875].
        (1.000_000_000_059_431_8 - h * (6.191_144_987_969_411E-1 - h * (2.218_084_473_657_601_4E-1 - h * (4.565_090_035_135_299E-2 - h * (5.545_521_007_735_379_5E-3 - h * (3.071_739_227_491_390_3E-4 - h * (4.276_659_783_590_871_4E-8 + 8.459_243_640_658_06E-10 * h))))))) / (1.0 - h * (1.872_428_636_958_916_3 - h * (1.568_549_723_607_765_2 - h * (7.657_648_983_658_903E-1 - h * (2.367_770_140_309_464E-1 - h * (4.676_254_890_319_496E-2 - h * (5.529_045_357_693_659E-3 - 3.082_202_041_792_714_7E-4 * h)))))))
    } else {
        1.0 + h * SQRT_PI_OVER_TWO * erfcx_cody(-FRAC_1_SQRT_2 * h)
    }
}

#[cfg(not(feature = "use_original_taylor_expansion"))]
#[inline(always)]
fn small_t_expansion_of_scaled_normalised_black(h: f64, t: f64) -> f64 {
    let a = y_prime(h);
    let h2 = h * h;
    let t2 = t * t;
    #[inline(always)]
    fn b0(a: f64) -> f64 { 2.0 * a }
    #[inline(always)]
    fn b1(a: f64, h2: f64) -> f64 {(-1.0+a*(3.0+h2))/3.0}
    #[inline(always)]
    fn b2(a: f64, h2: f64) -> f64 {(-7.0-h2+a*(15.0+h2*(10.0+h2)))/60.0}
    #[inline(always)]
    fn b3(a: f64, h2: f64) -> f64 {(-57.0+(-18.0-h2)*h2+a*(105.0+h2*(105.0+h2*(21.0+h2))))/2520.0}
    #[inline(always)]
    fn b4(a: f64, h2: f64) -> f64 {(-561.0+h2*(-285.0+(-33.0-h2)*h2)+a*(945.0+h2*(1260.0+h2*(378.0+h2*(36.0+h2)))))/181440.0}
    #[inline(always)]
    fn b5(a: f64, h2: f64) -> f64 {(-6555.0+h2*(-4680.0+h2*(-840.0+(-52.0-h2)*h2))+a*(10395.0+h2*(17325.0+h2*(6930.0+h2*(990.0+h2*(55.0+h2))))))/19958400.0}
    #[inline(always)]
    fn b6(a: f64, h2: f64) -> f64 {(-89055.0+h2*(-82845.0+h2*(-20370.0+h2*(-1926.0+(-75.0-h2)*h2)))+a*(135135.0+h2*(270270.0+h2*(135135.0+h2*(25740.0+h2*(2145.0+h2*(78.0+h2)))))))/3113510400.0}
    t * (b0(a) + t2 * (b1(a, h2) + t2 * (b2(a, h2) + t2 * (b3(a, h2) + t2 * (b4(a, h2) + t2 * (b5(a, h2) + b6(a, h2) * t2))))))
}

#[cfg(feature = "use_original_taylor_expansion")]
#[inline(always)]
fn small_t_expansion_of_scaled_normalised_black(h: f64, t: f64) -> f64 {
    let zeta = t * t;
    let g = h * h;
    let a = 1_f64 + h * SQRT_PI_OVER_TWO * erfcx_cody(-FRAC_1_SQRT_2 * h);
    let b_over_vega = 2.0 * t * (a + zeta * ((-1.0 + 3.0 * a + a * g) / 6.0 + zeta * ((-7.0 + 15.0 * a + g * (-1.0 + 10.0 * a + a * g)) / 120.0 + zeta * ((-57.0 + 105.0 * a + g * (-18.0 + 105.0 * a + g * (-1.0 + 21.0 * a + a * g))) / 5040.0 + zeta * ((-561.0 + 945.0 * a + g * (-285.0 + 1260.0 * a + g * (-33.0 + 378.0 * a + g * (-1.0 + 36.0 * a + a * g)))) / 362880.0 + zeta * ((-6555.0 + 10395.0 * a + g * (-4680.0 + 17325.0 * a + g * (-840.0 + 6930.0 * a + g * (-52.0 + 990.0 * a + g * (-1.0 + 55.0 * a + a * g))))) / 39916800.0 + ((-89055.0 + 135135.0 * a + g * (-82845.0 + 270270.0 * a + g * (-20370.0 + 135135.0 * a + g * (-1926.0 + 25740.0 * a + g * (-75.0 + 2145.0 * a + g * (-1.0 + 78.0 * a + a * g)))))) * zeta) / 6227020800.0))))));
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
                0.5 * ((0.5 * x).exp() * erfc_cody(q1) - (-0.5 * (h * h + t * t)).exp() * erfcx_cody(q2))
            }
        } else if q2 < CODYS_THRESHOLD {
            0.5 * ((-0.5 * (h * h + t * t)).exp() * erfcx_cody(q1) - (-0.5 * x).exp() * erfc_cody(q2))
        } else {
            0.5 * ((-0.5 * (h * h + t * t)).exp() * (erfcx_cody(q1) - erfcx_cody(q2)))
        };
    two_b.abs().max(0.0)
}

#[inline(always)]
fn normalised_vega(x: f64, s: f64) -> f64 {
    assert!(s > 0.0, "s must be positive, got: {s}");
    let h = x / s;
    let t = 0.5 * s;
    SQRT_TWO_PI.recip() * (-0.5 * (h *  h + t * t)).exp()
}

#[inline(always)]
fn inv_normalised_vega(x: f64, s: f64) -> f64 {
    assert!(s > 0.0, "s must be positive, got: {s}");
    let h = x / s;
    let t = 0.5 * s;
    SQRT_TWO_PI * (0.5 * (h * h + t * t)).exp()
}

#[inline(always)]
fn ln_normalised_vega(x: f64, s: f64) -> f64 {
    let ax = x.abs();
    if ax <= 0.0 {
        -HALF_OF_LN_TWO_PI - 0.125 * s * s
    } else if s <= 0.0 {
        f64::MIN
    } else {
        -HALF_OF_LN_TWO_PI - 0.5 * ((x / s).powi(2) + 0.25 * s * s)
    }
}

#[inline(always)]
fn b_u_over_b_max(s_c: f64) -> f64{
    if s_c >= 2.449_489_742_783_178 {
        let y = s_c.recip();
        let g = (-4.605_394_817_212_609E-2 + y * (1.375_163_082_077_259_1 + y * (1.455_319_886_249_397_7E1 + y * (8.488_089_220_080_239E1 + y * (2.983_680_162_805_663E2 + y * (6.169_692_835_129_17E2 + y * (6.589_280_957_677_407E2 - 1.229_189_712_271_654_4 * y))))))) / (1.0 + y * (9.327_034_903_790_405 + y * (5.436_378_146_588_073E1 + y * (2.030_420_459_952_177_3E2 + y * (5.079_647_179_123_228E2 + y * (8.698_830_313_690_185E2 + y * (8.881_238_333_960_678E2 + 5.206_084_752_279_256E2 * y)))))));
        // f = Φ(√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · g )
        return 0.894_954_297_278_031_3 + 0.113_984_531_941_499_06 * y * (-1.253_314_137_315_500_3 + y * g);
    }
    let g = (-6.063_099_880_334_851E-2 + s_c * (-1.344_864_378_589_371E-1 + s_c * (-1.416_368_116_424_721E-1 + s_c * (-8.976_383_086_137_545E-2 + s_c * (-3.512_133_741_041_69E-2 + s_c * (-8.143_812_878_548_491E-3 + s_c * (-8.733_991_026_156_887E-4 - 3.386_756_817_001_176_5E-9 * s_c))))))) / (1.0 + s_c * (2.722_003_340_655_505_5 + s_c * (3.650_335_036_015_884_6 + s_c * (3.018_638_953_766_389_6 + s_c * (1.652_734_794_196_848_7 + s_c * (5.959_161_649_351_221E-1 + s_c * (1.324_801_623_892_073E-1 + 1.421_206_743_529_177_8E-2 * s_c)))))));
    // f = c₀ + y²·(c₂+y·g)
    0.789_908_594_556_062_8 + (s_c * s_c) * (0.061_461_680_580_514_74 + s_c * g)
}

#[inline(always)]
fn b_l_over_b_max(s_c: f64) -> f64{
    if s_c < 2.6267851073127395 {
        if s_c < 0.7099295739719539 {
            // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x).
            // For small |x|, i.e., small y, we have f ≈ (exp(-1/π)/4-Φ(-√(π/2))/2)·y² + exp(-1/π)/(3·√(2π))·y³ + O(y⁴).
            //   c₂ := (exp(-1/π)/4-Φ(-√(π/2))/2) =  0.07560996640296361767172
            //   c₃ := exp(-1/π)/(3·√(2π))        = -0.09672719281339436290858
            // Nonlinear-Remez optimized minimax rational function of order (5,4) for g(y) := ((f/y²-c₂)/y-c₃)/y .   f = y²·(c₂+y·(c₃+y·g)) .
            let g = (8.074_107_237_288_286E-2 + s_c * (9.807_891_178_635_89E-2 + s_c * (3.976_063_144_567_705_5E-2 + s_c * (5.971_692_845_958_919E-3 + s_c * (-6.403_639_934_147_98E-6 + 4.542_510_209_361_606_4E-7 * s_c))))) / (1.0 + s_c * (1.859_497_767_228_766_5 + s_c * (1.365_880_147_571_179 + s_c * (4.613_270_710_865_565E-1 + 6.125_459_704_983_172E-2 * s_c))));
            // Branch I. Accuracy better than 7.43E-17 in perfect arithmetic.
            return (s_c * s_c) * (0.075_609_966_402_963_62 + s_c * (s_c * g - 0.096_727_192_813_394_37));
        }
        // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
        // Branch II. Accuracy better than 8.77E-17 in perfect arithmetic.
        return (1.979_573_792_759_858E-9 + s_c * (-2.708_128_856_468_558_7E-8 + s_c * (7.561_014_227_254_904E-2 + s_c * (6.917_130_174_466_835E-2 + s_c * (2.953_705_895_096_301_8E-2 + s_c * (6.584_925_270_230_231E-3 + 6.971_140_063_983_471E-4 * s_c)))))) / (1.0 + s_c * (2.194_144_852_558_658 + s_c * (2.129_710_354_999_518 + s_c * (1.157_148_318_717_978_3 + s_c * (3.783_162_225_306_046E-1 + s_c * (7.171_486_244_882_935E-2 + 6.636_197_582_786_12E-3 * s_c))))));
    }
    if s_c < 7.348469228349534 {
        // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bₗ(-y²/2)/bₘₐₓ(-y²/2).
        // Branch III. Accuracy better than 7.49E-17 in perfect arithmetic.
        return (-9.332_511_535_483_788E-5 + s_c * (5.311_803_397_279_465E-4 + s_c * (7.411_485_544_834_501E-2 + s_c * (7.403_965_818_682_282E-2 + s_c * (3.922_517_740_768_760_6E-2 + s_c * (1.002_291_337_825_409E-2 + 1.701_257_940_724_605_5E-3 * s_c)))))) / (1.0 + s_c * (2.221_723_813_222_813_4 + s_c * (2.344_181_670_708_740_4 + s_c * (1.391_232_364_627_114 + s_c * (5.323_125_844_350_184E-1 + s_c * (1.174_400_591_971_610_1E-1 + 1.619_540_589_593_093_7E-2 * s_c))))));
    }
    (1.450_007_229_724_060_4E-3 + s_c * (-1.511_669_248_501_119_6E-3 + s_c * (7.168_217_831_093_633E-2 + s_c * (3.921_610_857_820_463_6E-2 + s_c * (2.934_240_565_862_844_5E-2 + s_c * (5.183_252_617_163_152E-3 + 1.693_020_807_842_147_5E-3 * s_c)))))) / (1.0 + s_c * (1.617_631_350_230_541_5 + s_c * (1.682_315_917_528_153_2 + s_c * (8.487_830_756_737_222E-1 + s_c * (3.754_374_213_737_579E-1 + s_c * (7.126_137_099_644_303E-2 + 1.611_699_254_678_867_7E-2 * s_c))))))
}

#[inline(always)]
fn normalised_black(x: f64, s: f64) -> f64 {
    assert!(x < 0.0, "x: {x}");
    assert!(s > 0.0, "s: {s}");
    if is_region1(x, s) {
        asymptotic_expansion_of_scaled_normalised_black(x / s, 0.5 * s) * normalised_vega(x, s)
    }else if is_region2(x, s) {
        small_t_expansion_of_scaled_normalised_black(x / s, 0.5 * s) * normalised_vega(x, s)
    }else {
        normalised_black_call_with_optimal_use_of_codys_functions(x, s)
    }
}

const ETA: f64 = -13.0;
#[inline(always)]
fn is_region1(x: f64, s: f64) -> bool {
    x < s * ETA && s * (0.5 * s - (TAU + 0.5 + ETA)) + x < 0.0
}
#[inline(always)]
fn is_region2(x: f64, s: f64) -> bool {
    s * (s - (2.0 * TAU)) - x / ETA < 0.0
}

#[inline(always)]
fn scaled_normalised_black_and_ln_vega(x: f64, s: f64) -> (f64, f64) {
    assert!(x < 0.0, "x must be negative, got: {x}");
    assert!(s > 0.0, "s must be positive, got: {s}");
    let ln_vega = ln_normalised_vega(x, s);
    if is_region1(x, s) {
        (asymptotic_expansion_of_scaled_normalised_black(x / s, 0.5 * s), ln_vega)
    }
    else if is_region2(x, s) {
        (small_t_expansion_of_scaled_normalised_black(x / s, 0.5 * s), ln_vega)
    } else {
        (normalised_black_call_with_optimal_use_of_codys_functions(x, s) * (-ln_vega).exp(), ln_vega)
    }
}

#[inline(always)]
pub(crate) fn black(f: f64, k: f64, sigma: f64, t: f64, q: bool) -> f64 {
    let s = sigma * t.sqrt();
    if k == f {
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
    let y = z * z;
    let s2 = s * s;
    let phi_m = norm_cdf(-z);
    let phi = norm_pdf(z);
    let fpp = std::f64::consts::FRAC_PI_6 * y / (s2 * s) * phi_m * (8.0 * SQRT_THREE * s * ax + (3.0 * s2 * (s2 - 8.0) - 8.0 * x * x) * phi_m / phi) * (2.0 * y + 0.25 * s2).exp();
    let (fp, f) = if s.is_subnormal() {
        (1.0, 0.0)
    } else {
        let phi2 = phi_m * phi_m;
        let fp_val = std::f64::consts::TAU * y * phi2 * (y + 0.125 * s * s).exp();
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
        SQRT_PI_OVER_TWO * ((w + 0.125 * s * s).exp()) * w / s)
    };

    (f, fp, fpp)
}


#[inline(always)]
fn inverse_f_upper_map(f: f64) -> f64 {
    -2.0 * inverse_norm_cdf(f)
}

fn one_minus_erfcx(x: f64) -> f64 {
    if !(-1.0 / 5.0..=1.0 / 3.0).contains(&x) {
        1.0 - erfcx_cody(x)
    } else {
        x * (std::f64::consts::FRAC_2_SQRT_PI - x * (1.0000000000000002 + x * (1.1514967181784756 + x * (5.768_900_120_887_374E-1 + x * (1.406_918_874_460_965E-1 + 1.4069285713634565E-2 * x)))) / (1.0 + x * (1.9037494962421563 + x * (1.5089908593742723 + x * (6.248_608_165_864_026E-1 + x * (1.358008134514386E-1 + 1.2463320728346347E-2 * x))))))
    }
}

#[inline(always)]
fn implied_normalised_volatility_atm(beta: f64) -> f64 {
    2.0 * SQRT_2 * erfinv(beta)
}

#[inline(always)]
fn lets_be_rational(
    beta: f64, x: f64, n: u8,
) -> f64 {
    assert!(x  <= 0.0, "x must be non-positive, but got {x}");
    if beta <= 0. {
        return if beta == 0.0 {0.0} else {VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC};
    }
    let b_max = (0.5 * x).exp();
    if beta >= b_max {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    if x == 0.0 {
        return implied_normalised_volatility_atm(beta);
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
        assert!(x < 0.0, "x must be negative, but got {x}");
        let s_l = s_c - SQRT_PI_OVER_TWO * ome;
        debug_assert!(s_l > 0.0, "s_l must be positive, but got {s_l}");
        // let b_l = normalised_black(x, s_l);
        let b_l = b_l_over_b_max(s_c) * b_max;
        if beta < b_l {
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) = compute_f_lower_map_and_first_two_derivatives(x, s_l);
            let r2 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
            f = rational_cubic_interpolation(beta, 0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, r2);
            if !(f > 0.0) {
                let t = beta / b_l;
                f = (f_lower_map_l * t + b_l * (1.0 - t)) * t;
            }
            s = inverse_f_lower_map(x, f);
            assert!(s > 0.0, "s must be positive, but got {s}");
            let ln_beta = beta.ln();

            ds = 1.0_f64;
            while iterations < n && ds.abs() > f64::EPSILON * s {
                let (bx, ln_vega) = scaled_normalised_black_and_ln_vega(x, s);
                let ln_b = bx.ln() + ln_vega;
                let bpob = bx.recip();
                let h = x / s;
                let b_h2 = ((h * h) / s) - s / 4.0;
                let nu = (ln_beta - ln_b) * ln_b / ln_beta * bx;
                let lambda = ln_b.recip();
                let otlambda = lambda.mul_add2(2.0, 1.0);
                let h2 = b_h2 - bpob * otlambda;
                let c = 3.0 * (h / s).powi(2);
                let b_h3 = b_h2 * b_h2 - c - 0.25;
                let sq_bpob = bpob * bpob;
                let mu = 6.0 * lambda * (1.0 + lambda);
                let h3 = b_h3 + sq_bpob * (2.0 + mu) - (b_h2 * bpob * 3.0 * otlambda);
                ds = if x < -190.0 {
                    nu * householder4_factor(nu, h2, h3, ((b_h2 * (b_h3 - 0.5)) - ((b_h2 - 2.0 / s) * 2.0 * c)) - (bpob * (sq_bpob * (6.0 + lambda * (22.0 + lambda * (36.0 + lambda * 24.0))) - (b_h2 * bpob * (12.0 + 6.0 * mu))) - (b_h2 * bpob * 3.0 * otlambda) - (b_h3 * bpob * 4.0 * otlambda)))
                } else {
                    nu * householder3_factor(nu, h2, h3)
                };
                s += ds;
                assert!(s > 0.0, "s must be positive, but got {s}");
                iterations += 1;
            }
            return s;
        } else {
            let inv_v_c = SQRT_TWO_PI / b_max;
            let inv_v_l = inv_normalised_vega(x, s_l);
            let r_im = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l, b_c, s_l, s_c, inv_v_l, inv_v_c, 0.0, false);
            s = rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, inv_v_l, inv_v_c, r_im);
            assert!(s > 0.0, "s must be positive, but got {s}");
        }
    } else {
        let s_u = s_c + SQRT_PI_OVER_TWO * (2.0 - ome);
        assert!(s_u > 0.0, "s_u must be positive, but got {s_u}");
        let b_u = b_u_over_b_max(s_c) * b_max;
        if beta <= b_u {
            let inv_v_c = SQRT_TWO_PI / b_max;

            let inv_v_u = inv_normalised_vega(x, s_u);
            let r_u_m = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
                b_c, b_u, s_c, s_u, inv_v_c, inv_v_u, 0.0, false);
            s = rational_cubic_interpolation(beta, b_c, b_u, s_c, s_u, inv_v_c, inv_v_u, r_u_m);
            assert!(s > 0.0, "s must be positive, but got {s}");
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
                    let x_over_s_square = h * h / s;
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
        let h2 = h * h / s - s / 4.0;
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
            assert!((price - reprice).abs() <= 1.5 * f64::EPSILON, "{f},{k},{t},{sigma},{price},{reprice},{}", (price - reprice).abs() / f64::EPSILON);
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
