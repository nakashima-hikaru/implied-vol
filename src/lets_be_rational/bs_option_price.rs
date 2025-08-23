use crate::fused_multiply_add::MulAdd;
use crate::lets_be_rational::constants::{
    FRAC_1_SQRT_2_PI, HALF_OF_LN_2_PI, SIXTEENTH_ROOT_DBL_EPSILON, SQRT_2_PI, SQRT_PI_OVER_2,
};
use crate::lets_be_rational::special_function::SpecialFn;
use std::f64::consts::FRAC_1_SQRT_2;
use std::ops::Neg;

#[inline(always)]
pub fn normalised_black<SpFn: SpecialFn>(half_theta_x: f64, h: f64, t: f64) -> f64 {
    debug_assert!(t > 0.0);
    debug_assert!(half_theta_x < 0.0);
    if is_region1(h, t) {
        asymptotic_expansion_of_scaled_normalised_black(h, t) * normalised_vega(h, t)
    } else if is_region2(h, t) {
        small_t_expansion_of_scaled_normalised_black::<SpFn>(h, t) * normalised_vega(h, t)
    } else {
        normalised_black_with_optimal_use_of_codys_functions::<SpFn>(half_theta_x, h, t)
    }
}

pub const ETA: f64 = -13.0;

#[inline(always)]
fn is_region1(h: f64, t: f64) -> bool {
    h < ETA && t + h < TAU + 0.5 + ETA
}

#[inline(always)]
fn is_region2(h: f64, t: f64) -> bool {
    h.mul_add2(-0.5 / ETA, -TAU + t) < 0.0
}

#[inline(always)]
pub fn scaled_normalised_black_and_ln_vega<SpFn: SpecialFn>(
    half_theta_x: f64,
    h: f64,
    t: f64,
) -> (f64, f64) {
    debug_assert!(half_theta_x < 0.0);
    debug_assert!(t > 0.0);
    let ln_vega = ln_normalised_vega(h, t);
    if is_region1(h, t) {
        (
            asymptotic_expansion_of_scaled_normalised_black(h, t),
            ln_vega,
        )
    } else if is_region2(h, t) {
        (
            small_t_expansion_of_scaled_normalised_black::<SpFn>(h, t),
            ln_vega,
        )
    } else {
        (
            normalised_black_with_optimal_use_of_codys_functions::<SpFn>(half_theta_x, h, t)
                * (-ln_vega).exp(),
            ln_vega,
        )
    }
}

#[inline(always)]
fn asymptotic_expansion_of_scaled_normalised_black(h: f64, t: f64) -> f64 {
    #[inline(always)]
    const fn a0() -> f64 {
        2.0
    }
    #[inline(always)]
    fn a1(e: f64) -> f64 {
        e.mul_add2(-2.0, -6.0)
    }
    #[inline(always)]
    fn a2(e: f64) -> f64 {
        e.mul_add2(6.0, 60.0).mul_add2(e, 30.0)
    }
    #[inline(always)]
    fn a3(e: f64) -> f64 {
        e.mul_add2(-30.0, -6.3E2)
            .mul_add2(e, -1.05E3)
            .mul_add2(e, -2.1E2)
    }
    #[inline(always)]
    fn a4(e: f64) -> f64 {
        e.mul_add2(2.1E2, 7.56E3)
            .mul_add2(e, 2.646E4)
            .mul_add2(e, 1.764E4)
            .mul_add2(e, 1.89E3)
    }
    #[inline(always)]
    fn a5(e: f64) -> f64 {
        e.mul_add2(-1.89E3, -1.0395E5)
            .mul_add2(e, -6.237E5)
            .mul_add2(e, -8.7318E5)
            .mul_add2(e, -3.1185E5)
            .mul_add2(e, -2.079E4)
    }
    #[inline(always)]
    fn a6(e: f64) -> f64 {
        e.mul_add2(2.079E4, 1.62162E6)
            .mul_add2(e, 1.486_485E7)
            .mul_add2(e, 3.567_564E7)
            .mul_add2(e, 2.675_673E7)
            .mul_add2(e, 5.94594E6)
            .mul_add2(e, 2.7027E5)
    }
    #[inline(always)]
    fn a7(e: f64) -> f64 {
        e.mul_add2(-2.7027E5, -2.837_835E7)
            .mul_add2(e, -3.689_185_5E8)
            .mul_add2(e, -1.352_701_35E9)
            .mul_add2(e, -1.739_187_45E9)
            .mul_add2(e, -8.116_208_1E8)
            .mul_add2(e, -1.229_728_5E8)
            .mul_add2(e, -4.05405E6)
    }
    #[inline(always)]
    fn a8(e: f64) -> f64 {
        e.mul_add2(4.05405E6, 5.513_508E8)
            .mul_add2(e, 9.648_639E9)
            .mul_add2(e, 5.017_292_28E10)
            .mul_add2(e, 9.855_395_55E10)
            .mul_add2(e, 7.884_316_44E10)
            .mul_add2(e, 2.508_646_14E10)
            .mul_add2(e, 2.756_754E9)
            .mul_add2(e, 6.891_885E7)
    }
    #[inline(always)]
    fn a9(e: f64) -> f64 {
        e.mul_add2(-6.891_885E7, -1.178_512_335E10)
            .mul_add2(e, -2.671_294_626E11)
            .mul_add2(e, -1.869_906_238_2E12)
            .mul_add2(e, -5.209_024_520_7E12)
            .mul_add2(e, -6.366_585_525_3E12)
            .mul_add2(e, -3.472_683_013_8E12)
            .mul_add2(e, -8.013_883_878E11)
            .mul_add2(e, -6.678_236_565E10)
            .mul_add2(e, -1.309_458_15E9)
    }
    #[inline(always)]
    fn a10(e: f64) -> f64 {
        e.mul_add2(1.309_458_15E9, 2.749_862_115E11)
            .mul_add2(e, 7.837_107_027_75E12)
            .mul_add2(e, 7.105_643_705_16E13)
            .mul_add2(e, 2.664_616_389_435E14)
            .mul_add2(e, 4.618_668_408_354E14)
            .mul_add2(e, 3.848_890_340_295E14)
            .mul_add2(e, 1.522_637_936_82E14)
            .mul_add2(e, 2.664_616_389_435E13)
            .mul_add2(e, 1.741_579_339_5E12)
            .mul_add2(e, 2.749_862_115E10)
    }
    #[inline(always)]
    fn a11(e: f64) -> f64 {
        e.mul_add2(-2.749_862_115E10, -6.957_151_150_95E12)
            .mul_add2(e, -2.435_002_902_832_5E14)
            .mul_add2(e, -2.775_903_309_229_05E15)
            .mul_add2(e, -1.348_295_893_054_11E16)
            .mul_add2(e, -3.146_023_750_459_59E16)
            .mul_add2(e, -3.718_028_068_724_97E16)
            .mul_add2(e, -2.247_159_821_756_85E16)
            .mul_add2(e, -6.741_479_465_270_55E15)
            .mul_add2(e, -9.253_011_030_763_5E14)
            .mul_add2(e, -4.870_005_805_665E13)
            .mul_add2(e, -6.324_682_864_5E11)
    }
    #[inline(always)]
    fn a12(e: f64) -> f64 {
        e.mul_add2(6.324_682_864_5E11, 1.897_404_859_35E14)
            .mul_add2(e, 8.000_723_823_592_5E15)
            .mul_add2(e, 1.120_101_335_302_95E17)
            .mul_add2(e, 6.840_618_869_171_588E17)
            .mul_add2(e, 2.067_387_036_016_302E18)
            .mul_add2(e, 3.289_024_830_025_935E18)
            .mul_add2(e, 2.819_164_140_022_23E18)
            .mul_add2(e, 1.292_116_897_510_188_8E18)
            .mul_add2(e, 3.040_275_052_965_15E17)
            .mul_add2(e, 3.360_304_005_908_85E16)
            .mul_add2(e, 1.454_677_058_835E15)
            .mul_add2(e, 1.581_170_716_125E13)
    }
    #[inline(always)]
    fn a13(e: f64) -> f64 {
        e.mul_add2(-1.581_170_716_125E13, -5.549_909_213_598_75E15)
            .mul_add2(e, -2.774_954_606_799_375E17)
            .mul_add2(e, -4.680_423_436_801_613E18)
            .mul_add2(e, -3.510_317_577_601_209_5E19)
            .mul_add2(e, -1.333_920_679_488_459_6E20)
            .mul_add2(e, -2.748_685_036_521_674_2E20)
            .mul_add2(e, -3.171_559_657_525_009E20)
            .mul_add2(e, -2.061_513_777_391_255_6E20)
            .mul_add2(e, -7.410_670_441_602_553E19)
            .mul_add2(e, -1.404_127_031_040_483_7E19)
            .mul_add2(e, -1.276_479_119_127_712_5E18)
            .mul_add2(e, -4.624_924_344_665_625E16)
            .mul_add2(e, -4.269_160_933_537_5E14)
    }
    #[inline(always)]
    fn a14(e: f64) -> f64 {
        e.mul_add2(4.269_160_933_537_5E14, 1.733_279_339_016_225E17)
            .mul_add2(e, 1.013_968_413_324_491_6E19)
            .mul_add2(e, 2.027_936_826_648_983_3E20)
            .mul_add2(e, 1.832_385_775_507_831_3E21)
            .mul_add2(e, 8.551_133_619_036_546E21)
            .mul_add2(e, 2.215_520_983_114_014E22)
            .mul_add2(e, 3.311_108_282_456_109E22)
            .mul_add2(e, 2.897_219_747_149_095_4E22)
            .mul_add2(e, 1.477_013_988_742_676E22)
            .mul_add2(e, 4.275_566_809_518_273E21)
            .mul_add2(e, 6.663_221_001_846_66E20)
            .mul_add2(e, 5.069_842_066_622_458E19)
            .mul_add2(e, 1.559_951_405_114_602_5E18)
            .mul_add2(e, 1.238_056_670_725_875E16)
    }
    #[inline(always)]
    fn a15(e: f64) -> f64 {
        e.mul_add2(-1.238_056_670_725_875E16, -5.756_963_518_875_318E18)
            .mul_add2(e, -3.895_545_314_438_965_5E20)
            .mul_add2(e, -9.115_576_035_787_18E21)
            .mul_add2(e, -9.766_688_609_771_979E22)
            .mul_add2(e, -5.491_049_373_938_467_5E23)
            .mul_add2(e, -1.747_152_073_525_876E24)
            .mul_add2(e, -3.283_109_940_361_811E24)
            .mul_add2(e, -3.720_857_932_410_053E24)
            .mul_add2(e, -2.553_529_953_614_742E24)
            .mul_add2(e, -1.048_291_244_115_525_7E24)
            .mul_add2(e, -2.495_931_533_608_394_5E23)
            .mul_add2(e, -3.255_562_869_923_992_8E22)
            .mul_add2(e, -2.103_594_469_797_041_5E21)
            .mul_add2(e, -5.565_064_734_912_808E19)
            .mul_add2(e, -3.837_975_679_250_212_5E17)
    }
    #[inline(always)]
    fn a16(e: f64) -> f64 {
        e.mul_add2(3.837_975_679_250_212_5E17, 2.026_451_158_644_112E20)
            .mul_add2(e, 1.570_499_647_949_187E22)
            .mul_add2(e, 4.250_819_047_115_799E23)
            .mul_add2(e, 5.328_705_305_491_592E24)
            .mul_add2(e, 3.552_470_203_661_060_7E25)
            .mul_add2(e, 1.361_780_244_736_74E26)
            .mul_add2(e, 3.142_569_795_546_323E26)
            .mul_add2(e, 4.478_161_958_653_511E26)
            .mul_add2(e, 3.980_588_407_692_009_4E26)
            .mul_add2(e, 2.199_798_856_882_426E26)
            .mul_add2(e, 7.427_892_244_018_582E25)
            .mul_add2(e, 1.480_195_918_192_108_6E25)
            .mul_add2(e, 1.639_601_632_458_951_2E24)
            .mul_add2(e, 9.108_897_958_105_285E22)
            .mul_add2(e, 2.093_999_530_598_916E21)
            .mul_add2(e, 1.266_531_974_152_57E19)
    }
    const THRESHOLDS: [f64; 12] = [
        12.347, 12.958, 13.729, 14.718, 16.016, 17.769, 20.221, 23.816, 29.419, 38.93, 57.171,
        99.347,
    ];

    debug_assert!(h < ETA.abs().neg() && h < TAU + 0.5 - h + ETA);
    let e = (t / h).powi(2);
    let r = (h + t) * (h - t);
    let q = (h / r).powi(2);

    let idx = THRESHOLDS
        .iter()
        .position(|&threshold| -h - t + TAU + 0.5 < threshold)
        .unwrap_or(THRESHOLDS.len());
    let omega = if idx == THRESHOLDS.len() {
        q.mul_add2(a4(e), a3(e))
            .mul_add2(q, a2(e))
            .mul_add2(q, a1(e))
            .mul_add2(q, a0())
    } else {
        let mut omega_over_q = a16(e);
        if idx != 0 {
            omega_over_q = omega_over_q.mul_add2(q, a15(e));
            if idx != 1 {
                omega_over_q = omega_over_q.mul_add2(q, a14(e));
                if idx != 2 {
                    omega_over_q = omega_over_q.mul_add2(q, a13(e));
                    if idx != 3 {
                        omega_over_q = omega_over_q.mul_add2(q, a12(e));
                        if idx != 4 {
                            omega_over_q = omega_over_q.mul_add2(q, a11(e));
                            if idx != 5 {
                                omega_over_q = omega_over_q.mul_add2(q, a10(e));
                                if idx != 6 {
                                    omega_over_q = omega_over_q.mul_add2(q, a9(e));
                                    if idx != 7 {
                                        omega_over_q = omega_over_q.mul_add2(q, a8(e));
                                        if idx != 8 {
                                            omega_over_q = omega_over_q.mul_add2(q, a7(e));
                                            if idx != 9 {
                                                omega_over_q = omega_over_q.mul_add2(q, a6(e));
                                                if idx != 10 {
                                                    omega_over_q = omega_over_q.mul_add2(q, a5(e));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        q * omega_over_q
    };
    (t / r) * omega
}

#[inline(always)]
fn small_t_expansion_of_scaled_normalised_black<SpFn: SpecialFn>(h: f64, t: f64) -> f64 {
    #[inline(always)]
    fn b0(a: f64) -> f64 {
        2.0 * a
    }
    #[inline(always)]
    fn b1(a: f64, h2: f64) -> f64 {
        a.mul_add2(3.0 + h2, -1.0) / 3.0
    }
    #[inline(always)]
    fn b2(a: f64, h2: f64) -> f64 {
        h2.mul_add2(10.0 + h2, 15.0).mul_add2(a, -7.0 - h2) / 60.0
    }
    #[inline(always)]
    fn b3(a: f64, h2: f64) -> f64 {
        h2.mul_add2(21.0 + h2, 105.0)
            .mul_add2(h2, 105.0)
            .mul_add2(a, h2.mul_add2(-18.0 - h2, -57.0))
            / 2520.0
    }
    #[inline(always)]
    fn b4(a: f64, h2: f64) -> f64 {
        h2.mul_add2(36.0 + h2, 378.0)
            .mul_add2(h2, 1260.0)
            .mul_add2(h2, 945.0)
            .mul_add2(a, h2.mul_add2(-33.0 - h2, -285.0).mul_add2(h2, -561.0))
            / 181_440.0
    }
    #[inline(always)]
    fn b5(a: f64, h2: f64) -> f64 {
        h2.mul_add2(55.0 + h2, 990.0)
            .mul_add2(h2, 6930.0)
            .mul_add2(h2, 17325.0)
            .mul_add2(h2, 10395.0)
            .mul_add2(
                a,
                h2.mul_add2(-52.0 - h2, -840.0)
                    .mul_add2(h2, -4680.0)
                    .mul_add2(h2, -6555.0),
            )
            / 19_958_400.0
    }
    #[inline(always)]
    fn b6(a: f64, h2: f64) -> f64 {
        h2.mul_add2(78.0 + h2, 2145.0)
            .mul_add2(h2, 25740.0)
            .mul_add2(h2, 135_135.0)
            .mul_add2(h2, 270_270.0)
            .mul_add2(h2, 135_135.0)
            .mul_add2(
                a,
                h2.mul_add2(-75.0 - h2, -1926.0)
                    .mul_add2(h2, -20370.0)
                    .mul_add2(h2, -82845.0)
                    .mul_add2(h2, -89055.0),
            )
            / 3_113_510_400.0
    }
    let a = y_prime::<SpFn>(h);
    let h2 = h * h;
    let t2 = t * t;
    t2.mul_add2(b6(a, h2), b5(a, h2))
        .mul_add2(t2, b4(a, h2))
        .mul_add2(t2, b3(a, h2))
        .mul_add2(t2, b2(a, h2))
        .mul_add2(t2, b1(a, h2))
        .mul_add2(t2, b0(a))
        * t
}

#[inline(always)]
fn y_prime<SpFn: SpecialFn>(h: f64) -> f64 {
    // We copied the thresholds of -0.46875 and -4 from Cody.
    if h < -4.0 {
        // Nonlinear-Remez optimized minimax rational function of order (5,6) for g(w) := (Y'(h)/h²-1)/h² with w:=1/h².
        // The relative accuracy of Y'(h) ≈ w·(1+w·g(w)) is better than 9.8E-17 (in perfect arithmetic) on h in [-∞,-4] (i.e., on w in [0,1/16]).
        let w = (h * h).recip();
        w * (1.0 + y_prime_tail_expansion_rational_function_part(w))
    } else if h <= -0.46875 {
        // Remez-optimized minimax rational function of order (7,7) of relative accuracy better than 1.6E-16 (in perfect arithmetic) on h in [-4,-0.46875].
        h.mul_add2(8.459_243_640_658_06E-10, 4.276_659_783_590_871_4E-8)
            .mul_add2(-h, 3.071_739_227_491_390_3E-4)
            .mul_add2(-h, 5.545_521_007_735_379_5E-3)
            .mul_add2(-h, 4.565_090_035_135_299E-2)
            .mul_add2(-h, 2.218_084_473_657_601_4E-1)
            .mul_add2(-h, 6.191_144_987_969_411E-1)
            .mul_add2(-h, 1.000_000_000_059_431_8)
            / h.mul_add2(-3.082_202_041_792_714_7E-4, 5.529_045_357_693_659E-3)
                .mul_add2(-h, 4.676_254_890_319_496E-2)
                .mul_add2(-h, 2.367_770_140_309_464E-1)
                .mul_add2(-h, 7.657_648_983_658_903E-1)
                .mul_add2(-h, 1.568_549_723_607_765_2)
                .mul_add2(-h, 1.872_428_636_958_916_3)
                .mul_add2(-h, 1.0)
    } else {
        SpFn::erfcx(-FRAC_1_SQRT_2 * h).mul_add2(h * SQRT_PI_OVER_2, 1.0)
    }
}

#[inline(always)]
fn y_prime_tail_expansion_rational_function_part(w: f64) -> f64 {
    w.mul_add2(-6.681_824_903_261_685E4, -8.383_602_146_074_198E4)
        .mul_add2(w, -2.780_574_569_386_430_8E4)
        .mul_add2(w, -3.473_503_544_549_563_2E3)
        .mul_add2(w, -1.755_626_332_354_220_6E2)
        .mul_add2(w, -2.999_999_999_999_466)
        * w
        / w.mul_add2(6.928_651_867_980_375E4, 1.256_997_038_092_390_9E5)
            .mul_add2(w, 6.688_679_416_565_168E4)
            .mul_add2(w, 1.456_254_563_850_703_4E4)
            .mul_add2(w, 1.440_438_903_760_433_7E3)
            .mul_add2(w, 6.352_087_774_483_173_6E1)
            .mul_add2(w, 1.0)
}

const TAU: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

#[inline(always)]
fn normalised_black_with_optimal_use_of_codys_functions<SpFn: SpecialFn>(
    half_theta_x: f64,
    h: f64,
    t: f64,
) -> f64 {
    const CODYS_THRESHOLD: f64 = 0.46875;
    let q1 = -FRAC_1_SQRT_2 * (h + t);
    let q2 = -FRAC_1_SQRT_2 * (h - t);
    let two_b = if q1 < CODYS_THRESHOLD {
        if q2 < CODYS_THRESHOLD {
            half_theta_x
                .exp()
                .mul_add2(SpFn::erfc(q1), -((-half_theta_x).exp() * SpFn::erfc(q2)))
        } else {
            half_theta_x.exp().mul_add2(
                SpFn::erfc(q1),
                -((-0.5 * h.mul_add2(h, t * t)).exp() * SpFn::erfcx(q2)),
            )
        }
    } else if q2 < CODYS_THRESHOLD {
        (-0.5 * h.mul_add2(h, t * t))
            .exp()
            .mul_add2(SpFn::erfcx(q1), -((-half_theta_x).exp() * SpFn::erfc(q2)))
    } else {
        (-0.5 * h.mul_add2(h, t * t)).exp() * (SpFn::erfcx(q1) - SpFn::erfcx(q2))
    };
    (0.5 * two_b).max(0.0)
}

#[inline(always)]
pub fn normalised_vega(h: f64, t: f64) -> f64 {
    debug_assert!(t > 0.0);
    FRAC_1_SQRT_2_PI * (-0.5 * h.mul_add2(h, t * t)).exp()
}

#[inline(always)]
pub fn inv_normalised_vega(h: f64, t: f64) -> f64 {
    debug_assert!(t > 0.0);
    SQRT_2_PI * (0.5 * h.mul_add2(h, t * t)).exp()
}
#[inline(always)]
fn ln_normalised_vega(h: f64, t: f64) -> f64 {
    debug_assert!(h.abs() > 0.0);
    0.5f64.mul_add2(-h.mul_add2(h, t * t), -HALF_OF_LN_2_PI)
}

#[inline(always)]
pub fn black_input_unchecked<SpFn: SpecialFn, const IS_CALL: bool>(
    f: f64,
    k: f64,
    sigma: f64,
    t: f64,
) -> f64 {
    let s = sigma * t.sqrt();
    debug_assert!(s >= 0.0);
    if k == f {
        f * SpFn::erf((0.5 * FRAC_1_SQRT_2) * s)
    } else {
        (if IS_CALL { f - k } else { k - f }).max(0.0)
            + (if s <= 0.0 {
                0.0
            } else {
                debug_assert!(s > 0.0);
                let theta_x = (f / k).ln().abs().neg();
                debug_assert!(theta_x < 0.0);
                f.sqrt() * k.sqrt() * normalised_black::<SpFn>(0.5 * theta_x, theta_x / s, 0.5 * s)
            })
    }
}
