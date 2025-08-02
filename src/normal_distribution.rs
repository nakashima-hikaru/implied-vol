use crate::MulAdd;
use std::f64::consts::FRAC_1_SQRT_2;
use std::ops::Neg;

const FRAC_SQRT_2_PI: f64 = f64::from_bits(0x3fd9884533d43651);

#[inline(always)]
pub(crate) fn norm_pdf(x: f64) -> f64 {
    FRAC_SQRT_2_PI * (-0.5 * x * x).exp()
}

#[cfg(feature = "normal-distribution")]
#[inline(always)]
pub(crate) fn norm_cdf(z: f64) -> f64 {
    use crate::erf_cody::erfc_cody;
    const NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD: f64 = -10.0;
    const NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD: f64 = -67108864.0;
    if z <= NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD {
        let mut sum = 1.0;
        if z >= NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD {
            let zsqr = z * z;
            let mut i = 4.0_f64;
            let mut g = 1.0;
            let mut a = f64::MAX;
            loop {
                let lasta = a;
                let x = (i - 3.0) / zsqr;
                let y = x * ((i - 1.0) / zsqr);
                a = g * (x - y);
                sum -= a;
                g *= y;
                i += 4.0;
                a = a.abs();
                if !(lasta > a && a >= (sum * f64::EPSILON).abs()) {
                    break;
                }
            }
        }
        return -norm_pdf(z) * sum / z;
    }
    0.5 * erfc_cody(-z * FRAC_1_SQRT_2)
}

const U_MAX: f64 = 0.3413447460685429;
const U_MAX2: f64 = U_MAX * U_MAX;
#[inline(always)]
fn inverse_norm_cdfm_half_for_midrange_probabilities(u: f64) -> f64 {
    assert!(u.abs() <= U_MAX);
    let s = U_MAX2 - u * u;
    u * (s.mul_add2(
        s.mul_add2(
            s.mul_add2(
                s.mul_add2(
                    s.mul_add2(
                        s.mul_add2(-7.589_398_814_012_592_5, 1.342_332_435_026_538_6E2),
                        6.904_892_420_614_086E2,
                    ),
                    7.499_778_145_665_792E2,
                ),
                3.018_705_419_229_339E2,
            ),
            5.026_057_216_730_31E1,
        ),
        2.929_589_546_983_088,
    ) / s.mul_add2(
        s.mul_add2(
            s.mul_add2(
                s.mul_add2(
                    s.mul_add2(1.792_270_085_081_026E2, 4.791_239_145_097_567_3E2),
                    3.868_212_085_404_174_4E2,
                ),
                1.294_041_204_487_553E2,
            ),
            1.891_853_807_457_46E1,
        ),
        1.0,
    ))
}

#[inline(always)]
fn inverse_norm_cdf_for_low_probabilities(p: f64) -> f64 {
    assert!(p <= 0.15865525393146);
    let r = p.ln().neg().sqrt();
    if r < 2.05 {
        // Branch I: Accuracy better than 7.6E-17
        r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(-13.054_072_340_494_093, -83.383_894_003_636_97),
                        -74.594_687_726_045_93,
                    ),
                    65.451_292_110_261_45,
                ),
                47.170_590_600_740_69,
            ),
            3.691_562_302_945_566,
        ) / r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(0.000_182_951_748_520_535_3, 9.221_688_797_873_743),
                        59.270_122_556_046_076,
                    ),
                    71.813_812_182_579_26,
                ),
                20.837_211_328_697_755,
            ),
            1.0,
        )
    } else if r < 3.41 {
        // Branch II: Accuracy better than 9.4E-17
        r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(-1.201_314_787_943_552_6, -10.059_163_395_686_461),
                        -18.125_442_779_178_92,
                    ),
                    0.683_973_702_565_915_3,
                ),
                14.491_778_286_891_22,
            ),
            3.234_017_911_631_797,
        ) / r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(0.000_010_957_576_098_829_594, 0.848_848_921_991_492_5),
                        7.136_981_105_610_977,
                    ),
                    14.656_370_665_176_8,
                ),
                8.882_093_177_330_434,
            ),
            1.0,
        )
    } else if r < 6.7 {
        // Branch III: Accuracy better than 9.1E-17
        r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(-1.541_431_949_401_359_8E-1, -2.869_906_133_588_252_8),
                        -1.107_053_468_930_936_7E1,
                    ),
                    -5.163_392_911_552_553,
                ),
                9.948_372_431_703_657,
            ),
            3.125_223_578_008_758_3,
        ) / r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(1.356_598_356_444_129_7E-7, 1.089_797_223_413_183E-1),
                        2.030_707_606_430_904,
                    ),
                    8.108_634_112_236_153,
                ),
                7.076_769_154_309_171,
            ),
            1.,
        )
    } else if r < 12.9 {
        // Branch IV: Accuracy better than 9E-17
        r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(-0.016_123_033_183_901_45, -0.475_951_695_467_832_17),
                            -2.964_425_135_315_060_4,
                        ),
                        -0.065_127_593_753_781_67,
                    ),
                    -3.688_196_041_019_692,
                ),
                2.250_881_388_987_032,
            ),
            2.616_126_495_089_728_3,
        ) / r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            0.000_000_003_084_809_357_096_678_6,
                            0.011_400_087_282_177_594,
                        ),
                        0.336_637_464_056_264,
                    ),
                    2.128_203_027_215_319,
                ),
                3.251_745_516_903_592,
            ),
            1.0,
        )
    } else {
        // Branch V: Accuracy better than 9.5E-17
        r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(-0.001_056_635_772_720_258_4, -0.065_127_593_753_781_67),
                        -0.863_851_812_192_137_6,
                    ),
                    -2.589_445_156_846_573,
                ),
                -0.042_799_650_734_502_09,
            ),
            2.322_684_904_787_23,
        ) / r.mul_add2(
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(0.000000000023135343206304888, 0.000_747_144_799_216_722_6),
                        0.046_054_974_512_474_44,
                    ),
                    0.613_208_413_291_975,
                ),
                1.936_131_611_925_441_3,
            ),
            1.0,
        )
    }
}

#[inline(always)]
pub(crate) fn inverse_norm_cdf(p: f64) -> f64 {
    let u = p - 0.5;
    if u.abs() < U_MAX {
        return inverse_norm_cdfm_half_for_midrange_probabilities(u);
    }
    if u > 0.0 {
        -inverse_norm_cdf_for_low_probabilities(1.0 - p)
    } else {
        inverse_norm_cdf_for_low_probabilities(p)
    }
}

#[inline(always)]
pub(crate) fn erfinv(e: f64) -> f64 {
    if e.abs() < 2.0 * U_MAX {
        inverse_norm_cdfm_half_for_midrange_probabilities(0.5 * e) * FRAC_1_SQRT_2
    } else {
        (if e < 0.0 {
            inverse_norm_cdf_for_low_probabilities(e.mul_add2(0.5, 0.5))
        } else {
            -inverse_norm_cdf_for_low_probabilities(e.mul_add2(-0.5, 0.5))
        }) * FRAC_1_SQRT_2
    }
}
