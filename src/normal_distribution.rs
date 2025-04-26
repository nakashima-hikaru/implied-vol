use crate::erf_cody::erfc_cody;
use std::f64::consts::FRAC_1_SQRT_2;
use std::ops::Neg;
use crate::MulAdd;

const NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD: f64 = -10.0;
const NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD: f64 = -67108864.0;
// 1.0 / f64::sqrt(f64::EPSILON);

const FRAC_SQRT_2_PI: f64 = f64::from_bits(0x3fd9884533d43651);

#[inline(always)]
pub(crate) fn norm_pdf(x: f64) -> f64 {
    FRAC_SQRT_2_PI * (-0.5 * x.powi(2)).exp()
}

#[inline(always)]
pub(crate) fn norm_cdf(z: f64) -> f64 {
    if z <= NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD {
        let mut sum = 1.0;
        if z >= NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD {
            let zsqr = z.powi(2);
            let mut i: usize = 4;
            let mut g = 1.0;
            let mut a = f64::MAX;
            loop {
                let lasta = a;
                let x = (i - 3) as f64 / zsqr;
                let y = x * ((i - 1) as f64 / zsqr);
                a = g * (x - y);
                sum -= a;
                g *= y;
                i += 4;
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
    assert!(!(u.abs() > U_MAX));
    let s = U_MAX2 - u.powi(2);
    u * (s.mul_add2(
        s.mul_add2(
            s.mul_add2(
                s.mul_add2(
                    s.mul_add2(
                        s.mul_add2(-7.58939881401259242, 1.34233243502653864E2),
                        6.90489242061408612E2,
                    ),
                    7.4997781456657924E2,
                ),
                3.01870541922933937E2,
            ),
            5.0260572167303103E1,
        ),
        2.92958954698308805,
    ) / s.mul_add2(
        s.mul_add2(
            s.mul_add2(
                s.mul_add2(
                    s.mul_add2(1.79227008508102628E2, 4.79123914509756757E2),
                    3.86821208540417453E2,
                ),
                1.29404120448755281E2,
            ),
            1.8918538074574598E1,
        ),
        1.0,
    ))
}

#[inline(always)]
fn inverse_norm_cdf_for_low_probabilities(p: f64) -> f64 {
    assert!(!(p > 0.15865525393146));
    let r = p.ln().neg().sqrt();
    if r < 6.7 {
        if r < 3.41 {
            if r < 2.05 {
                // Branch I: Accuracy better than 7.6E-17
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(
                                r.mul_add2(-13.054072340494093704, -83.383894003636969722),
                                -74.594687726045926821,
                            ),
                            65.451292110261454609,
                        ),
                        47.170590600740689449,
                    ),
                    3.691562302945566191,
                ) / r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(
                                r.mul_add2(0.00018295174852053530579, 9.2216887978737432303),
                                59.270122556046077717,
                            ),
                            71.813812182579255459,
                        ),
                        20.837211328697753726,
                    ),
                    1.0,
                )
            } else {
                // Branch II: Accuracy better than 9.4E-17
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(
                                r.mul_add2(-1.2013147879435525574, -10.05916339568646151),
                                -18.1254427791789183,
                            ),
                            0.68397370256591532878,
                        ),
                        14.49177828689122096,
                    ),
                    3.2340179116317970288,
                ) / r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(
                                r.mul_add2(0.000010957576098829595323, 0.84884892199149255469),
                                7.1369811056109768745,
                            ),
                            14.656370665176799712,
                        ),
                        8.8820931773304337525,
                    ),
                    1.0,
                )
            }
        } else {
            // Branch III: Accuracy better than 9.1E-17
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(-0.15414319494013597492, -2.8699061335882526744),
                            -11.070534689309368061,
                        ),
                        -5.1633929115525534628,
                    ),
                    9.9483724317036560676,
                ),
                3.1252235780087584807,
            ) / r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(
                                r.mul_add2(0.00000013565983564441297634, 0.10897972234131828901),
                                0.0000000030848093570966787291,
                            ),
                            2.0307076064309043613,
                        ),
                        8.1086341122361532407,
                    ),
                    7.076769154309171622,
                ),
                1.0,
            )
        }
    } else {
        if r < 12.9 {
            // Branch IV: Accuracy better than 9E-17
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(
                                r.mul_add2(-0.01612303318390145052, -0.47595169546783216436),
                                -2.9644251353150605663,
                            ),
                            -0.065127593753781672404,
                        ),
                        -3.688196041019692267,
                    ),
                    2.250881388987032271,
                ),
                2.6161264950897283681,
            ) / r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(0.0000000030848093570966787291, 0.011400087282177594359),
                            0.33663746405626400164,
                        ),
                        2.1282030272153188194,
                    ),
                    3.2517455169035921495,
                ),
                1.0,
            )
        } else {
            // Branch V: Accuracy better than 9.5E-17
            r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(-0.0010566357727202585402, -0.065127593753781672404),
                            -0.86385181219213758847,
                        ),
                        -2.5894451568465728432,
                    ),
                    -0.042799650734502094297,
                ),
                2.3226849047872302955,
            ) / r.mul_add2(
                r.mul_add2(
                    r.mul_add2(
                        r.mul_add2(
                            r.mul_add2(0.000000000023135343206304888, 0.0007471447992167225483),
                            0.046054974512474443189,
                        ),
                        0.61320841329197493341,
                    ),
                    1.9361316119254412206,
                ),
                1.0,
            )
        }
    }
}

#[inline(always)]
pub(crate) fn inverse_norm_cdf(p: f64) -> f64 {
    let u = p - 0.5;
    if u.abs() < U_MAX {
        return inverse_norm_cdfm_half_for_midrange_probabilities(u);
    }
    if u.is_sign_positive() {
        -inverse_norm_cdf_for_low_probabilities(1.0 - p)
    } else {
        inverse_norm_cdf_for_low_probabilities(p)
    }
}
