use crate::MulAdd;
use std::f64::consts::FRAC_1_SQRT_2;
use std::ops::Neg;

#[inline(always)]
fn ab(z: f64) -> f64 {
    const A: [f64; 5] = [
        3.1611237438705656,
        113.864_154_151_050_16,
        377.485_237_685_302,
        3_209.377_589_138_469_4,
        0.185_777_706_184_603_15,
    ];
    const B: [f64; 4] = [
        23.601_290_952_344_122,
        244.024_637_934_444_17,
        1_282.616_526_077_372_3,
        2_844.236_833_439_171,
    ];
    z.mul_add2(A[4], A[0])
        .mul_add2(z, A[1])
        .mul_add2(z, A[2])
        .mul_add2(z, A[3])
        / z.mul_add2(1.0, B[0])
            .mul_add2(z, B[1])
            .mul_add2(z, B[2])
            .mul_add2(z, B[3])
}

#[inline(always)]
fn cd(y: f64) -> f64 {
    const C: [f64; 9] = [
        0.564_188_496_988_670_1,
        8.883_149_794_388_377,
        66.119_190_637_141_63,
        298.635_138_197_400_1,
        881.952_221_241_769,
        1_712.047_612_634_070_7,
        2_051.078_377_826_071_6,
        1_230.339_354_797_997_2,
        2.153_115_354_744_038_3E-8,
    ];
    const D: [f64; 8] = [
        15.744_926_110_709_835,
        117.693_950_891_312_5,
        537.181_101_862_009_9,
        1_621.389_574_566_690_3,
        3_290.799_235_733_459_7,
        4_362.619_090_143_247,
        3_439.367_674_143_721_6,
        1_230.339_354_803_749_5,
    ];
    C[8].mul_add2(y, C[0])
        .mul_add2(y, C[1])
        .mul_add2(y, C[2])
        .mul_add2(y, C[3])
        .mul_add2(y, C[4])
        .mul_add2(y, C[5])
        .mul_add2(y, C[6])
        .mul_add2(y, C[7])
        / (y + D[0])
            .mul_add2(y, D[1])
            .mul_add2(y, D[2])
            .mul_add2(y, D[3])
            .mul_add2(y, D[4])
            .mul_add2(y, D[5])
            .mul_add2(y, D[6])
            .mul_add2(y, D[7])
}

#[inline(always)]
fn pq(z: f64) -> f64 {
    const P: [f64; 6] = [
        0.305_326_634_961_232_36,
        0.360_344_899_949_804_45,
        0.125_781_726_111_229_26,
        0.016_083_785_148_742_275,
        6.587_491_615_298_378e-4,
        0.016_315_387_137_302_097,
    ];
    const Q: [f64; 5] = [
        2.568_520_192_289_822,
        1.872_952_849_923_460_4,
        0.527_905_102_951_428_5,
        0.060_518_341_312_441_32,
        0.002_335_204_976_268_691_8,
    ];
    z * (P[5]
        .mul_add2(z, P[0])
        .mul_add2(z, P[1])
        .mul_add2(z, P[2])
        .mul_add2(z, P[3])
        .mul_add2(z, P[4]))
        / ((z + Q[0])
            .mul_add2(z, Q[1])
            .mul_add2(z, Q[2])
            .mul_add2(z, Q[3])
            .mul_add2(z, Q[4]))
}

#[inline(always)]
fn smoothened_exponential_of_negative_square(y: f64) -> f64 {
    let y_tilde = (y * 16.0).trunc() / 16.0;
    (y_tilde * y_tilde).neg().exp() * (-(y - y_tilde) * (y + y_tilde)).exp()
}

#[inline(always)]
fn smoothened_exponential_of_positive_square(x: f64) -> f64 {
    let x_tilde = (x * 16.0).trunc() / 16.0;
    (x_tilde * x_tilde).exp() * ((x - x_tilde) * (x + x_tilde)).exp()
}

const THRESHOLD: f64 = 0.46875;
const XNEG: f64 = -26.628;
const XBIG: f64 = 26.543;

const ONE_OVER_SQRT_PI: f64 = 0.564_189_583_547_756_3;

#[inline(always)]
pub(crate) fn erf_cody(x: f64) -> f64 {
    let y = x.abs();
    if y <= THRESHOLD {
        //       |x| <= 0.46875
        return x * ab(y * y);
    }
    // Compute erfc(|x|)
    let erfc_abs_x = if y >= XBIG {
        0.0 // when |x| ≥ 26.543
    } else if y <= 4.0 {
        cd(y) // when 0.46875 < |x| <= 4.0
    } else {
        (ONE_OVER_SQRT_PI - pq((y * y).recip())) / y // when 4.0 < |x| < 26.543
    } * smoothened_exponential_of_negative_square(y);
    if x < 0.0 {
        erfc_abs_x - 1.0
    } else {
        1.0 - erfc_abs_x
    }
}

#[inline(always)]
pub(crate) fn erfc_cody(x: f64) -> f64 {
    let y = x.abs();
    if y <= THRESHOLD {
        return ab(y * y).neg().mul_add2(x, 1.0);
    }
    let erfc_abs_x = if y >= XBIG {
        0.0
    } else {
        (if y <= 4.0 {
            cd(y)
        } else {
            (FRAC_1_SQRT_2 - pq((y * y).recip())) / y
        }) * smoothened_exponential_of_negative_square(y)
    };

    if x < 0.0 {
        2.0 - erfc_abs_x
    } else {
        erfc_abs_x
    }
}

#[inline(always)]
fn erfcx_cody_above_threshold(y: f64) -> f64 {
    assert!(y > THRESHOLD, "{y} exceeds threshold");
    if y <= 4.0 {
        cd(y)
    } else {
        (FRAC_1_SQRT_2 - pq((y * y).recip())) / y
    }
}

#[inline(always)]
pub(crate) fn erfcx_cody(x: f64) -> f64 {
    let y = x.abs();
    if y <= THRESHOLD {
        let z = y * y;
        return z.exp() * (ab(z).neg().mul_add2(x, 1.0));
    }
    if x < XNEG {
        return f64::MAX;
    }
    let result = erfcx_cody_above_threshold(y);
    if x < 0.0 {
        let expx2 = smoothened_exponential_of_positive_square(x);
        return (expx2 + expx2) - result;
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::erf_cody::{THRESHOLD, XBIG, XNEG, erfc_cody, erfcx_cody};

    #[test]
    fn calerf_1() {
        let x = erfc_cody(THRESHOLD + f64::EPSILON);
        assert_eq!(x, 0.5073865267820618);
        let x = erfc_cody(THRESHOLD - f64::EPSILON);
        assert_eq!(x, 0.5073865267820623);
        let x = erfc_cody(-THRESHOLD - f64::EPSILON);
        assert_eq!(x, 1.4926134732179381);
        let x = erfc_cody(-THRESHOLD + f64::EPSILON);
        assert_eq!(x, 1.4926134732179377);

        let x = erfc_cody(4.0 + f64::EPSILON);
        assert_eq!(x, 1.541725790028002e-8);
        let x = erfc_cody(4.0 - f64::EPSILON);
        assert_eq!(x, 1.541725790028002e-8);
        let x = erfc_cody(-4.0 - f64::EPSILON);
        assert_eq!(x, 1.999999984582742);
        let x = erfc_cody(-4.0 + f64::EPSILON);
        assert_eq!(x, 1.999999984582742);

        let x = erfc_cody(XBIG + f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfc_cody(XBIG - f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfc_cody(-XBIG - f64::EPSILON);
        assert_eq!(x, 2.0);
        let x = erfc_cody(-XBIG + f64::EPSILON);
        assert_eq!(x, 2.0);

        let x = erfc_cody(0.0 + f64::EPSILON);
        assert_eq!(x, 0.9999999999999998);
        let x = erfc_cody(0.0 - f64::EPSILON);
        assert_eq!(x, 1.0000000000000002);

        let x = erfc_cody(XNEG + f64::EPSILON);
        assert_eq!(x, 2.0);
        let x = erfc_cody(XNEG - f64::EPSILON);
        assert_eq!(x, 2.0);
    }

    #[test]
    fn calerf_2() {
        let x = erfcx_cody(THRESHOLD + f64::EPSILON);
        assert_eq!(x, 0.6320696892495559);
        let x = erfcx_cody(THRESHOLD - f64::EPSILON);
        assert_eq!(x, 0.6320696892495563);
        let x = erfcx_cody(-THRESHOLD - f64::EPSILON);
        assert_eq!(x, 1.8594024168714227);
        let x = erfcx_cody(-THRESHOLD + f64::EPSILON);
        assert_eq!(x, 1.8594024168714214);

        let x = erfcx_cody(4.0 + f64::EPSILON);
        assert_eq!(x, 0.1369994576250614);
        let x = erfcx_cody(4.0 - f64::EPSILON);
        assert_eq!(x, 0.1369994576250614);
        let x = erfcx_cody(-4.0 - f64::EPSILON);
        assert_eq!(x, 17772220.904016286);
        let x = erfcx_cody(-4.0 + f64::EPSILON);
        assert_eq!(x, 17772220.904016286);

        let x = erfcx_cody(XBIG + f64::EPSILON);
        assert_eq!(x, 0.026624994527838793);
        let x = erfcx_cody(XBIG - f64::EPSILON);
        assert_eq!(x, 0.026624994527838793);
        let x = erfcx_cody(-XBIG - f64::EPSILON);
        assert_eq!(x, 1.8831722547514706e306);
        let x = erfcx_cody(-XBIG + f64::EPSILON);
        assert_eq!(x, 1.8831722547514706e306);

        let x = erfcx_cody(0.0 + f64::EPSILON);
        assert_eq!(x, 0.9999999999999998);
        let x = erfcx_cody(0.0 - f64::EPSILON);
        assert_eq!(x, 1.0000000000000002);

        let x = erfcx_cody(XNEG + f64::EPSILON);
        assert_eq!(x, 1.728618506590026e308);
        let x = erfcx_cody(XNEG - f64::EPSILON);
        assert_eq!(x, 1.728618506590026e308);
    }
}
