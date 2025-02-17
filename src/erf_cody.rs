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
const C: [f64; 9] = [
    0.564_188_496_988_670_1,
    8.883_149_794_388_377,
    66.119_190_637_141_63,
    298.635_138_197_400_1,
    881.952_221_241_769,
    1_712.047_612_634_070_7,
    2_051.078_377_826_071_6,
    1_230.339_354_797_997_2,
    2.153_115_354_744_038_3e-8,
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

const SQRPI: f64 = 0.564_189_583_547_756_3;
const THRESH: f64 = 0.46875;
const XINF: f64 = f64::MAX;
const XNEG: f64 = -26.628;
const XSMALL: f64 = 1.11e-16;
const XBIG: f64 = 26.543;
const XHUGE: f64 = 6.71e7;
const XMAX: f64 = 2.53e307;

#[inline(always)]
pub(crate) fn erfc_cody(x: f64) -> f64 {
    /* -------------------------------------------------------------------- */
    /* This subprogram computes approximate values for erfc(x). */
    /*   (see comments heading CALERF). */
    /*   Author/date: W. J. Cody, January 8, 1985 */
    /* -------------------------------------------------------------------- */
    #[inline(always)]
    fn finalize(y: f64) -> f64 {
        let ysq = (y * 16.0).trunc() / 16.0;
        let del = (y - ysq) * (y + ysq);
        (-ysq.powi(2)).exp() * (-del).exp()
    }
    let y = x.abs();
    let mut xden;
    let mut xnum;
    let mut result = 0.0;

    if y <= THRESH {
        let ysq = if y > XSMALL { y.powi(2) } else { 0.0 };
        xnum = A[4] * ysq;
        xden = ysq;

        for i in 0..3 {
            xnum = (xnum + A[i]) * ysq;
            xden = (xden + B[i]) * ysq;
        }
        return 1.0 - x * (xnum + A[3]) / (xden + B[3]);
    } else if y <= 4.0 {
        xnum = C[8] * y;
        xden = y;

        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        result = (xnum + C[7]) / (xden + D[7]);

        result *= finalize(y);
    } else if y >= XBIG {
        if x.is_sign_negative() {
            result = 2.0 - result;
        }
        return result;
    } else {
        let ysq = y.powi(2).recip();
        xnum = P[5] * ysq;
        xden = ysq;

        for i in 0..4 {
            xnum = (xnum + P[i]) * ysq;
            xden = (xden + Q[i]) * ysq;
        }
        result = ysq * (xnum + P[4]) / (xden + Q[4]);
        result = (SQRPI - result) / y;

        result *= finalize(y);
    }
    if x.is_sign_negative() {
        result = 2.0 - result;
    }
    result
}

#[inline(always)]
pub(crate) fn erfcx_cody(x: f64) -> f64 {
    /* ------------------------------------------------------------------ */
    /* This subprogram computes approximate values for exp(x*x) * erfc(x). */
    /*   (see comments heading CALERF). */
    /*   Author/date: W. J. Cody, March 30, 1987 */
    /* ------------------------------------------------------------------ */
    #[inline(always)]
    fn finalize(x: f64) -> f64 {
        let ysq = (x * 16.0).trunc() / 16.0;
        let del = (x - ysq) * (x + ysq);
        2.0_f64 * ysq.powi(2).exp() * del.exp()
    }
    if x < XNEG {
        return XINF;
    }

    let y = x.abs();
    let mut result = 0.0;

    if y <= THRESH {
        let ysq = if y > XSMALL { y.powi(2) } else { 0.0 };
        let mut xnum = A[4] * ysq;
        let mut xden = ysq;

        for i in 0..3 {
            xnum = (xnum + A[i]) * ysq;
            xden = (xden + B[i]) * ysq;
        }
        result = x * (xnum + A[3]) / (xden + B[3]);

        result = 1.0 - result;

        result *= ysq.exp();
        return result;
    } else if y <= 4.0 {
        let mut xnum = C[8] * y;
        let mut xden = y;

        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        result = (xnum + C[7]) / (xden + D[7]);
    } else if y >= XBIG {
        if y >= XMAX {
            if x.is_sign_negative() {
                if x < XNEG {
                    result = XINF;
                } else {
                    result = finalize(x) - result;
                }
            }
            return result;
        } else if y >= XHUGE {
            result = SQRPI / y;
            if x.is_sign_negative() {
                if x < XNEG {
                    result = XINF;
                } else {
                    result = finalize(x) - result;
                }
            }
            return result;
        }
    } else {
        let ysq = y.powi(2).recip();
        let mut xnum = P[5] * ysq;
        let mut xden = ysq;

        for i in 0..4 {
            xnum = (xnum + P[i]) * ysq;
            xden = (xden + Q[i]) * ysq;
        }
        result = ysq * (xnum + P[4]) / (xden + Q[4]);
        result = (SQRPI - result) / y;
    }
    if x.is_sign_negative() {
        result = finalize(x) - result;
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::erf_cody::{erfc_cody, erfcx_cody, THRESH, XBIG, XHUGE, XMAX, XNEG};

    #[test]
    fn calerf_1() {
        let x = erfc_cody(THRESH + f64::EPSILON);
        assert_eq!(x, 0.5073865267820618);
        let x = erfc_cody(THRESH - f64::EPSILON);
        assert_eq!(x, 0.5073865267820623);
        let x = erfc_cody(-THRESH - f64::EPSILON);
        assert_eq!(x, 1.4926134732179381);
        let x = erfc_cody(-THRESH + f64::EPSILON);
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

        let x = erfc_cody(XMAX + f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfc_cody(XMAX - f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfc_cody(-XMAX - f64::EPSILON);
        assert_eq!(x, 2.0);
        let x = erfc_cody(-XMAX + f64::EPSILON);
        assert_eq!(x, 2.0);

        let x = erfc_cody(XHUGE + f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfc_cody(XHUGE - f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfc_cody(-XHUGE - f64::EPSILON);
        assert_eq!(x, 2.0);
        let x = erfc_cody(-XHUGE + f64::EPSILON);
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
        let x = erfcx_cody(THRESH + f64::EPSILON);
        assert_eq!(x, 0.6320696892495559);
        let x = erfcx_cody(THRESH - f64::EPSILON);
        assert_eq!(x, 0.6320696892495563);
        let x = erfcx_cody(-THRESH - f64::EPSILON);
        assert_eq!(x, 1.8594024168714227);
        let x = erfcx_cody(-THRESH + f64::EPSILON);
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
        assert_eq!(x, 0.0);
        let x = erfcx_cody(XBIG - f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfcx_cody(-XBIG - f64::EPSILON);
        assert_eq!(x, 1.8831722547514706e306);
        let x = erfcx_cody(-XBIG + f64::EPSILON);
        assert_eq!(x, 1.8831722547514706e306);

        let x = erfcx_cody(XMAX + f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfcx_cody(XMAX - f64::EPSILON);
        assert_eq!(x, 0.0);
        let x = erfcx_cody(-XMAX - f64::EPSILON);
        assert_eq!(x, 1.7976931348623157e308);
        let x = erfcx_cody(-XMAX + f64::EPSILON);
        assert_eq!(x, 1.7976931348623157e308);

        let x = erfcx_cody(XHUGE + f64::EPSILON);
        assert_eq!(x, 8.408190514869691e-9);
        let x = erfcx_cody(XHUGE - f64::EPSILON);
        assert_eq!(x, 8.408190514869691e-9);
        let x = erfcx_cody(-XHUGE - f64::EPSILON);
        assert_eq!(x, 1.7976931348623157e308);
        let x = erfcx_cody(-XHUGE + f64::EPSILON);
        assert_eq!(x, 1.7976931348623157e308);

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
