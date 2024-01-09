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

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;
const SQRPI: f64 = 0.564_189_583_547_756_3;
const THRESH: f64 = 0.46875;
const SIXTEN: f64 = 16.0;
const XINF: f64 = 1.79e308;
const XNEG: f64 = -26.628;
const XSMALL: f64 = 1.11e-16;
const XBIG: f64 = 26.543;
const XHUGE: f64 = 6.71e7;
const XMAX: f64 = 2.53e307;

/// This function calculates the complementary error function (erfc) using the algorithm outlined in the paper:
/// "Approximating the erfc and related functions" by W.J. Cody, Math. Comp., 1969, Volume 23, No 107, Pages 631-637.
///
/// # Arguments
/// * `x` - The value for which to calculate the complementary error function.
/// * `jint` - The indicator for the type of calculation to perform.
///     - 0: Calculate the complementary error function.
///     - 1: Calculate the scaled complementary error function.
///     - 2: Calculate the scaled complementary error function multiplied by `exp(x^2)`.
///
/// # Returns
/// The calculated complementary error function value.
pub fn calerf(x: f64, jint: i32) -> f64 {
    let y = x.abs();
    let mut ysq = ZERO;
    let mut xden;
    let mut xnum;
    let mut result = ZERO;

    if y <= THRESH {
        if y > XSMALL {
            ysq = y * y;
        }
        xnum = A[4] * ysq;
        xden = ysq;

        for i in 0..3 {
            xnum = (xnum + A[i]) * ysq;
            xden = (xden + B[i]) * ysq;
        }
        result = x * (xnum + A[3]) / (xden + B[3]);

        if jint != 0 {
            result = ONE - result;
        }

        if jint == 2 {
            result *= ysq.exp();
        }
        return result;
    } else if y <= FOUR {
        xnum = C[8] * y;
        xden = y;

        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        result = (xnum + C[7]) / (xden + D[7]);

        if jint != 2 {
            ysq = (y * SIXTEN).trunc() / SIXTEN;
            let del = (y - ysq) * (y + ysq);
            result *= (-ysq * ysq).exp() * (-del).exp();
        }
    } else if y >= XBIG {
        if jint != 2 || y >= XMAX {
            return fix_up(x, jint, result);
        } else if y >= XHUGE {
            result = SQRPI / y;
            return fix_up(x, jint, result);
        }
    } else {
        ysq = ONE / (y * y);
        xnum = P[5] * ysq;
        xden = ysq;

        for i in 0..4 {
            xnum = (xnum + P[i]) * ysq;
            xden = (xden + Q[i]) * ysq;
        }
        result = ysq * (xnum + P[4]) / (xden + Q[4]);
        result = (SQRPI - result) / y;

        if jint != 2 {
            ysq = (y * SIXTEN).trunc() / SIXTEN;
            let del = (y - ysq) * (y + ysq);
            result *= (-ysq * ysq).exp() * (-del).exp();
        }
    }
    fix_up(x, jint, result)
}

/// Calculates the result based on the given parameters.
///
/// # Arguments
///
/// * `x` - A `f64` value representing a number.
/// * `jint` - An `i32` value representing the integer part of a number.
/// * `result` - A `f64` value representing the current result.
///
/// # Returns
///
/// Returns a `f64` value representing the calculated result.
///
/// # Remarks
///
/// This function applies different calculations to `result` based on the value of `jint`
/// and the sign of `x`. The calculated result is then returned.
///
/// - If `jint` is 0:
///   - `result` is adjusted to be in the range [0.5, 1) by using the formula `(HALF - result) + HALF`.
///   - If `x` is less than 0, the result is negated.
///
/// - If `jint` is 1:
///   - If `x` is less than 0, the result is adjusted to be in the range [1, 2) by using the formula `TWO - result`.
///
/// - Otherwise:
///   - If `x` is less than 0:
///     - If `x` is less than XNEG, the result is set to XINF.
///     - Otherwise, the result is calculated based on the given formula and stored in `result`.
///
/// The final `result` value is then returned.
fn fix_up(x: f64, jint: i32, result: f64) -> f64 {
    let mut result = result;
    if jint == 0 {
        result = (HALF - result) + HALF;
        if x < ZERO {
            result = -result;
        }
    } else if jint == 1 {
        if x < ZERO {
            result = TWO - result;
        }
    } else if x < ZERO {
        if x < XNEG {
            result = XINF;
        } else {
            let ysq = (x * SIXTEN).trunc() / SIXTEN;
            let del = (x - ysq) * (x + ysq);
            let y = (ysq * ysq).exp() * del.exp();
            result = (y + y) - result;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::erf_cody::{calerf, FOUR, THRESH, XBIG, XHUGE, XMAX, XNEG, ZERO};

    #[test]
    fn calerf_0() {
        let x = calerf(THRESH + f64::EPSILON, 0);
        assert_eq!(x, 0.49261347321793825);
        let x = calerf(THRESH - f64::EPSILON, 0);
        assert_eq!(x, 0.49261347321793775);
        let x = calerf(-THRESH - f64::EPSILON, 0);
        assert_eq!(x, -0.49261347321793825);
        let x = calerf(-THRESH + f64::EPSILON, 0);
        assert_eq!(x, -0.49261347321793775);

        let x = calerf(FOUR + f64::EPSILON, 0);
        assert_eq!(x, 0.9999999845827421);
        let x = calerf(FOUR - f64::EPSILON, 0);
        assert_eq!(x, 0.9999999845827421);
        let x = calerf(-FOUR - f64::EPSILON, 0);
        assert_eq!(x, -0.9999999845827421);
        let x = calerf(-FOUR + f64::EPSILON, 0);
        assert_eq!(x, -0.9999999845827421);

        let x = calerf(XBIG + f64::EPSILON, 0);
        assert_eq!(x, 1.0);
        let x = calerf(XBIG - f64::EPSILON, 0);
        assert_eq!(x, 1.0);
        let x = calerf(-XBIG - f64::EPSILON, 0);
        assert_eq!(x, -1.0);
        let x = calerf(-XBIG + f64::EPSILON, 0);
        assert_eq!(x, -1.0);

        let x = calerf(XMAX + f64::EPSILON, 0);
        assert_eq!(x, 1.0);
        let x = calerf(XMAX - f64::EPSILON, 0);
        assert_eq!(x, 1.0);
        let x = calerf(-XMAX - f64::EPSILON, 0);
        assert_eq!(x, -1.0);
        let x = calerf(-XMAX + f64::EPSILON, 0);
        assert_eq!(x, -1.0);

        let x = calerf(XHUGE + f64::EPSILON, 0);
        assert_eq!(x, 1.0);
        let x = calerf(XHUGE - f64::EPSILON, 0);
        assert_eq!(x, 1.0);
        let x = calerf(-XHUGE - f64::EPSILON, 0);
        assert_eq!(x, -1.0);
        let x = calerf(-XHUGE + f64::EPSILON, 0);
        assert_eq!(x, -1.0);

        let x = calerf(ZERO + f64::EPSILON, 0);
        assert_eq!(x, 2.5055050636335897e-16);
        let x = calerf(ZERO - f64::EPSILON, 0);
        assert_eq!(x, -2.5055050636335897e-16);

        let x = calerf(XNEG + f64::EPSILON, 0);
        assert_eq!(x, -1.0);
        let x = calerf(XNEG - f64::EPSILON, 0);
        assert_eq!(x, -1.0);
    }

    #[test]
    fn calerf_1() {
        let x = calerf(THRESH + f64::EPSILON, 1);
        assert_eq!(x, 0.5073865267820618);
        let x = calerf(THRESH - f64::EPSILON, 1);
        assert_eq!(x, 0.5073865267820623);
        let x = calerf(-THRESH - f64::EPSILON, 1);
        assert_eq!(x, 1.4926134732179381);
        let x = calerf(-THRESH + f64::EPSILON, 1);
        assert_eq!(x, 1.4926134732179377);

        let x = calerf(FOUR + f64::EPSILON, 1);
        assert_eq!(x, 1.541725790028002e-8);
        let x = calerf(FOUR - f64::EPSILON, 1);
        assert_eq!(x, 1.541725790028002e-8);
        let x = calerf(-FOUR - f64::EPSILON, 1);
        assert_eq!(x, 1.999999984582742);
        let x = calerf(-FOUR + f64::EPSILON, 1);
        assert_eq!(x, 1.999999984582742);

        let x = calerf(XBIG + f64::EPSILON, 1);
        assert_eq!(x, 0.0);
        let x = calerf(XBIG - f64::EPSILON, 1);
        assert_eq!(x, 0.0);
        let x = calerf(-XBIG - f64::EPSILON, 1);
        assert_eq!(x, 2.0);
        let x = calerf(-XBIG + f64::EPSILON, 1);
        assert_eq!(x, 2.0);

        let x = calerf(XMAX + f64::EPSILON, 1);
        assert_eq!(x, 0.0);
        let x = calerf(XMAX - f64::EPSILON, 1);
        assert_eq!(x, 0.0);
        let x = calerf(-XMAX - f64::EPSILON, 1);
        assert_eq!(x, 2.0);
        let x = calerf(-XMAX + f64::EPSILON, 1);
        assert_eq!(x, 2.0);

        let x = calerf(XHUGE + f64::EPSILON, 1);
        assert_eq!(x, 0.0);
        let x = calerf(XHUGE - f64::EPSILON, 1);
        assert_eq!(x, 0.0);
        let x = calerf(-XHUGE - f64::EPSILON, 1);
        assert_eq!(x, 2.0);
        let x = calerf(-XHUGE + f64::EPSILON, 1);
        assert_eq!(x, 2.0);

        let x = calerf(ZERO + f64::EPSILON, 1);
        assert_eq!(x, 0.9999999999999998);
        let x = calerf(ZERO - f64::EPSILON, 1);
        assert_eq!(x, 1.0000000000000002);

        let x = calerf(XNEG + f64::EPSILON, 1);
        assert_eq!(x, 2.0);
        let x = calerf(XNEG - f64::EPSILON, 1);
        assert_eq!(x, 2.0);
    }

    #[test]
    fn calerf_2() {
        let x = calerf(THRESH + f64::EPSILON, 2);
        assert_eq!(x, 0.6320696892495559);
        let x = calerf(THRESH - f64::EPSILON, 2);
        assert_eq!(x, 0.6320696892495563);
        let x = calerf(-THRESH - f64::EPSILON, 2);
        assert_eq!(x, 1.8594024168714227);
        let x = calerf(-THRESH + f64::EPSILON, 2);
        assert_eq!(x, 1.8594024168714214);

        let x = calerf(FOUR + f64::EPSILON, 2);
        assert_eq!(x, 0.1369994576250614);
        let x = calerf(FOUR - f64::EPSILON, 2);
        assert_eq!(x, 0.1369994576250614);
        let x = calerf(-FOUR - f64::EPSILON, 2);
        assert_eq!(x, 17772220.904016286);
        let x = calerf(-FOUR + f64::EPSILON, 2);
        assert_eq!(x, 17772220.904016286);

        let x = calerf(XBIG + f64::EPSILON, 2);
        assert_eq!(x, 0.0);
        let x = calerf(XBIG - f64::EPSILON, 2);
        assert_eq!(x, 0.0);
        let x = calerf(-XBIG - f64::EPSILON, 2);
        assert_eq!(x, 1.8831722547514706e306);
        let x = calerf(-XBIG + f64::EPSILON, 2);
        assert_eq!(x, 1.8831722547514706e306);

        let x = calerf(XMAX + f64::EPSILON, 2);
        assert_eq!(x, 0.0);
        let x = calerf(XMAX - f64::EPSILON, 2);
        assert_eq!(x, 0.0);
        let x = calerf(-XMAX - f64::EPSILON, 2);
        assert_eq!(x, 1.79e308);
        let x = calerf(-XMAX + f64::EPSILON, 2);
        assert_eq!(x, 1.79e308);

        let x = calerf(XHUGE + f64::EPSILON, 2);
        assert_eq!(x, 8.408190514869691e-9);
        let x = calerf(XHUGE - f64::EPSILON, 2);
        assert_eq!(x, 8.408190514869691e-9);
        let x = calerf(-XHUGE - f64::EPSILON, 2);
        assert_eq!(x, 1.79e308);
        let x = calerf(-XHUGE + f64::EPSILON, 2);
        assert_eq!(x, 1.79e308);

        let x = calerf(ZERO + f64::EPSILON, 2);
        assert_eq!(x, 0.9999999999999998);
        let x = calerf(ZERO - f64::EPSILON, 2);
        assert_eq!(x, 1.0000000000000002);

        let x = calerf(XNEG + f64::EPSILON, 2);
        assert_eq!(x, 1.728618506590026e308);
        let x = calerf(XNEG - f64::EPSILON, 2);
        assert_eq!(x, 1.728618506590026e308);
    }
}
