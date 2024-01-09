const A: [f64; 5] = [
    3.1611237438705656,
    113.864154151050156,
    377.485237685302021,
    3209.37758913846947,
    0.185777706184603153,
];
const B: [f64; 4] = [
    23.6012909523441209,
    244.024637934444173,
    1282.61652607737228,
    2844.23683343917062,
];
const C: [f64; 9] = [
    0.564188496988670089,
    8.88314979438837594,
    66.1191906371416295,
    298.635138197400131,
    881.95222124176909,
    1712.04761263407058,
    2051.07837782607147,
    1230.33935479799725,
    2.15311535474403846e-8,
];
const D: [f64; 8] = [
    15.7449261107098347,
    117.693950891312499,
    537.181101862009858,
    1621.38957456669019,
    3290.79923573345963,
    4362.61909014324716,
    3439.36767414372164,
    1230.33935480374942,
];
const P: [f64; 6] = [
    0.305326634961232344,
    0.360344899949804439,
    0.125781726111229246,
    0.0160837851487422766,
    6.58749161529837803e-4,
    0.0163153871373020978,
];
const Q: [f64; 5] = [
    2.56852019228982242,
    1.87295284992346047,
    0.527905102951428412,
    0.0605183413124413191,
    0.00233520497626869185,
];

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;
const SQRPI: f64 = 0.56418958354775628695;
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
pub(crate) fn calerf(x: f64, jint: i32) -> f64 {
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
    } else {
        if y >= XBIG {
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
                result = ((-ysq * ysq).exp() * (-del).exp()) * result;
            }
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
    } else {
        if x < ZERO {
            if x < XNEG {
                result = XINF;
            } else {
                let ysq = (x * SIXTEN).trunc() / SIXTEN;
                let del = (x - ysq) * (x + ysq);
                let y = (ysq * ysq).exp() * del.exp();
                result = (y + y) - result;
            }
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
