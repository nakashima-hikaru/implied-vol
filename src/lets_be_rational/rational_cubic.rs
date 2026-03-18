use crate::fused_multiply_add::MulAdd;

const MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 =
    -(1f64 - 0.000_000_014_901_161_193_847_656);
const MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = 2f64 / (f64::EPSILON * f64::EPSILON);

#[inline]
pub fn rational_cubic_interpolation(
    x_minus_x_l: f64,
    h: f64,
    y: (f64, f64),
    d: (f64, f64),
    r: f64,
) -> f64 {
    if h == 0.0_f64 {
        return 0.5 * (y.0 + y.1);
    }
    let t = x_minus_x_l / h;
    if r < MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE {
        let omt = 1.0 - t;
        let omt2 = omt * omt;
        return y.0.mul_add2(r, h * d.0).mul_add2(t, y.0 * omt).mul_add2(
            omt2,
            t * t * y.1.mul_add2(r, -h * d.1).mul_add2(omt, y.1 * t),
        ) / (r - 3.0).mul_add2(t * omt, 1.0);
    }
    (1.0 - t).mul_add2(y.0, y.1 * t)
}

#[inline]
fn rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
    h: f64,
    y_diff: f64,
    d: (f64, f64),
    second_derivative_l: f64,
) -> f64 {
    let numerator = second_derivative_l.mul_add2(0.5 * h, d.1 - d.0);
    if numerator == 0.0_f64 {
        return 0.0;
    }
    let denominator = y_diff / h - d.0;
    if denominator == 0.0_f64 {
        if numerator > 0.0 {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        }
    } else {
        numerator / denominator
    }
}

#[inline]
fn rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
    h: f64,
    y_diff: f64,
    d_r: f64,
    d_diff: f64,
    second_derivative_r: f64,
) -> f64 {
    let numerator = second_derivative_r.mul_add2(0.5 * h, d_diff);
    if numerator == 0.0_f64 {
        return 0.;
    }
    let denominator = d_r - y_diff / h;
    if denominator == 0.0_f64 {
        return if numerator > 0.0 {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        };
    }
    numerator / denominator
}

#[inline]
fn is_zero(x: f64) -> bool {
    x.abs() < f64::MIN_POSITIVE
}

#[inline]
fn minimum_rational_cubic_control_parameter<
    const PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS: bool,
>(
    d: (f64, f64),
    s: f64,
) -> f64 {
    let monotonic = d.0 * s >= 0.0 && d.1 * s >= 0.0;
    let convex = d.0 <= s && s <= d.1;
    let concave = d.0 >= s && s >= d.1;
    if !monotonic && !convex && !concave {
        return MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
    }
    let d_r_m_d_l = d.1 - d.0;
    let d_r_m_s = d.1 - s;
    let s_m_d_l = s - d.0;

    let mut r1 = -f64::MAX;
    let mut r2 = r1;

    if monotonic {
        if !is_zero(s) {
            r1 = (d.1 + d.0) / s;
        } else if PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS {
            r1 = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
        }
    }

    if convex || concave {
        if !(is_zero(s_m_d_l) || is_zero(d_r_m_s)) {
            r2 = (d_r_m_d_l / d_r_m_s).abs().max((d_r_m_d_l / s_m_d_l).abs());
        } else if PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS {
            r2 = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
        }
    } else if monotonic && PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS {
        r2 = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
    }

    MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE.max(r1.max(r2))
}

#[inline]
pub fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side<
    const PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS: bool,
>(
    h: f64,
    y_diff: f64,
    d: (f64, f64),
    second_derivative_l: f64,
) -> f64 {
    let r = rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
        h,
        y_diff,
        d,
        second_derivative_l,
    );
    let s = y_diff / h;
    let r_min =
        minimum_rational_cubic_control_parameter::<PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS>(d, s);
    r.max(r_min)
}

#[inline]
pub fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side<
    const PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS: bool,
>(
    h: f64,
    y_diff: f64,
    d: (f64, f64),
    second_derivative_r: f64,
) -> f64 {
    let r = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
        h,
        y_diff,
        d.1,
        d.1 - d.0,
        second_derivative_r,
    );
    let s = y_diff / h;
    let r_min =
        minimum_rational_cubic_control_parameter::<PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS>(d, s);
    r.max(r_min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rational_cubic_interpolation() {
        let x_l = 0.0;
        assert_eq!(
            rational_cubic_interpolation(
                0.5 - x_l,
                1.0 - x_l,
                (1.0, 2.0),
                (0.0, 0.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.5
        );
        let x_l = 0.0;
        assert_eq!(
            rational_cubic_interpolation(
                0.0 - x_l,
                1.0 - x_l,
                (1.0, 3.0),
                (0.0, 1.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.0
        );
        let x_l = 0.0;
        assert_eq!(
            rational_cubic_interpolation(
                0.2 - x_l,
                1.0 - x_l,
                (1.0, 2.0),
                (0.0, 0.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.200_000_000_000_000_2
        );
        let x_l = 0.0;
        assert_eq!(
            rational_cubic_interpolation(
                0.0 - x_l,
                1.0 - x_l,
                (1.0, 2.0),
                (0.0, 2.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.0
        );
        let x_l = 0.5;
        assert_eq!(
            rational_cubic_interpolation(
                0.5 - x_l,
                0.5 - x_l,
                (1.0, 2.0),
                (0.0, 1.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.5
        );

        let x_l = 0.0;
        let r = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - f64::EPSILON;
        assert_eq!(
            rational_cubic_interpolation(0.5 - x_l, 1.0 - x_l, (1.0, 2.0), (0.0, 0.0), r),
            1.5
        );
        let x_l = 0.0;
        let r = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - f64::EPSILON;
        assert_eq!(
            rational_cubic_interpolation(0.0 - x_l, 1.0 - x_l, (1.0, 3.0), (0.0, 1.0), r),
            1.0
        );
        let x_l = 0.0;
        let r = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - f64::EPSILON;
        assert_eq!(
            rational_cubic_interpolation(0.2 - x_l, 1.0 - x_l, (1.0, 2.0), (0.0, 0.0), r),
            1.200_000_000_000_000_2
        );
        let x_l = 0.0;
        assert_eq!(
            rational_cubic_interpolation(
                0.0 - x_l,
                1.0 - x_l,
                (1.0, 2.0),
                (0.0, 2.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) - f64::EPSILON,
            0.999_999_999_999_999_8
        );
        let x_l = 0.5;
        assert_eq!(
            rational_cubic_interpolation(
                0.5 - x_l,
                0.5 - x_l,
                (1.0, 2.0),
                (0.0, 1.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) - f64::EPSILON,
            1.499_999_999_999_999_8
        );

        let x_l = 0.0;
        let r = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE + f64::EPSILON;
        assert_eq!(
            rational_cubic_interpolation(0.5 - x_l, 1.0 - x_l, (1.0, 2.0), (0.0, 0.0), r),
            1.5
        );
        let x_l = 0.0;
        let r = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE + f64::EPSILON;
        assert_eq!(
            rational_cubic_interpolation(0.0 - x_l, 1.0 - x_l, (1.0, 3.0), (0.0, 1.0), r),
            1.0
        );
        let x_l = 0.0;
        let r = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE + f64::EPSILON;
        assert_eq!(
            rational_cubic_interpolation(0.2 - x_l, 1.0 - x_l, (1.0, 2.0), (0.0, 0.0), r),
            1.200_000_000_000_000_2
        );
        let x_l = 0.0;
        assert_eq!(
            rational_cubic_interpolation(
                0.0 - x_l,
                1.0 - x_l,
                (1.0, 2.0),
                (0.0, 2.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) + f64::EPSILON,
            1.000_000_000_000_000_2
        );
        let x_l = 0.5;
        assert_eq!(
            rational_cubic_interpolation(
                0.5 - x_l,
                0.5 - x_l,
                (1.0, 2.0),
                (0.0, 1.0),
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) + f64::EPSILON,
            1.500_000_000_000_000_2
        );
    }

    #[test]
    fn test_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side() {
        let x_l = 1.0;
        let x_r = 2.0;
        let y_l = 3.0;
        let y_r = 4.0;
        let d_l = 1.0;
        let d_r = 2.0;
        let second_derivative_r = 0.5;

        let output = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
            x_r - x_l,
            y_r - y_l,
            d_r,
            d_r - d_l,
            second_derivative_r,
        );

        assert_eq!(output, 1.25);
    }
}
