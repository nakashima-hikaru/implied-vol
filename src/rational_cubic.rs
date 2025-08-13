const MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = -(1f64 - 0.000_000_014_901_161_193_847_656);
const MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = 2f64 / (f64::EPSILON * f64::EPSILON);

#[inline(always)]
pub const fn rational_cubic_interpolation(
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
        let t2 = t * t;
        let omt2 = omt * omt;
        return (y.1 * t2 * t
            + (r * y.1 - h * d.1) * t2 * omt
            + (r * y.0 + h * d.0) * t * omt2
            + y.0 * omt2 * omt)
            / (1.0 + (r - 3.0) * t * omt);
    }
    y.1 * t + y.0 * (1.0 - t)
}

#[inline(always)]
const fn rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
    h: f64,
    y_diff: f64,
    d: (f64, f64),
    second_derivative_l: f64,
) -> f64 {
    let numerator = 0.5 * h * second_derivative_l + (d.1 - d.0);
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

#[inline(always)]
const fn rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
    h: f64,
    y_diff: f64,
    d_r: f64,
    d_diff: f64,
    second_derivative_r: f64,
) -> f64 {
    let numerator = 0.5 * h * second_derivative_r + d_diff;
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

#[inline(always)]
const fn minimum_rational_cubic_control_parameter<
    const PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS: bool,
>(
    d: (f64, f64),
    s: f64,
) -> f64 {
    let monotonic = d.0 * s >= 0.0 && d.1 * s >= 0.0;
    let convex_or_concave = (d.0 <= s && s <= d.1) || (d.0 >= s && s >= d.1);
    if !monotonic && !convex_or_concave {
        return MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
    }
    let r1 = if monotonic && s != 0.0 {
        (d.1 + d.0) / s
    } else if monotonic && PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS {
        MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
    } else {
        f64::MIN
    };
    let r2 = if convex_or_concave {
        let s_m_d_l = s - d.0;
        let d_r_m_s = d.1 - s;
        let d_r_m_d_l = d.1 - d.0;
        if s_m_d_l != 0.0 || d_r_m_s == 0.0 {
            (d_r_m_d_l / d_r_m_s.min(s_m_d_l)).abs()
        } else if PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            r1
        }
    } else if monotonic && PREFER_SHAPE_PRESERVATION_OVER_SMOOTHNESS {
        MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
    } else {
        r1
    };
    r1.max(r2)
        .max(MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE)
}

#[inline(always)]
pub const fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side<
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

#[inline(always)]
pub const fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side<
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
            1.2000000000000002
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
            1.2000000000000002
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
            0.9999999999999998
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
            1.4999999999999998
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
            1.2000000000000002
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
            1.0000000000000002
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
            1.5000000000000002
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
