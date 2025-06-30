const MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = -(1f64 - 0.000000014901161193847656);
const MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = 2f64 / (f64::EPSILON * f64::EPSILON);

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(crate) fn rational_cubic_interpolation(
    x: f64,
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    r: f64,
) -> f64 {
    let h = x_r - x_l;
    if h == 0.0_f64 {
        return 0.5 * (y_l + y_r);
    }
    let t = (x - x_l) / h;
    if r < MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE {
        let omt = 1.0 - t;
        let t2 = t.powi(2);
        let omt2 = omt.powi(2);
        return (y_r * t2 * t
            + (r * y_r - h * d_r) * t2 * omt
            + (r * y_l + h * d_l) * t * omt2
            + y_l * omt2 * omt)
            / (1.0 + (r - 3.0) * t * omt);
    }
    y_r * t + y_l * (1.0 - t)
}

#[inline(always)]
pub(crate) fn rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_l: f64,
) -> f64 {
    let h = x_r - x_l;
    let numerator = 0.5 * h * second_derivative_l + (d_r - d_l);
    if numerator == 0.0_f64 {
        return 0.0;
    }
    let denominator = (y_r - y_l) / h - d_l;
    if denominator == 0.0_f64 {
        if numerator.is_sign_positive() {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        }
    } else {
        numerator / denominator
    }
}

#[inline(always)]
pub(crate) fn rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_r: f64,
) -> f64 {
    let h = x_r - x_l;
    let numerator = 0.5 * h * second_derivative_r + (d_r - d_l);
    if numerator == 0.0_f64 {
        return 0.;
    }
    let denominator = d_r - (y_r - y_l) / h;
    if denominator == 0.0_f64 {
        return if numerator.is_sign_positive() {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        };
    }
    numerator / denominator
}

#[inline(always)]
pub(crate) fn minimum_rational_cubic_control_parameter(
    d_l: f64,
    d_r: f64,
    s: f64,
    prefer_shape_preservation_over_smoothness: bool,
) -> f64 {
    let monotonic = d_l * s >= 0.0 && d_r * s >= 0.0;
    let convex_or_concave = (d_l <= s && s <= d_r) || (d_l >= s && s >= d_r);
    if !monotonic && !convex_or_concave {
        return MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
    }
    let d_r_m_d_l = d_r - d_l;
    let d_r_m_s = d_r - s;
    let s_m_d_l = s - d_l;
    let r1 = if monotonic && s != 0.0 {
        (d_r + d_l) / s
    } else if monotonic && s == 0.0 && prefer_shape_preservation_over_smoothness {
        MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
    } else {
        f64::MIN
    };
    let r2 = if convex_or_concave {
        if s_m_d_l != 0.0 && d_r_m_s != 0.0 {
            (d_r_m_d_l / d_r_m_s.min(s_m_d_l)).abs()
        } else if prefer_shape_preservation_over_smoothness {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            f64::MIN
        }
    } else if monotonic && prefer_shape_preservation_over_smoothness {
        MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
    } else {
        f64::MIN
    };
    r1.max(r2)
        .max(MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE)
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(crate) fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_l: f64,
    prefer_shape_preservation_over_smoothness: bool,
) -> f64 {
    let r = rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
        x_l,
        x_r,
        y_l,
        y_r,
        d_l,
        d_r,
        second_derivative_l,
    );
    let r_min = minimum_rational_cubic_control_parameter(
        d_l,
        d_r,
        (y_r - y_l) / (x_r - x_l),
        prefer_shape_preservation_over_smoothness,
    );
    r.max(r_min)
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(crate) fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_r: f64,
    prefer_shape_preservation_over_smoothness: bool,
) -> f64 {
    let r = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
        x_l,
        x_r,
        y_l,
        y_r,
        d_l,
        d_r,
        second_derivative_r,
    );
    let r_min = minimum_rational_cubic_control_parameter(
        d_l,
        d_r,
        (y_r - y_l) / (x_r - x_l),
        prefer_shape_preservation_over_smoothness,
    );
    r.max(r_min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rational_cubic_interpolation() {
        assert_eq!(
            rational_cubic_interpolation(
                0.5,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                0.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.5
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.0,
                0.0,
                1.0,
                1.0,
                3.0,
                0.0,
                1.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.0
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.2,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                0.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.2000000000000002
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.0,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                2.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.0
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.5,
                0.5,
                0.5,
                1.0,
                2.0,
                0.0,
                1.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ),
            1.5
        );

        assert_eq!(
            rational_cubic_interpolation(
                0.5,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                0.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - f64::EPSILON
            ),
            1.5
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.0,
                0.0,
                1.0,
                1.0,
                3.0,
                0.0,
                1.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - f64::EPSILON
            ),
            1.0
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.2,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                0.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - f64::EPSILON
            ),
            1.2000000000000002
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.0,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                2.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) - f64::EPSILON,
            0.9999999999999998
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.5,
                0.5,
                0.5,
                1.0,
                2.0,
                0.0,
                1.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) - f64::EPSILON,
            1.4999999999999998
        );

        assert_eq!(
            rational_cubic_interpolation(
                0.5,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                0.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE + f64::EPSILON
            ),
            1.5
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.0,
                0.0,
                1.0,
                1.0,
                3.0,
                0.0,
                1.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE + f64::EPSILON
            ),
            1.0
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.2,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                0.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE + f64::EPSILON
            ),
            1.2000000000000002
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.0,
                0.0,
                1.0,
                1.0,
                2.0,
                0.0,
                2.0,
                MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
            ) + f64::EPSILON,
            1.0000000000000002
        );
        assert_eq!(
            rational_cubic_interpolation(
                0.5,
                0.5,
                0.5,
                1.0,
                2.0,
                0.0,
                1.0,
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
            x_l,
            x_r,
            y_l,
            y_r,
            d_l,
            d_r,
            second_derivative_r,
        );

        assert_eq!(output, 1.25);
    }
}
