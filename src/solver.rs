//! Black inverse solvers selected independently for each calculation.
//!
//! [`Hybrid`] and [`Jaeckel`] are always available. The `flashiv` and
//! `experimental` features add their respective solver types without changing
//! the default. Select a type with the builders' `calculate_with` method.

use crate::{DefaultSpecialFn, SpecialFn, lets_be_rational};
use std::marker::PhantomData;

mod sealed {
    pub trait Sealed {
        const MINIMUM_TIME_VALUE: f64 = f64::MIN_POSITIVE;
        const CHECK_ROUNDED_INTERIOR_CAP: bool = true;
    }
}

/// A supported Black inverse solver.
///
/// Implementations preserve their own ATM, price-boundary, and convergence
/// contracts. This trait is sealed; customize special functions through the
/// solver's type parameter rather than implementing an additional solver.
pub trait BlackSolver: sealed::Sealed {
    /// Invert an OTM price divided by `sqrt(F*K)` to total volatility.
    ///
    /// `log_moneyness` is `ln(F/K)`. Invalid inputs return `None`.
    /// Zero price returns zero volatility. Upper-cap classification and
    /// numerical convergence behavior follow the selected solver's contract.
    fn implied_total_volatility(log_moneyness: f64, normalised_price: f64) -> Option<f64>;
}

/// The default precision-qualified mixture of Jäckel and restricted `FlashIV`.
pub struct Hybrid<SpFn = DefaultSpecialFn>(PhantomData<SpFn>);

/// Jäckel's Let's Be Rational inverse, with no cross-method fallback.
pub struct Jaeckel<SpFn = DefaultSpecialFn>(PhantomData<SpFn>);

/// The fixed-count `FlashIV` inverse following the paper's Algorithm 1.
///
/// Available with the `flashiv` feature. It never falls back to Jäckel or adds
/// the optional `FlashIV+` final price correction. Its accuracy follows the
/// paper method, rather than the default hybrid's attainable-precision target.
#[cfg(feature = "flashiv")]
pub struct FlashIv<SpFn = DefaultSpecialFn>(PhantomData<SpFn>);

/// The author's own Black inverse solver, using fixed special functions.
///
/// Available with the `experimental` feature. Its explicit FMA arithmetic,
/// original internal dispatcher, and exact normalized-cap evaluation are
/// independent of the optional `fma` feature and custom `SpecialFn` providers.
#[cfg(feature = "experimental")]
pub struct Experimental;

impl<SpFn: SpecialFn> sealed::Sealed for Hybrid<SpFn> {}
impl<SpFn: SpecialFn> sealed::Sealed for Jaeckel<SpFn> {}
#[cfg(feature = "flashiv")]
impl<SpFn: SpecialFn> sealed::Sealed for FlashIv<SpFn> {}
#[cfg(feature = "experimental")]
impl sealed::Sealed for Experimental {
    const MINIMUM_TIME_VALUE: f64 = 0.0;
    const CHECK_ROUNDED_INTERIOR_CAP: bool = false;
}

#[inline]
fn valid_input(x: f64, beta: f64) -> bool {
    x.is_finite() && beta.is_finite() && beta >= 0.0
}

#[inline]
fn standard_inverse<SpFn: SpecialFn, const HYBRID: bool>(x: f64, beta: f64) -> Option<f64> {
    if !valid_input(x, beta) {
        return None;
    }
    if beta == 0.0 {
        return Some(0.0);
    }
    let theta_x = -x.abs();
    let b_max = (0.5 * theta_x).exp();
    if beta >= b_max {
        return (beta == b_max).then_some(f64::INFINITY);
    }
    if x == 0.0 {
        return Some(lets_be_rational::implied_normalised_volatility_atm::<SpFn>(
            beta,
        ));
    }
    Some(if HYBRID {
        lets_be_rational::hybrid_unchecked::<SpFn>(beta, theta_x, b_max)
    } else {
        lets_be_rational::lets_be_rational_unchecked::<SpFn>(beta, theta_x, b_max)
    })
}

impl<SpFn: SpecialFn> BlackSolver for Hybrid<SpFn> {
    #[inline]
    fn implied_total_volatility(x: f64, beta: f64) -> Option<f64> {
        standard_inverse::<SpFn, true>(x, beta)
    }
}

impl<SpFn: SpecialFn> BlackSolver for Jaeckel<SpFn> {
    #[inline]
    fn implied_total_volatility(x: f64, beta: f64) -> Option<f64> {
        standard_inverse::<SpFn, false>(x, beta)
    }
}

#[cfg(feature = "flashiv")]
impl<SpFn: SpecialFn> BlackSolver for FlashIv<SpFn> {
    #[inline]
    fn implied_total_volatility(x: f64, beta: f64) -> Option<f64> {
        if !valid_input(x, beta) {
            return None;
        }
        if beta == 0.0 {
            return Some(0.0);
        }
        let theta_x = -x.abs();
        let b_max = (0.5 * theta_x).exp();
        if beta >= b_max {
            return (beta == b_max).then_some(f64::INFINITY);
        }
        crate::flashiv::normalised_volatility::<SpFn>(beta, theta_x, b_max)
    }
}

#[cfg(feature = "experimental")]
impl BlackSolver for Experimental {
    #[inline]
    fn implied_total_volatility(x: f64, beta: f64) -> Option<f64> {
        if !valid_input(x, beta) {
            return None;
        }
        if beta == 0.0 {
            return Some(0.0);
        }
        if x == 0.0 && beta == 1.0 {
            return Some(f64::INFINITY);
        }
        let s = crate::experimental::implied_total_volatility(x.abs(), beta);
        (s.is_finite() && s > 0.0).then_some(s)
    }
}

#[inline(always)]
pub(crate) fn implied_black_volatility<S: BlackSolver, const IS_CALL: bool>(
    price: f64,
    f: f64,
    k: f64,
    t: f64,
) -> Option<f64> {
    let intrinsic = (if IS_CALL { f - k } else { k - f }).max(0.0);
    if t == 0.0 {
        return (price == intrinsic).then_some(0.0);
    }
    let cap = if IS_CALL { f } else { k };
    if price >= cap {
        return (price == cap).then_some(f64::INFINITY);
    }
    let beta = if intrinsic > 0.0 {
        price - intrinsic
    } else {
        price
    } / (f.sqrt() * k.sqrt());
    if beta <= S::MINIMUM_TIME_VALUE {
        return (beta >= 0.0).then_some(0.0);
    }
    let x = if f == k {
        0.0
    } else {
        lets_be_rational::bs_option_price::negative_abs_log_moneyness(f, k)
    };
    let total_volatility = S::implied_total_volatility(x, beta)?;
    // The market cap was handled above. Only an infinite normalized result
    // needs the interior-cap check; finite inversions already checked that cap.
    if S::CHECK_ROUNDED_INTERIOR_CAP
        && f != k
        && total_volatility.is_infinite()
        && beta >= (0.5 * -x.abs()).exp()
    {
        return None;
    }
    Some(total_volatility / t.sqrt())
}
