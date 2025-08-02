mod erf_cody;
pub mod normal_distribution;

use crate::special_function::erf_cody::{erf_cody, erfc_cody, erfcx_cody};
use crate::special_function::normal_distribution::{erfinv, inverse_norm_cdf, norm_pdf};

/// The `SpecialFn` trait provides a collection of special mathematical functions that are commonly
/// used in statistics, probability, and scientific computations. These methods are implemented
/// to calculate specific special functions such as the error function (`erf`) and its complement,
/// inverse error functions, the normal probability density function, and the inverse cumulative
/// distribution function for a standard normal distribution.
///
/// # Required Methods
///
/// ### `erf`
///
/// Computes the error function `erf(x)`. The error function gives the probability of a random variable
/// from a normal distribution falling within a certain range of standard deviations.
///
/// ### `erfc`
///
/// Computes the complementary error function `erfc(x)` which is defined as `1 - erf(x)`.
/// It is often used to simplify formulas for survival or tail probability distributions.
///
/// ### `erfcx`
///
/// Computes the scaled complementary error function `erfcx(x)` defined as `exp(x^2) * erfc(x)`.
/// It is often used to avoid numerical instability in calculations involving the complementary error function.
///
/// ### `erfinv`
///
/// Computes the inverse error function `erfinv(x)`. This function determines the value of the argument
/// that produces the given value for the error function, i.e., it satisfies `erf(erfinv(x)) = x`.
///
/// ### `inverse_norm_cdf`
///
/// Computes the inverse of the cumulative distribution function (quantile function) for a standard
/// normal distribution. That is, for a given probability `x` in the range `[0, 1]`, it returns
/// the value `z` such that the probability of a standard normal random variable `N(0, 1)` being
/// less than `z` equals `x`.
///
/// ### `norm_pdf`
///
/// Computes the probability density function of the standard normal distribution `N(0, 1)`.
/// This function evaluates the PDF at a given input `x`.
pub trait SpecialFn {
    fn erf(x: f64) -> f64;
    fn erfc(x: f64) -> f64;
    fn erfcx(x: f64) -> f64;
    fn erfinv(x: f64) -> f64;
    fn inverse_norm_cdf(x: f64) -> f64;
    fn norm_pdf(x: f64) -> f64;
}

/// A struct representing the default implementation of a special function.
pub struct DefaultSpecialFn;
impl SpecialFn for DefaultSpecialFn {
    #[inline(always)]
    fn erf(x: f64) -> f64 {
        erf_cody(x)
    }

    #[inline(always)]
    fn erfc(x: f64) -> f64 {
        erfc_cody(x)
    }
    #[inline(always)]
    fn erfcx(x: f64) -> f64 {
        erfcx_cody(x)
    }
    #[inline(always)]
    fn erfinv(x: f64) -> f64 {
        erfinv(x)
    }
    #[inline(always)]
    fn inverse_norm_cdf(x: f64) -> f64 {
        inverse_norm_cdf(x)
    }
    #[inline(always)]
    fn norm_pdf(x: f64) -> f64 {
        norm_pdf(x)
    }
}

#[cfg(feature = "normal-distribution")]
#[inline(always)]
pub(super) fn norm_cdf(z: f64) -> f64 {
    normal_distribution::norm_cdf(z)
}
