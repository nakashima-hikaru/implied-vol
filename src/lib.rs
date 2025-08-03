//! This module provides a Rust implementation of Peter Jäckel's original C++ code for calculating
//! implied volatilities in financial derivatives.
//! To learn more about the algorithms, please refer to Peter Jäckel's papers [Let's Be Rational](http://www.jaeckel.org/LetsBeRational.pdf) and [Implied Normal Volatility](http://www.jaeckel.org/ImpliedNormalVolatility.pdf).
//!
//! # Features
//!
//! - Calculation of implied Black volatility
//! - Calculation of the price of a European option using the Black-Scholes model
//! - Calculation of implied normal volatility
//! - Calculation of the price of an option using Bachelier's model
//!
//! All models support both call and put options.
//!
//! # Examples
//!
//! Check out the documentation for each function for practical examples.
//!
//! # Usage
//!
//! Import the crate in your Rust project by adding the following to your `Cargo.toml`
//!
//! ```toml
//! [dependencies]
//! implied-vol = "1.2.3"
//! ```
//!
//! Then, in your code, bring the functions you need into scope with:
//!
//! ```rust
//! let black_vol = implied_vol::implied_black_volatility(20.0, 100.0, 90.0, 30.0, true);
//! assert_eq!(black_vol, 0.07011701801482094);
//!
//! let price = implied_vol::calculate_european_option_price_by_black_scholes(100.0, 90.0, 0.07011701801482094, 30.0, true);
//! assert!(((price - 20.0) / price).abs() <= 2.0 * f64::EPSILON);
//!
//! let normal_vol = implied_vol::implied_normal_volatility(20.0, 100.0, 90.0, 30.0, true);
//! assert_eq!(normal_vol, 6.614292466299764);
//!
//! let price = implied_vol::calculate_european_option_price_by_bachelier(100.0, 90.0, 6.614292466299764, 30.0, true);
//! assert!(((price - 20.0) / price).abs()<= 2.0 * f64::EPSILON);
//! ```
//!
//! Moreover, you can use some internal functions if you specify feature flags:
//!
//! ```toml
//! [dependencies]
//! implied-vol = { versions = "1.0.0", features = ["normal-distribution", "error-function"] }
//! ```
//!
//! For detailed explanations of each feature, please refer to the README.md file.
pub use crate::special_function::{DefaultSpecialFn, SpecialFn};

mod bachelier;
mod constants;
mod lets_be_rational;
mod rational_cubic;
pub mod special_function;

trait MulAdd {
    fn mul_add2(self, a: f64, b: f64) -> f64;
}
#[cfg(feature = "fma")]
#[inline(always)]
fn fma_function(r: f64, a: f64, b: f64) -> f64 {
    r.mul_add(a, b)
}

#[cfg(not(feature = "fma"))]
#[inline(always)]
fn fma_function(r: f64, a: f64, b: f64) -> f64 {
    a * r + b
}
impl MulAdd for f64 {
    #[inline(always)]
    fn mul_add2(self, a: f64, b: f64) -> f64 {
        fma_function(self, a, b)
    }
}

/// Calculates the implied black volatility using a transformed rational guess with limited iterations.
///
/// # Arguments
///
/// * `option_price` - The current price of the option.
/// * `forward` - The current forward price of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `expiry` - The time to expiration in years.
/// * `is_call` - A boolean flag indicating whether the option is a call (true) or put (false).
///
/// # Returns
///
/// The implied black volatility.
///
/// # Examples
///
/// ```
/// let black_vol = implied_vol::implied_black_volatility(20.0, 100.0, 90.0, 30.0, true);
/// assert_eq!(black_vol, 0.07011701801482094);
/// ```
#[inline]
pub fn implied_black_volatility(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    lets_be_rational::implied_black_volatility::<DefaultSpecialFn>(
        option_price,
        forward,
        strike,
        expiry,
        is_call,
    )
}

/// Computes the implied Black-Scholes volatility of an option given its price.
///
/// # Type Parameters
/// - `SpFn`: A type implementing the `SpecialFn` trait. This type is used internally
///   by the `lets_be_rational` library for specialized mathematical computations.
///
/// # Arguments
/// - `option_price`: The observed market price of the option.
/// - `forward`: The forward price of the underlying asset.
/// - `strike`: The strike price of the option.
/// - `expiry`: The time to expiry (in years) of the option.
/// - `is_call`: A boolean flag indicating the type of the option:
///   - `true` if the option is a call.
///   - `false` if the option is a put.
///
/// # Returns
/// - Returns a `f64` representing the implied volatility of the option.
///
/// # Examples
/// ```
/// use implied_vol::implied_black_volatility_custom;
/// use implied_vol::special_function::DefaultSpecialFn;
///
/// let option_price = 10.0;
/// let forward = 100.0;
/// let strike = 105.0;
/// let expiry = 1.0; // 1 year
/// let is_call = true;
///
/// let implied_vol = implied_black_volatility_custom::<DefaultSpecialFn>(
///     option_price,
///     forward,
///     strike,
///     expiry,
///     is_call,
/// );
///
/// println!("Implied Volatility: {}", implied_vol);
/// ```
#[inline]
pub fn implied_black_volatility_custom<SpFn: SpecialFn>(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    lets_be_rational::implied_black_volatility::<SpFn>(
        option_price,
        forward,
        strike,
        expiry,
        is_call,
    )
}

/// Calculates the price of a European option using the Black-Scholes formula.
///
/// # Arguments
///
/// * `forward` - The current value of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `volatility` - The volatility of the underlying asset.
/// * `expiry` - The time to expiration of the option.
/// * `is_call` - A boolean flag indicating whether the option is a call (true) or put (false).
///
/// # Returns
///
/// The price of the European option.
///
/// # Examples
///
/// ```
/// let price = implied_vol::calculate_european_option_price_by_black_scholes(100.0, 90.0, 0.07011701801482094, 30.0, true);
/// assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
/// ```
#[inline]
pub fn calculate_european_option_price_by_black_scholes(
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    lets_be_rational::black::<DefaultSpecialFn>(forward, strike, volatility, expiry, is_call)
}

/// Calculates the price of a European option using the Black-Scholes model with a custom special function.
///
/// # Type Parameters
/// - `SpFn`: A custom type that implements the `SpecialFn` trait to handle special mathematical functions used in the computation.
///
/// # Parameters
/// - `forward` (`f64`): The forward price of the underlying asset.
/// - `strike` (`f64`): The strike price of the option.
/// - `volatility` (`f64`): The implied volatility of the underlying asset, expressed as a decimal.
/// - `expiry` (`f64`): The time to expiry of the option, measured in years.
/// - `is_call` (`bool`): A flag to specify the option type:
///   - `true` for a call option.
///   - `false` for a put option.
///
/// # Returns
/// - `f64`: The Black-Scholes price of the European option, calculated using the custom special function.
///
/// # Example
/// ```
/// use implied_vol::calculate_european_option_price_by_black_scholes_custom;
/// use implied_vol::special_function::DefaultSpecialFn;
/// let forward_price = 100.0; // Forward price of the underlying asset
/// let strike_price = 105.0; // Strike price of the option
/// let volatility = 0.2; // Implied volatility (20%)
/// let time_to_expiry = 1.0; // Time to expiry in years
/// let is_call = true; // Call option
///
/// let price = calculate_european_option_price_by_black_scholes_custom::<DefaultSpecialFn>(
///     forward_price,
///     strike_price,
///     volatility,
///     time_to_expiry,
///     is_call,
/// );
///
/// println!("European option price: {}", price);
/// ```
#[inline]
pub fn calculate_european_option_price_by_black_scholes_custom<SpFn: SpecialFn>(
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    lets_be_rational::black::<SpFn>(forward, strike, volatility, expiry, is_call)
}

/// Calculates the implied normal volatility.
///
/// # Arguments
///
/// * `price` - The market price of the option.
/// * `forward` - The forward price of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `expiry` - The time to expiration in years.
/// * `is_call` - A boolean flag indicating whether the option is a call (true) or put (false).
///
/// # Returns
///
/// The implied normal volatility as a `f64` value.
///
/// # Examples
///
/// ```
/// let normal_vol = implied_vol::implied_normal_volatility(20.0, 100.0, 90.0, 30.0, true);
/// assert_eq!(normal_vol, 6.614292466299764);
/// ```
pub fn implied_normal_volatility(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    bachelier::implied_normal_volatility::<DefaultSpecialFn>(
        option_price,
        forward,
        strike,
        expiry,
        is_call,
    )
}

/// Calculates the price of an option using Bachelier's model.
///
/// # Arguments
///
/// * `forward` - The forward price of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `volatility` - The volatility of the underlying asset.
/// * `expiry` - The time to expiration in years.
/// * `is_call` - A boolean flag indicating whether the option is a call (true) or a put (false).
///
/// # Returns
///
/// The price of the European option.
///
/// # Examples
///
/// ```
/// let price = implied_vol::calculate_european_option_price_by_bachelier(100.0, 90.0, 6.614292466299764, 30.0, true);
/// assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
/// ```
#[inline]
pub fn calculate_european_option_price_by_bachelier(
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    bachelier::bachelier(forward, strike, volatility, expiry, is_call)
}

#[cfg(feature = "error-function")]
/// Calculates the scaled complementary error function of `x`.
///
/// The scaled complementary error function is defined as: `erfcx(x) = exp(x^2) * erfc(x)`,
/// where `erfc(x)` is the complementary error function.
///
/// # Arguments
///
/// * `x` - The input value to calculate the scaled complementary error function for.
///
/// # Returns
///
/// The result of calculating the scaled complementary error function of `x`.
///
/// # Example
///
/// ```
/// let result = implied_vol::erfcx(0.5);
/// assert!((result - 0.6156903441929259) / result <= f64::EPSILON);
/// ```
#[inline]
pub fn erfcx(x: f64) -> f64 {
    DefaultSpecialFn::erfcx(x)
}

#[cfg(feature = "error-function")]
/// Calculates the complementary error function.
///
/// # Arguments
///
/// * `x` - The input number for which the complementary error function needs to be calculated.
///
/// # Returns
///
/// The result of the complementary error function calculation.
///
/// # Example
///
/// ```
/// let result = implied_vol::erfc(0.5);
/// assert!((result - 0.4795001221869535) / result <= f64::EPSILON);
/// ```
#[inline]
pub fn erfc(x: f64) -> f64 {
    DefaultSpecialFn::erfc(x)
}

/// Calculates the probability density function of a standard normal distribution.
///
/// # Arguments
///
/// * `x` - The value at which to calculate the probability density function.
///
/// # Returns
///
/// The probability density function value at the given `x` value.
///
/// # Examples
///
/// ```
/// let pdf = implied_vol::norm_pdf(0.0);
/// assert!((pdf - 0.3989422804014327) / pdf <= f64::EPSILON);
/// ```
#[cfg(feature = "normal-distribution")]
#[inline]
pub fn norm_pdf(x: f64) -> f64 {
    DefaultSpecialFn::norm_pdf(x)
}
/// Calculates the cumulative distribution function (CDF) of the standard normal distribution.
///
/// # Arguments
///
/// * `x` - The value at which to calculate the CDF.
///
/// # Returns
///
/// The CDF value for `x` in the standard normal distribution, ranging from 0 to 1.
///
/// # Examples
///
/// ```
/// let cdf = implied_vol::norm_cdf(1.5);
/// assert!((cdf - 0.9331927987311419) / cdf <= f64::EPSILON);
/// ```
#[cfg(feature = "normal-distribution")]
#[inline]
pub fn norm_cdf(x: f64) -> f64 {
    DefaultSpecialFn::norm_cdf(x)
}

#[cfg(feature = "normal-distribution")]
/// Calculates the inverse cumulative distribution function (CDF).
///
/// The inverse CDF is also known as the quantile function or percent-point function.
/// It returns the value x such that P(X < x) = probability, where X follows a standard normal distribution.
///
/// # Arguments
///
/// * `x` - The probability value between 0 and 1.
///
/// # Examples
///
/// ```
/// let probability = 0.8;
/// let inverse_cdf = implied_vol::inverse_norm_cdf(probability);
/// assert!((inverse_cdf - 0.8416212335729144) / inverse_cdf <= f64::EPSILON);
/// ```
///
/// # Panics
///
/// This function will panic if the given probability value is outside the range [0, 1].
#[inline]
pub fn inverse_norm_cdf(x: f64) -> f64 {
    DefaultSpecialFn::inverse_norm_cdf(x)
}
