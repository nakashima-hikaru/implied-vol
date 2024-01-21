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
//! implied-vol = "0.2"
//! ```
//!
//! Then, in your code, bring the functions you need into scope with:
//!
//! ```rust
//! use implied_vol::{implied_black_volatility, calculate_european_option_price_by_black_scholes, implied_normal_volatility, calculate_european_option_price_by_bachelier};
//!
//! let black_vol = implied_black_volatility(20.0, 100.0, 90.0, 30.0, true);
//! assert_eq!(black_vol, 0.07011701801482094);
//!
//! let price = calculate_european_option_price_by_black_scholes(100.0, 90.0, 0.07011701801482094, 30.0, true);
//! assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
//!
//! let normal_vol = implied_normal_volatility(20.0, 100.0, 90.0, 30.0, true);
//! assert_eq!(normal_vol, 6.614292466299764);
//!
//! let price = calculate_european_option_price_by_bachelier(100.0, 90.0, 6.614292466299764, 30.0, true);
//! assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
//! ```
mod erf_cody;
mod normal_distribution;
mod rational_cubic;
pub(crate) mod lets_be_rational;
pub(crate) mod bachelier;
pub(crate) mod constants;

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
/// use implied_vol::implied_black_volatility;
/// let black_vol = implied_black_volatility(20.0, 100.0, 90.0, 30.0, true);
/// assert_eq!(black_vol, 0.07011701801482094);
/// ```
#[inline]
pub fn implied_black_volatility(option_price: f64, forward: f64, strike: f64, expiry: f64, is_call: bool) -> f64 {
    lets_be_rational::implied_black_volatility(option_price, forward, strike, expiry, is_call)
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
/// use implied_vol::calculate_european_option_price_by_black_scholes;
/// let price = calculate_european_option_price_by_black_scholes(100.0, 90.0, 0.07011701801482094, 30.0, true);
/// assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
/// ```
#[inline]
pub fn calculate_european_option_price_by_black_scholes(forward: f64, strike: f64, volatility: f64, expiry: f64, is_call: bool) -> f64 {
    lets_be_rational::black(forward, strike, volatility, expiry, is_call)
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
/// use implied_vol::implied_normal_volatility;
/// let normal_vol = implied_normal_volatility(20.0, 100.0, 90.0, 30.0, true);
/// assert_eq!(normal_vol, 6.614292466299764);
/// ```
pub fn implied_normal_volatility(option_price: f64, forward: f64, strike: f64, expiry: f64, is_call: bool) -> f64 {
    bachelier::implied_normal_volatility(option_price, forward, strike, expiry, is_call)
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
/// use implied_vol::calculate_european_option_price_by_bachelier;
/// let price = calculate_european_option_price_by_bachelier(100.0, 90.0, 6.614292466299764, 30.0, true);
/// assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
/// ```
#[inline]
pub fn calculate_european_option_price_by_bachelier(forward: f64, strike: f64, volatility: f64, expiry: f64, is_call: bool) -> f64 {
    bachelier::bachelier(forward, strike, volatility, expiry, is_call)
}