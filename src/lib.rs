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
pub use crate::special_function::{DefaultSpecialFn, SpecialFn};

mod bachelier;
mod constants;
mod fused_multiply_add;
mod lets_be_rational;
mod rational_cubic;
pub mod special_function;

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
/// The implied Black volatility.
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
/// - The implied Black volatility of the option.
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
/// The price of the European option based on the Black-Scholes model.
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
