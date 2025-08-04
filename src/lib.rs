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
//! implied-vol = "2"
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

mod bachelier_impl;
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
#[inline]
pub fn implied_normal_volatility(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    bachelier_impl::implied_normal_volatility::<DefaultSpecialFn>(
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
    bachelier_impl::bachelier_price(forward, strike, volatility, expiry, is_call)
}

/// A module implementing pricing and implied volatility computation using the Black-Scholes model.
pub mod black {
    use crate::{SpecialFn, lets_be_rational};

    /// Computes the price of a European option using the Black-Scholes model.
    ///
    /// # Type Parameters
    /// - `SpFn`: A type that implements the `SpecialFn` trait.
    ///
    /// # Arguments
    /// - `forward`: The forward price of the underlying asset.
    /// - `strike`: The strike price of the option.
    /// - `volatility`: The Black volatility of the underlying asset (in annualized terms).
    /// - `expiry`: Time to expiry of the option, in years.
    /// - `is_call`: A boolean flag indicating whether the option is a call (`true`)
    ///   or a put (`false`).
    ///
    /// # Returns
    /// - The price of the option, calculated under the Black-Scholes model.
    ///
    /// # Example
    /// ```rust
    /// use implied_vol::black::price;
    /// use implied_vol::DefaultSpecialFn;
    ///
    /// let forward = 100.0;
    /// let strike = 95.0;
    /// let volatility = 0.2; // 20% annualized volatility
    /// let expiry = 1.0; // 1 year
    /// let is_call = true; // Call option
    ///
    /// let option_price = price::<DefaultSpecialFn>(forward, strike, volatility, expiry, is_call);
    /// println!("The option price is: {}", option_price);
    /// ```
    #[inline]
    pub fn price<SpFn: SpecialFn>(
        forward: f64,
        strike: f64,
        volatility: f64,
        expiry: f64,
        is_call: bool,
    ) -> f64 {
        lets_be_rational::black::<SpFn>(forward, strike, volatility, expiry, is_call)
    }

    /// Computes the implied Black volatility from the European option price.
    ///
    /// # Parameters
    ///
    /// - `option_price`: The observed price of the option (call or put).
    /// - `forward`: The forward price of the underlying asset.
    /// - `strike`: The strike price of the option.
    /// - `expiry`: The time to expiration of the option, typically expressed in years.
    /// - `is_call`: A boolean flag indicating whether the option is a call (`true`) or a put (`false`).
    ///
    /// # Type Parameters
    ///
    /// - `SpFn`: A special function type implementing numerical methods, such as approximations for
    ///   probability distributions, required for calculating the implied volatility. This is a generic
    ///   placeholder for specialization traits used in numerical methods.
    ///
    /// # Returns
    ///
    /// The implied Black volatility.
    ///
    /// # Examples
    ///
    /// ```
    /// use implied_vol::black::implied_volatility;
    /// use implied_vol::DefaultSpecialFn;
    ///
    /// let option_price = 10.0;
    /// let forward = 100.0;
    /// let strike = 95.0;
    /// let expiry = 1.0; // 1 year
    /// let is_call = true;
    ///
    /// let iv = implied_volatility::<DefaultSpecialFn>(option_price, forward, strike, expiry, is_call);
    /// println!("Implied Volatility: {}", iv);
    /// ```
    #[inline]
    pub fn implied_volatility<SpFn: SpecialFn>(
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
}

pub mod bachelier {
    use crate::{SpecialFn, bachelier_impl};

    /// Calculates the European option price under the Bachelier model.
    ///
    /// # Parameters
    /// - `forward`: The forward price of the underlying asset.
    /// - `strike`: The strike price of the option.
    /// - `volatility`: The normal volatility of the underlying asset (constant for this model).
    /// - `expiry`: The time to expiry of the option, expressed in years.
    /// - `is_call`: A boolean indicating whether the option is a call option (`true`) or a put option (`false`).
    ///
    /// # Returns
    /// The computed price of the option.
    ///
    /// # type parameters
    /// - `SpFn`: A trait bound representing a special function type, used internally by the computation.
    ///
    /// # Example
    /// ```
    /// use implied_vol::bachelier::price;
    /// use implied_vol::DefaultSpecialFn;
    ///
    /// let forward = 100.0;
    /// let strike = 95.0;
    /// let volatility = 10.0;
    /// let expiry = 1.0;
    /// let is_call = true;
    ///
    /// let option_price = price::<DefaultSpecialFn>(forward, strike, volatility, expiry, is_call);
    /// println!("Option Price: {}", option_price);
    /// ```
    #[inline]
    pub fn price<SpFn: SpecialFn>(
        forward: f64,
        strike: f64,
        volatility: f64,
        expiry: f64,
        is_call: bool,
    ) -> f64 {
        bachelier_impl::bachelier_price(forward, strike, volatility, expiry, is_call)
    }

    /// Computes the implied normal volatility for a European option price.
    ///
    /// # Type Parameters
    /// * `SpFn`: A trait bound for objects that implement the `SpecialFn` trait, used
    ///   internally by the computation of the implied normal volatility.
    ///
    /// # Arguments
    /// * `option_price` - The observed market price of the option.
    /// * `forward` - The forward price of the underlying asset.
    /// * `strike` - The strike price of the option.
    /// * `expiry` - The time to expiry of the option (in years).
    /// * `is_call` - A boolean indicating whether the option is a call option (`true`)
    ///   or a put option (`false`).
    ///
    /// # Returns
    /// The implied normal volatility.
    ///
    /// # Example
    /// ```
    /// use implied_vol::bachelier::implied_volatility;
    /// use implied_vol::DefaultSpecialFn;
    ///
    /// let option_price = 1.50;
    /// let forward = 100.0;
    /// let strike = 102.0;
    /// let expiry = 0.5;
    /// let is_call = true;
    ///
    /// let vol = implied_volatility::<DefaultSpecialFn>(
    ///     option_price,
    ///     forward,
    ///     strike,
    ///     expiry,
    ///     is_call,
    /// );
    ///
    /// println!("Implied Volatility: {}", vol);
    /// ```
    #[inline]
    pub fn implied_volatility<SpFn: SpecialFn>(
        option_price: f64,
        forward: f64,
        strike: f64,
        expiry: f64,
        is_call: bool,
    ) -> f64 {
        bachelier_impl::implied_normal_volatility::<SpFn>(
            option_price,
            forward,
            strike,
            expiry,
            is_call,
        )
    }
}
