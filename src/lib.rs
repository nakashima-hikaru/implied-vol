//! This module provides a pure Rust implementation of calculating implied volatilities.
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
//! let black_vol = implied_vol::implied_black_volatility(20.0, 100.0, 90.0, 30.0, true).unwrap();
//! assert_eq!(black_vol, 0.07011701801482094);
//!
//! let price = implied_vol::black_scholes_option_price(100.0, 90.0, 0.07011701801482094, 30.0, true);
//! assert!(((price - 20.0) / price).abs() <= 2.0 * f64::EPSILON);
//!
//! let normal_vol = implied_vol::implied_normal_volatility(20.0, 100.0, 90.0, 30.0, true).unwrap();
//! assert_eq!(normal_vol, 6.614292466299764);
//!
//! let price = implied_vol::bachelier_option_price(100.0, 90.0, 6.614292466299764, 30.0, true);
//! assert!(((price - 20.0) / price).abs()<= 2.0 * f64::EPSILON);
//! ```

use bon::Builder;
#[cfg(feature = "bench")]
pub mod cxx;

mod fused_multiply_add;
mod lets_be_rational;

pub use crate::lets_be_rational::special_function::DefaultSpecialFn;
pub use crate::lets_be_rational::special_function::SpecialFn;

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
/// The implied Black volatility if it exists, or `None`.
///
/// # Examples
///
/// ```
/// let black_vol = implied_vol::implied_black_volatility(20.0, 100.0, 90.0, 30.0, true).unwrap();
/// assert_eq!(black_vol, 0.07011701801482094);
/// ```
///
/// # Note
/// You must check that the inputs are valid on your own, as this function does not perform input validation.
/// Otherwise, this function can panic.
/// - `option_price` is non-negative and finite.
/// - `forward` is positive and finite.
/// - `strike` is positive and finite.
/// - `expiry` is positive (but can be positive infinity).
#[inline(always)]
#[must_use]
pub fn implied_black_volatility(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> Option<f64> {
    if is_call {
        lets_be_rational::implied_black_volatility_input_unchecked::<DefaultSpecialFn, true>(
            option_price,
            forward,
            strike,
            expiry,
        )
    } else {
        lets_be_rational::implied_black_volatility_input_unchecked::<DefaultSpecialFn, false>(
            option_price,
            forward,
            strike,
            expiry,
        )
    }
}

/// Calculates the price of a European option using the Black-Scholes formula.
///
/// # Arguments
///
/// * `forward` - The current value of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `volatility` - The Black volatility of the underlying asset.
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
/// let price = implied_vol::black_scholes_option_price(100.0, 90.0, 0.07011701801482094, 30.0, true);
/// assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
/// ```
/// # Note
/// You must check that the inputs are valid on your own, as this function does not perform input validation.
#[inline(always)]
#[must_use]
pub fn black_scholes_option_price(
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    if is_call {
        lets_be_rational::bs_option_price::black_input_unchecked::<DefaultSpecialFn, true>(
            forward, strike, volatility, expiry,
        )
    } else {
        lets_be_rational::bs_option_price::black_input_unchecked::<DefaultSpecialFn, false>(
            forward, strike, volatility, expiry,
        )
    }
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
/// let normal_vol = implied_vol::implied_normal_volatility(20.0, 100.0, 90.0, 30.0, true).unwrap();
/// assert_eq!(normal_vol, 6.614292466299764);
/// ```
/// # Note
/// You must check that the inputs are valid on your own, as this function does not perform input validation.
#[inline(always)]
#[must_use]
pub fn implied_normal_volatility(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> Option<f64> {
    if is_call {
        lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<
            DefaultSpecialFn,
            true,
        >(option_price, forward, strike, expiry)
    } else {
        lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<
            DefaultSpecialFn,
            false,
        >(option_price, forward, strike, expiry)
    }
}

/// Calculates the price of an option using Bachelier's model.
///
/// # Arguments
///
/// * `forward` - The forward price of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `volatility` - The normal volatility of the underlying asset.
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
/// let price = implied_vol::bachelier_option_price(100.0, 90.0, 6.614292466299764, 30.0, true);
/// assert!((price - 20.0).abs()<= 2.0 * f64::EPSILON * 20.0);
/// ```
/// # Note
/// You must check that the inputs are valid on your own, as this function does not perform input validation.
#[inline(always)]
#[must_use]
pub fn bachelier_option_price(
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
) -> f64 {
    if is_call {
        lets_be_rational::bachelier_impl::bachelier_price::<true>(
            forward, strike, volatility, expiry,
        )
    } else {
        lets_be_rational::bachelier_impl::bachelier_price::<false>(
            forward, strike, volatility, expiry,
        )
    }
}

/// A struct representing the parameters required for calculating the price of an option
/// using the Black-Scholes model.
///
/// This struct is marked with the `Builder` derive macro, which provides a builder pattern
/// for constructing instances of `PriceBlackScholes`. The builder derives `Clone` and `Debug`
/// traits for convenient usage.
///
/// # Fields
/// - `forward` (f64): The forward price of the underlying asset.
/// - `strike` (f64): The strike price of the option.
/// - `volatility` (f64): The volatility of the underlying asset.
/// - `expiry` (f64): The time to expiration of the option, expressed in years.
/// - `is_call` (bool): A flag indicating whether the option is a call option (`true`)
///   or a put option (`false`).
///
/// # Builder
/// - The `builder` provides a convenient way to create an instance of `PriceBlackScholes`.
/// - The custom `finish_fn` is named `build_internal` and has private visibility to
///   restrict its direct usage, ensuring encapsulation.
///
/// # Example
/// ```rust
/// use implied_vol::PriceBlackScholes;
///
/// let option_params = PriceBlackScholes::builder()
///     .forward(100.0)
///     .strike(95.0)
///     .volatility(0.2)
///     .expiry(1.0)
///     .is_call(true)
///     .build()
///     .unwrap();
/// ```
#[derive(Builder)]
#[builder(derive(Clone, Debug))]
#[builder(finish_fn(vis = "", name = build_internal))]
pub struct PriceBlackScholes {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
}

impl<S: price_black_scholes_builder::IsComplete> PriceBlackScholesBuilder<S> {
    pub fn build(self) -> Option<PriceBlackScholes> {
        let price_black_scholes = self.build_internal();
        if !price_black_scholes.forward.is_finite() {
            return None;
        }
        if !price_black_scholes.strike.is_finite() {
            return None;
        }
        if matches!(
            price_black_scholes.volatility.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        if matches!(
            price_black_scholes.expiry.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        Some(price_black_scholes)
    }
}

impl PriceBlackScholes {
    #[must_use]
    pub fn calculate<SpFn: SpecialFn>(&self) -> f64 {
        if self.is_call {
            lets_be_rational::bs_option_price::black_input_unchecked::<SpFn, true>(
                self.forward,
                self.strike,
                self.volatility,
                self.expiry,
            )
        } else {
            lets_be_rational::bs_option_price::black_input_unchecked::<SpFn, false>(
                self.forward,
                self.strike,
                self.volatility,
                self.expiry,
            )
        }
    }
}

#[derive(Builder)]
#[builder(derive(Clone, Debug))]
#[builder(finish_fn(vis = "", name = build_internal))]
pub struct ImpliedBlackVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

impl<S: implied_black_volatility_builder::IsComplete> ImpliedBlackVolatilityBuilder<S> {
    pub fn build(self) -> Option<ImpliedBlackVolatility> {
        let implied_black_volatility = self.build_internal();

        if matches!(
            implied_black_volatility.forward.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_black_volatility.forward.is_infinite()
        {
            return None;
        }
        if matches!(
            implied_black_volatility.strike.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_black_volatility.strike.is_infinite()
        {
            return None;
        }
        if matches!(
            implied_black_volatility.expiry.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        if matches!(
            implied_black_volatility.option_price.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_black_volatility.option_price.is_infinite()
        {
            return None;
        }
        Some(implied_black_volatility)
    }
}

impl ImpliedBlackVolatility {
    #[must_use]
    pub fn calculate<SpFn: SpecialFn>(&self) -> Option<f64> {
        if self.is_call {
            lets_be_rational::implied_black_volatility_input_unchecked::<SpFn, true>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        } else {
            lets_be_rational::implied_black_volatility_input_unchecked::<SpFn, false>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        }
    }
}

#[derive(Builder)]
#[builder(derive(Clone, Debug))]
#[builder(finish_fn(vis = "", name = build_internal))]
pub struct PriceBachelier {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
}

impl<S: price_bachelier_builder::IsComplete> PriceBachelierBuilder<S> {
    pub fn build(self) -> Option<PriceBachelier> {
        let price_bachelier = self.build_internal();
        if !price_bachelier.forward.is_finite() {
            return None;
        }
        if !price_bachelier.strike.is_finite() {
            return None;
        }
        if matches!(
            price_bachelier.volatility.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        if matches!(
            price_bachelier.expiry.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        Some(price_bachelier)
    }
}

impl PriceBachelier {
    #[must_use]
    pub fn calculate<SpFn: SpecialFn>(&self) -> f64 {
        if self.is_call {
            lets_be_rational::bachelier_impl::bachelier_price::<true>(
                self.forward,
                self.strike,
                self.volatility,
                self.expiry,
            )
        } else {
            lets_be_rational::bachelier_impl::bachelier_price::<false>(
                self.forward,
                self.strike,
                self.volatility,
                self.expiry,
            )
        }
    }
}

#[derive(Builder)]
#[builder(derive(Clone, Debug))]
#[builder(finish_fn(vis = "", name = build_internal))]
pub struct ImpliedNormalVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

impl<S: implied_normal_volatility_builder::IsComplete> ImpliedNormalVolatilityBuilder<S> {
    pub fn build(self) -> Option<ImpliedNormalVolatility> {
        let implied_normal_volatility = self.build_internal();
        if !implied_normal_volatility.forward.is_finite() {
            return None;
        }
        if !implied_normal_volatility.strike.is_finite() {
            return None;
        }
        if matches!(
            implied_normal_volatility.expiry.partial_cmp(&0.0_f64),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        if matches!(
            implied_normal_volatility.option_price.partial_cmp(&0.0_f64),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_normal_volatility.option_price.is_infinite()
        {
            return None;
        }
        Some(implied_normal_volatility)
    }
}

impl ImpliedNormalVolatility {
    #[must_use]
    #[inline]
    pub fn calculate<SpFn: SpecialFn>(&self) -> Option<f64> {
        assert!(
            self.option_price >= 0.0_f64
                && self.option_price.is_finite()
                && self.forward.is_finite()
                && self.strike.is_finite()
                && self.expiry >= 0.0_f64
        ); // never panics because this is validated in the builder
        if self.is_call {
            lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<SpFn, true>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        } else {
            lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<SpFn, false>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strike_anomaly() {
        for k in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let price = 100.0;
            let f = 100.0;
            let t = 1.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn forward_anomaly() {
        for f in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let price = 100.0;
            let k = 100.0;
            let t = 1.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn price_anomaly() {
        for price in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let f = 100.0;
            let t = 1.0;
            let k = 100.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn time_anomaly() {
        for t in [f64::NAN, f64::NEG_INFINITY] {
            let price = 10.0;
            let f = 100.0;
            let k = 100.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn time_inf() {
        let price = 10.0;
        let f = 100.0;
        let k = 100.0;
        let t = f64::INFINITY;
        const Q: bool = true;

        let vol = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(Q)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>()
            .unwrap();
        assert_eq!(vol, 0.0);
    }
}
