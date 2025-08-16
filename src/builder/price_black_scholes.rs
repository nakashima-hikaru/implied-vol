use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

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
    pub fn build_unchecked(self) -> PriceBlackScholes {
        self.build_internal()
    }
}

impl PriceBlackScholes {
    #[must_use]
    #[inline(always)]
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
