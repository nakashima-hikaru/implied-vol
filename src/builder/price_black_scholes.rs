use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

/// Builder-backed container for computing undiscounted European option prices
/// under the Black–Scholes model.
///
/// Construct instances with `PriceBlackScholes::builder()`. The builder performs
/// basic domain validation when you call `.build()`. Use `.build_unchecked()` to
/// skip validation when you know inputs are already valid.
///
/// Fields:
/// - `forward`: forward price of the underlying (F). Must be finite.
/// - `strike`: strike price (K). Must be finite.
/// - `volatility`: volatility (σ). Must be finite and `σ >= 0`.
/// - `expiry`: time to expiry (T). Must be finite and `T >= 0`.
/// - `is_call`: `true` to price a call option, `false` to price a put option.
///
/// The `calculate::<SpFn>()` method performs the numerical evaluation and uses a
/// `SpecialFn` implementation for any special-function approximations required
/// by the numerical routines.
#[derive(Builder)]
#[builder(derive(Clone, Debug), finish_fn(vis = "", name = build_internal))]
pub struct PriceBlackScholes {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
}

impl<S: price_black_scholes_builder::IsComplete> PriceBlackScholesBuilder<S> {
    /// Validate builder inputs and construct `PriceBlackScholes`.
    ///
    /// Validation performed:
    /// - `forward` must be finite.
    /// - `strike` must be finite.
    /// - `volatility` must be non-negative (`σ >= 0`) and finite.
    /// - `expiry` must be non-negative (`T >= 0`) and not NaN.
    ///
    /// Returns `Some(PriceBlackScholes)` when all checks pass; otherwise returns
    /// `None`.
    pub fn build(self) -> Option<PriceBlackScholes> {
        let price_black_scholes = self.build_internal();
        if !price_black_scholes.forward.is_finite() || price_black_scholes.forward <= 0.0 {
            return None;
        }
        if !price_black_scholes.strike.is_finite() || price_black_scholes.strike <= 0.0 {
            return None;
        }
        if matches!(
            price_black_scholes.volatility.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || !price_black_scholes.volatility.is_finite()
        {
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

    /// Construct `PriceBlackScholes` without performing validation checks.
    ///
    /// Use this when you have externally guaranteed that the inputs are valid
    /// or when you want to avoid the runtime cost of validation.
    pub fn build_unchecked(self) -> PriceBlackScholes {
        self.build_internal()
    }
}

impl PriceBlackScholes {
    /// Compute the undiscounted Black–Scholes option price for the stored inputs.
    ///
    /// # Type parameter
    /// - `SpFn: SpecialFn` — implementation used for internal special-function
    ///   evaluations (e.g., approximations used by the algorithm). Use the crate's
    ///   `DefaultSpecialFn` for the default behavior or provide a custom
    ///   implementation to change numerical characteristics.
    ///
    /// # Returns
    /// The computed undiscounted option price as `f64`.
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
