use super::{Set, Unset};
use crate::{SpecialFn, lets_be_rational};
use std::marker::PhantomData;

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
/// - `volatility`: volatility (sigma). Must be finite and `sigma >= 0`.
/// - `expiry`: time to expiry (T). Must be finite and `T >= 0`.
/// - `is_call`: `true` to price a call option, `false` to price a put option.
///
/// The `calculate::<SpFn>()` method performs the numerical evaluation and uses a
/// `SpecialFn` implementation for any special-function approximations required
/// by the numerical routines.
pub struct PriceBlackScholes {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
}

#[derive(Clone, Debug)]
pub struct PriceBlackScholesBuilder<
    Forward = Unset,
    Strike = Unset,
    Volatility = Unset,
    Expiry = Unset,
    IsCall = Unset,
> {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
    _marker: PhantomData<(Forward, Strike, Volatility, Expiry, IsCall)>,
}

impl PriceBlackScholes {
    #[must_use]
    #[inline(always)]
    pub const fn builder() -> PriceBlackScholesBuilder {
        PriceBlackScholesBuilder::new()
    }

    /// Compute the undiscounted Black–Scholes option price for the stored inputs.
    ///
    /// # Type parameter
    /// - `SpFn: SpecialFn` — implementation used for internal special-function
    ///   evaluations. Use the crate's `DefaultSpecialFn` for the default behavior or provide a custom
    ///   implementation to change numerical characteristics.
    ///
    /// # Returns
    /// The computed undiscounted European option price.
    #[must_use]
    #[inline]
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

impl PriceBlackScholesBuilder {
    #[must_use]
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            forward: 0.0,
            strike: 0.0,
            volatility: 0.0,
            expiry: 0.0,
            is_call: false,
            _marker: PhantomData,
        }
    }
}

impl<Forward, Strike, Volatility, Expiry, IsCall>
    PriceBlackScholesBuilder<Forward, Strike, Volatility, Expiry, IsCall>
{
    #[must_use]
    #[inline(always)]
    pub const fn forward(
        self,
        forward: f64,
    ) -> PriceBlackScholesBuilder<Set, Strike, Volatility, Expiry, IsCall> {
        PriceBlackScholesBuilder {
            forward,
            strike: self.strike,
            volatility: self.volatility,
            expiry: self.expiry,
            is_call: self.is_call,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn strike(
        self,
        strike: f64,
    ) -> PriceBlackScholesBuilder<Forward, Set, Volatility, Expiry, IsCall> {
        PriceBlackScholesBuilder {
            forward: self.forward,
            strike,
            volatility: self.volatility,
            expiry: self.expiry,
            is_call: self.is_call,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn volatility(
        self,
        volatility: f64,
    ) -> PriceBlackScholesBuilder<Forward, Strike, Set, Expiry, IsCall> {
        PriceBlackScholesBuilder {
            forward: self.forward,
            strike: self.strike,
            volatility,
            expiry: self.expiry,
            is_call: self.is_call,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn expiry(
        self,
        expiry: f64,
    ) -> PriceBlackScholesBuilder<Forward, Strike, Volatility, Set, IsCall> {
        PriceBlackScholesBuilder {
            forward: self.forward,
            strike: self.strike,
            volatility: self.volatility,
            expiry,
            is_call: self.is_call,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    #[allow(clippy::wrong_self_convention)]
    pub const fn is_call(
        self,
        is_call: bool,
    ) -> PriceBlackScholesBuilder<Forward, Strike, Volatility, Expiry, Set> {
        PriceBlackScholesBuilder {
            forward: self.forward,
            strike: self.strike,
            volatility: self.volatility,
            expiry: self.expiry,
            is_call,
            _marker: PhantomData,
        }
    }
}

impl PriceBlackScholesBuilder<Set, Set, Set, Set, Set> {
    /// Build without performing any validation.
    ///
    /// This constructor constructs the `PriceBlackScholes` directly from
    /// the builder's fields and does **not** check for NaNs, infinities, or
    /// sign constraints. Use only when you are certain the inputs are valid
    /// or when you want to avoid the cost of runtime validation.
    #[must_use]
    #[inline(always)]
    pub const fn build_unchecked(self) -> PriceBlackScholes {
        PriceBlackScholes {
            forward: self.forward,
            strike: self.strike,
            volatility: self.volatility,
            expiry: self.expiry,
            is_call: self.is_call,
        }
    }

    /// Validate builder inputs and construct `PriceBlackScholes`.
    ///
    /// Validation performed:
    /// - `forward` must be positive and finite.
    /// - `strike` must be positive and finite.
    /// - `volatility` must be non-negative (`sigma >= 0`) but can be positive infinite.
    /// - `expiry` must be non-negative (`T >= 0`) but can be positive infinite.
    ///
    /// Returns `Some(PriceBlackScholes)` when all checks pass; otherwise returns
    /// `None`.
    #[must_use]
    #[inline(always)]
    pub const fn build(self) -> Option<PriceBlackScholes> {
        let price_black_scholes = self.build_unchecked();
        if !price_black_scholes.forward.is_finite() || !(price_black_scholes.forward > 0.0) {
            return None;
        }
        if !price_black_scholes.strike.is_finite() || !(price_black_scholes.strike > 0.0) {
            return None;
        }
        if !(price_black_scholes.volatility >= 0.0) {
            return None;
        }
        if !(price_black_scholes.expiry >= 0.0) {
            return None;
        }
        Some(price_black_scholes)
    }
}
