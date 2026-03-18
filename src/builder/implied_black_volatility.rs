use super::{Set, Unset};
use crate::{SpecialFn, lets_be_rational};
use std::marker::PhantomData;

/// Builder-backed container that represents the inputs required to compute
/// the **implied Black volatility** for an undiscounted European option.
///
/// Use `ImpliedBlackVolatility::builder()` to construct this
/// type. The builder performs input validation when you call `.build()`.
/// Or use `.build_unchecked()` to skip validation.
///
/// Fields:
/// - `forward`: forward price of the underlying (F). Must be finite and `> 0`.
/// - `strike`: strike price (K). Must be finite and `> 0`.
/// - `expiry`: time to expiry (T). Must be finite and `>= 0`.
/// - `is_call`: `true` for a call option, `false` for a put option.
/// - `option_price`: observed undiscounted option price (P). Must be finite and `>= 0`.
///
/// This struct is consumed by `calculate::<SpFn>()` which performs the numerical
/// inversion `BS(F, K, T, sigma) = P` to find the implied volatility `sigma`.
pub struct ImpliedBlackVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

#[derive(Clone, Debug)]
pub struct ImpliedBlackVolatilityBuilder<
    Forward = Unset,
    Strike = Unset,
    Expiry = Unset,
    IsCall = Unset,
    OptionPrice = Unset,
> {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
    _marker: PhantomData<(Forward, Strike, Expiry, IsCall, OptionPrice)>,
}

impl ImpliedBlackVolatility {
    #[must_use]
    #[inline(always)]
    pub const fn builder() -> ImpliedBlackVolatilityBuilder {
        ImpliedBlackVolatilityBuilder::new()
    }

    /// Compute the implied Black volatility `sigma` for the stored inputs.
    ///
    /// Returns:
    /// - `Some(sigma)` if an implied volatility consistent with `BS(F, K, T, sigma) = P`
    ///   exists, and the numerical routine converges to a finite value.
    /// - `None` if the given `option_price` is outside the attainable range for
    ///   the supplied model parameters.
    ///
    /// # Type parameter
    /// - `SpFn: SpecialFn` — implementation of the `SpecialFn` trait used for
    ///   internal special-function computations. Use `DefaultSpecialFn` for the
    ///   crate-provided default behavior, or supply your own implementation
    ///   to change numerical characteristics.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};
    ///
    /// let iv = ImpliedBlackVolatility::builder()
    ///     .option_price(10.0)
    ///     .forward(100.0)
    ///     .strike(100.0)
    ///     .expiry(1.0)
    ///     .is_call(true)
    ///     .build().unwrap();
    ///
    /// let sigma = iv.calculate::<DefaultSpecialFn>().unwrap();
    /// assert!(sigma.is_finite());
    /// ```
    #[must_use]
    #[inline(always)]
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

impl ImpliedBlackVolatilityBuilder {
    #[must_use]
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            forward: 0.0,
            strike: 0.0,
            expiry: 0.0,
            is_call: false,
            option_price: 0.0,
            _marker: PhantomData,
        }
    }
}

impl<Forward, Strike, Expiry, IsCall, OptionPrice>
    ImpliedBlackVolatilityBuilder<Forward, Strike, Expiry, IsCall, OptionPrice>
{
    #[must_use]
    #[inline(always)]
    pub const fn forward(
        self,
        forward: f64,
    ) -> ImpliedBlackVolatilityBuilder<Set, Strike, Expiry, IsCall, OptionPrice> {
        ImpliedBlackVolatilityBuilder {
            forward,
            strike: self.strike,
            expiry: self.expiry,
            is_call: self.is_call,
            option_price: self.option_price,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn strike(
        self,
        strike: f64,
    ) -> ImpliedBlackVolatilityBuilder<Forward, Set, Expiry, IsCall, OptionPrice> {
        ImpliedBlackVolatilityBuilder {
            forward: self.forward,
            strike,
            expiry: self.expiry,
            is_call: self.is_call,
            option_price: self.option_price,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn expiry(
        self,
        expiry: f64,
    ) -> ImpliedBlackVolatilityBuilder<Forward, Strike, Set, IsCall, OptionPrice> {
        ImpliedBlackVolatilityBuilder {
            forward: self.forward,
            strike: self.strike,
            expiry,
            is_call: self.is_call,
            option_price: self.option_price,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    #[allow(clippy::wrong_self_convention)]
    pub const fn is_call(
        self,
        is_call: bool,
    ) -> ImpliedBlackVolatilityBuilder<Forward, Strike, Expiry, Set, OptionPrice> {
        ImpliedBlackVolatilityBuilder {
            forward: self.forward,
            strike: self.strike,
            expiry: self.expiry,
            is_call,
            option_price: self.option_price,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn option_price(
        self,
        option_price: f64,
    ) -> ImpliedBlackVolatilityBuilder<Forward, Strike, Expiry, IsCall, Set> {
        ImpliedBlackVolatilityBuilder {
            forward: self.forward,
            strike: self.strike,
            expiry: self.expiry,
            is_call: self.is_call,
            option_price,
            _marker: PhantomData,
        }
    }
}

impl ImpliedBlackVolatilityBuilder<Set, Set, Set, Set, Set> {
    /// Build without performing any validation.
    ///
    /// This constructor constructs the `ImpliedBlackVolatility` directly from
    /// the builder's fields and does **not** check for NaNs, infinities, or
    /// sign constraints. Use only when you are certain the inputs are valid
    /// or when you want to avoid the cost of runtime validation.
    #[must_use]
    #[inline(always)]
    pub const fn build_unchecked(self) -> ImpliedBlackVolatility {
        ImpliedBlackVolatility {
            forward: self.forward,
            strike: self.strike,
            expiry: self.expiry,
            is_call: self.is_call,
            option_price: self.option_price,
        }
    }

    /// Validate inputs and build an `ImpliedBlackVolatility`.
    ///
    /// Performs the following validation checks:
    /// - `forward` must be positive and finite.
    /// - `strike` must be positive and finite.
    /// - `expiry` must be non-negative (`T >= 0`) but can be positive infinite.
    /// - `option_price` must be a finite, non-negative number.
    ///
    /// Returns `Some(ImpliedBlackVolatility)` when all checks pass, otherwise `None`.
    ///
    /// # Rationale
    /// These checks ensure the constructed object lies within the mathematical domain
    /// required by the Black–Scholes pricing function. Use `build_unchecked()` to skip validation.
    #[must_use]
    #[inline(always)]
    pub const fn build(self) -> Option<ImpliedBlackVolatility> {
        let implied_black_volatility = self.build_unchecked();

        if !(implied_black_volatility.forward > 0.0)
            || implied_black_volatility.forward.is_infinite()
        {
            return None;
        }
        if !(implied_black_volatility.strike > 0.0) || implied_black_volatility.strike.is_infinite()
        {
            return None;
        }
        if !(implied_black_volatility.expiry >= 0.0) {
            return None;
        }
        if !(implied_black_volatility.option_price >= 0.0)
            || implied_black_volatility.option_price.is_infinite()
        {
            return None;
        }
        Some(implied_black_volatility)
    }
}

#[cfg(test)]
mod tests {
    use crate::DefaultSpecialFn;
    use crate::builder::implied_black_volatility::ImpliedBlackVolatility;

    #[test]
    const fn normal_const() {
        const PRICE: f64 = 10.0;
        const F: f64 = 100.0;
        const K: f64 = 100.0;
        const T: f64 = 1.0;
        const Q: bool = true;

        const IV_BUILDER: Option<ImpliedBlackVolatility> = ImpliedBlackVolatility::builder()
            .option_price(PRICE)
            .forward(F)
            .strike(K)
            .expiry(T)
            .is_call(Q)
            .build();
        assert!(IV_BUILDER.is_some());
    }

    #[test]
    fn strike_anomaly() {
        const Q: bool = true;
        for k in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0] {
            let price = 100.0;
            let f = 100.0;
            let t = 1.0;
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
    fn strike_boundary() {
        const Q: bool = true;
        let price = 100.0;
        let f = 100.0;
        let k = f64::MIN_POSITIVE;
        let t = 1.0;
        assert!(
            ImpliedBlackVolatility::builder()
                .option_price(price)
                .forward(f)
                .strike(k)
                .expiry(t)
                .is_call(Q)
                .build()
                .is_some()
        );
    }

    #[test]
    fn forward_anomaly() {
        const Q: bool = true;
        for f in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0] {
            let price = 100.0;
            let k = 100.0;
            let t = 1.0;
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
        const Q: bool = true;
        for price in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let f = 100.0;
            let t = 1.0;
            let k = 100.0;
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
    fn price_below_intrinsic() {
        const Q: bool = true;
        let price = 10.0;
        let f = 120.0;
        let t = 1.0;
        let k = 100.0;
        let iv_builder = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(Q)
            .build();
        assert!(iv_builder.is_some());
        let iv = iv_builder.unwrap();
        assert!(iv.calculate::<DefaultSpecialFn>().is_none());
    }

    #[test]
    fn time_anomaly() {
        const Q: bool = true;
        for t in [f64::NAN, f64::NEG_INFINITY] {
            let price = 10.0;
            let f = 100.0;
            let k = 100.0;
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
    fn time_zero() {
        // the price is below intrinsic
        let price = 10.0;
        let f = 120.0;
        let k = 100.0;
        let t = 0.0;
        let q = true;

        let vol = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        assert!(vol.is_none());

        let price = 20.0;
        let f = 120.0;
        let k = 100.0;
        let t = 0.0;

        let vol = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(q)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        assert_eq!(vol.unwrap(), 0.0);
    }

    #[test]
    fn time_inf() {
        const Q: bool = true;
        let price = 10.0;
        let f = 100.0;
        let k = 100.0;
        let t = f64::INFINITY;

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
