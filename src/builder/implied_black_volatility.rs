use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

/// Builder-backed container that represents the inputs required to compute
/// the **implied Black volatility** for an undiscounted European option.
///
/// Use the generated `ImpliedBlackVolatility::builder()` to construct this
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
/// inversion `BS(F, K, T, σ) = P` to find the implied volatility `σ`.
#[derive(Builder)]
#[builder(derive(Clone, Debug), finish_fn(vis = "", name = build_internal), const)]
pub struct ImpliedBlackVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

impl<S: implied_black_volatility_builder::IsComplete> ImpliedBlackVolatilityBuilder<S> {
    /// Validate inputs and build an `ImpliedBlackVolatility`.
    ///
    /// Performs the following validation checks:
    /// - `forward` must be a finite number and greater than 0.
    /// - `strike` must be a finite number and greater than 0.
    /// - `expiry` must be finite and non-negative (`T >= 0`).
    /// - `option_price` must be a finite, non-negative number.
    ///
    /// Returns `Some(ImpliedBlackVolatility)` when all checks pass, otherwise `None`.
    ///
    /// # Rationale
    /// These checks ensure the constructed object lies within the mathematical domain
    /// required by the Black–Scholes pricing function. Use `build_unchecked()` to skip validation.
    pub fn build(self) -> Option<ImpliedBlackVolatility> {
        let implied_black_volatility = self.build_internal();

        if implied_black_volatility.forward.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
            || implied_black_volatility.forward.is_infinite()
        {
            return None;
        }
        if implied_black_volatility.strike.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
            || implied_black_volatility.strike.is_infinite()
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

    /// Build without performing any validation.
    ///
    /// This constructor constructs the `ImpliedBlackVolatility` directly from
    /// the builder's fields and does **not** check for NaNs, infinities, or
    /// sign constraints. Use only when you are certain the inputs are valid
    /// or when you want to avoid the cost of runtime validation.
    pub fn build_unchecked(self) -> ImpliedBlackVolatility {
        self.build_internal()
    }
}

impl ImpliedBlackVolatility {
    /// Compute the implied Black volatility `σ` for the stored inputs.
    ///
    /// Returns:
    /// - `Some(σ)` if an implied volatility consistent with `BS(F, K, T, σ) = P`
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

#[cfg(test)]
mod tests {
    use crate::DefaultSpecialFn;
    use crate::builder::implied_black_volatility::ImpliedBlackVolatility;

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
