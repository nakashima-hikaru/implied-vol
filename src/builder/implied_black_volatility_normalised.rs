use super::{Set, Unset};
use crate::{
    SpecialFn, explicit_black,
    solver::{BlackSolver, Hybrid},
};
use std::marker::PhantomData;

/// Builder-backed container for computing the **normalised implied Black volatility**.
///
/// In this context:
/// - `log_moneyness` ($x$) is $\ln(F/K)$.
/// - `normalised_price` ($b$) is the time value divided by $\sqrt{FK}$.
/// - The result is the **total volatility** ($v = \sigma \sqrt{T}$).
pub struct ImpliedBlackVolatilityNormalised {
    log_moneyness: f64,
    normalised_price: f64,
}

#[derive(Clone, Debug)]
pub struct ImpliedBlackVolatilityNormalisedBuilder<LogMoneyness = Unset, NormalisedPrice = Unset> {
    log_moneyness: f64,
    normalised_price: f64,
    _marker: PhantomData<(LogMoneyness, NormalisedPrice)>,
}

impl ImpliedBlackVolatilityNormalised {
    #[must_use]
    #[inline(always)]
    pub const fn builder() -> ImpliedBlackVolatilityNormalisedBuilder {
        ImpliedBlackVolatilityNormalisedBuilder::new()
    }

    /// Compute total implied volatility with the default hybrid solver.
    ///
    /// Solver features do not change this selection. Use [`Self::calculate_with`] for
    /// an explicit solver. Returns `None` for unattainable prices.
    #[must_use]
    #[inline]
    pub fn calculate<SpFn: SpecialFn>(&self) -> Option<f64> {
        self.calculate_with::<Hybrid<SpFn>>()
    }

    /// Compute total implied volatility with the selected solver type.
    ///
    /// Each solver owns its ATM, price-boundary, and convergence behavior.
    /// `Experimental` uses a compensated exact-cap comparison and its own
    /// special functions. Solver features can be enabled together.
    ///
    /// ```
    /// use implied_vol::{ImpliedBlackVolatilityNormalised, solver::Jaeckel};
    /// let iv = ImpliedBlackVolatilityNormalised::builder()
    ///     .log_moneyness(-0.1).normalised_price(0.05)
    ///     .build().unwrap();
    /// assert!(iv.calculate_with::<Jaeckel>().unwrap().is_finite());
    /// ```
    #[must_use]
    #[inline]
    pub fn calculate_with<S: BlackSolver>(&self) -> Option<f64> {
        S::implied_total_volatility(self.log_moneyness, self.normalised_price)
    }

    /// Compute the total implied volatility ($v$) via the inverse-Gaussian explicit formula
    /// from arXiv:2604.24480.
    ///
    /// Returns `None` if the intermediate inverse-Gaussian quantile is unrepresentable.
    #[must_use]
    #[inline]
    pub fn calculate_explicit<SpFn: SpecialFn>(&self) -> Option<f64> {
        explicit_black::implied_black_volatility_normalised::<SpFn>(
            self.log_moneyness,
            self.normalised_price,
        )
    }
}

impl ImpliedBlackVolatilityNormalisedBuilder {
    #[must_use]
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            log_moneyness: 0.0,
            normalised_price: 0.0,
            _marker: PhantomData,
        }
    }
}

impl<LogMoneyness, NormalisedPrice>
    ImpliedBlackVolatilityNormalisedBuilder<LogMoneyness, NormalisedPrice>
{
    #[must_use]
    #[inline(always)]
    pub const fn log_moneyness(
        self,
        log_moneyness: f64,
    ) -> ImpliedBlackVolatilityNormalisedBuilder<Set, NormalisedPrice> {
        ImpliedBlackVolatilityNormalisedBuilder {
            log_moneyness,
            normalised_price: self.normalised_price,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn normalised_price(
        self,
        normalised_price: f64,
    ) -> ImpliedBlackVolatilityNormalisedBuilder<LogMoneyness, Set> {
        ImpliedBlackVolatilityNormalisedBuilder {
            log_moneyness: self.log_moneyness,
            normalised_price,
            _marker: PhantomData,
        }
    }
}

impl ImpliedBlackVolatilityNormalisedBuilder<Set, Set> {
    /// Build without performing any validation.
    #[must_use]
    #[inline(always)]
    pub const fn build_unchecked(self) -> ImpliedBlackVolatilityNormalised {
        ImpliedBlackVolatilityNormalised {
            log_moneyness: self.log_moneyness,
            normalised_price: self.normalised_price,
        }
    }

    /// Validate builder inputs and construct `ImpliedBlackVolatilityNormalised`.
    #[must_use]
    #[inline(always)]
    pub const fn build(self) -> Option<ImpliedBlackVolatilityNormalised> {
        let iv = self.build_unchecked();
        if !iv.log_moneyness.is_finite() {
            return None;
        }
        if !iv.normalised_price.is_finite() || !(iv.normalised_price >= 0.0) {
            return None;
        }
        Some(iv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DefaultSpecialFn;

    #[test]
    fn atm_price_upper_bound_is_checked_before_inversion() {
        let upper = ImpliedBlackVolatilityNormalised::builder()
            .log_moneyness(0.0)
            .normalised_price(1.0)
            .build()
            .unwrap();
        assert_eq!(upper.calculate::<DefaultSpecialFn>(), Some(f64::INFINITY));
        assert_eq!(
            upper.calculate_explicit::<DefaultSpecialFn>(),
            Some(f64::INFINITY)
        );

        let above_upper = ImpliedBlackVolatilityNormalised::builder()
            .log_moneyness(0.0)
            .normalised_price(1.01)
            .build()
            .unwrap();
        assert!(above_upper.calculate::<DefaultSpecialFn>().is_none());
        assert!(
            above_upper
                .calculate_explicit::<DefaultSpecialFn>()
                .is_none()
        );
    }

    #[test]
    fn test_normalised_iv_roundtrip() {
        let x = 0.1;
        let v = 0.2;
        let b = crate::PriceBlackScholesNormalised::builder()
            .log_moneyness(x)
            .total_volatility(v)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();

        let v2 = ImpliedBlackVolatilityNormalised::builder()
            .log_moneyness(x)
            .normalised_price(b)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>()
            .unwrap();

        assert!((v - v2).abs() < 1e-12);
    }

    #[test]
    fn test_normalised_iv_roundtrip_explicit() {
        let x = -0.35;
        let v = 0.6;
        let b = crate::PriceBlackScholesNormalised::builder()
            .log_moneyness(x)
            .total_volatility(v)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();

        let v2 = ImpliedBlackVolatilityNormalised::builder()
            .log_moneyness(x)
            .normalised_price(b)
            .build()
            .unwrap()
            .calculate_explicit::<DefaultSpecialFn>()
            .unwrap();

        assert!((v - v2).abs() < 1e-12);
    }
}
