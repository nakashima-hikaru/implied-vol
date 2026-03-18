use super::{Set, Unset};
use crate::{SpecialFn, lets_be_rational};
use std::marker::PhantomData;

/// Builder-backed container for computing the **normalised** Black–Scholes option price.
///
/// In this context, "normalised" refers to the formulation used in Jäckel's "LetsBeRational":
/// - `log_moneyness` ($x$) is $\ln(F/K)$.
/// - `total_volatility` ($v$) is $\sigma \sqrt{T}$.
/// - The result is the **normalised price** ($b$), which represents the time value
///   divided by $\sqrt{FK}$.
///
/// To get the actual undiscounted price:
/// $Price = \sqrt{FK} \cdot (\max(e^{x/2} - e^{-x/2}, 0) + b)$ if it is a call.
/// 
/// However, note that $b$ is symmetric for calls and puts of the same strike.
pub struct PriceBlackScholesNormalised {
    log_moneyness: f64,
    total_volatility: f64,
}

#[derive(Clone, Debug)]
pub struct PriceBlackScholesNormalisedBuilder<LogMoneyness = Unset, TotalVolatility = Unset> {
    log_moneyness: f64,
    total_volatility: f64,
    _marker: PhantomData<(LogMoneyness, TotalVolatility)>,
}

impl PriceBlackScholesNormalised {
    #[must_use]
    #[inline(always)]
    pub const fn builder() -> PriceBlackScholesNormalisedBuilder {
        PriceBlackScholesNormalisedBuilder::new()
    }

    /// Compute the normalised Black–Scholes price $b$ for the stored inputs.
    ///
    /// # Returns
    /// The computed normalised price $b$ (time value / $\sqrt{FK}$).
    #[must_use]
    #[inline]
    pub fn calculate<SpFn: SpecialFn>(&self) -> f64 {
        let x = self.log_moneyness;
        let v = self.total_volatility;

        if v <= 0.0 {
            return 0.0;
        }

        if x == 0.0 {
            return SpFn::erf((0.5 * std::f64::consts::FRAC_1_SQRT_2) * v);
        }

        // Internal implementation uses half_theta_x = 0.5 * x, h = x / v, t = 0.5 * v
        // where theta_x is always non-positive in the core routine for symmetry.
        let theta_x = -x.abs();
        lets_be_rational::bs_option_price::normalised_black::<SpFn>(
            0.5 * theta_x,
            theta_x / v,
            0.5 * v,
        )
    }
}

impl PriceBlackScholesNormalisedBuilder {
    #[must_use]
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            log_moneyness: 0.0,
            total_volatility: 0.0,
            _marker: PhantomData,
        }
    }
}

impl<LogMoneyness, TotalVolatility> PriceBlackScholesNormalisedBuilder<LogMoneyness, TotalVolatility> {
    #[must_use]
    #[inline(always)]
    pub const fn log_moneyness(
        self,
        log_moneyness: f64,
    ) -> PriceBlackScholesNormalisedBuilder<Set, TotalVolatility> {
        PriceBlackScholesNormalisedBuilder {
            log_moneyness,
            total_volatility: self.total_volatility,
            _marker: PhantomData,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn total_volatility(
        self,
        total_volatility: f64,
    ) -> PriceBlackScholesNormalisedBuilder<LogMoneyness, Set> {
        PriceBlackScholesNormalisedBuilder {
            log_moneyness: self.log_moneyness,
            total_volatility,
            _marker: PhantomData,
        }
    }
}

impl PriceBlackScholesNormalisedBuilder<Set, Set> {
    /// Build without performing any validation.
    #[must_use]
    #[inline(always)]
    pub const fn build_unchecked(self) -> PriceBlackScholesNormalised {
        PriceBlackScholesNormalised {
            log_moneyness: self.log_moneyness,
            total_volatility: self.total_volatility,
        }
    }

    /// Validate builder inputs and construct `PriceBlackScholesNormalised`.
    ///
    /// Validation:
    /// - `log_moneyness` must be finite.
    /// - `total_volatility` must be non-negative and finite.
    #[must_use]
    #[inline(always)]
    pub const fn build(self) -> Option<PriceBlackScholesNormalised> {
        let price = self.build_unchecked();
        if !price.log_moneyness.is_finite() {
            return None;
        }
        if !price.total_volatility.is_finite() || !(price.total_volatility >= 0.0) {
            return None;
        }
        Some(price)
    }
}

#[cfg(test)]
mod tests {
    use crate::DefaultSpecialFn;
    use super::*;

    #[test]
    fn test_normalised_price_atm() {
        let p = PriceBlackScholesNormalised::builder()
            .log_moneyness(0.0)
            .total_volatility(0.2)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>();
        
        // ATM price (b) should be roughly erf(v/(2*sqrt(2)))
        assert!(p > 0.0);
        assert!(p < 1.0);
    }
}
