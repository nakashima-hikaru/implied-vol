use super::{Set, Unset};
use crate::{SpecialFn, lets_be_rational};
use std::marker::PhantomData;
use std::ops::Neg;

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

    /// Compute the total implied volatility ($v$) for the stored inputs.
    ///
    /// # Returns
    /// `Some(v)` if it converges, otherwise `None`.
    #[must_use]
    #[inline]
    pub fn calculate<SpFn: SpecialFn>(&self) -> Option<f64> {
        let x = self.log_moneyness;
        let b = self.normalised_price;

        if b <= 0.0 {
            return (b == 0.0).then_some(0.0);
        }

        if x == 0.0 {
            return Some(lets_be_rational::implied_normalised_volatility_atm::<SpFn>(b));
        }

        let theta_x = x.abs().neg();
        let b_max = (0.5 * theta_x).exp();
        
        if b >= b_max {
             return (b == b_max).then_some(f64::INFINITY);
        }

        Some(lets_be_rational::lets_be_rational_unchecked::<SpFn>(b, theta_x, b_max))
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

impl<LogMoneyness, NormalisedPrice> ImpliedBlackVolatilityNormalisedBuilder<LogMoneyness, NormalisedPrice> {
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
    use crate::DefaultSpecialFn;
    use super::*;

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
}
