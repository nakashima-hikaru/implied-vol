use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

#[derive(Builder)]
#[builder(const, derive(Clone, Debug),
finish_fn(name = build_unchecked,
doc{
/// Build without performing any validation.
///
/// This constructor constructs the `PriceBachelier` directly from
/// the builder's fields and does **not** check for NaNs, infinities, or
/// sign constraints. Use only when you are certain the inputs are valid
/// or when you want to avoid the cost of runtime validation.
})
)]
pub struct PriceBachelier {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
}

impl<S: price_bachelier_builder::IsComplete> PriceBachelierBuilder<S> {
    pub const fn build(self) -> Option<PriceBachelier> {
        let price_bachelier = self.build_unchecked();
        if !price_bachelier.forward.is_finite() {
            return None;
        }
        if !price_bachelier.strike.is_finite() {
            return None;
        }
        if !(price_bachelier.volatility >= 0.0) {
            return None;
        }
        if !(price_bachelier.expiry >= 0.0) {
            return None;
        }
        Some(price_bachelier)
    }
}

impl PriceBachelier {
    #[must_use]
    #[inline(always)]
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
