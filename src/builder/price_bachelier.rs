use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

#[derive(Builder)]
#[builder(derive(Clone, Debug), finish_fn(vis = "", name = build_internal), const)]
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
    pub fn build_unchecked(self) -> PriceBachelier {
        self.build_internal()
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
