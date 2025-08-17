use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

#[derive(Builder)]
#[builder(const, derive(Clone, Debug),
finish_fn(name = build_unchecked,
doc{
/// Build without performing any validation.
///
/// This constructor constructs the `ImpliedNormalVolatility` directly from
/// the builder's fields and does **not** check for NaNs, infinities, or
/// sign constraints. Use only when you are certain the inputs are valid
/// or when you want to avoid the cost of runtime validation.
})
)]
pub struct ImpliedNormalVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

impl<S: implied_normal_volatility_builder::IsComplete> ImpliedNormalVolatilityBuilder<S> {
    pub const fn build(self) -> Option<ImpliedNormalVolatility> {
        let implied_normal_volatility = self.build_unchecked();
        if !implied_normal_volatility.forward.is_finite() {
            return None;
        }
        if !implied_normal_volatility.strike.is_finite() {
            return None;
        }
        if !(implied_normal_volatility.expiry >= 0.0_f64) {
            return None;
        }
        if !(implied_normal_volatility.option_price >= 0.0_f64)
            || implied_normal_volatility.option_price.is_infinite()
        {
            return None;
        }
        Some(implied_normal_volatility)
    }
}

impl ImpliedNormalVolatility {
    #[must_use]
    #[inline(always)]
    pub fn calculate<SpFn: SpecialFn>(&self) -> Option<f64> {
        if self.is_call {
            lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<SpFn, true>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        } else {
            lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<SpFn, false>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        }
    }
}
