use super::{Set, Unset};
use crate::{SpecialFn, lets_be_rational};
use std::marker::PhantomData;

pub struct ImpliedNormalVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

#[derive(Clone, Debug)]
pub struct ImpliedNormalVolatilityBuilder<
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

impl ImpliedNormalVolatility {
    #[must_use]
    #[inline(always)]
    pub const fn builder() -> ImpliedNormalVolatilityBuilder {
        ImpliedNormalVolatilityBuilder::new()
    }

    #[must_use]
    #[inline]
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

impl ImpliedNormalVolatilityBuilder {
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
    ImpliedNormalVolatilityBuilder<Forward, Strike, Expiry, IsCall, OptionPrice>
{
    #[must_use]
    #[inline(always)]
    pub const fn forward(
        self,
        forward: f64,
    ) -> ImpliedNormalVolatilityBuilder<Set, Strike, Expiry, IsCall, OptionPrice> {
        ImpliedNormalVolatilityBuilder {
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
    ) -> ImpliedNormalVolatilityBuilder<Forward, Set, Expiry, IsCall, OptionPrice> {
        ImpliedNormalVolatilityBuilder {
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
    ) -> ImpliedNormalVolatilityBuilder<Forward, Strike, Set, IsCall, OptionPrice> {
        ImpliedNormalVolatilityBuilder {
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
    ) -> ImpliedNormalVolatilityBuilder<Forward, Strike, Expiry, Set, OptionPrice> {
        ImpliedNormalVolatilityBuilder {
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
    ) -> ImpliedNormalVolatilityBuilder<Forward, Strike, Expiry, IsCall, Set> {
        ImpliedNormalVolatilityBuilder {
            forward: self.forward,
            strike: self.strike,
            expiry: self.expiry,
            is_call: self.is_call,
            option_price,
            _marker: PhantomData,
        }
    }
}

impl ImpliedNormalVolatilityBuilder<Set, Set, Set, Set, Set> {
    #[must_use]
    #[inline(always)]
    pub const fn build_unchecked(self) -> ImpliedNormalVolatility {
        ImpliedNormalVolatility {
            forward: self.forward,
            strike: self.strike,
            expiry: self.expiry,
            is_call: self.is_call,
            option_price: self.option_price,
        }
    }

    #[must_use]
    #[inline(always)]
    pub const fn build(self) -> Option<ImpliedNormalVolatility> {
        let implied_normal_volatility = self.build_unchecked();
        if !implied_normal_volatility.forward.is_finite() {
            return None;
        }
        if !implied_normal_volatility.strike.is_finite() {
            return None;
        }
        if !(implied_normal_volatility.expiry >= 0.0) {
            return None;
        }
        if !(implied_normal_volatility.option_price >= 0.0)
            || implied_normal_volatility.option_price.is_infinite()
        {
            return None;
        }
        Some(implied_normal_volatility)
    }
}
