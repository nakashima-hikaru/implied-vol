use super::{Set, Unset};
use crate::{SpecialFn, lets_be_rational};
use std::marker::PhantomData;

pub struct PriceBachelier {
    forward: f64,
    strike: f64,
    volatility: f64,
    expiry: f64,
    is_call: bool,
}

#[derive(Clone, Debug)]
pub struct PriceBachelierBuilder<
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

impl PriceBachelier {
    #[must_use]
    #[inline(always)]
    pub const fn builder() -> PriceBachelierBuilder {
        PriceBachelierBuilder::new()
    }

    #[must_use]
    #[inline]
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

impl PriceBachelierBuilder {
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
    PriceBachelierBuilder<Forward, Strike, Volatility, Expiry, IsCall>
{
    #[must_use]
    #[inline(always)]
    pub const fn forward(
        self,
        forward: f64,
    ) -> PriceBachelierBuilder<Set, Strike, Volatility, Expiry, IsCall> {
        PriceBachelierBuilder {
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
    ) -> PriceBachelierBuilder<Forward, Set, Volatility, Expiry, IsCall> {
        PriceBachelierBuilder {
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
    ) -> PriceBachelierBuilder<Forward, Strike, Set, Expiry, IsCall> {
        PriceBachelierBuilder {
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
    ) -> PriceBachelierBuilder<Forward, Strike, Volatility, Set, IsCall> {
        PriceBachelierBuilder {
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
    ) -> PriceBachelierBuilder<Forward, Strike, Volatility, Expiry, Set> {
        PriceBachelierBuilder {
            forward: self.forward,
            strike: self.strike,
            volatility: self.volatility,
            expiry: self.expiry,
            is_call,
            _marker: PhantomData,
        }
    }
}

impl PriceBachelierBuilder<Set, Set, Set, Set, Set> {
    #[must_use]
    #[inline(always)]
    pub const fn build_unchecked(self) -> PriceBachelier {
        PriceBachelier {
            forward: self.forward,
            strike: self.strike,
            volatility: self.volatility,
            expiry: self.expiry,
            is_call: self.is_call,
        }
    }

    #[must_use]
    #[inline(always)]
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
