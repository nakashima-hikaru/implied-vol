#![allow(dead_code)]

use implied_vol::{DefaultSpecialFn, PriceBachelier, PriceBlackScholes};

#[derive(Clone, Copy, Debug)]
pub struct PriceCase {
    pub forward: f64,
    pub strike: f64,
    pub volatility: f64,
    pub expiry: f64,
    pub is_call: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct ImpliedCase {
    pub option_price: f64,
    pub forward: f64,
    pub strike: f64,
    pub expiry: f64,
    pub is_call: bool,
}

pub const ATM_CALL: PriceCase = PriceCase {
    forward: 100.0,
    strike: 100.0,
    volatility: 0.2,
    expiry: 1.0,
    is_call: true,
};

pub const DEEP_ITM_PUT_LONG: PriceCase = PriceCase {
    forward: 0.25,
    strike: 2.0,
    volatility: 0.5,
    expiry: 10.0,
    is_call: false,
};

pub const DEEP_OTM_CALL_LONG: PriceCase = PriceCase {
    forward: 0.25,
    strike: 2.0,
    volatility: 0.5,
    expiry: 10.0,
    is_call: true,
};

pub const NEAR_ATM_CALL_SHORT: PriceCase = PriceCase {
    forward: 100.0,
    strike: 100.25,
    volatility: 0.2,
    expiry: 1.0e-3,
    is_call: true,
};

pub const EDGE_ZERO_EXPIRY_CALL: PriceCase = PriceCase {
    forward: 120.0,
    strike: 100.0,
    volatility: 0.3,
    expiry: 0.0,
    is_call: true,
};

#[must_use]
pub const fn cxx_q(is_call: bool) -> f64 {
    if is_call { 1.0 } else { -1.0 }
}

#[must_use]
pub fn black_price(case: PriceCase) -> f64 {
    PriceBlackScholes::builder()
        .forward(case.forward)
        .strike(case.strike)
        .volatility(case.volatility)
        .expiry(case.expiry)
        .is_call(case.is_call)
        .build_unchecked()
        .calculate::<DefaultSpecialFn>()
}

#[must_use]
pub fn bachelier_price(case: PriceCase) -> f64 {
    PriceBachelier::builder()
        .forward(case.forward)
        .strike(case.strike)
        .volatility(case.volatility)
        .expiry(case.expiry)
        .is_call(case.is_call)
        .build_unchecked()
        .calculate::<DefaultSpecialFn>()
}

#[must_use]
pub fn black_implied_case(case: PriceCase) -> ImpliedCase {
    ImpliedCase {
        option_price: black_price(case),
        forward: case.forward,
        strike: case.strike,
        expiry: case.expiry,
        is_call: case.is_call,
    }
}

#[must_use]
pub fn normal_implied_case(case: PriceCase) -> ImpliedCase {
    ImpliedCase {
        option_price: bachelier_price(case),
        forward: case.forward,
        strike: case.strike,
        expiry: case.expiry,
        is_call: case.is_call,
    }
}
