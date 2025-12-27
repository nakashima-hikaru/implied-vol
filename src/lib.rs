//! Safe, ergonomic builders for Black–Scholes / Bachelier pricing and implied volatility.
//!
//! This crate exposes builders for computing:
//! - the implied **Black** volatility,
//! - the implied **Normal** (Bachelier) volatility,
//! - undiscounted European option prices under the **Black–Scholes** model,
//! - undiscounted European option prices under the **Bachelier** (Normal) model.
//!
//! Additionally, the crate provides a `SpecialFn` trait for implementing special functions such as `erf`
//! and `normal_cdf` used by the underlying algorithms.
//!
//! # Getting started
//!
//! Add the crate to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! implied-vol = "2.0.0"
//! ```
//!
//! To enable aggressive fused-multiply-add optimizations (when available), enable the
//! `fma` feature in `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! implied-vol = { version = "2.0.0", features = ["fma"] }
//! ```
//!
//! # Models and notation
//!
//! For clarity, we restrict the documentation examples to *undiscounted* European call options.
//! The same APIs apply to put options and to the Bachelier model (normal model).
//!
//! Let `BS(F, K, T, σ)` denote the undiscounted European call price under the Black–Scholes model,
//! where:
//! - `F`: forward price (underlying), `0 < F < ∞`,
//! - `K`: strike price, `0 < K < ∞`,
//! - `T`: time to expiry, `0 ≤ T < ∞`,
//! - `σ`: volatility, `0 ≤ σ < ∞`.
//!
//! Given `F`, `K`, `T` and an observed option price `P`, the **implied Black volatility**
//! is the value `σ` such that `BS(F, K, T, σ) = P`. In other words, it is the inverse
//! of `BS` with respect to `σ`.
//!
//! # Usage and error handling
//!
//! The crate provides two ergonomic entry points:
//! - **Free functions** (`implied_black_volatility`, `implied_normal_volatility`, etc.) that
//!   validate inputs and return the requested output directly.
//! - **Builders** (enabled by the default `builders` feature) that allow you to validate once
//!   and call `.calculate()` repeatedly.
//!
//! Validation behavior for builders:
//! - `build()` performs input validation and returns `Option<...>`; it yields `None` when
//!   the provided parameters are outside the mathematical domain of the target function.
//! - `build_unchecked()` constructs the object without validation (for callers who prefer
//!   to do their own checks or avoid the runtime cost).
//!
//! All heavy numerical work is done in `calculate()` (or inside the free functions). For
//! implied-volatility entry points, return type is `Option<f64>`; for price entry points,
//! return type is `f64`.
//!
//! ## Special functions
//!
//! Some algorithms require special mathematical functions.
//! Those are abstracted behind the `SpecialFn` trait. The crate provides
//! a default implementation named `DefaultSpecialFn` (based on the original author's code).
//! If you need to swap in a different implementation (for testing or higher-precision math),
//! implement the `SpecialFn` trait and call `calculate::<YourSpecialFn>()`.
//!
//! ## Feature flags
//!
//! - `builders` *(default)*: enables the builder structs and methods. Disable it if you want
//!   a minimal API surface and prefer the free functions.
//! - `fma`: enables aggressive fused-multiply-add optimizations where supported.
//!
//! ## Free-function example
//!
//! ```rust
//! use implied_vol::{DefaultSpecialFn, implied_black_volatility};
//!
//! let sigma = implied_black_volatility::<DefaultSpecialFn>(
//!     10.0,   // option_price
//!     100.0,  // forward
//!     100.0,  // strike
//!     1.0,    // expiry
//!     true,   // is_call
//! ).unwrap();
//! assert!(sigma.is_finite());
//! ```
//!
//! ## `PriceBlackScholes` (builder example, requires `builders` feature)
//!
//! `PriceBlackScholes::builder()` validates inputs (via `build()`), returning `None` when
//! the inputs are not in the domain of `BS`. Use `build_unchecked()` to skip validation.
//!
//! ```rust
//! use implied_vol::{DefaultSpecialFn, PriceBlackScholes};
//!
//! // Valid inputs -> builder.build() returns Some(...)
//! let builder = PriceBlackScholes::builder()
//!     .forward(100.0)
//!     .strike(100.0)
//!     .volatility(0.2)
//!     .expiry(1.0)
//!     .is_call(true);
//!
//! let price_builder = builder.build();
//! assert!(price_builder.is_some());
//! let price = price_builder.unwrap().calculate::<DefaultSpecialFn>();
//! assert!(price.is_finite());
//!
//! // Skip validation:
//! let price_builder = PriceBlackScholes::builder()
//!     .forward(100.0)
//!     .strike(100.0)
//!     .volatility(0.2)
//!     .expiry(1.0)
//!     .is_call(true)
//!     .build_unchecked();
//! let price = price_builder.calculate::<DefaultSpecialFn>();
//! assert!(price.is_finite());
//!
//! // Invalid inputs -> build() returns None
//! let invalid = PriceBlackScholes::builder()
//!     .forward(f64::INFINITY) // invalid forward
//!     .strike(100.0)
//!     .volatility(0.2)
//!     .expiry(1.0)
//!     .is_call(true)
//!     .build();
//! assert!(invalid.is_none());
//! ```
//!
//! ## `ImpliedBlackVolatility` (example)
//!
//! `ImpliedBlackVolatility::builder()` validates that `BS(F, K, T, ·)` is well-defined for
//! the supplied `F`, `K`, `T` and that the provided `option_price` is a finite, non-negative number.
//! - `build()` returns `None` if the inputs fall outside the model domain.
//! - After a successful `build()`, calling `calculate()` returns `Some(σ)` when the given
//!   `option_price` lies in the image of `BS(F, K, T, ·)`. If the price is outside that image,
//!   `calculate()` returns `None`.
//!
//! ```rust
//! use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};
//!
//! // Valid inputs -> build() returns Some(...), calculate() may return Some(σ).
//! let builder = ImpliedBlackVolatility::builder()
//!     .option_price(10.0)
//!     .forward(100.0)
//!     .strike(100.0)
//!     .expiry(1.0)
//!     .is_call(true);
//!
//! let iv_builder = builder.build();
//! assert!(iv_builder.is_some());
//! let sigma_opt = iv_builder.unwrap().calculate::<DefaultSpecialFn>();
//! assert!(sigma_opt.is_some()); // implied vol found
//!
//! // Skip validation:
//! let sigma = ImpliedBlackVolatility::builder()
//!     .option_price(10.0)
//!     .forward(100.0)
//!     .strike(100.0)
//!     .expiry(1.0)
//!     .is_call(true)
//!     .build_unchecked()
//!     .calculate::<DefaultSpecialFn>();
//! assert!(sigma.is_some());
//!
//! // If model parameters are invalid -> build() returns None
//! let invalid_builder = ImpliedBlackVolatility::builder()
//!     .option_price(10.0)
//!     .forward(f64::INFINITY) // invalid forward
//!     .strike(100.0)
//!     .expiry(1.0)
//!     .is_call(true)
//!     .build();
//! assert!(invalid_builder.is_none());
//!
//! // If the option price is outside the attainable range, calculate() returns None.
//! let out_of_range = ImpliedBlackVolatility::builder()
//!     .option_price(110.0) // too large for F=100,K=100
//!     .forward(100.0)
//!     .strike(100.0)
//!     .expiry(1.0)
//!     .is_call(true)
//!     .build()
//!     .unwrap()
//!     .calculate::<DefaultSpecialFn>();
//! assert!(out_of_range.is_none());
//! ```
#[cfg(feature = "builders")]
mod builder;
#[cfg(feature = "cxx_bench")]
pub mod cxx;
mod fused_multiply_add;
mod lets_be_rational;

pub use crate::lets_be_rational::special_function::DefaultSpecialFn;
pub use crate::lets_be_rational::special_function::SpecialFn;
#[cfg(feature = "builders")]
pub use builder::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn free_function_black_matches_normalised() {
        let forward = 120.0;
        let strike = 100.0;
        let expiry = 1.25;
        let option_price = lets_be_rational::bs_option_price::black_input_unchecked::<
            DefaultSpecialFn,
            true,
        >(forward, strike, 0.35, expiry);
        let sigma = implied_black_volatility::<DefaultSpecialFn>(
            option_price,
            forward,
            strike,
            expiry,
            true,
        )
        .unwrap();

        let intrinsic = (forward - strike).max(0.0);
        let normalised_time_value = if intrinsic > 0.0 {
            option_price - intrinsic
        } else {
            option_price
        } / (forward.sqrt() * strike.sqrt());
        let log_moneyness = forward.ln() - strike.ln();
        let sigma_normalised = implied_black_volatility_normalised::<DefaultSpecialFn>(
            normalised_time_value,
            log_moneyness,
            expiry,
        )
        .unwrap();

        assert!((sigma - sigma_normalised).abs() < 1e-12);
    }

    #[test]
    fn free_function_normal_matches_normalised() {
        let forward = 101.0;
        let strike = 100.0;
        let expiry = 1.0;
        let sigma_true = 0.25;
        let price = lets_be_rational::bachelier_impl::bachelier_price::<true>(
            forward, strike, sigma_true, expiry,
        );

        let sigma =
            implied_normal_volatility::<DefaultSpecialFn>(price, forward, strike, expiry, true)
                .unwrap();
        let intrinsic = (forward - strike).max(0.0);
        let absolute_moneyness = (forward - strike).abs();
        let sigma_normalised = implied_normal_volatility_normalised::<DefaultSpecialFn>(
            price,
            intrinsic,
            absolute_moneyness,
            expiry,
        )
        .unwrap();

        let error = (sigma - sigma_true).abs();
        let error_normalised = (sigma - sigma_normalised).abs();
        assert!(error < 1e-4, "error: {error}");
        assert!(
            error_normalised < 1e-4,
            "normalised error: {error_normalised}"
        );
    }
}

/// Compute the implied **Black** volatility from raw (non-normalised) inputs.
///
/// This function performs the same validation as the `ImpliedBlackVolatility` builder:
/// - `forward` and `strike` must be positive and finite.
/// - `expiry` must be non-negative.
/// - `option_price` must be finite and non-negative.
///
/// If the inputs are in-domain and the provided price lies inside the image of the
/// Black–Scholes pricing function, returns `Some(σ)`. Otherwise returns `None`.
#[must_use]
#[inline(always)]
pub fn implied_black_volatility<SpFn: SpecialFn>(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> Option<f64> {
    if !(forward > 0.0) || forward.is_infinite() {
        return None;
    }
    if !(strike > 0.0) || strike.is_infinite() {
        return None;
    }
    if !(expiry >= 0.0) {
        return None;
    }
    if !(option_price >= 0.0) || option_price.is_infinite() {
        return None;
    }

    if is_call {
        lets_be_rational::implied_black_volatility_input_unchecked::<SpFn, true>(
            option_price,
            forward,
            strike,
            expiry,
        )
    } else {
        lets_be_rational::implied_black_volatility_input_unchecked::<SpFn, false>(
            option_price,
            forward,
            strike,
            expiry,
        )
    }
}

/// Compute implied **Black** volatility from *normalised* inputs.
///
/// Supply the already normalised time value and log-moneyness to bypass per-call
/// normalisation overhead. Inputs are validated as follows:
/// - `normalised_time_value` must be finite and non-negative.
/// - `log_moneyness` must be finite.
/// - `expiry` must be non-negative.
///
/// Returns `Some(σ)` when the inputs are valid and within the attainable range,
/// otherwise returns `None`.
#[must_use]
#[inline(always)]
pub fn implied_black_volatility_normalised<SpFn: SpecialFn>(
    normalised_time_value: f64,
    log_moneyness: f64,
    expiry: f64,
) -> Option<f64> {
    lets_be_rational::implied_black_volatility_normalised_input::<SpFn>(
        normalised_time_value,
        log_moneyness,
        expiry,
    )
}

/// Compute the implied **Normal (Bachelier)** volatility from raw inputs.
///
/// Validation mirrors the `ImpliedNormalVolatility` builder:
/// - `forward` and `strike` must be finite.
/// - `expiry` must be non-negative.
/// - `option_price` must be finite and non-negative.
///
/// Returns `Some(σ)` when inputs are in-domain and the price is attainable,
/// otherwise returns `None`.
#[must_use]
#[inline(always)]
pub fn implied_normal_volatility<SpFn: SpecialFn>(
    option_price: f64,
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
) -> Option<f64> {
    if !forward.is_finite() || !strike.is_finite() {
        return None;
    }
    if !(expiry >= 0.0) {
        return None;
    }
    if !(option_price >= 0.0) || option_price.is_infinite() {
        return None;
    }

    if is_call {
        lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<SpFn, true>(
            option_price,
            forward,
            strike,
            expiry,
        )
    } else {
        lets_be_rational::bachelier_impl::implied_normal_volatility_input_unchecked::<SpFn, false>(
            option_price,
            forward,
            strike,
            expiry,
        )
    }
}

/// Compute the implied **Normal (Bachelier)** volatility from already normalised
/// inputs: the intrinsic price, absolute moneyness, and total option price.
///
/// Use this when you have pre-computed intrinsic value and absolute moneyness and
/// wish to avoid repeating that work per call. Validation rules:
/// - `absolute_moneyness` must be finite and non-negative.
/// - `intrinsic_price` and `option_price` must be finite, with `option_price >= intrinsic_price`.
/// - `expiry` must be non-negative.
///
/// Returns `Some(σ)` when the inputs are valid and attainable; otherwise returns `None`.
#[must_use]
#[inline(always)]
pub fn implied_normal_volatility_normalised<SpFn: SpecialFn>(
    option_price: f64,
    intrinsic_price: f64,
    absolute_moneyness: f64,
    expiry: f64,
) -> Option<f64> {
    lets_be_rational::bachelier_impl::implied_normal_volatility_normalised_input::<SpFn>(
        option_price,
        intrinsic_price,
        absolute_moneyness,
        expiry,
    )
}
