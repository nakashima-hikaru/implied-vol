//! Safe, ergonomic builders for Black–Scholes / Bachelier pricing and implied volatility.
//!
//! This crate exposes builders for computing:
//! - the implied **Black** volatility,
//! - the implied **Normal** (Bachelier) volatility,
//! - undiscounted European option prices under the **Black–Scholes** model,
//! - undiscounted European option prices under the **Bachelier** (Normal) model.
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
//! The crate provides builder types that validate inputs at construction time and then
//!  exposes a `calculate()` method for performing the computation. Validation behavior is:
//!
//! - `build()` performs input validation and returns `Option<...>`; it yields `None` when
//!   the provided parameters are outside the mathematical domain of the target function.
//! - `build_unchecked()` constructs the object without validation (for callers who prefer
//!   to do their own checks or avoid the runtime cost).
//!
//! All heavy numerical work is done in `calculate()`; for implied-volatility builders
//! `calculate()` returns `Option<f64>` (or `None` when the given option price is not in the
//! function's image), and for price builders `calculate()` returns `f64`.
//!
//! ## Special functions
//!
//! Some algorithms require special mathematical functions.
//! Those are abstracted behind the `SpecialFn` trait. The crate provides
//! a default implementation named `DefaultSpecialFn` (based on the original author's code).
//! If you need to swap in a different implementation (for testing or higher-precision math),
//! implement the `SpecialFn` trait and call `calculate::<YourSpecialFn>()`.
//!
//! ## `PriceBlackScholes` (example)
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
mod builder;
#[cfg(feature = "cxx_bench")]
pub mod cxx;
mod fused_multiply_add;
mod lets_be_rational;

pub use crate::lets_be_rational::special_function::DefaultSpecialFn;
pub use crate::lets_be_rational::special_function::SpecialFn;
pub use builder::*;
