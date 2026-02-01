# Implied Vol

[![Crates.io](https://img.shields.io/crates/v/implied-vol)](https://crates.io/crates/implied-vol)
[![Actions status](https://github.com/nakashima-hikaru/implied-vol/actions/workflows/ci.yaml/badge.svg)](https://github.com/nakashima-hikaru/implied-vol/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`implied-vol` is a high-performance, pure Rust library for calculating implied volatility,
implemented based on the methods described in Peter Jäckel's seminal papers.

## Usage

This crate exposes builders for computing:

- the implied **Black** volatility,
- the implied **Normal** (Bachelier) volatility,
- (undiscounted) European option prices under the **Black–Scholes** model,
- (undiscounted) European option prices under the **Bachelier** model.

Additionally, the crate provides a trait for implementing custom special functions, which can be used to customize
the calculation of implied volatilities and option prices.

Add the following to your `Cargo.toml`:

```toml
[dependencies]
implied-vol = "2.0"
```

The calculations are performed via builders that allow you to handle errors.

### Example

```rust
use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};
let iv_builder = ImpliedBlackVolatility::builder()
.option_price(10.0)
.forward(100.0)
.strike(100.0)
.expiry(1.0)
.is_call(true)
.build().unwrap();

let iv = iv_builder.calculate::<DefaultSpecialFn>().unwrap();
assert!(iv.is_finite());
```

More details can be found in the [crate documentation](https://docs.rs/implied-vol/2.0/implied_vol/).

## Source References

This crate implements algorithms from two key papers by Peter Jäckel:

1. [Let's Be Rational](http://www.jaeckel.org/LetsBeRational.pdf) — A method for accurately and efficiently
   extracting Black implied volatility from option prices.

2. [Implied Normal Volatility](http://www.jaeckel.org/ImpliedNormalVolatility.pdf) — An analytical formula for
   computing implied normal (Bachelier) volatility from vanilla option prices.

Both papers and related materials are available on [Peter Jäckel's website](http://www.jaeckel.org/).

## Key differences from the original Authors' Implementation

- Enables the use of custom special functions for calculating implied volatilities and option prices.
- Modified to return results earlier for ATM cases.
- Utilizes explicit FMA instructions when available.
    - Results are undefined when inputs contain NaNs. (This is because this crate provides builders that handle errors.)
    - Some implementations of asymptotic expansions differ from the original authors' implementation for FMA
      utilization.

## Performance

Benchmark results, available via our [GitHub Actions](https://github.com/nakashima-hikaru/implied-vol/actions),
compare the execution speed against FFI to Jäckel’s original reference C++ implementation.
With aggressive compiler optimizations applied to both implementations, this Rust crate often outperforms the C++ FFI
version.

## Precision

Random tests have confirmed that the normalized implied volatility function under the Black-Scholes model is accurate to within 3 ULPs (Units in the Last Place) of the true implied volatility.

## Cargo Feature Flags

* `fma`: Enables Fused Multiply-Add (FMA) instructions when supported by the target CPU, providing a slight performance
  boost over the default implementation.
