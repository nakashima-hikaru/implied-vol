# Implied Vol

[![Crates.io](https://img.shields.io/crates/v/implied-vol)](https://crates.io/crates/implied-vol)
[![Actions status](https://github.com/nakashima-hikaru/implied-vol/actions/workflows/ci.yaml/badge.svg)](https://github.com/nakashima-hikaru/implied-vol/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

More information about this crate can be found in
the [crate documentation](https://docs.rs/implied-vol/1.0/implied_vol/).

## About

`implied-vol` is a high-performance, pure Rust implementation of Peter Jäckel's implied volatility calculations. This
library serves as a robust Rust reimplementation of the methodologies presented in Jäckel's works.

## Source Works

Our library follows the methods presented in two pivotal papers by Peter Jäckel:

1. [Let's Be Rational](http://www.jaeckel.org/LetsBeRational.pdf): This work presents an approach to deduce Black’s
   volatility from option prices with high precision.

2. [Implied Normal Volatility](http://www.jaeckel.org/ImpliedNormalVolatility.pdf): Here, Jäckel provides an analytical
   formula to calculate implied normal volatility (also known as Bachelier volatility) from vanilla option prices.

Both resources can be accessed at [Peter Jäckel's homepage](http://www.jaeckel.org/).

## Performance

Peter Jäckel, the author of the original paper, asserts that "the calculation of a single implied volatility is now down
to just under 270 nanoseconds" based on his machine's benchmark measurements. By examining the benchmark measurements
performed on this crate's [GitHub Actions](https://github.com/nakashima-hikaru/implied-vol/actions), it becomes clear
that comparable performance is being achieved.

## Precision

On our machine, the relative error for both implied Black volatility and implied normal
volatility calculations is confirmed to be less than twice the machine epsilon in random tests.

Community contributions are always welcome!

## Cargo Feature Flags

- `normal-distribution`: Provide functions related to standard normal distribution used in calculation of implied
  volatility
- `error-function`: Provide functions related to error function used in calculation of implied volatility

## License

This project is licensed under the [MIT license](https://github.com/nakashima-hikaru/implied-vol/blob/main/LICENSE).