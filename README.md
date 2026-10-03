# Implied Vol

[![Crates.io](https://img.shields.io/crates/v/implied-vol)](https://crates.io/crates/implied-vol)
[![Actions status](https://github.com/nakashima-hikaru/implied-vol/actions/workflows/ci.yaml/badge.svg)](https://github.com/nakashima-hikaru/implied-vol/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Black and Bachelier implied volatility and European option prices in pure Rust.
The default build has no required dependencies. Builders validate inputs, and
calculations support custom special functions through the `SpecialFn` trait.

## Getting started

```toml
[dependencies]
implied-vol = "2.1"
```

```rust
use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};

let option = ImpliedBlackVolatility::builder()
    .option_price(10.0)
    .forward(100.0)
    .strike(100.0)
    .expiry(1.0)
    .is_call(true)
    .build()
    .unwrap();

let volatility = option.calculate::<DefaultSpecialFn>().unwrap();
assert!(volatility.is_finite());
```

Prices are undiscounted; apply discounting outside the library. The Black
implied-volatility builder returns annualized volatility. The normalized builder
takes log-moneyness and normalized price, and returns total volatility `sigma * sqrt(T)`.

The crate also provides builders for Bachelier implied volatility, Black and
Bachelier option prices, and normalized Black prices. See the
[API documentation](https://docs.rs/implied-vol/) for their inputs and return values.

## Choose a Black solver

`calculate::<SpFn>()` always uses `solver::Hybrid`. Select a solver explicitly with
`calculate_with::<S>()` on either Black implied-volatility builder:

```rust
use implied_vol::{ImpliedBlackVolatilityNormalised, solver::Jaeckel};

let option = ImpliedBlackVolatilityNormalised::builder()
    .log_moneyness(-0.1)
    .normalised_price(0.05)
    .build()
    .unwrap();
let volatility = option.calculate_with::<Jaeckel>().unwrap();
```

| Solver type | Availability | Method |
|---|---|---|
| `solver::Hybrid` | Always | Let's Be Rational + FlashIV for selected low-price inputs |
| `solver::Jaeckel` | Always | Let's Be Rational |
| `solver::FlashIv` | `flashiv` feature | FlashIV paper's fixed-count Algorithm 1 |
| `solver::Experimental` | `experimental` feature | Author's own implementation |

Solver features add types and leave the default solver selection unchanged. Both
optional solvers can be enabled together:

```toml
implied-vol = { version = "2.1", features = ["flashiv", "experimental"] }
```

Then use `calculate_with::<solver::FlashIv>()` or
`calculate_with::<solver::Experimental>()`. Hybrid, Jaeckel, and FlashIv default to
`DefaultSpecialFn`; use `solver::Jaeckel<MySpecialFn>` to choose your own provider.
Experimental preserves its own special functions and explicit FMA arithmetic.

`calculate_explicit::<SpFn>()` remains a separate inverse-Gaussian Black formula.
Black solver selection does not change pricing or Bachelier inversion.

## Features

| Feature | Effect |
|---|---|
| `flashiv` | Makes `solver::FlashIv` available |
| `experimental` | Makes `solver::Experimental` available |
| `fma` | Enables optional FMA arithmetic in shared numerical kernels |
| `cxx_bench` | Builds the bundled C++ comparison implementation; requires a C++ compiler |

## Performance and precision

Median time per calculation on Apple M1, using native CPU optimization and LTO,
the optional `fma` feature disabled, and 15 samples per solver. The mixed
workloads contain 4,096 seeded synthetic OTM calls; prices and builders are
prepared before timing.

| Solver | Mixed normalized | Mixed full API | Accuracy |
|---|---:|---:|---|
| Hybrid (default) | 181.4 ns | 202.6 ns | Targets Jäckel's maximum attainable precision; some edge cases fall short |
| Jaeckel | 198.4 ns | 218.9 ns | Targets Jäckel's maximum attainable precision; some edge cases fall short |
| FlashIv | 178.4 ns | 200.3 ns | Paper method; can lose accuracy near ATM |
| Experimental | 129.4 ns | 145.1 ns | Met Jäckel's precision target on all 114,961 reference cases |

Experimental took 28.7%/28.4% less time than Hybrid on the mixed normalized/full
workloads, but was slower on several middle and near-ATM paths. These are local
measurements; no solver is fastest in every region.

Jäckel's maximum attainable precision accounts for how option-price rounding
affects implied volatility: more sensitive inputs allow a larger error. The
accuracy column describes current results, with different test coverage for each
solver. Experimental passed all available reference cases; accuracy over its
entire input domain remains unproven. Detailed error definitions and test results
are in the numerical notes below.

FlashIv follows the paper's fixed iteration count and has a different accuracy
trade-off from Hybrid. Its direct erfcx subtraction can lose digits near ATM;
returning a finite result does not certify convergence.

- [Performance results and reproducible benchmarks](docs/performance.md)
- [Algorithms, accuracy contracts, and numerical limits](docs/numerics.md)

## Algorithm references

| Implementation / API | Algorithm and reference |
|---|---|
| `solver::Hybrid` | [Let's Be Rational](https://www.jaeckel.org/LetsBeRational.pdf), with [FlashIV](https://arxiv.org/abs/2605.29102v1) for selected low-price inputs |
| `solver::Jaeckel` | Peter Jäckel's [Let's Be Rational](https://www.jaeckel.org/LetsBeRational.pdf) |
| `solver::FlashIv` | [FlashIV, Algorithm 1](https://arxiv.org/abs/2605.29102v1) |
| `solver::Experimental` | Author's own implementation |
| `ImpliedNormalVolatility` | Peter Jäckel's [Implied Normal Volatility](https://www.jaeckel.org/ImpliedNormalVolatility.pdf), for Bachelier inversion |
| `calculate_explicit()` on the Black implied-volatility builders | [An Explicit Solution to Black-Scholes Implied Volatility](https://arxiv.org/abs/2604.24480), using the inverse-Gaussian representation |

The FlashIV paths in Hybrid and FlashIv use initial estimates from
[Li and Lee, equation 42 and seed coefficients](https://mpra.ub.uni-muenchen.de/6867/1/MPRA_paper_6867.pdf).
