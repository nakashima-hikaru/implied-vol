# Performance

## Local comparison

Median calculation times measured on 2026-10-06 with Apple M1, macOS 27.0.1,
Rust 1.98.1, native CPU optimization, and LTO. All four solvers use the same
executable with `flashiv,experimental` enabled and the optional `fma` feature
disabled. Prices and builders are prepared before timing.

Each input has 24 samples per solver, with 200,000 calls per round and balanced
solver order. The mixed workloads contain 4,096 seeded synthetic OTM calls.
These warm single-threaded measurements are specific to this workload and host;
cold caches, parallel batches, tail latency, and native Windows are not covered.

| Input | Hybrid | Jaeckel | FlashIv | Experimental |
|---|---:|---:|---:|---:|
| ATM | 7.6 ns | 7.6 ns | 162.1 ns | 10.4 ns |
| Lowest region | 181.7 ns | 262.3 ns | 146.4 ns | 105.4 ns |
| Lower middle | 111.2 ns | 109.2 ns | 158.5 ns | 107.5 ns |
| Upper middle | 112.3 ns | 110.5 ns | 160.3 ns | 130.5 ns |
| Highest region | 187.5 ns | 183.5 ns | 175.0 ns | 90.8 ns |
| Near ATM | 104.8 ns | 102.4 ns | 165.3 ns | 74.9 ns |
| Near ATM, wider | 95.3 ns | 93.7 ns | 164.2 ns | 111.5 ns |
| Mixed normalized | 176.9 ns | 198.1 ns | 178.0 ns | 123.9 ns |
| Mixed full API | 197.1 ns | 219.2 ns | 200.2 ns | 140.6 ns |
| Deep OTM, full API | 127.9 ns | 127.7 ns | 186.4 ns | 126.7 ns |
| Near-ATM short expiry, full API | 169.9 ns | 171.1 ns | 190.6 ns | 92.7 ns |
| OTM long expiry, full API | 133.9 ns | 135.7 ns | 191.1 ns | 150.7 ns |

Experimental uses 30.0%/28.7% less time than Hybrid on the mixed normalized/full
workloads. Upper-middle, wider near-ATM, and OTM long-expiry full-API inputs favor
Hybrid or Jaeckel. Choose a solver using both speed and accuracy for your inputs.

### Accuracy

The accuracy comparison uses the same numerical source and FMA-off policy on
81,873 shared exact-binary64 `(a,b)` reference inputs. This separate population
includes dense small-domain and seed-boundary inputs, extreme moneyness, and
cap neighbors; it is not the timing workload. Subnormal-price and noninterior inputs are excluded.

`rho_J` is relative volatility error divided by Jäckel's attainable-precision
estimate, evaluated against independent high-precision roots. Smaller values
are better; `rho_J<1` meets that comparison target. Maxima and exceedance counts
use finite positive outputs, with non-finite outputs reported separately.

| Solver | Maximum `rho_J` | Finite outputs with `rho_J>=1` | Non-finite outputs |
|---|---:|---:|---:|
| Hybrid | 8.019177 | 3,416 | 3,075 |
| Jaeckel | 11.430756 | 3,599 | 3,075 |
| FlashIv | 1.572e+08 | 17,763 | 3,075 |
| Experimental | 0.495459 | 0 | 0 |

The 3,075 non-finite outputs are infinities at prices equal to a rounded-down
cap, although the mathematical cap is above the input price and its exact root
is finite. This boundary handling is separate from finite-output error.
FlashIv has its paper method's accuracy trade-off; `rho_J<1` is a common
comparison yardstick, not its promised precision contract. Finite checks do
not establish a whole-domain guarantee. See
[accuracy contracts and limits](numerics.md#accuracy-contracts).

[Current measurement data](results.json) provides the table values, measurement
settings, source fingerprints, and worst-case accuracy inputs.

## Run the benchmarks

```sh
RUSTFLAGS="-C target-cpu=native" \
  cargo bench --bench implied_black_vol_paths --features flashiv,experimental -- \
  200000 3 --solver hybrid
```

Repeat with `--solver jaeckel`, `--solver flashiv`, and `--solver experimental`.
The features make optional solver types available; `--solver` selects the
measured solver. Add `fma` consistently to compare that arithmetic policy.
Use the same features, flags, inputs, and call counts, rotate solver order,
and finish builds before timing.

The benchmark reports ns/call and a checksum. It includes seven fixed
normalized paths, two mixed workloads, and three full-API inputs. Optional
`--experimental-paths`, `--as1-paths`, and `--dpoly-paths` add focused
Experimental workloads; their inputs are also prepared outside timing.

The reported Rust measurements permit FMA only at explicit `mul_add` sites.
The `RUSTFLAGS` override above disables the repository's configured implicit
contraction. Adding `-C llvm-args=-fp-contract=fast` changes the arithmetic
policy and should be reported separately. Experimental uses explicit FMA
independently of the optional crate feature.

## Compare with C++

The optional `cxx_bench` feature builds the bundled C++ Let's Be Rational
implementation and requires a C++ compiler. To compare Hybrid, Jaeckel, and C++
on the same benchmark inputs:

```sh
RUSTFLAGS="-C target-cpu=native" \
  cargo bench --bench implied_black_vol_paths --features cxx_bench --no-run
python3 benches/compare_solvers.py <benchmark-executable> \
  --features cxx_bench --flags="-C target-cpu=native" \
  --output solvers.json --summary solvers.md
```

Use the executable path printed by Cargo. Add `fma` consistently to the build
and `--features` argument for the explicit-FMA comparison. The report rotates
solver order, checks within-solver checksum stability, and reports medians.
Cross-solver checksums may differ; validity preflight is not an accuracy oracle.

The C++ build uses `-Ofast` and `-ffp-contract=fast`, while the Rust command
above uses explicit FMA only. Interpret timing differences alongside these
arithmetic policies and the solvers' numerical limits.
