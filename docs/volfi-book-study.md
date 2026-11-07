# Volfi whole-book direct charts in Rust: scalar baseline and numerical limits

**Research date:** 2026-11-07. **Status:** implemented as opt-in, **not recommended
for production or as an Experimental replacement**.

## Method and provenance

- Original: Wolfgang Schadner, [volfi](https://github.com/wol-fi/volfi),
  `volfi` revision `7347f68025fea2463b30f67c714a7448a64d1c39`,
  `include/volfi/volfi_wb.hpp` and generated
  `include/volfi/volfi_wb_tables.hpp`.
- Original's method/paper: *Implied Volatility in One Straight Line:
  Machine Precision at Vector-Hardware Throughput*, see the original
  repository's `docs/volfi_v0.3.1_paper.pdf`.
- The coefficient literals of `NW_G0` (23), `NW_F1` (40), `NW_P`
  (137) and `NW_D` (171) are transcribed without numeric re-fitting.
  `src/volfi_book/tables.rs` holds the generated data.
- Original copyright (c) 2026 Wolfgang Schadner, BSD-3-Clause;
  license copied to `src/volfi_book/LICENSE-BSD-3-CLAUSE`.

The Rust port includes **only** the whole-book A/B direct charts, not
the original volfi's separately routed four-chart inverse, vector intrinsics,
or CUDA implementation. The `solver::VolfiBook` adapter falls back to
the library's pre-existing `solver::Experimental` outside the A/B chart.
This makes its end-to-end routing benchmark a comparison of **the chart
with Experimental fallback**, **not the official complete C++ volfi**.

Rust standard-library `exp_m1` and `ln_1p` replace the original's
custom ~1-ULP vectorizable kernels; arithmetic bits, portability, and
the author's reported native SIMD throughput are **not transferable**.

## Coordinates, exact transformations and acceptance

Let `h=abs(log(F/K))` and `b=OTM price/sqrt(F*K)` be the library's
normalized inputs. The exact Black price cap is `exp(-h/2)`.
Define `c = b exp(h/2)` and intrinsic coordinate

```text
a = log1p(expm1(h)/c)
  = log1p(2 sinh(h/2)/b).
```

The second equality is **mathematical**, not bitwise: the corresponding
`VolfiBookSinh` implementation is provided as a separate research
variant and tested on the same roots and benchmark inputs.

Direct **region A** requires `a < 2*pi` and `h <= 0.35*a`.
It evaluates a one-variable Chebyshev approximation to `G(a)` and
an analytic polynomial in the scaled square `(h/(a G(a)))^2`.

Direct **region B** requires `2*pi <= a <= 700` and
`h <= min(4, 0.16*a+0.85)`. Its expansion uses

```text
z = h^2 / (4*pi^2)
q = z / (sqrt(1+z)+1)^2,
```

the conformal coordinate of the original whole-book chart. Both charts
evaluate 15/16 rows of supplied polynomial coefficients; input/outside
conditions call Experimental. These conditions jointly imply `h<=4`,
so the optimized variant can cheaply reject larger h before coordinate work.

## Independent accuracy

Reference source: `tests/experimental_adaptive_seeds.rs` (975 exact
binary64 pairs with independent Black roots and `eta_J` values verified
at 160 and 240 decimal digits in the preceding project).

```text
eta_J = 2^-52 * (1 + b/(s_star * V_star))
rho_J = abs(s_hat/s_star - 1)/eta_J
V_star = exp(-h^2/(2*s_star^2)-s_star^2/8)/sqrt(2*pi).
```

Results on GitHub Actions, both `volfi_book` and `volfi_book,fma`:
- **495 of 975** reference inputs entered the direct A/B charts.
- Maximum direct-chart `rho_J=2.17876517`, attained near
  `h=0.999999999999999889`, `b=6.14421235332821065e-6`.
- Hyperbolic-sine variant has the **same** observed maximum `rho_J`;
  its maximum relative output difference from the original-coordinate
  port is `4.46724986e-16`.
- Tests assert `rho_J<64` only as a *smoke guard*. A test pass
  **does not assert the strict `rho_J<1` precision target**.
- Existing Experimental and default solver paths were not changed.

Therefore the sampled direct-chart implementation **fails** the
existing Experimental target on at least one accepted input. This is
a finite counterexample to asserting `rho_J<1` for the port; it is not
a proof that the original C++ author's different numerical implementation
has the same error at that input.

## Same-host scalar benchmark (original-coordinate port)

**Runner:** Linux x86-64, Azure AMD EPYC; `rustc 1.99.1` /
LLVM 23.1.1, `RUSTFLAGS="-C target-cpu=native -C target-feature=+fma"`,
bench profile with LTO. The `volfi_book` feature was enabled, but not
the optional crate `fma` feature in this table.

**Measurement:** one executable, shared prepared `(h,b)` inputs, warmed
up, rotating ABBA/BAAB 4 blocks × 2 rounds, 150,000 calls per
invocation; 16 samples per solver and case. Positive changes mean
the port is slower.

| Case | Experimental ns | VolfiBook ns | Port change | Direct accepted |
|---|---:|---:|---:|---:|
| book_a | 84.10 | 166.07 | +97.48% | 1/1 |
| book_b | 189.98 | 229.57 | +20.84% | 1/1 |
| book_seam | 190.89 | 229.88 | +20.43% | 1/1 |
| fallback | 187.84 | 191.10 | +1.74% | 0/1 |
| mixed_book | 189.34 | 221.22 | +16.84% | 4096/4096 |
| mixed_market | 196.71 | 252.95 | +28.59% | 1041/4096 |
| mixed_full | 218.68 | 269.40 | +23.19% | 1041/4096 |

The separately benchmarked `sinh` coordinate is 1.48–2.32%
slower on three directly routed single-quote cases and 0.20%
slower on `mixed_book`, but 2.85% faster on the mixed market
population (which includes fallback). Differences in routing/code layout
and CPU environment may contribute. **No evidence supports using it as a
faster direct solver.**

FMA-on and FMA-off are separate CI matrix hosts; do not compare absolute
nanoseconds across those jobs. Original run logs, CPU details,
round-level samples, and hashed executables are attached to:

- [FMA-off original-coordinate and sinh paired samples](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42374038581/artifacts/13574118862)
- [FMA-on original-coordinate and sinh paired samples](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42374038581/artifacts/13574398436)

Earlier no-sinh baseline from the initial integration
[also showed slower scalar performance](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42373276103/artifacts/13574616912).
The populations and timing arrangements are comparable within
each experiment, not across different CI hosts.

All cases here are synthetic. The upstream published SPX market quote
data are licensed and **not distributed**. The `mixed_book` collection
is constructed to lie in the chart, while `mixed_market` samples
option prices from a deliberately broader region. Direct coverage of the
latter is a workload property and **not** a universal market coverage rate.
This is a scalar **latency** test, not an AVX2/AVX-512 batch throughput or
end-to-end official volfi port benchmark.

## Theory-first follow-up research

1. **Weighted approximation / attainable precision.** Construct a direct
   approximation targeting the residual weight
   `|B(s)-b|/(s*vega+b)`, corresponding locally to `rho_J`.
   Replace the generic uniform chart polynomial with conditional,
   cell-specific degree budgets; establish the inverse derivative bound
   across each cell before deleting coefficients.
2. **Source-rounding certificates.** Bound the upstream `G` and
   conformal-polynomial representation errors, coordinate transformation
   errors (`exp_m1`, `ln_1p`, `sinh`) and finite precision source order.
   Include exact price-cap and seam behavior. A sampled error below
   one would not itself be an analytic certificate.
3. **Selective direct routing.** Introduce a cheap, mathematically
   certified gate accepting the direct result only if a derived
   upper bound on `rho_J` is below one. Otherwise use Experimental.
   Evaluating a full price residual to decide acceptance may erase
   all speed benefits: this gate must be cheaper than the saved finish.
4. **Cross-quote SIMD.** The imported source has scalar independent
   two-lane approximation machinery in the older Experimental
   implementation. volfi's reported throughput advantage relies on
   batched same-coordinate vector execution. Separate the batch
   throughput comparison from the single-quote latency comparison,
   including mask/blend penalties when inputs straddle A/B regions.
5. **Original C++ reference comparison.** To compare against the
   *complete* published engine, build the original `volfi` scalar
   and AVX2/AVX-512 benchmark on the same CPU and same normalized
   workload; it must be labelled separately from this port.

## Reproduce

```sh
cargo test --locked --features volfi_book --test experimental_adaptive_seeds volfi_book_ -- --nocapture
cargo test --locked --features volfi_book,fma
cargo clippy --locked --features volfi_book -- -D warnings
cargo fmt --all --check
python3 -m unittest discover -s benches -p test_compare_volfi.py
cargo bench --bench volfi_compare --features volfi_book --no-run
# Use the single executable path emitted by Cargo's --message-format=json.
python3 benches/compare_volfi.py /path/to/volfi_compare_binary \
  --iterations 150000 --rounds 2 --blocks 4 \
  --reference experimental --candidate volfi-book \
  --features volfi_book --flags '-C target-cpu=native -C target-feature=+fma' \
  --output /tmp/volfi-book-comparison.json
```

The comparison reports its compiler, execution-platform and executable
SHA-256 provenance. No timings are portable guarantees. The optional
research feature neither changes the default `Hybrid` nor the
existing `Experimental` implementation.
