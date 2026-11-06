# Adaptive physical predictor and central error-budget allocation

Implemented on 2026-11-07 against `59c5ec6`. The
[machine-readable record](adaptive-physical-improvements.json) binds the
sources, conditional certificates, independent reference corpus and local
measurements. The preceding [cost reductions](further-theory-improvements.md)
form the baseline.

## Adaptive upper predictor

The upper large-moneyness route now chooses reciprocal-derivative index
**3, 4 or 5** for its initial physical-price predictor. Let `n=M(y)*y/s` be
its physical Newton proposal and `r=abs(x*n)`, the magnitude of the leading change in the Gaussian
exponent. Index 3 is used at `r<=1/8`, index 4 at `r<=3/8`, and index 5
otherwise. The proposal is the exact Newton ratio only with an exact
inverse-normal seed and exact Mills evaluation. Higher recurrence terms are
evaluated only when needed.

For the normalized local series
`f(x+n*t)/f(x)=1-t-q2*t²-q3*t³-...`, let `c_j` denote the reciprocal-series
coefficients. The updates are `n*c2/c3`, `n*c3/c4` and `n*c4/c5`.
Exact polynomial identities verify these forms. With exact derivatives and
Newton ratio, their individual local Householder orders are 4, 5 and 6.
These orders do not describe a global convergence rate of the adaptive solver.

The thresholds are cost choices, not residual-error certificates. The source
checker binds each retained recurrence expression to the reviewed original
and binds both branches. Finite/in-bracket checks, the precise residual gate,
reliable sign updates and forced bisection remain unchanged. The conditional
large-route final-error bound remains **`rho_J<0.873762`**. There is no claim
that every input passes the first precise gate.

In the 52,617-input replay, all 21,003 upper large inputs still use exactly
one precise evaluation and no fast residual evaluation or bisection:

| Predictor index | Inputs | Largest first precise residual magnitude |
|---|---:|---:|
| 3 | 7,165 | 5.412949e-5 |
| 4 | 11,428 | 3.302644e-5 |
| 5 | 2,410 | 8.902104e-6 |

All three are below `2^-14`. The 26,067 lower large inputs retain one fast
and one precise evaluation each. Instrumentation reproduces the corresponding
uninstrumented outputs. These are finite observations; the safeguards cover
missed initial handoffs.

## Central error budget

The four normal-shape intervals now use degrees **16/18/18/24**, previously
18/20/20/25. This removes another **2/2/2/1 two-lane Horner stages**. Stored
coefficients, the 48-coefficient Black conversion, both compensated final
stages, coordinate-low correction and runtime guards are retained.

The previous deleted-tail reserve of `0.01u` was stricter than necessary.
Recomputing the whole error transport permits a `0.2u` reserve, with actual
whole-interval tail bounds below `0.177446u`, `0.091231u`, `0.105882u` and
`0.109107u`. The checker includes derivative-times-low terms, coordinate
errors, source rounding, analytic conversion error and final rounding.
Its maximum conditional bound rises from `0.570599` to **`rho_J<0.668277`**,
still below the target of one. This trades surplus approximation precision
for less arithmetic; individual outputs can become less accurate. See the
[central derivation](central-expansion.md).

## Accuracy and validation

Both FMA feature policies pass **52,617 unique exact binary64 normal-price
inputs**, with no invalid output or `rho_J>=1`. Their output bytes agree.
Maximum `rho_J` remains `0.6172411593798071`. Of 45 changed outputs, 13 improve
and 32 worsen; 26 changes are in the small route and 19 in the large route.
The maximum among changed outputs is `0.4631368812`.

The new predictor fixture contains 95 exact-price inputs near ideal switching
thresholds, neighboring represented prices and interior points. Independent
Gaussian-erfc roots and conditioning values agree at 160 and 240 decimal
digits. Ideal thresholds need not coincide with the rounded runtime switch;
instrumentation records actual branch coverage separately.

All **2,988 additional reference fixtures** pass both FMA policies: 1,026
central, 625 rank-compression, 975 adaptive-seed, 156 common-finish, 111 large
iteration and 95 new predictor inputs. The central set changes 48 outputs
and the new predictor set changes one; the other sets retain their bits.
Maximum `rho_J` is `0.4000845225` in the central set and `0.2866641568` in the
new predictor set.

One archived C++ bit anchor, `(a,b)=(binary64(.1),binary64(.05))`, changes by
one ulp. Its independently computed exact-input root and conditioning agree
at 160 and 240 digits; `rho_J` changes from `0.1628260976` to `0.1765337164`.
That input now has an independent-root precision test, and all twelve other
archived anchors remain. The source bridge binds this exact test migration.

CXX-enabled combined suites pass **124 tests per FMA policy**, including
doctests. Warning-denying Clippy passes both policies, formatting passes, and
all 23 Python benchmark-contract tests pass. Library compilation passes for
Windows MSVC x86/x64 and WebAssembly under both policies; these six checks
are cross-compilations, not native execution. The central checker passes its
positive control and rejects 20 mutations; the large checker rejects 18.
The direct-D and AS1 consumers also accept the exact source bridges, with
their ten seed defects rejected by both consumers.

## Timing

The frozen baseline and candidate use an identical benchmark harness,
native CPU optimization and bench-profile LTO. Six rotating ABBA/BAAB blocks
use 300,000 calls and two rounds per invocation. No compiler, proof or
validation jobs from this task run concurrently. The host is macOS AArch64
with rustc 1.98.1 / LLVM 22.1.8.

| Workload | Baseline ns | New ns | Observed change |
|---|---:|---:|---:|
| central_deferred | 61.87 | 60.33 | -2.50% |
| central_shape1 | 63.50 | 62.34 | -1.83% |
| central_shape2 | 63.64 | 62.48 | -1.82% |
| central_extended | 66.87 | 66.72 | -0.21% |
| mixed_near_atm | 71.01 | 69.60 | -1.98% |
| large_upper | 221.23 | 216.11 | -2.31% |
| large_upper_small_step | 221.00 | 212.26 | -3.96% |
| large_upper_full_step | 292.09 | 291.95 | -0.05% |
| large_near_cap | 289.34 | 284.13 | -1.80% |
| large_lower | 180.46 | 181.86 | +0.78% |
| large_asymptotic | 182.99 | 184.48 | +0.82% |
| mixed_large | 193.49 | 193.41 | -0.04% |
| mixed_normalised | 124.82 | 124.85 | +0.03% |
| mixed_full | 142.16 | 142.08 | -0.06% |

All six paired blocks improve for the first three central cases,
`mixed_near_atm`, `large_upper`, `large_upper_small_step` and `large_near_cap`.
The two new cases exercise predictor indices 3 and 5, respectively. The
unchanged lower large path is about 0.8% slower in all six blocks; changes
to the function's compiled layout may contribute, but the cause is not
isolated. General mixed workloads show no material gain. Unchanged Hybrid
control medians range from -0.19% to +1.40%, and individual control blocks
range from -1.34% to +1.93%. Every block is retained.

The affected-case improvements are local observations, not a portable
speedup guarantee. The JSON record preserves all workloads and block changes.
To reproduce the comparison, build the baseline revision with the final
benchmark harness copied into it, build both executables with the stated
flags and run `benches/compare_experimental.py` with `--iterations 300000
--rounds 2 --blocks 6`.

## Reproduction and limits

Use the same hash-bound `implied-black-volatility` reference checkout as the
previous reports:

```sh
REFERENCE=/path/to/implied-black-volatility
.venv/bin/python scripts/verify_central_expansion.py --reference-root "$REFERENCE"
.venv/bin/python scripts/test_central_expansion.py --reference-root "$REFERENCE"
.venv/bin/python scripts/verify_large_fast.py --reference-root "$REFERENCE" --out target/proof/large-fast-proof.json
.venv/bin/python scripts/verify_large_direct.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
.venv/bin/python scripts/test_large_reduction_binding.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
.venv/bin/python scripts/generate_large_predictor_seams.py --output /tmp/large-predictor.rs
.venv/bin/python scripts/verify_as1_atan.py --reference-root "$REFERENCE"
cargo test --locked --offline --features flashiv,experimental,cxx_bench
cargo test --locked --offline --features flashiv,experimental,fma,cxx_bench
cargo clippy --locked --offline --features flashiv,experimental,cxx_bench -- -D warnings
cargo clippy --locked --offline --features flashiv,experimental,fma,cxx_bench -- -D warnings
cargo fmt --all --check
python3 -m unittest discover -s benches -p 'test_compare*.py'
```

The component bounds retain their archived-model, IEEE arithmetic and libm
premises. They do not establish a whole-domain or correct-rounding theorem.
Raw local artifacts and replay helpers are in
`target/adaptive-physical-20261107/`, which is ignored by Git.
