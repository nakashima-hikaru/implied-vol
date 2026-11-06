# Experimental: matching approximation work to the precision budget

Implemented on 2026-11-06 against the existing working tree on `d51a263`.
The preceding rank compression, degree-16 Mills evaluation and widened central
domain were already present and are part of this comparison's baseline.
The [machine-readable report](theoretical-cost-reduction.json) records source
hashes, conditional bounds, independent-reference provenance and measurements.

## Adopted changes

The improvement is to spend only the approximation work needed by the final
precision budget. It changes two kernels, without widening an acceptance gate.

The central normal-shape polynomial now uses degrees **18/20/20/25** across
its four existing intervals, instead of 25 throughout. Exact represented
coefficient-tail bounds, including padded coordinates and the derivative
times the coordinate low part, justify removing **7/5/5 SIMD Horner stages**.
Both compensated final stages and the low-coordinate correction remain.
The maximum conditional bound is unchanged at `rho_J<0.543694`, set by the
last interval. Tables retain their original bits and storage size. See the
[central derivation](central-expansion.md).

The large solver now proposes a Householder update with reciprocal derivative
index **d=3**, instead of d=5. For the log-price residual `f`, write

```text
n = -f/f', A = f''/f', B = f'''/f',
step = n*(1+A*n/2)/(1+A*n+B*n*n/6).
```

The checker verifies this identity from `3*(1/f)''/(1/f)'''` and checks its
inverse series through cubic order using exact rational polynomial arithmetic.
It has fourth-order local convergence, versus sixth-order for the previous
proposal. It removes the calculations of `H2`, `H3`, `C` and `D` from the old
derivative recurrence. This is useful because the proposal need only reach
the compensated finalizer's `2^-14` residual gate.

The existing precise gate, fast-sign margin, finite/in-bracket checks and
forced midpoint every fourth update are preserved. The conditional stopping
argument relies on those safeguards, with no universal one-step convergence
premise. The final cubic correction and its conditional original-input bound
`rho_J<0.873762` remain unchanged. See the
[large derivation](large-direct-proof.md#shortened-logarithmic-iteration).
For general background on local convergence and bisection, see
[NIST DLMF §3.8](https://dlmf.nist.gov/3.8).

## Accuracy and evaluator counts

Both feature policies freshly replayed **52,617 unique exact binary64
normal-price inputs** against the archived independent high-precision roots.
There are no invalid results or `rho_J>=1`. The maximum remains
`0.6172411593798071`; FMA-on and FMA-off produce identical output bytes.
All 46 changed outputs are in the large route: 23 improve and 23 worsen.
Their maximum `rho_J` is `0.4590707024760558`. The objective is lower work
within the same precision contract, not a pointwise error improvement.

Additional physical-input fixtures use independent Gaussian-erfc Black roots
whose represented root-high, root-low and conditioning values agree at 160
and 240 decimal digits. Both feature policies pass:

| Fixture set | Inputs | Baseline maximum rho | New maximum rho | Changed outputs |
|---|---:|---:|---:|---:|
| Central guards, polynomial seams and interiors | 1,026 | 0.400085 | 0.400085 | 1 |
| Large dispatch, handoff and deep tails | 111 | 0.471943 | 0.471943 | 1 |
| Existing rank-compression fixtures | 625 | 0.455327 | 0.455327 | 0 |
| Existing adaptive-seed fixtures | 975 | 0.476817 | 0.476817 | 0 |
| Existing common-finish fixtures | 156 | 0.446674 | 0.446674 | 0 |

Instrumented replay reproduces the uninstrumented output bytes on the full
corpus. All 26,067 lower large inputs still use one fast and one precise
evaluation; all 21,003 upper inputs still use one precise evaluation and no
fast evaluation. Neither version uses a bisection on these inputs.
The largest first precise residual below the cap-side split increases from
`5.003285e-8` to `4.125034e-5`, still below `2^-14 ≈ 6.103516e-5`.
The upper maximum remains `8.902104e-6`. Thus this corpus shows cheaper
proposals with no added evaluations. It is finite performance evidence,
separate from the safeguarded stopping argument.

## Timing

Six rotating ABBA/BAAB blocks use 300,000 calls and two rounds per invocation.
The two immutable executables use identical benchmark source, native CPU
optimization and bench-profile LTO with rustc 1.98.1 / LLVM 22.1.8 on
`aarch64-apple-darwin`. No compiler or proof jobs from this task ran during
measurement. Inputs are prepared before timing.

| Workload | Baseline ns | New ns | Observed change |
|---|---:|---:|---:|
| central_deferred (degree 18) | 77.24 | 70.35 | -8.92% |
| central_shape1 (degree 20) | 77.55 | 72.23 | -6.86% |
| central_shape2 (degree 20) | 81.88 | 72.16 | -11.88% |
| central_extended (degree 25) | 81.82 | 77.70 | -5.03% |
| mixed_near_atm | 85.92 | 79.76 | -7.17% |
| large_lower | 194.24 | 190.23 | -2.06% |
| large_asymptotic | 193.39 | 189.49 | -2.02% |
| mixed_large | 207.98 | 201.40 | -3.16% |
| mixed_normalised | 128.44 | 127.44 | -0.78% |
| mixed_full | 144.85 | 144.61 | -0.17% |

All six paired blocks improve for the first three central cases and the
three listed large cases. Unchanged Hybrid controls have overall median
changes from -0.59% to +0.21%, with individual paired blocks from -2.41% to
+1.36%. Some unchanged Experimental paths also shift; notably, the degree-25
central interval evaluates the same polynomial but shares the changed
dispatch/code layout. Its measured improvement is not a deleted-term benefit.
These observations support the intended cost reduction on these workloads,
but do not establish a portable or universal latency improvement, especially
for changes comparable to control variation. No sample blocks were removed.

The JSON includes every workload summary and paired-block change. Raw samples,
executable hashes and local replay helpers are retained under
`target/theory-20261106/`. These local artifacts are ignored by Git.

## Reproduction and scope

```sh
REFERENCE=/path/to/implied-black-volatility
.venv/bin/python scripts/verify_central_expansion.py --reference-root "$REFERENCE"
.venv/bin/python scripts/test_central_expansion.py --reference-root "$REFERENCE"
.venv/bin/python scripts/verify_large_fast.py --reference-root "$REFERENCE" --out target/proof/large-fast-proof.json
.venv/bin/python scripts/verify_large_direct.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
.venv/bin/python scripts/test_large_reduction_binding.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
.venv/bin/python scripts/verify_direct_d.py --reference-root "$REFERENCE"
.venv/bin/python scripts/verify_as1_atan.py --reference-root "$REFERENCE"
.venv/bin/python scripts/generate_central_expansion_seams.py --output tests/experimental_central_expansion.rs
.venv/bin/python scripts/generate_large_iteration_cases.py --output tests/experimental_large_iteration.rs
cargo test --locked --offline --features flashiv,experimental
cargo test --locked --offline --features flashiv,experimental,fma
cargo fmt --all --check
cargo clippy --locked --offline --features flashiv,experimental --lib --bench implied_black_vol_paths -- -D warnings
python3 -m unittest discover -s benches -p test_compare_experimental.py
```

The new component certificates recompute the changed approximation/arithmetic
and preserve exact source bridges. Thirteen central and ten large deliberate
mutations are rejected alongside positive controls. Complete Rust tests pass
**120 tests per FMA policy**, including six doctests, as do formatting, the
library/path-benchmark Clippy check and nine Python benchmark-contract tests.

The bounds remain conditional on the hash-bound archived analytic models,
seed/Mills certificates and the documented floating-point/libm assumptions.
They cover the changed components, not an end-to-end Lean theorem or every
input in the full public API. Native execution here is on macOS AArch64.
