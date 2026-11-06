# Further reductions: central conversion and fast Mills evaluation

Implemented on 2026-11-06 after the
[normal-shape and Householder changes](theoretical-cost-reduction.md).
Those changes form this report's frozen baseline. The
[machine-readable record](further-theory-improvements.json) binds the current
sources, reference corpus, proof outputs, executables and measurements.

## Central conversion

The Black-conversion polynomial previously evaluated all 98 stored
coefficients in ten rows. It now retains the lowest powers of each row with
lengths **12/11/9/7/5/3/1/0/0/0**, evaluating **48 coefficients**. The final
three rows vanish and row 6 becomes a constant. Source coefficient bits and
normal-shape degrees are unchanged. Constant-degree Horner evaluation saves
**22 two-lane FMAs plus six scalar FMAs** per accepted central call.

The choice is supported by the continuous model rectangle
`A=a²<=.203, w=v²<=.391`. An exact coefficient norm bounds the deleted
monomials; that diagnostic alone is not the precision proof. The checker
instead rebuilds the nonlinear conversion residual `K(C)-exp(w/8)T` using
the shortened polynomial `C`, with exact rational interval arithmetic,
analytic tails and the derivative bound needed to transfer residual to root
error. It also recomputes the retained Horner and coordinate-rounding ledger.

The resulting relative conversion-model bound is below `6.913492e-18`.
The complete helper's conditional bound is **`rho_J<0.570599`**, compared
with `0.543694` before this reduction. The larger bound remains comfortably
inside the target `rho_J<1`. This trades unused approximation accuracy for
fewer operations. The [derivation](central-expansion.md) records the exact
row selection and error transport; acceptance guards are unchanged.

## Fast Mills evaluation

The inexpensive large-route Mills table now uses degree **13 instead of 16**.
It therefore matches the existing degree-13 asymptotic approximation. Table
evaluations and mixed table/asymptotic pairs omit three Horner stages; mixed
pairs no longer need leading zero padding. The compensated final evaluator
still uses table degree 16 and asymptotic degree 18.

This distinction follows the accuracy needed by each calculation: the fast
evaluation supplies reliable bracket signs and a handoff decision, while the
compensated evaluation determines the final result. Recomputing all 34 table
cells and 114 continuous residual boxes gives the same worst fast residual
bound **`7.752447e-13 < 2^-38`**. The sign and stopping margins, precise gate
and conditional final bound `rho_J<0.873762` remain valid. See the
[large-route derivation](large-direct-proof.md#degree-13-fast-mills-evaluation).

Shortening the upper physical predictor uniformly to reciprocal derivative
index 3 or 4 was also probed. Some saved inputs missed the first precise gate,
which could add evaluator calls. That alternative was not adopted.

## Accuracy and validation

Both FMA feature policies pass the same **52,617 unique exact binary64
normal-price reference inputs**, with no invalid result or `rho_J>=1`.
Maximum `rho_J` remains `0.6172411593798071`, and the policies produce
identical output bytes. Three outputs change, all in the large route: one
improves and two worsen, with changed-point maximum `rho_J=0.3971683194`.

All **2,893 additional reference fixtures** retain their output bits and pass:
1,026 central seams/interiors, 111 large handoff/tail cases, 625 rank-compression
cases, 975 adaptive-seed cases and 156 common-finish cases. Their exact-input
root references are independent of the candidate algorithms.

Instrumented replay matches the corresponding uninstrumented bytes. All
26,067 lower large inputs still use one fast and one precise evaluation; all
21,003 upper inputs use one precise evaluation and no fast evaluation. Neither
version uses bisection on this corpus. The largest first precise residual is
unchanged at `4.125034e-5` below and `8.902104e-6` above the cap-side split.
These counts describe this finite corpus, not universal one-step convergence.

The final CXX-enabled combined test suites pass **122 tests per FMA policy**,
including doctests. Unscoped warning-denying Clippy passes both policies.
Formatting and all 23 Python benchmark-contract tests pass. Library checks
also pass for Windows MSVC x86/x64 and WebAssembly with both policies; those
six checks are cross-compilations, not native execution.

The central checker passes one positive control and rejects 19 deliberate
source/proof mutations. The large checker rejects 13 faults, with its three
new fast-evaluator faults independently rejected by both the fast component
and finalizer integration checkers. The existing rank/seed generator reproduces
its represented coefficients, including its optional SVD provenance check;
the seed bridge rejects ten faults through both of its consumers.

## Timing

Paired measurements compare this report's frozen baseline and final candidate
using identical benchmark source, native CPU optimization, bench-profile LTO,
and six rotating ABBA/BAAB blocks of 300,000 calls and two rounds. No compiler
or proof jobs from this task ran during measurement. The toolchain is rustc
1.98.1 / LLVM 22.1.8 on macOS AArch64.

| Workload | Baseline ns | New ns | Observed change |
|---|---:|---:|---:|
| central_deferred | 72.83 | 65.23 | -10.44% |
| central_shape1 | 74.82 | 67.14 | -10.27% |
| central_shape2 | 74.72 | 67.23 | -10.02% |
| central_extended | 80.48 | 70.77 | -12.08% |
| mixed_near_atm | 82.96 | 75.64 | -8.83% |
| large_lower | 196.77 | 193.04 | -1.89% |
| large_upper | 240.39 | 236.84 | -1.48% |
| large_asymptotic | 196.15 | 196.26 | +0.05% |
| mixed_large | 208.61 | 206.72 | -0.91% |
| mixed_normalised | 132.27 | 132.51 | +0.18% |
| mixed_full | 150.45 | 150.05 | -0.27% |

All six paired blocks improve for the four central cases, `large_lower`,
`large_upper` and `mixed_large`. Unchanged Hybrid control medians range from
-0.86% to +0.44%, while individual control blocks range from -7.71% to +3.33%.
All blocks are retained. The central reduction is clear in this experiment;
the smaller large-route timing changes do not establish a firm latency gain
beyond host variation. Reduced arithmetic is established independently of
timing. The linked JSON preserves every workload summary and raw artifact
hash; these synthetic observations are not a portable latency guarantee.

## Reproduction and limits

Use the same hash-bound `implied-black-volatility` reference checkout required
by the previous reports:

```sh
REFERENCE=/path/to/implied-black-volatility
.venv/bin/python scripts/verify_central_expansion.py --reference-root "$REFERENCE"
.venv/bin/python scripts/test_central_expansion.py --reference-root "$REFERENCE"
.venv/bin/python scripts/verify_large_fast.py --reference-root "$REFERENCE" --out target/proof/large-fast-proof.json
.venv/bin/python scripts/verify_large_direct.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
.venv/bin/python scripts/test_large_reduction_binding.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
.venv/bin/python scripts/verify_as1_atan.py --reference-root "$REFERENCE"
.venv/bin/python scripts/verify_seed_bindings.py --reference-root "$REFERENCE"
.venv/bin/python scripts/generate_seed_reduction.py --reference-root "$REFERENCE" --check --check-svd
cargo test --locked --offline --features flashiv,experimental,cxx_bench
cargo test --locked --offline --features flashiv,experimental,fma,cxx_bench
cargo clippy --locked --offline --features flashiv,experimental,cxx_bench -- -D warnings
cargo clippy --locked --offline --features flashiv,experimental,fma,cxx_bench -- -D warnings
cargo fmt --all --check
python3 -m unittest discover -s benches -p 'test_compare*.py'
```

The continuous bounds are conditional component results under the documented
archived-model, floating-point and libm premises. They are distinct from the
finite reference and timing evidence and do not establish a whole-domain or
correct-rounding theorem. Local raw artifacts and replay helpers live under
`target/further-theory-20261106/`, which is ignored by Git.
