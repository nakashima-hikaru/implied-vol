# Experimental: per-rectangle seeds and localized handoff

Implemented on 2026-10-23 against revision `0bd4312`. The
[machine-readable report](adaptive-theory-improvements.json) preserves the
source/toolchain hashes, independent-reference provenance, conditional bounds,
adverse accuracy cases, evaluation counts, and timing summaries. This follows
the [uniform-degree and direct-cubic changes](theory-improvements.md).

## Adopted changes

The 26 small/core rank tables now select degrees independently: twelve use
6/5, twelve use 6/6, one uses 7/5, and one uses 7/6. Exact rational Chebyshev
reduction generates the coefficients from the original represented factors.
The continuous final-error budget determines which degree pairs are admissible;
the selection is not claimed globally optimal.

Coefficient payload falls from **14,400 to 11,328 bytes** (21.33%). Explicit
degree metadata adds 52 bytes, giving a **net saving of 3,020 bytes** (20.97%)
for these tables and metadata. This is static data payload, not RSS or binary
size. Padded ranks are retained. Rank-three/four evaluations use 22, 24 or 26
two-lane Horner stages instead of 28; rank-two uses 11 instead of 14. Fixed
specializations preserve lane and rank-reduction order. The cubic common
finish remains unchanged. See the [seed derivation](seed-reduction.md).

The large-route cubic finalizer now accepts a precise residual high part up
to `2^-14`, previously `2^-16`. The stronger proof bounds the actual accepted
inverse-Taylor segment: the lower log target implies `x<.35`, and both routes
then satisfy `abs(lambda)<.4`, `abs(alpha)<1`, `beta<5`, and `abs(gamma)<8`.
This reduces the cubic remainder bound sufficiently for the wider gate.
The source's `2^-34` fast prefilter and cubic arithmetic remain unchanged.
Target low parts and the entire segment are included in the bound; the
original-price transfer retains its wider prior allowance. See the
[large derivation](large-direct-proof.md).

## Accuracy and proof scope

Both FMA feature settings passed **52,617 unique exact binary64 normal-price
inputs**, with no invalid result or `rho_J>=1`. The maximum remains
`0.6172411593798071`. FMA-on and FMA-off output bytes are identical.

The isolated large change preserves every output bit. The seed and combined
candidates each change 222 outputs; their maximum `rho_J` among changed cases
is `0.4095085658923388`. Five existing exact-output seam fixtures change by
one ulp. Independent roots at 160 and 240 digits qualify each change:

| Input `(a bits, b bits)` | Previous rho | Current rho |
|---|---:|---:|
| `3fb999999999999a, 3fc3333333333334` | 0.138140 | 0.142204 |
| `3ff0000000000000, 3fd0000000000001` | 0.198262 | 0.157299 |
| `3ff0000000000001, 3fbfffffffffffff` | 0.242627 | 0.331754 |
| `4014000000000000, 3f8fffffffffffff` | 0.347366 | 0.275547 |
| `4014000000000000, 3f90000000000000` | 0.260552 | 0.362361 |

All 15 existing seam cases remain. The goal is cheaper seeds within the
accuracy contract; individual errors do not all decrease.

A new physical-price corpus contains 975 boundary-neighborhood and interior
inputs. Independent root high/low and conditioning values agree at 160 and
240 digits. Both final binaries pass, with maximum `rho_J=0.4768161275734432`
and 65 changed outputs. Large component fixtures expand from 306 to 564,
covering the new gate, its neighbors, Mills seams and the lower-routing seam
with both signs of the target low part. These component targets are separate
from the physical-price corpus. Finite tests are not continuous proofs.

| Component | Conditional rho upper bound, rounded upward |
|---|---:|
| Large, including original-price transfer | 0.874279206 |
| Wing common finish | 0.593445984 |
| Small-rank common finish | 0.884890543 |
| Finite-rank common finish | 0.785460849 |

The integrated checker regenerates every degree and coefficient entry and
recomputes the direct-D, perturbed-seed and final cubic budgets. It audits
5,474 seed leaves and transports 5,457 after the existing guard excludes 17
wing leaves. The large checker consumes a freshly recomputed fast-component
certificate and binds the actual target, gate and source operations. The AS1
atan-helper bound remains `1.3586386736442848e-16`.

Eight deliberate seed defects are rejected by both the shared direct-D bridge
and the AS1 consumer (16 rejections); eleven large/source/certificate controls
are also rejected. These checks test binding integrity, separately from the
mathematical inequalities. All proofs remain conditional on the archived
seed/Mills/analytic certificates and primitive rounding/libm assumptions.
They do not establish an end-to-end Lean theorem or vendor-libm guarantee.

## Work counts and timing

A separate instrumented replay preserves the uninstrumented outputs for all
52,617 inputs. On its 47,070 large inputs, both versions call the precise
evaluator 47,070 times and the fast evaluator 26,067 times. **The wider gate
saves zero evaluations in this corpus.** It broadens the certified handoff
region without demonstrating fewer iterations on these inputs.

Timing used immutable baseline and combined binaries, the same harness,
rustc 1.98.1 / LLVM 22.1.8 on `aarch64-apple-darwin`, native target CPU and the
LTO bench profile. Six ABBA/BAAB blocks used 150,000 calls and two rounds per
invocation. No validation/compiler processes from this task ran concurrently.

| Workload | Baseline ns | Combined ns | Observed change |
|---|---:|---:|---:|
| central_deferred | 77.49 | 77.50 | +0.01% |
| finite_seed | 114.06 | 113.42 | -0.56% |
| wing_seed | 109.77 | 109.79 | +0.02% |
| mixed_near_atm | 89.32 | 86.62 | -3.02% |
| mixed_normalised | 126.42 | 126.53 | +0.08% |
| mixed_full | 143.58 | 143.77 | +0.13% |
| large_lower | 203.10 | 202.99 | -0.05% |
| mixed_large | 208.96 | 209.05 | +0.04% |

These are observations, **excluded from firm latency claims**. Unchanged
Hybrid controls have small overall median differences but a paired-block
outlier of 13.21%; candidate blocks also contain large outliers. Host load was
high. All samples and checksums are retained under
`target/adaptive-theory-20261023/combined.timing.json`; the report includes their
hash and per-block summaries. No blocks were silently removed. The supported
benefits are the smaller static tables, fewer Horner stages and wider
conditionally certified handoff region.

## Validation and reproduction

- `cargo test --locked --offline --features flashiv,experimental` and the
  same command with `,fma`: 117 tests each, including six doctests.
- `cargo fmt --all --check` and Clippy for the library and path benchmark
  with `-D warnings`: passed.
- Library compilation for Windows MSVC x86/x64 and WebAssembly: passed;
  these are compilation checks, not native execution.
- Integrated direct-D, AS1, large-fast and large-direct certificates, plus
  coefficient regeneration: passed against the current source.

The repository generators and derivations document the independent fixture
and certificate commands. Use the reviewed `implied-black-volatility` archive
as `--reference-root`; keep its source hashes fixed. Benchmark reproduction
requires separate targets, `cargo bench --no-run` with the stated features and
`RUSTFLAGS='-C target-cpu=native'`, frozen executables, and identical input
populations. A successful `cargo build --release --bench` does not use the
same bench-profile LTO settings and was not used for the reported comparison.
