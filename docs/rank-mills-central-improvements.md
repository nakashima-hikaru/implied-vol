# Rank compression, cheaper Mills evaluation and wider central expansion

Applied against `d51a2631019695887468d23d9b6c988ba0a3bc10` on 2026-10-23.
The [machine-readable report](rank-mills-central-improvements.json) records
source and executable hashes, independent-reference provenance, conditional
bounds, accuracy changes, counters and paired timings.

## Adopted changes

Small rectangle 2 now uses rank two instead of rank three, with the existing
degree-6/5 specialization. Its seed uses 11 rather than 22 two-lane Horner
stages. This removes 224 bytes of coefficient payload; all 26 small/core tables
now total 11,104 bytes, plus the unchanged 52 bytes of degree metadata.
The represented factors are frozen in the generator. Verification uses exact
bivariate Chebyshev coefficient differences, independently of how SVD found
the factors. It includes rounded coordinates, old/new source evaluation,
positive seed floors, transported root geometry and the cubic finish.
The conditional final bound on all 139 rectangle-2 leaves is below `0.590475`.
See [the seed derivation](seed-reduction.md).

The compensated large Mills table evaluator uses degree 16 instead of 24.
All three compensated final Horner stages and low-argument correction remain;
asymptotic lanes retain degree 18. Two table lanes save eight Horner stages;
a mixed table/asymptotic pair saves six. The exact deleted coefficient tail
satisfies `(|tail1|+|tail2|)/(s*u)<0.010460408` on the accepted region.
The checker also recomputes approximation and source-rounding errors against
the true Mills function, then transports them to the original price input.
Its new conditional final bound is below `0.873762`.

The upper large predictor uses the same physical recurrence stopped at order
five instead of seven. Its output still passes the existing bracket and
precise residual checks. The stopping proof uses the preserved safeguarded
iteration and forced bisections; it does not assume that the predictor reaches
the finalizer gate in one step. See [the large derivation](large-direct-proof.md).

The deferred central expansion extends its guards from `a<=0.36, q.high<=0.5`
to `a<=0.45, q.high<=0.625`. Its coefficients and arithmetic are unchanged.
The continuous conversion remainder is freshly recomputed on
`a²<=0.203, v²<=0.391`; its relative contribution is below `9.395e-19`.
Recomputing native source-order rounding gives conditional `rho_J<0.543694`.
The normal-shape certificate retains its existing padded z covers. This is
a modest extension certified through the existing expansion; no conformal
wing expansion is introduced. See [the central derivation](central-expansion.md).

## Accuracy and work counts

Both FMA feature policies pass the same **52,617 unique exact binary64
normal-price reference inputs**, with no invalid output or `rho_J>=1`.
The maximum remains `0.6172411593798071`. There are 29 changed outputs; their
maximum `rho_J` is `0.4294627846756894`. Both feature policies produce identical
output bytes. Saved independent roots are finite evidence, not a global proof.

New physical input fixtures independently agree at 160 and 240 decimal digits:

| Corpus | Inputs | Baseline max rho | New max rho | Changed outputs |
| --- | ---: | ---: | ---: | ---: |
| Rectangle-2 interior | 625 | 0.427383 | 0.455327 | 44 |
| Central guard seams/interiors | 234 | 0.364741 | 0.319375 | 19 |
| Existing adaptive seed seams | 975 | 0.476817 | 0.476817 | 11 |

Individual errors can increase while remaining inside the contract. The
unchanged 156 common-finish fixtures also pass both policies, with identical
baseline and new output bits.

An isolated instrumented replay matches the corresponding frozen output bytes.
For all 47,070 large inputs, precise calls remain 47,070 and fast calls remain
26,067. All 21,003 upper inputs enter the finalizer after the initial physical
predictor. Their maximum first residual increases from `1.123001e-7` to
`8.902104e-6`, below `2^-14`. This is empirical one-step evidence only.
The reductions save arithmetic within evaluations, rather than evaluator calls.

## Timing observations

Immutable baseline and combined executables used identical benchmark source,
separate source/target directories, rustc 1.98.1 / LLVM 22.1.8,
`aarch64-apple-darwin`, native target CPU and bench-profile LTO. Six rotating
ABBA/BAAB blocks used 300,000 calls and two rounds per invocation; case and
solver order also alternate. No compilation or proof work from this task ran
during the timing window. The raw samples and executable SHA-256 hashes are
retained in `target/proposal-20261023/combined.timing.json`.

| Workload | Baseline ns | New ns | Observed change |
| --- | ---: | ---: | ---: |
| central_deferred | 77.41 | 76.48 | -1.19% |
| central_extended | 108.02 | 81.51 | -24.54% |
| rank_small2 | 105.91 | 104.23 | -1.59% |
| finite_seed | 113.04 | 112.40 | -0.56% |
| wing_seed | 109.09 | 109.56 | +0.43% |
| mixed_near_atm | 85.97 | 84.57 | -1.63% |
| mixed_normalised | 126.11 | 126.41 | +0.23% |
| mixed_full | 145.33 | 143.84 | -1.03% |
| large_lower | 201.87 | 192.25 | -4.77% |
| large_upper | 260.22 | 230.57 | -11.39% |
| mixed_large | 207.78 | 204.48 | -1.59% |

These are **observations excluded from firm latency claims**. Unchanged Hybrid
controls have paired-block variation up to 11.61%, and one mixed-full candidate
block changes by +87.83%. Every block is retained. The fixed central extension
and large changes shorten the arithmetic as intended, but noisy timing does
not establish a universal or independently qualified speedup.

## Other proposal candidates

`scripts/probe_quartic_seeds.py` reproduces the proposed joint degree/correction
experiment at 90 and 130 digits, on both mathematical and represented-input
grids. For small rectangle 1, degree 6/5 gives maximum ideal `rho_J` about
`9.0608` with the cubic and `0.004247` with the quartic. Degree 5/5 still fails
with the quartic (`4.8293`). Reducing 7/5 to 6/5 saves two two-lane Horner
stages; a reused-geometry quartic adds roughly two scalar multiplies and three
scalar FMAs. A net cost improvement is not established, and the quartic still
needs a continuous fifth-inverse-derivative and source-rounding proof.
The runtime therefore retains the certified cubic. Optional rank-two probes
identify other potential rectangles, but these finite ideal grids do not
qualify additional runtime changes.

The [weighted direct-D probe](weighted-direct-d-research.md) makes the proposed
`abs(Dhat-D)/(1+D)` weight explicit by fitting `g=D/(1+D)` and reconstructing
`D=g/(1-g)`. It retains represented coefficient bits and exact positive
denominator bounds. Some cells improve and a near-ATM candidate worsens.
Its finite ideal-arithmetic grids do not certify a continuous approximation
remainder, FMA evaluation, final solver precision or lower cost. It is retained
as a reproducible research probe, with no runtime replacement.

## Proof scope and reproduction

The continuous component bounds are conditional on hash-bound archived
seed/Mills/normal-shape certificates and primitive rounding assumptions.
They establish neither a whole-domain solver theorem nor completed Lean or
vendor-libm verification. The changed domains, approximation errors and
source-order rounding are recomputed; old final-error reports are not used
as bounds for the changed kernels.

The source bridges reverse only explicitly reviewed guard/representation
edits, then require the original source hashes. Deliberate coefficient,
rank, loop-shape, compensation, predictor and domain defects must be rejected.

Both complete feature-policy test runs pass **119 tests each**, including six
doctests. Formatting, the library/path-benchmark Clippy check with `-D warnings`,
and nine Python benchmark-contract tests pass. Windows MSVC x86/x64 and
WebAssembly library compilation with all three solver features passes; these
are compilation checks, with no native target execution. Ten seed defects are
rejected by both bridge consumers, and seven large plus seven central defects
are rejected by their respective source checks.

```sh
python3 scripts/generate_seed_reduction.py --reference-root "$REFERENCE" --check --check-svd
python3 scripts/verify_direct_d.py --reference-root "$REFERENCE"
python3 scripts/verify_as1_atan.py --reference-root "$REFERENCE"
python3 scripts/verify_large_fast.py --reference-root "$REFERENCE" --out target/proof/large-fast-proof.json
python3 scripts/verify_large_direct.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
python3 scripts/verify_central_expansion.py --reference-root "$REFERENCE"
python3 scripts/verify_seed_bindings.py --reference-root "$REFERENCE"
python3 scripts/test_large_reduction_binding.py --reference-root "$REFERENCE" --fast-certificate target/proof/large-fast-proof.json
python3 scripts/test_central_expansion.py --reference-root "$REFERENCE"
```

Use the project's Python environment for commands requiring mpmath. The
optional SVD provenance check is separate from the exact represented-factor
certificate. The fixture generators reproduce the committed binary64 root
high/low and eta values at two precision settings.
