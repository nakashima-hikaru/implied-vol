# Experimental: direct inverse correction and reduced seeds

Implemented on 2026-10-23 against revision `1f703b3`. The
[machine-readable report](theory-improvements.json) records source and executable
hashes, exact-input reference provenance, conditional bounds, adverse cases,
and all timing summaries. This follows the
[compact-table changes](experimental-improvements.md).

## Adopted changes

The large route now maps the compensated log-price residual directly to a
relative volatility correction. It retains the compensated current volatility
and removes the old coordinate correction and its square-root fallback.
The conditional source-order error budget permits precise handoff at `2^-16`,
instead of `2^-22`. The fast residual, target construction and brackets retain
their existing arithmetic. See the [derivation](large-direct-proof.md).

The small/core rank seeds now use degree eight in the first coordinate and
degree six in the second. Exact Chebyshev-mode reduction replaces naive
truncation of power coefficients. The 26 tables shrink from **17,600 to 14,400
bytes**, saving **3,200 bytes of static coefficient payload**. Rank-two
evaluation removes six two-lane Horner stages, and rank-four evaluation removes
eight. This does not claim an equivalent reduction in process RSS or file size.

A direct fourth-derivative bound for the inverse price proves that the reduced
seeds still work with the existing cubic finish. The initially tested quartic
finish was therefore unnecessary and was not adopted. Lane order, scalar rank
reduction, price residual and final FMA remain unchanged. See the
[seed derivation and generator](seed-reduction.md).

## Accuracy and conditional certificates

Both FMA feature settings passed **52,617 unique exact binary64 normal-price
inputs**, with no invalid outputs or `rho_J>=1`. Maximum `rho_J` is
`0.6172411593798071`. The two policies produce identical output bytes.

The isolated large change preserves every output bit in this corpus. The
integrated implementation changes 173 outputs because of the reduced seeds;
the maximum `rho_J` among changed cases is `0.4077440718963207`.

Two existing seam fixtures change by one ulp. Independent Black inversion at
100 and 160 digits confirms both: the first improves from `rho_J=0.227569` to
`0.052776`; the second increases from `0.275547` to `0.347366`. All 15 seam
cases and their exact output regression checks remain. This is not a claim
that every individual error decreased.

Another 306 component fixtures target the large handoff gates and Mills seams.
Their target/root/conditioning bits agree at 160 and 240 digits. They check
the direct correction even where ordinary solver iterations would skip the
gate boundary; they are distinct from the physical-price input corpus.

| Component | Conditional `rho_J` upper bound |
|---|---:|
| Large, including original-price target transfer | 0.896240684 |
| Wing common finish | 0.593445984 |
| Reduced small-rank seeds | 0.797816752 |
| Reduced finite-rank seeds | 0.785245451 |

The small checker audits all 5,474 seed leaves, with the existing guard
excluding 17 wing leaves. It regenerates all reduced coefficients and expands
the inherited geometry/residual bounds before evaluating the cubic remainder
and rounding ledger. The large checker separately verifies the new finalizer
and its integration with the fast component. The unchanged AS1 checker uses a
narrow, explicitly approximation-changing source bridge for the seed edits.

Six deliberate large defects, five seed defects and four AS1 defects were
rejected. The large controls were rerun against the final integrated checker.
Coefficient regeneration and source bridges preserve reviewed original hashes;
changing an expected hash alone does not qualify new arithmetic.

These are conditional component certificates. The archived seed cover, Mills
analysis, primitive error contracts and external cap-separation result remain
premises. Finite references and these checks do not establish an end-to-end
Lean theorem or a whole-domain/vendor-libm guarantee.

## Timing observations and limits

Baseline and candidates used separate targets and frozen executables, Rust
1.98.1 / LLVM 22.1.8 on `aarch64-apple-darwin`, native CPU flags and LTO.
The benchmark now accepts `--case NAME`, allowing each population to be paired
closely in time. The final comparison used six alternating ABBA/BAAB blocks,
150,000 calls and two rounds, with unchanged Hybrid controls.

| Final integrated observation | Baseline ns | Candidate ns | Change |
|---|---:|---:|---:|
| `large_lower` | 213.13 | 205.39 | -3.63% |
| `large_upper` | 275.71 | 266.17 | -3.46% |
| `large_asymptotic` | 200.60 | 192.74 | -3.92% |
| `mixed_large` | 216.47 | 208.34 | -3.76% |
| `finite_seed` | 113.67 | 113.71 | +0.04% |
| `wing_seed` | 110.37 | 109.44 | -0.85% |
| `mixed_normalised` | 130.30 | 127.14 | -2.42% |
| `mixed_full` | 143.43 | 142.00 | -1.00% |

**These observations are excluded from firm latency claims.** Earlier runs had
unchanged-control median changes up to about 28% and 15%. In the final run,
control medians were within 1.74%, but individual paired control blocks still
varied by up to 73%, and one `mixed_large` candidate block regressed by 26.9%.
All blocks and unfavorable observations are retained. The established benefits
are fewer operations and smaller tables; a stable latency comparison remains
unqualified on this busy host.

A positive continued-fraction alternative received a separate real-arithmetic
feasibility probe. On the tested depth grid, enclosure widths below one tenth
of the conditioning budget first appeared at depths 33–96 for the sampled
`h=2..8`, `t=0.5,0.75,1`. Depths below 33 were not tested. Binary64 rounding and
runtime were not qualified, so that alternative was not adopted.

## Reproduction and checks

Pass the archived `implied-black-volatility` checkout explicitly. For example:

```sh
rtk proxy .venv/bin/python scripts/generate_seed_reduction.py --reference-root /path/to/implied-black-volatility --check
rtk proxy .venv/bin/python scripts/verify_direct_d.py --reference-root /path/to/implied-black-volatility
rtk proxy .venv/bin/python scripts/verify_large_fast.py --reference-root /path/to/implied-black-volatility --out /tmp/large-fast.json
rtk proxy .venv/bin/python scripts/verify_large_direct.py --reference-root /path/to/implied-black-volatility --fast-certificate /tmp/large-fast.json
rtk proxy .venv/bin/python scripts/verify_as1_atan.py --reference-root /path/to/implied-black-volatility
rtk cargo test --locked --offline --features flashiv,experimental
rtk cargo test --locked --offline --features flashiv,experimental,fma
```

Both complete feature-scoped suites passed 116 test functions, including
doctests. Formatting and warning-denying Clippy for the library and path
benchmark passed. Library compilation passed for Windows x86/x64 and
WebAssembly; those checks do not constitute execution on those targets.

Frozen binaries, complete raw samples and proof outputs remain locally under
`target/theory-improvements-20261023/`, which is ignored by Git. The repository
scripts, generators, fixtures and JSON summary retain the reproducible methods
and provenance; the raw local working directory is not a published artifact.
