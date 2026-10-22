# Experimental: compact tables and shorter evaluations

Implemented and checked on 2026-10-22/23 against baseline
`c0ecb7c45d144ccebf5b411fef6568030cebd04a`.
[Machine-readable results](experimental-improvements.json) retain source and
executable hashes, reference provenance, conditional bounds, all workload
summaries, and excluded timing runs. The earlier `results.json` remains a
record of its original source revision.

## Changes

- Small and large solvers share 44 Mills coefficient rows, storing 25 high
  parts and the three low parts actually read by compensated Horner evaluation.
  Every retained bit matches the archived table. The precise evaluation's
  coefficient order, lane order, reduction order and explicit FMA semantics
  are preserved.
- The 253 degree entries use `u8`, converted to `usize` at indexing sites.
- Large-region fast Mills evaluation uses degree 16 instead of 24. The
  compensated final evaluation retains degree 24; the asymptotic degrees
  are unchanged.
- Direct-D selects fewer recurrence terms at larger `h.high`:

| `t` | `h.high` | Terms |
|---|---|---:|
| `(0.5, 0.75]` | `<6`, `[6,7)`, `>=7` | 14, 12, 10 |
| `(0.75, 1]` | `<5`, `[5,6)`, `>=6` | 16, 14, 12 |

The dispatch domain and the counts for `t<=0.5` are unchanged.

## Storage

On a 64-bit target, coefficient payload falls from 17,600 to 9,856 bytes,
and degree payload decreases by 1,771 bytes: **9,515 bytes less static data**.
The baseline linker already merged the two identical original Mills arrays;
this calculation does not count their source duplication as another saving.

In the frozen native benchmark, Mach-O `__TEXT,__const` shrank from 78,680
to 69,016 bytes. Code size grew by 560 bytes, and page-aligned mapped segment
sizes stayed unchanged. The executable file grew by 208 bytes. These are
payload/section measurements, not a claim of reduced process RSS or file size.

## Accuracy and conditional bounds

- **52,617 unique exact-input normal-price references**: no invalid outputs,
  no `rho_J>=1`, and no output-bit changes versus baseline for each isolated
  candidate and the integrated implementation. Both FMA feature settings of
  the integrated build have maximum `rho_J=0.6172411593798071`.
- **123 new final-root cases** have independently inverted references whose
  high/low root and conditioning-estimate bits agree at 160 and 240 decimal
  digits. Actual seed tracing confirms all four shortened direct-D paths are
  exercised; 105 cases shorten the series and 18 retain their old count.
- **144 direct-D component cases** cover both sides of `h.high=5,6,7`, signed
  nonzero `h.low`, and both sides of the relevant `t` cutoffs. Both FMA
  settings pass the final-root and component checks.
- The exact rational direct-D checker compares complete truncation-plus-rounding
  bounds on all 77 cells. Twelve cells use fewer terms, and no cell's bound
  increases. It audits 5,474 seed rows and transports 5,457 accepted rows;
  the existing `h<=8` guard excludes 17 wing rows. Maximum D error is
  `0.8211800596659627u`, with conditional final bounds `0.604976705238934`
  (wing), `0.8155064548897512` (small), and `0.7976690502946979` (finite).
- The new large-fast checker covers 34 continuous table cells and 114 residual
  boxes. Its residual error bound is `7.752446906743081e-13 < 2^-38`, including
  the omitted represented coefficients, arithmetic rounding and geometry.
  It checks sign, precise-handoff and forced-bisection stopping margins.
- The AS1 checker passes with unchanged arithmetic and non-Mills constants.
  The source bindings preserve the earlier arithmetic hashes through an exact
  reversal of reviewed load/type notation, alongside independent coefficient
  bit comparisons. Deliberately changed high, low and degree entries were
  rejected even after refreshing the superficial compact-layout hashes.

These remain conditional component certificates. They retain the inherited
seed/Mills, elementary-function, FMA/sqrt, target/bracket and final-correction
assumptions described by their checkers. Finite references and these component
bounds do not establish an end-to-end Lean theorem or a whole-domain guarantee.
The reference corpus is the recorded 52,617-case set, not an asserted replay of
an older report's differently sized corpus.

## Timing

The baseline and each candidate used separate targets and frozen executables,
`flashiv,experimental`, FMA feature off, native CPU flags and bench-profile LTO.
The main comparison used 200,000 calls, two rounds and four alternating
ABBA/BAAB blocks, including extended, AS1 and dpoly populations. Checksums
remained stable and unchanged between versions. Rebuilding the final source
produced exactly the same integrated benchmark executable bytes.

Selected medians from the first integrated run:

| Workload | Baseline ns | Integrated ns | Change |
|---|---:|---:|---:|
| `large_lower` | 223.35 | 211.50 | -5.31% |
| `large_upper` | 280.25 | 268.51 | -4.19% |
| `mixed_large` | 222.39 | 217.10 | -2.38% |
| `highest` | 94.33 | 89.76 | -4.84% |
| `mixed_normalised` | 128.62 | 129.79 | +0.91% |
| `mixed_full` | 145.60 | 146.46 | +0.59% |
| `dpoly_mixed` | 109.13 | 110.48 | +1.24% |

All four paired blocks improved the listed large-region cases. The ordinary
mixed cases have mixed block signs; `dpoly_mixed` includes a +26% outlier block.
The unchanged Hybrid controls' median changes stayed within 0.75% in this run.
This supports a local large-region improvement, not a general speedup.

Later runs encountered substantial concurrent host load. The isolated fast and
direct-D runs showed unchanged Hybrid median changes up to 29.54% and 14.73%.
A further eight-block targeted direct-D comparison also had unstable unchanged
controls while host load exceeded 300. **Those timing runs are excluded from
latency conclusions**, with their summaries preserved. Direct-D's reduced
operation counts and accuracy are established; its isolated latency benefit
still needs measurement on a quiet host.

## Reproduction and checks

The archived reference checkout is `implied-black-volatility`; pass its path
explicitly. Checkers bind the required source/dependency hashes.

```sh
rtk proxy .venv/bin/python scripts/verify_direct_d.py --reference-root /path/to/implied-black-volatility
rtk proxy .venv/bin/python scripts/verify_large_fast.py --reference-root /path/to/implied-black-volatility
rtk proxy .venv/bin/python scripts/verify_as1_atan.py
rtk proxy .venv/bin/python scripts/generate_direct_d_seams.py --output /tmp/direct_d_roots.rs --component-output /tmp/direct_d_components.rs
rtk cargo test --locked --offline --features flashiv,experimental
rtk cargo test --locked --offline --features flashiv,experimental,fma
```

Both feature-scoped test suites, final expanded regression tests, formatting,
and warning-denying Clippy checks for the library and new integration test
passed. Library compilation also passed for Windows x86/x64 and
`wasm32-unknown-unknown`; this does not constitute execution on those targets.

Local raw timing/proof outputs, binary corpus replay tools, trace harnesses and
immutable executables are in `target/experimental-improvements-20261022/`.
That ignored working directory is not part of the committed source; the JSON
report retains its provenance and summaries, and the checked-in generators and
proof checkers reproduce the new boundary references and component bounds.
