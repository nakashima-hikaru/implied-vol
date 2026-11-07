# Experimental rank-three factor evaluation: scalar-tail A/B experiment

2026-11-07. Baseline: `fcd44963881830e8a07596da1d7c117e0792e901`.

## Change selected

For the rank-three seed tables, the original evaluator executes four
`F64x2` Horner chains: two pairs of x-factors and two pairs of y-factors.
Only the first three factors of each coordinate participate in the final
rank reduction. The last SIMD lane of the second pair in each coordinate is
unused.

The selected evaluator keeps the **first two factors in each coordinate**
in the original two-lane Horner chain and evaluates the third x-factor and
third y-factor with independent scalar `f64::mul_add` chains. It avoids
evaluating the unused fourth factors entirely.

For degrees `(dx, dy)`, the original has `2dx + 2dy` two-lane FMA
stages; the candidate has `dx + dy` two-lane stages **plus**
`dx + dy` scalar FMA stages. This reduces evaluations of unused
lane elements but does not establish a reduction in the number of native
machine instructions. Live-factor coefficient bits and Horner order are
unchanged, as is the scalar rank reduction order. Other ranks and all
downstream corrections are unchanged.

No approximation is replaced or truncated.

## Rejected alternatives

Two other bit-preserving rank-three layouts were measured under the
same pinned-base, identical-harness comparison protocol:

| Alternative | rank_small3 | finite_seed | mixed_normalised | Outcome |
|---|---:|---:|---:|---|
| Three cross-coordinate SIMD pairs | +2.40% | +3.31% | +0.53% | Rejected |
| Two original pairs plus one cross-coordinate pair | +1.97% | +1.98% | +0.35% | Rejected |
| Two original SIMD pairs plus scalar third factors | **-2.32%** | **-2.76%** | **-1.86%** | Selected |

Signs are **candidate latency relative to pinned baseline** (negative is
faster). Each row comes from a distinct GitHub Actions host, so absolute
times must **not** be compared across rows. The comparisons within a row
use immutable executables built on the **same** host.

Raw samples, executable hashes, and per-block comparisons:
- [First cross-coordinate candidate](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42369781944/artifacts/13573088000)
- [Hybrid-pair candidate](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42370072121/artifacts/13573442543)
- [Selected scalar-tail candidate](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42370280616/artifacts/13573785337)

## Selected candidate measurement

GitHub Actions `ubuntu-24.04`, AMD EPYC 9V74, rustc 1.99.1, LLVM 23.1.1,
`--features experimental`, `RUSTFLAGS="-C target-cpu=native -C target-feature=+fma"`,
bench-profile LTO. Builds use identical benchmark source, distinct target
directories, pinned baseline, and the same compiler flags. The comparison
uses 100,000 timed calls per invocation, two rounds and four rotating
ABBA/BAAB blocks (16 samples per version and workload). All benchmark
checksums agree exactly across versions.

| Workload | Baseline ns | Candidate ns | Change | Four paired-block changes |
|---|---:|---:|---:|---|
| rank_small3 | 151.24 | 147.73 | -2.32% | -2.42%, -3.45%, -3.08%, -2.40% |
| finite_seed | 154.09 | 149.84 | -2.76% | -5.15%, -2.66%, -2.89%, -2.85% |
| lower_middle | 152.87 | 149.22 | -2.38% | -2.53%, -0.91%, -3.96%, -2.46% |
| upper_middle | 184.68 | 181.40 | -1.78% | -1.89%, -1.99%, -1.07%, -1.76% |
| mixed_normalised | 143.21 | 140.55 | -1.86% | -2.09%, -1.47%, -2.10%, -1.39% |
| mixed_full | 162.98 | 160.63 | -1.44% | -1.33%, -1.48%, -1.69%, -1.69% |
| rank_core4 (unmodified control) | 146.19 | 145.76 | -0.29% | +0.12%, -0.32%, -2.11%, -0.17% |
| Hybrid mixed_normalised (control) | 177.61 | 177.51 | -0.06% | -0.23%, +0.66%, -1.53%, +2.27% |

The separate Hybrid `mixed_full` control moves +0.87% and has higher
per-block variation. Broad mixed-case measurements remain host-specific,
and the exact distribution of affected rank-three inputs is not asserted
by timing alone.

## Correctness evidence and limitations

- `src/experimental/small/rank3_tests.rs` contains a frozen copy of the
  original four-pair evaluator and checks **binary64 output bits** across
  14 x 14 coordinate pairs for every rank-three small/core/upper-gap table.
  This includes both signed zeros and coordinate boundaries.
- CI passes unit/integration tests, Clippy `-D warnings` and rustfmt with
  `experimental`, `experimental,fma`, and combined feature matrices.
- The A/B benchmark's independent executable checksums agree across **all**
  tested workloads, not only the three targeted inputs.
- No exact root error or precision-budget improvement is claimed: the
  intended result is identical numerical output at lower cost.
- The source-hash-bound analytic component verifiers are **not** automatically
  re-certified merely because the implementation passes bitwise tests.
  They require source-bridge requalification; additionally, no full-domain,
  correct-rounding, or platform-independent performance theorem is claimed.
- Benchmarks currently run natively on x86_64 Linux only; the existing Apple
  M1 historical measurements are not comparable with these figures.

The temporary A/B workflow was used to obtain the linked artifacts; the
normal repository CI remains the merge gate.
