# Performance

The stable benchmark measures seven fixed Black solver paths and two mixed
workloads. It prepares prices and builders before timing and reports
nanoseconds per calculation plus a checksum. The mixed workloads contain 4,096
seeded synthetic OTM calls; they are not a market distribution.
With `--solver experimental --experimental-paths`, it also measures eight
additional fixed inputs and seeded mixed large-moneyness and near-ATM inputs.

## Reproduce a comparison

```sh
RUSTFLAGS="-C target-cpu=native" \
  cargo bench --bench implied_black_vol_paths --features flashiv,experimental -- \
  200000 15 --solver hybrid
```

Repeat with `--solver jaeckel`, `--solver flashiv`, and `--solver experimental`.
The solver is selected once before the timed loops. `flashiv` and `experimental`
make their solver types available; enabling them does not select the default solver.
Hybrid and Jaeckel need no solver feature.

Use the same feature set, compiler flags, workload, and iteration count for each
comparison. Add `fma` consistently when measuring that
arithmetic policy. For paired measurements, rotate solver order and retain every
sample; finish builds and correctness checks before starting timing.

The reference measurements use native CPU optimization and LTO, with explicit
FMA only. They do not enable global fast-math or implicit FMA contraction.
Measurements made with `-C llvm-args=-fp-contract=fast` use a different arithmetic
policy and should be kept separate.

## Local comparison

The 2026-10-03 comparison used Apple M1 and Rust 1.98.1, warm single-threaded
calls, 200,000 calculations per case, and 15 samples per solver in five rotating
batches of three rounds. All four solver types use one executable with both
optional solvers enabled and the optional `fma` feature disabled. A separately
built immutable pre-change executable provides paired FlashIv and Hybrid
controls. The benchmark has no CPU affinity on macOS;
scheduler and frequency variation remain possible. These timings exclude input
construction and are specific to this workload and machine.

| Input | Hybrid | Jaeckel | FlashIv | Experimental |
|---|---:|---:|---:|---:|
| ATM | 7.8 ns | 7.9 ns | 164.5 ns | 10.6 ns |
| Lowest region | 193.4 ns | 262.9 ns | 147.2 ns | 174.9 ns |
| Lower middle | 110.5 ns | 110.4 ns | 158.6 ns | 129.0 ns |
| Upper middle | 111.5 ns | 111.2 ns | 160.6 ns | 130.5 ns |
| Highest region | 183.5 ns | 183.7 ns | 178.5 ns | 92.1 ns |
| Near ATM | 104.1 ns | 103.9 ns | 165.7 ns | 107.7 ns |
| Near ATM, wider | 95.4 ns | 95.5 ns | 164.6 ns | 145.0 ns |
| Mixed normalized | 181.4 ns | 198.4 ns | 178.4 ns | 129.4 ns |
| Mixed full API | 202.6 ns | 218.9 ns | 200.3 ns | 145.1 ns |

Experimental reduced mixed normalized/full calculation time by 28.7%/28.4%
against Hybrid in this measurement. Several fixed middle and near-ATM paths
were slower. FlashIv follows the paper's fixed-count Algorithm 1. It reduced
mixed normalized/full time from 326.8/354.7 ns to 178.4/200.3 ns against the
pre-change control, or 45.4%/43.5%. Hybrid control checksums matched on every
benchmark case, and its mixed timings changed by less than 0.1%.

The pure FlashIv ATM case now follows the paper's common iteration instead of
the previous closed-form ATM inverse, increasing that case from 12.1 to
164.5 ns. Ordinary middle-price cases also remain slower than LBR. The mixed
speedup is a workload result, not evidence that FlashIv is faster everywhere.

The same paired comparison with the optional `fma` feature enabled gave:

| Workload | FlashIv before | FlashIv after | Jaeckel | Hybrid |
|---|---:|---:|---:|---:|
| Mixed normalized | 289.5 ns | 161.3 ns | 172.7 ns | 157.6 ns |
| Mixed full API | 315.2 ns | 182.2 ns | 192.0 ns | 178.1 ns |

Hybrid control checksums also matched under this arithmetic policy.

An earlier solver API comparison measured Hybrid at 181.5/202.7 ns against the paired earlier default's
180.0/200.0 ns: 0.8%/1.3% more time on these mixed workloads. Checksums matched,
and separate tests verified unchanged default outputs under both FMA policies.
This is a local measurement of the API change, not a universal overhead bound.

These measurements do not cover cold caches, parallel batches, tail latency,
or a live market workload.
The bundled C++ comparison checks are separate from the four-solver benchmark.

## Experimental polynomial optimization

The 2026-10-03 optimization pairs independent Horner chains in the wing/upper
tensor seeds and the near-ATM conversion polynomial. On this Apple M1, those
chains use two NEON lanes. Each chain retains its coefficients, degree, and
explicit FMA sequence; the final scalar reductions are unchanged. The shared
upper route also reuses the cap exponential already computed by its dispatcher.

Separate immutable baseline/candidate executables used Rust 1.98.1, native CPU
optimization, LTO, and explicit FMA without implicit contraction. Both used the
same extended benchmark source. Four rotating ABBA blocks with two rounds per
invocation supplied 16 samples per version, solver, and feature policy, at
200,000 calls per case. The medians were:

| Experimental input | FMA off, before → after | FMA on, before → after |
|---|---:|---:|
| Near ATM, existing fixed case | 107.4 → 74.8 ns (-30.4%) | 107.3 → 74.7 ns (-30.4%) |
| Mixed near ATM | 122.2 → 93.9 ns (-23.2%) | 123.8 → 92.9 ns (-25.0%) |
| Microscopic input | 133.4 → 118.0 ns (-11.5%) | 133.1 → 118.0 ns (-11.4%) |
| Wing tensor seed | 127.1 → 126.0 ns (-0.9%) | 127.1 → 125.8 ns (-1.0%) |
| Existing mixed normalized | 129.5 → 129.2 ns (-0.2%) | 129.9 → 129.1 ns (-0.6%) |
| Existing mixed full API | 145.1 → 145.3 ns (+0.1%) | 145.3 → 144.9 ns (-0.3%) |
| Mixed large moneyness | 218.3 → 218.1 ns (-0.1%) | 218.2 → 218.0 ns (-0.1%) |

The new near-ATM workload has 4,096 seeded OTM prices, with
`abs(x)` uniform in `[1e-6, 0.25)` and generating total volatility uniform in
`[0.1, 0.5)`. The large workload samples `a` in `[10, 1416)` and log relative
price parameter in `(-300, -0.01]`, retaining normal prices. These are synthetic
distributions. Existing broad mixed inputs improved little; the benefit is
concentrated in the near-ATM conversion path. Hybrid control mixed medians
changed by at most 0.2%, and every timed checksum matched. Individual fixed
cases also showed compiler/layout sensitivity across feature builds; the
complete samples are retained in the [measurement record](experimental-speed-2026-10-03.json).
No performance claim is made for other architectures or market distributions.

Moving the large-domain derivative divisions after convergence regressed its
mixed workload by 2.3%. Inspection found that the compiler stopped inlining
the paired Mills evaluator. Forcing that inline, or returning a reduced pair,
also failed to improve the large mixed workload in subsequent screening.
These candidates were rejected; the large-domain implementation is unchanged.

Final native release revalidation covered 117,044 independent inputs per FMA
policy: 46,688 archived core roots, 2,083 additional tiny/ATM/seam/cap roots,
and 68,273 paper inputs. All public outputs matched the baseline bit for bit,
with zero invalid outputs or `rho_J>=1`. Maximum core-coordinate rho was
0.6049, and maximum paper rho was 0.7733 (0.8690 after reannualization).
All these adapted prices are positive normal values. The finite validation
does not expand the mathematical or platform assumptions in
[the numerical contract](numerics.md).

The paired runner retains executable SHA-256 hashes, all samples, and checksum
checks. It records caller-supplied build flags/features and the current host
compiler; it does not infer build provenance from binaries. The saved record
also includes the source/build hashes verified for this comparison.
Build the baseline and candidate in separate target directories, copying the
same benchmark source into both checkouts before building, then run:

```sh
python3 benches/compare_experimental.py /path/to/baseline /path/to/candidate \
  --output comparison.json
```

Repeat with identically built FMA executables and
`--features flashiv,experimental,fma`. The existing no-affinity, warm-cache,
single-threaded measurement limitations apply.

## Actions C++ comparison

The [2026-10-03 Actions benchmark](https://github.com/nakashima-hikaru/implied-vol/actions/runs/37083979975/job/111090342937)
compares the default Hybrid with the bundled C++ LBR. It enables `cxx_bench`,
and optionally `fma`; it does not select pure FlashIv. That run predates the
near-ATM stability change below.

| Black implied-volatility case | Rust, FMA off | C++ | Rust, FMA on | C++ |
|---|---:|---:|---:|---:|
| Legacy OTM call, input preparation included | 210.65 ns | 181.78 ns | 192.23 ns | 181.65 ns |
| Deep OTM call, calculation only | 196.64 ns | 180.16 ns | 179.98 ns | 180.43 ns |
| Near-ATM short call, calculation only | 238.93 ns | 198.22 ns | 204.55 ns | 200.83 ns |

The C++ build uses `-Ofast` and `-ffp-contract=fast` in both runs. The Rust
benchmark command overrides the repository's configured compiler flags and
permits FMA only at explicit `mul_add` sites when the `fma` feature is enabled.
This differs from v2.0's Actions command, which also allowed implicit FMA
contraction in Rust. Consequently the FMA-off row compares different arithmetic
policies; even the FMA-on comparison retains C++ fast-math differences.

The legacy Rust benchmark now includes input preparation in its timed loop
and prevents constant propagation of inputs and outputs. Its v2.0 version
prepared the builder outside that loop. Use the split `calculate_only` and
`cpp_direct` measurements to distinguish calculation cost from API preparation.
The timing region and arithmetic policy changes prevent a direct regression
claim from older Actions numbers alone. Hosted-runner variation remains another
limit of these measurements.

## Near-ATM stability change

The 2026-10-03 comparison used the same Apple M1 and Rust 1.98.1, native CPU
optimization, LTO, and separate immutable release/candidate executables. Each
case used 300,000 calls and 16 samples per version, solver, and FMA policy,
with alternating paired order. Medians on the existing mixed workloads were:

| Explicit FMA | Solver | Normalized, before → after | Full API, before → after |
|---|---|---:|---:|
| Off | Hybrid | 181.8 → 181.4 ns (-0.2%) | 202.8 → 202.7 ns (0.0%) |
| Off | Jaeckel | 198.0 → 198.5 ns (+0.2%) | 219.0 → 219.2 ns (+0.1%) |
| On | Hybrid | 154.7 → 157.3 ns (+1.7%) | 174.5 → 178.3 ns (+2.2%) |
| On | Jaeckel | 171.7 → 173.3 ns (+1.0%) | 189.8 → 192.6 ns (+1.5%) |

Separating the near-ATM lowest-price correction from the ordinary kernel kept
Jaeckel's fixed lowest-region cost at 261.8 → 263.0 ns without FMA and
227.4 → 227.7 ns with FMA. All benchmark checksums matched the release baseline.
These timings measure ordinary inputs outside the new numerical path, so they
qualify its dispatch overhead rather than the cost of solving microscopic
inputs. Small increases remain in the FMA build; this is a numerical repair,
not a demonstrated speedup. The workload and machine limitations above apply.
