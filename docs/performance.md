# Performance

The stable benchmark measures seven fixed Black solver paths and two mixed
workloads. It prepares prices and builders before timing and reports
nanoseconds per calculation plus a checksum. The mixed workloads contain 4,096
seeded synthetic OTM calls; they are not a market distribution.

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

Apple M1, warm single-threaded calls, 200,000 calculations per case, 15 samples per
solver, five rotating batches of three rounds. All four solver types use one
executable with both optional solvers enabled. An immutable earlier default
executable provides a paired control. The benchmark has no CPU affinity on macOS;
scheduler and frequency variation remain possible. These timings exclude input
construction and are specific to this workload and machine.

| Input | Hybrid | Jaeckel | FlashIv | Experimental |
|---|---:|---:|---:|---:|
| ATM | 7.9 ns | 7.9 ns | 12.1 ns | 10.6 ns |
| Lowest region | 190.2 ns | 261.9 ns | 340.0 ns | 174.8 ns |
| Lower middle | 109.7 ns | 109.6 ns | 346.7 ns | 129.0 ns |
| Upper middle | 111.2 ns | 111.3 ns | 289.8 ns | 130.5 ns |
| Highest region | 183.5 ns | 183.5 ns | 258.3 ns | 92.3 ns |
| Near ATM | 103.3 ns | 103.8 ns | 349.1 ns | 107.6 ns |
| Near ATM, wider | 94.5 ns | 94.4 ns | 282.1 ns | 145.0 ns |
| Mixed normalized | 181.5 ns | 198.0 ns | 326.3 ns | 129.4 ns |
| Mixed full API | 202.7 ns | 218.9 ns | 354.8 ns | 145.0 ns |

Experimental reduced mixed normalized/full calculation time by 28.7%/28.5%
against Hybrid in this measurement. Several fixed middle and near-ATM paths
were slower. FlashIv denotes this crate's safeguarded numerical variant, which
differs from the author's fixed-count implementation.

The current Hybrid measured 181.5/202.7 ns against the paired earlier default's
180.0/200.0 ns: 0.8%/1.3% more time on these mixed workloads. Checksums matched,
and separate tests verified unchanged default outputs under both FMA policies.
This is a local measurement of the API change, not a universal overhead bound.

[Raw measurements](../performance/solver-api-2026-10-03.json) contain all samples,
checksums, compiler flags, source/executable hashes, and the retained earlier
layout trials. Solver selection occurs outside every timed loop.

## Other recorded measurements

- [Earlier solver comparison](../performance/solver-modes-2026-10-03.json) uses
  implicit FMA contraction and a different selection API.
- [Earlier research-solver comparison](../performance/experimental-2026-10-03.json)
  measures build-selected solvers before the explicit API.
- [Restricted hybrid qualification](../performance/flashiv-2026-10-02.json) records
  the lower-price route and its FMA checks.
- [Earlier baseline and exponential reuse](../performance/2026-10-02.json) records
  the earlier baseline, isolated candidates, and their retained samples.

These files preserve distinct source and arithmetic configurations. They do not
measure cold caches, parallel batches, tail latency, or a live market workload.
The bundled C++ comparison checks are separate from the four-solver benchmark.
