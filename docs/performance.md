# Performance

The stable benchmark measures seven fixed Black solver paths, two mixed
workloads, and three full-API inputs from the Actions benchmarks. It prepares
prices and builders before timing and reports
nanoseconds per calculation plus a checksum. The mixed workloads contain 4,096
seeded synthetic OTM calls; they are not a market distribution.
With `--solver experimental --experimental-paths`, it also measures eight
additional fixed inputs and seeded mixed large-moneyness and near-ATM inputs.

## Reproduce a comparison

```sh
RUSTFLAGS="-C target-cpu=native" \
  cargo bench --bench implied_black_vol_paths --features flashiv,experimental -- \
  200000 3 --solver hybrid
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

The 2026-10-04 comparison measures the current implementation at `19d5eab`,
including the Experimental polynomial degree reductions and Hybrid low-price
dispatch optimization below. It used Apple M1, macOS 27.0.1, Rust 1.98.1,
`-C target-cpu=native`, the bench profile with LTO, and
`flashiv,experimental`, with the optional `fma` feature disabled. All four
solvers use the same immutable executable and the same twelve workloads.
Prices and builders are constructed outside timing.

Each case has 24 samples per solver: 200,000 calculations per round, three
rounds per invocation, and eight blocks. Four cyclic solver rotations are
followed by their reversed rotations, so each solver occupies each position
twice. All builds finished before timing, and no proof or reference jobs ran
during measurement. The benchmark has no CPU affinity on macOS; scheduler
and frequency variation remain possible.

| Input | Hybrid | Jaeckel | FlashIv | Experimental |
|---|---:|---:|---:|---:|
| ATM | 7.6 ns | 7.6 ns | 162.7 ns | 10.5 ns |
| Lowest region | 182.6 ns | 263.8 ns | 147.2 ns | 162.3 ns |
| Lower middle | 111.6 ns | 109.8 ns | 159.5 ns | 125.6 ns |
| Upper middle | 112.9 ns | 111.0 ns | 161.0 ns | 131.0 ns |
| Highest region | 188.6 ns | 184.6 ns | 175.9 ns | 91.3 ns |
| Near ATM | 105.3 ns | 103.0 ns | 166.0 ns | 75.3 ns |
| Near ATM, wider | 96.0 ns | 94.2 ns | 165.0 ns | 128.7 ns |
| Mixed normalized | 178.4 ns | 199.2 ns | 179.1 ns | 128.1 ns |
| Mixed full API | 199.3 ns | 220.8 ns | 201.2 ns | 144.2 ns |
| Deep OTM, full API | 129.1 ns | 128.6 ns | 187.3 ns | 127.2 ns |
| Near-ATM short expiry, full API | 171.0 ns | 171.7 ns | 191.4 ns | 93.2 ns |
| Legacy OTM, full API | 134.6 ns | 136.5 ns | 191.9 ns | 150.3 ns |

Experimental used 28.2%/27.6% less time than Hybrid on the mixed
normalized/full workloads. The middle-price, wider near-ATM, and legacy full
cases favor Hybrid or Jaeckel. These medians describe the current solvers on
this workload; differences from older runs are not paired optimization results.
The separate before/after measurements below qualify each accepted change.

All within-solver checksums were stable, and the executable and numerical-source
hashes were unchanged after measurement. Cross-solver checksums may differ.
The finite-positive-output preflight is separate from the independent-root
accuracy checks described in [the numerical notes](numerics.md#accuracy-contracts).
Full samples, invocation order, block medians, build/source/executable hashes,
and the measurement runner are retained in the
[current comparison record](local-solvers-2026-10-04.json).

To reproduce this sampling schedule, use the command above with three rounds
and run the four solver selections in each of the eight recorded orders.
These are warm single-threaded synthetic inputs; the comparison does not cover
cold caches, parallel batches, tail latency, or native Windows performance.

## Earlier FlashIV and solver comparison (2026-10-03)

The 2026-10-03 comparison used Apple M1 and Rust 1.98.1, warm single-threaded
calls, 200,000 calculations per case, and 15 samples per solver in five rotating
batches of three rounds. All four solver types use one executable with both
optional solvers enabled and the optional `fma` feature disabled. A separately
built immutable pre-change executable provides paired FlashIv and Hybrid
controls. The benchmark has no CPU affinity on macOS;
scheduler and frequency variation remain possible. These timings exclude input
construction and are specific to this workload and machine.
These historical values precede the subsequent Experimental and Hybrid
optimizations documented below; the current four-solver comparison is above.

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
The C++ comparison added below uses the same benchmark adapter and inputs.

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
These candidates were rejected; the polynomial optimization retained the
large-domain implementation.

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

## Experimental dispatch and logarithm optimization

The subsequent 2026-10-03 comparison starts from `84a7f1f`, which already
contains the polynomial optimization above. It adopts three changes:

- A conservative dyadic lower bound avoids the cap exponential when the upper
  route cannot accept the price. In an untimed pass over the existing mixed
  workload, this skips 747 exponentials per 4,096 inputs.
- The deferred route rejects guaranteed failures before evaluating `sinhc`.
  Keeping these checks inside the route limits their cost to its callers.
  The [numerical notes](numerics.md) give the rejection margins.
- The initial large-domain logarithm uses the validated positive normal price
  to omit generic special-case checks and a zero-correction division. Both
  compensated components and all subsequent arithmetic are preserved.

The same separate immutable native/LTO builds and four rotating ABBA blocks
supplied 16 samples per version at 200,000 calls per case.
The feature columns toggle the crate's `fma` feature. Experimental uses explicit
FMA in both, while generated prices can differ between feature builds.
Medians on Apple M1 with Rust 1.98.1, explicit FMA, and no implicit contraction
were:

| Experimental input | FMA off, before → after | FMA on, before → after |
|---|---:|---:|
| Existing mixed normalized | 129.5 → 127.8 ns (-1.3%) | 129.2 → 127.7 ns (-1.2%) |
| Existing mixed full API | 145.4 → 143.5 ns (-1.3%) | 145.3 → 143.8 ns (-1.0%) |
| Mixed near ATM | 94.2 → 89.6 ns (-4.9%) | 93.2 → 90.0 ns (-3.4%) |
| Mixed large moneyness | 218.5 → 215.2 ns (-1.5%) | 218.5 → 215.4 ns (-1.4%) |
| Near ATM, wider fixed case | 145.1 → 128.4 ns (-11.5%) | 144.9 → 128.4 ns (-11.4%) |
| Finite seed | 123.4 → 118.4 ns (-4.0%) | 123.3 → 118.4 ns (-3.9%) |
| Near ATM, successful fixed case | 74.9 → 75.5 ns (+0.9%) | 74.9 → 75.7 ns (+1.1%) |

Successful fixed deferred cases pay about 0.6–0.9 ns for the new checks;
their mixed workload benefits from cheaper failed attempts. Hybrid mixed
control medians changed by at most 0.3%, and every timed checksum matched.
The workload definitions and local, warm-cache, single-threaded limitations
above apply. These synthetic inputs do not establish a market-workload or
other-platform speedup.

Degree-interleaved rank chains, paired quadrature atoms, and bit-based normal
`frexp` did not establish a broad mixed-input improvement in screening. An
additional deferred lower-price rejection improved the microscopic fixed case
by 5.8%, but regressed the wing fixed case by 8.1%; it was rejected.

Both feature policies again passed all 117,044 independent references with
zero bitwise differences, invalid outputs, or `rho_J>=1`. A separate unfiltered
660-input replay checked guard thresholds, dyadic integer jumps, logarithm
mantissa/exponent seams, and subnormal continuation. It preserved all six
output lanes and classifications, including five existing invalid cases per
lane. This latter replay checks parity and makes no new root-accuracy claim.
The two added regression tests check threshold neighbors and both compensated
log components across all 2,046 normal exponents.

Source/build hashes, all final timing samples, screening evidence, guard
derivations, and reference provenance are saved in the
[dispatch measurement record](experimental-dispatch-speed-2026-10-03.json).

## Experimental AS1 seed refinement

The additional AS1 inverse-logarithm seed term improves precision margin while
adding arithmetic. On the local Apple Silicon host, separate baseline/candidate
release builds used Rust 1.98.1, `-C target-cpu=native`, LTO, and
`experimental,cxx_bench`, with and without the `fma` feature. Each case has
24 samples per version: 200,000 calls, two rounds, six alternating ABBA blocks.
Other validation jobs were paused during measurement. The FMA columns select
the crate feature; AS1 uses explicit FMA in both columns.

| Input | FMA off, before → after | FMA on, before → after |
|---|---:|---:|
| AS1 central | 80.40 → 81.98 ns (+1.96%) | 82.96 → 84.57 ns (+1.94%) |
| AS1 high-a guard | 80.33 → 82.03 ns (+2.11%) | 82.99 → 84.62 ns (+1.97%) |
| AS1 deep tail | 81.84 → 83.39 ns (+1.90%) | 84.94 → 84.90 ns (-0.05%) |
| AS1 tiny scaling | 101.28 → 104.44 ns (+3.12%) | 104.60 → 107.60 ns (+2.87%) |
| AS1 mixed, 4,096 inputs | 81.62 → 83.02 ns (+1.71%) | 84.21 → 85.54 ns (+1.58%) |
| Ordinary mixed normalized | 128.09 → 128.16 ns (+0.06%) | 131.75 → 131.95 ns (+0.16%) |
| Ordinary mixed full | 144.40 → 143.64 ns (-0.52%) | 148.03 → 148.19 ns (+0.11%) |

The AS1 mixed inputs were traced through the actual Rust dispatcher: all 4,096
use AS1 in both versions. C++ timing controls and per-block variation are retained
in the [measurement record](experimental-as1-2026-10-03.json). C++ retains its
existing fast-math flags and serves as a timing control, not an accuracy oracle.
Small changes in the ordinary mixed rows and the FMA-on deep-tail row are within
the observed variation. These local measurements do not establish Windows or
universal performance.

Run the additional fixed and mixed cases with:

```sh
RUSTFLAGS="-C target-cpu=native" \
  cargo bench --bench implied_black_vol_paths --features experimental,cxx_bench -- \
  200000 2 --solver experimental --as1-paths
```

Repeat with `--features experimental,cxx_bench,fma` for the other feature policy.
Input construction is outside the timed loops, and the existing default workload
is unchanged. Accuracy checks permit changed outputs and compare them against
independent roots; checksum stability is checked within each binary.

## Experimental polynomial degrees under the precision contract

The 2026-10-04 comparison starts from `ca050c6`, including the AS1 seed
refinement above. It accepts the degree-five AS1 atan approximation (A1) and
the cell-dependent degree of the intermediate-`t` Mills divided difference
(A4). Acceptance requires `rho_J<1`, rather than matching the old output bits.
The [numerical discussion](numerics.md#accuracy-contracts) distinguishes the
conditional component bounds from finite independent-root checks.

Separate frozen baseline/candidate builds used Apple M1, macOS 27.0.1,
Rust 1.98.1, `-C target-cpu=native`, the bench profile with LTO, and
`flashiv,experimental`, with and without `fma`. C++ was not enabled in these
timings. The identical harness measured 29 Experimental cases and 12 Hybrid
controls, using 200,000 calls per round, two rounds and six alternating ABBA
blocks: 24 samples per version per case. Other build, proof, and reference
jobs were stopped during timing. The columns toggle the crate feature;
Experimental's affected arithmetic uses explicit FMA in both.

| Input | FMA off, before → after | FMA on, before → after |
|---|---:|---:|
| AS1 mixed, 4,096 inputs | 82.83 → 80.72 ns (-2.55%) | 82.87 → 80.75 ns (-2.55%) |
| Divided difference, fixed | 170.62 → 157.96 ns (-7.42%) | 170.76 → 157.84 ns (-7.56%) |
| Divided difference, 160 inputs | 169.49 → 158.85 ns (-6.27%) | 169.57 → 159.00 ns (-6.24%) |
| Lowest branch | 171.08 → 161.56 ns (-5.56%) | 170.92 → 161.57 ns (-5.47%) |
| Ordinary mixed normalized | 127.55 → 127.69 ns (+0.11%) | 127.75 → 127.84 ns (+0.07%) |
| Ordinary mixed full | 143.53 → 143.56 ns (+0.03%) | 144.12 → 144.20 ns (+0.05%) |
| Deep OTM full | 126.21 → 126.67 ns (+0.37%) | 126.18 → 126.71 ns (+0.41%) |

All six paired blocks improved for the AS1 mixed and divided-difference
workloads. Instrumented copies of the final combined source confirmed that all
4,096 AS1 mixed inputs use AS1 and all 161 divided-difference inputs (one fixed
plus 160 mixed) reach the shortened polynomial in the accepted final correction
under both feature policies. The ordinary mixed workloads were essentially
unchanged. Deep OTM full consistently cost about 0.5 ns more in both policies;
this trade-off is retained alongside the targeted gains. Hybrid controls remained within 0.13%
of their baseline medians. These warm single-threaded synthetic workloads do
not establish native Windows or universal performance.

The dynamic degree loop was faster than a separately qualified constant-degree
dispatch: the latter improved divided-difference mixed timing by only about
3.0–3.2%, versus 6.2% for the dynamic loop in the screening run. The proposed
large-route reduction in precise evaluations (A5) was also profiled: all 47,091
archived inputs with `a>10` already performed exactly one precise evaluation,
so that path was left unchanged. This observation applies to that population.

The [validation and measurement record](experimental-rho-2026-10-04.json)
retains per-block comparisons, source and executable hashes, candidate screens,
input classifications and proof assumptions. To repeat the comparison, build
both versions in separate target directories with this same benchmark harness,
then compare the frozen executables:

```sh
python3 benches/compare_experimental.py /path/to/baseline /path/to/candidate \
  --iterations 200000 --rounds 2 --blocks 6 \
  --experimental-paths --as1-paths --dpoly-paths \
  --features flashiv,experimental --flags="-C target-cpu=native" \
  --output /path/to/comparison.json
```

Repeat with binaries built using `flashiv,experimental,fma` and update
`--features` accordingly. These provenance arguments describe the build; they
do not recompile or change the supplied binaries. The additional AS1 and
divided-difference populations are opt-in and are constructed outside timing.

## Hybrid low-price dispatch optimization

The 2026-10-03 change starts the existing restricted FlashIV path before
constructing LBR's interpolation nodes for `abs(x)>=0.01`,
`b<=0.0005*b_max`, and a normal cap. The [numerical notes](numerics.md)
derive why this guard is strictly inside the existing lowest-price region.
If the attempt declines, the existing LBR path runs once. The special-function
arithmetic, iteration counts, and other numerical methods are unchanged.

Separate immutable baseline/candidate executables used Apple M1, Rust 1.98.1,
native optimization, LTO, and matching SDK 27.0. Each feature policy enabled
only `cxx_bench`, optionally with `fma`. Four alternating ABBA blocks with
two rounds per invocation supplied 16 samples per version, solver, and case
at 200,000 calls per sample. The same twelve-case benchmark source was copied
into both snapshots. Medians in ns/call were:

| Hybrid input | FMA off, before → after | FMA on, before → after |
|---|---:|---:|
| Lowest region | 193.15 → 184.32 (-4.58%) | 162.10 → 150.31 (-7.28%) |
| Mixed normalized | 181.29 → 177.08 (-2.32%) | 157.45 → 154.39 (-1.94%) |
| Mixed full API | 202.84 → 197.47 (-2.65%) | 177.84 → 174.50 (-1.88%) |
| Lower middle | 110.39 → 113.01 (+2.38%) | 89.99 → 92.03 (+2.27%) |
| Deep OTM, full API | 128.71 → 128.09 (-0.48%) | 109.07 → 110.58 (+1.38%) |
| Near-ATM short expiry, full API | 198.26 → 170.06 (-14.22%) | 150.80 → 152.45 (+1.09%) |

Some fixed ordinary paths regressed by up to 2.93% without FMA and 2.27%
with FMA. The deep-OTM and short near-ATM inputs do not enter the shortcut;
their changes reflect dispatch cost and generated-code layout. In particular,
the short-case FMA-off improvement is not a faster FlashIV calculation.
Jaeckel control mixed medians were within 0.1% without FMA but increased
about 2% with FMA, demonstrating Rust layout sensitivity. C++ controls
changed by at most 0.26%. Every timed checksum matched its baseline.

In the final FMA-on build, Hybrid mixed normalized/full medians were
154.39/174.50 ns against C++'s 170.29/187.81 ns. The short near-ATM full
case still took 152.45 ns against C++'s 141.16 ns, or 8.00% more time.
The C++ fast-math policy described below remains different from Rust's.
These local synthetic, warm-cache, single-threaded results do not establish
a universal speedup or a change to the Linux Actions result.

Moving the ATM branch before its cap exponential improved the fixed ATM
case, but the combined candidate increased FMA-on short-case time by 6.81%.
It was rejected. Normalization sorting, cap reuse, and outlining the middle
Black evaluator also failed to establish a suitable broad improvement;
the outline screen additionally had possible concurrent-build interference
and is excluded from adoption evidence.

Final revalidation passed 100 library/integration tests and six doctests per
FMA policy, plus formatting, Clippy for the library, stable benchmark and
changed integration targets, and 14 comparison-report tests. Independent
100/180-digit roots cover 18 neighboring dispatch inputs. An unfiltered
baseline/candidate replay covered 477,112 normalized/full call/put inputs
and 954,224 Hybrid/Jaeckel outputs per policy with zero bit differences.
It retains existing extreme-input failures and NaNs; this is preservation
evidence rather than a new accuracy guarantee. The
[measurement record](hybrid-speed-2026-10-03.json) retains all final samples,
source/build/executable hashes, runner and replay sources, classifications,
and rejected-candidate qualification.

## Actions C++ comparison

The [latest 2026-10-03 Actions benchmark](https://github.com/nakashima-hikaru/implied-vol/actions/runs/37110794388/job/111168142389)
compares the default Hybrid with the bundled C++ LBR. It enables `cxx_bench`,
and optionally `fma`. Hybrid retains LBR in the middle and upper regions and
uses the restricted FlashIV method in selected lowest-region cases. Improving
that lowest region does not make every LBR path faster.

| Black implied-volatility case | Rust, FMA off | C++ | Rust, FMA on | C++ |
|---|---:|---:|---:|---:|
| Legacy OTM call, input preparation included | 189.24 ns | 165.95 ns | 167.94 ns | 165.24 ns |
| Deep OTM call, calculation only | 181.69 ns | 165.37 ns | 170.05 ns | 166.02 ns |
| Near-ATM short call, calculation only | 224.09 ns | 207.46 ns | 198.68 ns | 207.18 ns |
| ATM call, calculation only | 21.21 ns | 23.53 ns | 19.29 ns | 23.60 ns |

The FMA-on Hybrid is 2.4% slower in the deep OTM calculation and 4.1% faster
in the near-ATM calculation. The FMA-off near-ATM C++ sample reported
`+/- 84.11 ns`, so that particular difference is noisy. The default solver
still has slower paths; these results do not show that Rust wins everywhere.

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
Legacy C++ loops black-box only the result, whereas the split C++ loops also
black-box the inputs.
The timing region and arithmetic policy changes prevent a direct regression
claim from older Actions numbers alone. Hosted-runner variation remains another
limit of these measurements.

The added **Compare Hybrid with LBR** jobs explicitly select Hybrid, Jaeckel,
and C++ in the same compiled benchmark, using the same
fixed and 4,096-input seeded workloads. Prices and builders are prepared before
timing. Every solver's outputs are checked for finite positive values before
timing; the report rotates solver order, checks within-solver checksum stability,
and saves all samples and the executable hash. Cross-solver checksums may differ
because the solvers have different numerical contracts. The C++ arithmetic
policy above remains in effect. These jobs supplement the existing benchmark
logs with a median comparison table in the Actions summary and downloadable
measurement artifacts; the new jobs have not yet run remotely.
The full-API rows include the deep-OTM, near-ATM short-expiry, and legacy OTM
inputs from the existing Actions benchmarks. The legacy row prepares its input
before timing here, so it measures calculation cost rather than the legacy
job's preparation-plus-calculation loop.

To run the same comparison locally with the crate's `fma` feature disabled:

```sh
RUSTFLAGS="-C target-cpu=native" cargo bench --bench implied_black_vol_paths \
  --features cxx_bench --no-run
python3 benches/compare_solvers.py <benchmark-executable> \
  --features cxx_bench --flags "-C target-cpu=native" \
  --output solvers.json --summary solvers.md
```

Use the executable path printed by Cargo. Add `fma` to the feature list for the
explicit-FMA comparison in shared kernels.

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
