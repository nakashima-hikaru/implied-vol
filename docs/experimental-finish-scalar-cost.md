# Experimental common finish: scalar/twofold cost reduction

**Study date:** 2026-11-08 (Asia/Tokyo).
**Scope:** Existing `solver::Experimental`, small/wing/finite seed
routes, whose common finish evaluates a compensated normalized Black
price residual and a cubic inverse. **Not** a volfi replacement.
**Baseline:** `main@c6bb0c1313ac4dc3a23e18ac4752ac21fc87e46f`
(`src/experimental/small.rs` Git blob
`6a6c5e55719c64de2134d8803482a892c02eb638`).

## Mathematical and implementation rationale

Write the existing common finish as

```text
h = a/s0,       t = s0/2,
V = exp(-(h^2+t^2)/2)/sqrt(2*pi),
n = (b - B(a,s0))/(s0*V),
H2 = h^2-t^2,   Q = 3h^2+t^2,
delta = n*(1 + n*(-H2/2+n*(2*H2^2+Q)/6)),
s = fma(s0, delta, s0).
```

No terms of the cubic inverse are removed and the existing fitted seed,
pricing, compensated Mills polynomial, domain preflights, and fallbacks
are unchanged. Optimization acts on **how the double-double geometry
and normalization factors are represented**, not on the Black model.

### A. Multiply a twofold number by one scalar

The original `scale(a,b)` called `mul(a, Pair::from(b))`.
For a finite scalar `b`, its second double-double component is exactly
zero. The specialization computes

```text
hi = RN(a_hi*b)
e0 = fma(a_hi,b,-hi)
e1 = fma(a_lo,b,e0)
normalize(hi,e1)
```

instead of also evaluating zero low-factor products. The two *live*
FMA operations keep their source order and the normalization is the
same `add` routine. Exceptional infinities/NaNs are not within the
accepted normal-price contract.

### B. Normal binary64 scalar square

`t=Pair::from(s0/2)` always has low part zero, so `t*t` can
use the error-free `TwoProd(t_hi,t_hi)` high/low pair:

```text
high = RN(t_hi*t_hi)
low  = fma(t_hi,t_hi,-high).
```

For a **normal rounded product**, RN of the represented exact
product returns `high`; no extra twofold normalization is needed.
If the rounded square is zero/subnormal, the optimized source falls
back to the original `mul(t,t)`. This handles extreme underflow
without silently applying the normal-product theorem.

### C. Halve a normalized twofold square sum

The original `scale(add(hh,tt),-0.5)` is replaced in the common
normal range by an exact power-of-two scaling of the already-normalized
high and low parts. An explicit guard `high>=2^-500` sends
extreme near-zero inputs through the old generic scaling.

The high component remains normal under the guard. An unusually tiny
low component can still be subnormal or signed zero; a universal,
bit-for-bit identity for every possible binary64 pair, including the
sign of an unused low zero, is **not asserted**. The observed final
solver outputs are checked independently.

### D. Exact power-of-two denominator

Inside `forward_D`, `t.lo==0` by construction and `2*t.hi`
is a radix-two exact multiplication throughout the accepted finite
geometry. The old twofold scalar scaling of `t` is replaced with
`Pair::from(2*t.hi)`. The compensated Mills residual numerator and
the final pair division are unchanged.

### E. Exponential times a unit-high pair

When `ev=exp(exponent.hi)`, the original expression

```text
mul(Pair(ev,0), Pair(1,exponent.lo))
```

has an exactly-one high factor and a zero low component. The
specialization evaluates `normalize(ev,ev*exponent.lo)` instead,
and retains the subsequent multiplication by the existing
represented `1/sqrt(2*pi)` double-double constant. This removes
redundant arithmetic without introducing a different price model.

## Numerical regression

- A self-contained frozen copy of the previous common-finish
  expression compares **exact binary64 output bits** across
  2,048 synthetic seed/price combinations (accepted cases only)
  and deliberately chosen price/geometry seam inputs.
- The scalar square/scaling tests exercise **64,000**
  deterministic geometries, including a separate underflow
  boundary grid. The exponential unit-high multiplication has
  **20,000** deterministic exact-bit tests.
- Existing feature-matrix unit, integration, doctest, and
  independent high-precision normal-price fixtures remain in
  the repository and run in CI.
- A source bridge performs **exact literal reconstruction**
  of the entire original `small.rs`; it verifies the resulting
  Git blob hash equals the pinned main blob. The bridge only
  permits the explicitly enumerated arithmetic transformations.
- For every tested A/B workload (including broad mixed cases),
  baseline and optimized executable output **checksums agree**.
  The checksum agreement is a useful full-population regression,
  **not a proof of pointwise byte equality**.

The `src/experimental/small.rs` file changes, so the previous
source-hash-bound **continuous `rho_J<1` certificate is not
automatically re-qualified** by these finite tests or the source
bridge. In particular, the theorem requires accounting for
normality, underflow and vendor `exp` behavior in the new
source operation order. The working assumption remains round-to-nearest
binary64, explicit FMA, and gradual underflow; fast-math reassociation
is not included in the numerical contract.

## Timings: pinned source, same host, rotating ABBA

The two immutable executables were built with identical benchmark
sources on each GitHub Actions runner, rustc 1.99.1 / LLVM 23.1.1,
`RUSTFLAGS="-C target-cpu=native -C target-feature=+fma"`,
bench-profile LTO, with synthetic inputs and builders prepared outside
the timed region. Each comparison has four ABBA/BAAB blocks × two rounds
and 120,000 calls/round (16 measurements per version and case).

**Trial 1 — pinned main vs geometry/scalar specializations**,
AMD EPYC 7763:

| Workload | Main ns | Candidate ns | Change |
|---|---:|---:|---:|
| rank_small2 | 190.44 | 177.28 | -6.91% |
| rank_small3 | 200.04 | 184.96 | -7.54% |
| rank_core4 | 182.32 | 168.94 | -7.34% |
| finite_seed | 187.07 | 171.67 | -8.23% |
| wing_seed | 186.11 | 171.08 | -8.08% |
| mixed_normalised | 191.29 | 185.26 | -3.15% |
| mixed_full | 211.50 | 204.54 | -3.29% |
| central_deferred (control) | 79.07 | 79.31 | +0.31% |
| large_lower (control) | 223.44 | 224.41 | +0.43% |

**Independent trial 2**, AMD EPYC 9V45:

| Workload | Main ns | Candidate ns | Change |
|---|---:|---:|---:|
| rank_small2 | 98.45 | 93.89 | -4.63% |
| rank_small3 | 104.19 | 100.52 | -3.52% |
| finite_seed | 96.07 | 90.70 | -5.58% |
| mixed_normalised | 99.31 | 96.99 | -2.34% |
| mixed_full | 117.79 | 115.63 | -1.83% |
| large_lower (control) | 155.48 | 155.61 | +0.08% |

**Incremental trial 3 — exponential unit-high specialization only**,
pinned previous optimized candidate vs new exponential factor,
AMD EPYC 7763:

| Workload | Before ns | After ns | Change |
|---|---:|---:|---:|
| rank_small2 | 177.19 | 171.95 | -2.96% |
| rank_small3 | 184.93 | 181.14 | -2.05% |
| rank_core4 | 168.73 | 165.34 | -2.00% |
| finite_seed | 171.53 | 168.64 | -1.69% |
| mixed_normalised | 185.12 | 184.15 | -0.53% |
| mixed_full | 204.94 | 203.55 | -0.68% |

Do not calculate a synthetic total speedup by combining **absolute ns**
from different CI instances. All four paired blocks in trial 1
improve for rank/small/finite paths; the non-finish control cases are
approximately flat. The same pattern recurs on the different AMD host.

Raw evidence:
- [First pinned-main A/B on AMD EPYC 7763](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42443067198/artifacts/13601922971)
- [Independent rerun and residual histogram on AMD EPYC 9V45](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42443298542/artifacts/13602811467)
- [Incremental exponential-factor experiment](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42443516439/artifacts/13602003389)

**Interpretation:** these are changes to the existing solver, not a
new algorithm. The demonstrated hotspot improvement is material (typically
3–8% in the target rank/finite/wing workloads); mixed workload gains
are smaller because some inputs use other branches. These are
host/workload-specific timings, not portable speed guarantees.

## Why not remove the cubic correction?

A deterministic sample of 4,096 existing `mixed_normalised`
input prices reached `finish` on **1,927** inputs.
The observed counts by `|n|` interval were:

| Residual interval | Count |
|---|---:|
| <2^-32 | 13 |
| [2^-32,2^-28) | 81 |
| [2^-28,2^-24) | 256 |
| [2^-24,2^-20) | 271 |
| [2^-20,2^-16) | 1,182 |
| >=2^-16 | 124 |

For all accepted geometry with `0<=h<=8`,
`0<t<5.5`, `H2=h²-t²`, `Q=3h²+t²`,
the cubic coefficient is bounded by

```text
|c3| <= (2*(64)^2 + (3*64+5.5^2))/6 < 1403.
```

Thus at `|n|<2^-24`, `|c3*n^3|<1403*2^-72`,
well below binary64 relative precision, *before source rounding
and price conditioning allowances*. But this applies to only
**350/1927=18.2%** of recorded calls. Adding a branch to
omit the cubic term for a small minority is unlikely to pay
for its own control-flow overhead; no modified cubic code is adopted.

The counts are a sampled workload property, not a uniform
seed-residual bound or proof that the quadratic inverse would satisfy
the existing source precision certificate.

## Reproduction and remaining proof obligation

```sh
python3 scripts/verify_finish_scalar_bridge.py
cargo fmt --all --check
cargo clippy --features experimental -- -D warnings
cargo test --features experimental
# Performance: see benches/study_finish_cost.py and the pinned executable
# recipe in the A/B job artifact. Use identical flags/harness for both revisions.
```

The next research direction is to certify the new source-operation
rounding ledger against the archived continuous seed covers,
then investigate branch-specific geometry simplifications only where
the dynamic A/B experiments show positive return. A direct
seed-only acceptance path is not presently certified or benchmarked;
the existing cubic finish therefore remains intact.
