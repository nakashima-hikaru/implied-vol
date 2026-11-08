# Transfer of the Experimental common-finish rounding certificate

**Date:** 2026-11-08. **Scope:** Conditional `small::finish`
small/core/wing branches of the existing `Experimental` solver.
This is a source-level IEEE-754 analysis of the accepted branches,
**not** a new unconditional precision certificate for all possible
binary64 inputs or all underlying libm implementations.

## Baseline and source binding

- Reviewed baseline: `main@c6bb0c1313ac4dc3a23e18ac4752ac21fc87e46f`.
- Exact original `src/experimental/small.rs` blob:
  `6a6c5e55719c64de2134d8803482a892c02eb638`.
- Candidate: [PR #15](https://github.com/nakashima-hikaru/implied-vol/pull/15).
- `scripts/verify_finish_scalar_bridge.py` recognizes only five
  enumerated helper/finish changes, reverses all of them, and checks
  the restored **complete source** against that exact Git blob.
- `scripts/verify_finish_transfer.py` checks exact rational
  bounds for the imported proof domain, then exercises a deterministic
  independent IEEE-754 implementation with adverse negative controls
  and exact **bit-pattern** comparisons.

Reconstructing the old source is an exact source **binding**, not
by itself a proof that modified floating-point computations are equal.
The transfer argument therefore requires the following separate lemmas.

## Near-ATM scope exclusion

The original `small.rs::near_atm` route also called the
shared `scale` helper when computing `w=(s^2)/8`.
That route has **different** precision premises from
the common small/core/wing finish, so the generic
`scale` optimization cannot simply import the
common-finish bound there.

The candidate now explicitly invokes
`mul(mul(s,s),Pair::from(0.125))` for this particular
call, **exactly the expression the original
`scale(mul(s,s),0.125)` expanded to**. Thus the
separately reviewed near-ATM calculation retains
its previous arithmetic and no new claim about
that route's full domain is necessary.
The source bridge checks the exact inverse.

## 1. Inherited proof-domain hypotheses

The existing continuous seed/common-finish analysis qualifies
conditional accepted leaves, not all `(a,b)`. For those leaves,
the [inherited seed cover and finish reserves](seed-reduction.md)
include

```text
s0 > 2^-350,       0 <= h < 8,      0 < t < 5.5,
t = RN(s0/2),      h = a/s0 (with compensated low component),
0 < D < 2,        small/core/wing certified-domain guards,
correctly rounded binary64 basic arithmetic and FMA,
gradual underflow, fixed source operation order.
```

The old proof additionally tracks compensated Mills-model
error, `0.01 < D < 2` for transported leaves, and final
cubic remainder/rounding budgets. Those parts are **unchanged**.

Since `s0 > 2^-350` and dividing a normal binary64
by 2 is exact in this range, `t>2^-351`, hence

```text
t² > 2^-702 > 2^-900.
```

With `h<8` and `t<5.5`, the Gaussian exponent has
`-(h²+t²)/2 > -377/8 = -47.125`; the resulting positive
`exp` factor is greater than `2^-70` in real arithmetic
(the elementary bound `ln 2>0.69` suffices).
It is normal and nowhere near underflow.

These are deductions from the inherited *accepted leaf*
conditions; they do not assert that every input received by
the public API takes such a leaf.

## 2. Scalar double-double multiplication (live FMA order)

The generic `mul((a_hi,a_lo),(b,0))` calculates

```text
h = RN(a_hi*b)
l = FMA(a_hi,b,-h)
l = FMA(a_hi,0,l)
l = FMA(a_lo,b,l)
add((h,0),(RN(l+RN(a_lo*0)),0)).
```

The specialization omits only the zero multiplicands,
retaining the first and last **live** FMA and the same
normalizing `add` function. Under finite accepted
operands, the skipped products are exactly zero.
The signed-zero case is handled by the final
normalization: the old and new normalized pair components
are bitwise equal in the tested contract.

Outside the finite accepted domain, the source makes no
promise to preserve NaN signaling or infinity behavior
of the removed operations. The public solver's input and
geometry guards already exclude those cases.

## 3. Scalar square / TwoProd and the underflow counterexample

The old `mul((t,0),(t,0))` forms

```text
p = RN(t*t)
e = FMA(t,t,-p)
add((p,0),(e,0)).
```

The new `square_scalar(t)` returns the `(p,e)` pair
without the final normalizer **only if**
`p >= 2^-900` and `p` is finite. Otherwise it calls
the **unchanged generic `mul`**.

For a binary64 `t=M*2^k` with at most 53 significant
bits, the exact square has at most 106 significant bits.
`RN(t*t)>=2^-900` implies `t>2^-451` even allowing
for a product just below `2^-900` rounded upward,
and thus `k>=-503`. The exact product and the
difference `t²-p` lie on a binary grid no finer than
`2^-1006`, which is **above the normal minimum**
`2^-1022`. The FMA error is exactly representable,
and `|e|<=0.5 ulp(p)`. Therefore normalizing the
exact `(p,e)` pair cannot change `p`, and it
preserves `e` as well. The expression `e+0.0`
normalizes an exact `-0.0` in the same way as the
generic multiplication's final zero product.

**Why the guard is necessary:** An earlier implementation
checked only `p.is_normal()`. A specific positive
binary64 input

```text
t.to_bits() = 0x200ab8bf40df0b0a
```

has normal `p` but a subnormal residual. The generic
renormalization can move the high component by one ULP.
That prior shortcut is therefore **not equivalent**.
The permanent regression test checks this exact witness.
For the inherited `s0>2^-350` leaves, the new guard
always accepts because `t²>2^-702`.

## 4. Halving a normalized twofold sum

For `A=add(hh,tt)`, the previous source evaluated
`mul(A,(-0.5,0))`; the new source directly
halves both components only when `A.hi>=2^-500`.
The high-component scaling is then an exact radix-two
exponent shift and is normal. The inherited normalized
twofold relation `|A.lo|<=0.5 ulp(A.hi)` means the
normalizer cannot change the high component after
the exact shift. If the low component is subnormal,
both paths use the same binary64 rounding of its
half. The explicit `+0.0` on the low component
matches the generic add's canonical signed zero.

The negative-control input

```text
A.hi.to_bits() = 0x2505fec697039bba
A.lo.to_bits() = 0x0000000000000001
```

demonstrates why an unqualified `(-A.lo)/2`
return is not bitwise identical: the shortcut would
return `-0.0` where the old routine returns `+0.0`.
The new source includes the zero canonicalization
and uses generic scalar scaling below the guard.

This lemma relies on `A` being normalized, as
provided by the existing `add` routine and its
original accepted-domain rounding contract.

## 5. Doubled denominator and unit-high exponential factor

The `t=Pair::from(s0/2)` value has **zero low
component**. For the finite accepted `t`,
`2t` is a representable exponent shift, and the
old `scale(t,2)` returns `Pair::from(2*t)`
(modulo normalized zeros). The new denominator
has exactly the same bits.

For `ev=exp(exponent.hi)` (positive, normal),
`mul((ev,0),(1,exponent.lo))` contains a unit-high
factor. Its nonzero product is
`RN(ev*exponent.lo)`, and its subsequent `add`
is exactly the same as the new
`add((ev,0),(RN(ev*exponent.lo),0))`.
The resulting pair feeds the **unchanged** subsequent
`1/sqrt(2pi)` factor and Black-price residual.

No new approximation to `exp`, `log`, Vega or
the Black price has been introduced.

## 6. Transfer of the normalized Black error bound

Let `F_old(a,b)` and `F_new(a,b)` denote the
**complete binary64 solvers** on a common inherited
accepted leaf, with identical external primitive
results and rounding mode. Under the preceding
operation-by-operation hypotheses:

```text
same seed, guards and finite input
   -> same compensated h/t geometry
   -> same normalized h², t² and exponent pairs
   -> same exp factor, compensated Vega and D
   -> same n, H2, Q and cubic inverse
   -> F_new(a,b) == F_old(a,b) as a binary64 output.
```

Consequently, **conditional on these exact IEEE
equalities and the inherited reviewed certificate,**
the old upper bound on

```text
rho_J = |s_computed/s_true - 1| /
        (2^-52 * (1 + b/(s_true*V_true)))
```

transfers unchanged to the new implementation on those
leaves. No new inverse-series remainder is incurred:
the cubic calculation itself is byte-for-byte
unchanged.

**The distinction is essential:** the archived
continuous proof for the old solver was source-hash
bound and has not been re-executed in this workspace.
The finite helper checks and manual IEEE lemmas
are not a separately mechanized all-binary64 theorem.
We therefore record a **conditional transfer argument**,
not an unconditional global recertification.

## 7. Executable checks and replay requirements

The current PR runs, in its normal CI matrix for
`experimental` and combined features:

```sh
python3 scripts/verify_finish_scalar_bridge.py
python3 scripts/verify_finish_transfer.py --samples 80000
cargo fmt --all --check
cargo test --features experimental
```

The independent Python binary64 ledger compares
the old/modified helper pairs bit-for-bit on
80,000 deterministic samples each for scalar
products, halving and exponential factors, and
87,680 squares including guard neighbors.
It also asserts that **both old unsafe
shortcuts fail** on their constructed negative
controls. Rust tests additionally cover **300,000
deterministic normalized pairs** with biased
exponents from 1 through 1533 (unbiased -1022
through +510, chosen so tested squares are finite),
varying high/low components and signed zeros.
They compare the actual old and new f64 helper
outputs bit-for-bit and protect the near-ATM
operation-order handoff. The existing
independent high-precision root regressions and
normal/fma CI variants also pass.

These are finite validation results, not
exhaustive verification of all binary64 states.

For a **fresh replay of the original source-hash-
bound continuous certificate**, the required
external reference archive is not in this
GitHub repository. When present on a local
machine, execute:

```sh
python3 scripts/replay_finish_certificate.py \
  --reference-root /path/to/implied-black-volatility \
  --out target/proof/finish-source-transfer.json
```

The replay tool first runs the current
`verify_finish_transfer.py` source and f64
equivalence checks, then reconstructs the exact
reviewed `small.rs` blob, creates a temporary
canonical repository and runs the existing
`scripts/verify_direct_d.py` against its
**actual** archived dependencies and SIMD backend.
Both stages are mandatory for a successful result.
The tool rejects missing archives, source mismatches,
or a recomputed error bound `rho_J>=1`; it
cannot fabricate archive inputs, substitute
cached JSON bounds, or silently pass on finite
samples alone. It never modifies the active
optimized source.

**Outstanding for an unconditional release:**
supply the reviewed reference archive, run
the full conditional source-bound check,
and obtain an independent formal review of
the IEEE helper lemmas (or mechanized FP
equivalence in a suitable tool). Then
benchmark on Apple Silicon as well.

## 8. Performance after underflow/zero correction

The corrected common finish was benchmarked
against the pinned main executable on the same
GitHub Actions x86-64 host with rustc 1.99.1,
native/FMA flags, LTO and 4 alternating ABBA
blocks (two rounds, 120,000 calls per run).

| Workload | Main ns | Corrected ns | Change |
|---|---:|---:|---:|
| rank_small3 | 200.46 | 181.81 | **-9.30%** |
| rank_core4 | 182.51 | 165.05 | **-9.56%** |
| finite_seed | 187.00 | 169.31 | **-9.46%** |
| wing_seed | 186.07 | 167.96 | **-9.73%** |
| mixed_normalised | 191.54 | 183.95 | **-3.96%** |
| mixed_full | 211.24 | 204.40 | **-3.24%** |
| central_shape1 (control) | 84.48 | 84.50 | +0.02% |

The gains in the four compared blocks persist
after adding the correctness guards. They are
not universal CPU/runtime guarantees.
[Exact inputs, all block-level samples and
compiler/CPU hashes](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42445003369/artifacts/13602129801).

The original `Experimental` seed, cubic inverse,
application-facing API and fallback routes
remain unchanged.

### Independent final repeat after near-ATM isolation

The `near_atm` route uses the original generic
twofold multiplication rather than inheriting the
small/core/wing common-finish precision conditions.
This change has its own exact inverse source bridge
and 300,000 wide-exponent helper regressions. The
last same-executable, pinned-main scalar A/B on an
AMD EPYC Linux runner gives:

| Workload | Main ns | Candidate ns | Change |
|---|---:|---:|---:|
| rank_small3 | 113.22 | 109.85 | **-2.98%** |
| rank_core4 | 100.03 | 96.76 | **-3.27%** |
| finite_seed | 103.74 | 96.64 | **-6.84%** |
| wing_seed | 108.21 | 102.33 | **-5.44%** |
| mixed_normalised | 106.77 | 104.71 | **-1.93%** |
| mixed_full | 126.57 | 123.76 | **-2.22%** |
| large_lower (control) | 166.69 | 166.72 | +0.01% |

All four paired blocks improve for each of the
six listed target workloads. Changes are runner
and population specific; this table must not be
mixed with the preceding runner's absolute
nanoseconds to form a before/after improvement.

[Complete last-run raw JSON, paired blocks,
binary SHA-256 hashes and CPU/compiler record](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42447206905/artifacts/13603594714).


## Adoption decision — 2026-11-08

**Accepted for the opt-in `solver::Experimental` implementation**, not
for the default Hybrid solver and not as a newly published,
unconditionally certified algorithm. The objective is to reduce
unnecessary double-double operations *without changing* the
underlying approximation or numerical model.

The decision uses the following independent gates:

1. **Bit-level end-to-end comparator.** Two separate native/FMA
   release-mode executables, built from pinned
   `main@c6bb0c1313ac4dc3a23e18ac4752ac21fc87e46f`
   and the final candidate, ran the same deterministic test harness
   over 16 constructed regimes (ATM, normal-price interiors,
   log-price tails, upper cap neighbors, tiny a/b, and
   small/large branch thresholds), 8,192 cases per regime.
   The entire **131,072-line output files were byte-for-byte equal**,
   not just a checksum. Of these outputs, **107,712** were
   positive finite binary64 results, **22,781** were zero,
   and **579** were `None`. The output SHA-256 in both
   builds is `efc126994a49ee117542da4fa854441a58fbc1d7a3c5581260e08558627b1fb1`.
   [Complete output and provenance artifact](https://github.com/nakashima-hikaru/implied-vol/actions/runs/42452983641/artifacts/13605534124).

2. **Source and rounding domain.** The exact inverse
   transformation reconstructs the reviewed old source
   (Git blob `6a6c5e55719c64de2134d8803482a892c02eb638`).
   The FP transfer ledger accounts for normality, the
   `2^-900` exact-product threshold, the `2^-500` halving
   guard, the counterexamples and signed zero. The
   near-ATM route explicitly retains the original generic
   operation because it has a different proof domain.

3. **No loss of existing tests.** The full feature-matrix
   `experimental` unit/integration/independent-root tests,
   source binding, 300,000 wide-exponent comparisons and
   negative controls pass. The exact-output comparator is
   finite evidence, not a universal all-input proof.

4. **Measured return.** Pinned-main paired Linux AMD x86-64
   scalar benchmarks improve the broader
   `mixed_normalised` case by **1.93%** and `mixed_full`
   by **2.22%** on the final near-ATM-preserving version;
   `finite_seed` improves **6.84%** and `wing_seed`
   **5.44%**. The result is not claimed for every CPU,
   and user-level total-portfolio impact depends on the
   trade population.

**Limit of acceptance:** the independently reviewed external
`implied-black-volatility` proof archive was not accessible
for a fresh rerun; a fully mechanized all-binary64
equivalence theorem is also unavailable. Existing source-bound
`rho_J<1` conclusions are therefore inherited **conditionally**
on the original accepted-leaf contracts and the explicitly
documented floating-point equivalence lemmas. This adoption
does **not** assert a fresh unconditional certificate. The
fail-closed archive replay tool remains available for later
independent requalification. Because this is confined to
opt-in Experimental and the complete output comparison is
bit-identical on a broad finite corpus, the modest performance
improvement is judged worth adopting without changing the
default or claiming stronger numerical guarantees.

