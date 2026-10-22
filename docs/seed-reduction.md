# Reduced rank seeds and the cubic common finish

The small/core rank factors use per-rectangle Chebyshev degrees. The 26 tables
occupy **11,328 bytes**, down from 14,400 bytes with uniform degree 8/6 and
17,600 bytes in the original archive. The explicit degree arrays contain
another 52 bytes. This counts coefficient and degree payload, not process
memory or binary size. The tables retain the existing padded rank layout.

| x/y degrees | Small rectangles | Core rectangles | Table count |
|---|---|---|---:|
| 6/5 | 0, 2, 3, 4, 5 | 4, 9, 14, 16, 17, 18, 19 | 12 |
| 6/6 | — | 0, 1, 2, 5, 6, 7, 8, 10, 11, 12, 13, 15 | 12 |
| 7/5 | 1 | — | 1 |
| 7/6 | — | 3 | 1 |

Each rank-three/four evaluation uses 22, 24 or 26 two-lane Horner stages,
instead of 28 for uniform degree 8/6; rank-two evaluation uses 11 instead of
14. These operation counts do not establish a runtime speedup. The compiler
sees fixed degree and array-size specializations selected by generated
metadata. Upper-gap factors retain degree 8/8. Lane pairing, scalar rank
reduction, price residual and the cubic common finish remain unchanged.

The generator works from the exact binary64 coefficients in the hash-bound
reference archive. It converts each power polynomial to Chebyshev form using
rational arithmetic, removes the high modes, converts back exactly, and rounds
each resulting coefficient once to binary64. Directly deleting high power
coefficients would give a different and much larger error.

```sh
python3 scripts/generate_seed_reduction.py --reference-root /path/to/implied-black-volatility --check
python3 scripts/verify_direct_d.py --reference-root /path/to/implied-black-volatility
```

`--check` regenerates and compares every coefficient and degree entry. The integrated checker
recomputes the direct-D cell bounds and the reduced-seed proof. It does not use
a previously generated JSON report as a premise. Its narrow source bridge
reconstructs the reviewed baseline and enforces the original whole-file hashes
before importing any unchanged kernel certificate. It binds the exact degree
dispatch, array lengths, Horner loops and rectangle-to-degree selection. A
changed selector cannot be admitted by updating a whole-file hash. The same
bridge preserves the unchanged AS1/FMA-provider checks; it does not make AS1
responsible for proving the changed seeds. The approximation change has its
own bound and is not treated as a representation equivalence.

The degree map is selected from finite candidate pairs by the final continuous
error budget. It is not claimed globally optimal. Uniform degree 6/6 fails the
current conservative bound on small rectangle 1 and one leaf of core rectangle
3, so those rectangles retain degree 7 in their first coordinate. Such a failed
bound is not an observed floating-point failure.

## Seed perturbation

For a normalized argument with magnitude at most `1 + 32u`, the omitted
Chebyshev modes are bounded by their coefficient 1-norm times
`T_10(1 + 32u)`, where `u = 2^-53`. Coefficient rounding, original and reduced
Horner evaluation, product rounding and rank reduction are bounded separately
in source order. Each continuous leaf supplies an interval for its actual
coordinate. This gives a finite relative bound `epsilon` between the original
computed seed `s_old` and the reduced computed seed `s0`.

The old root radius `r_old` transports to

```text
r = (r_old + epsilon) / (1 - epsilon).
```

Finite vega-ratio bounds transport the normalized residual and the geometry
`h=a/s0`, `t=s0/2`, `H=h²`, `T=t²`. The checker verifies the full enlarged
domain required by the imported compensated reserves on every leaf:
`0<=h<9`, `0<=t<3`, `h-t>-1.5`, `h+t<9`, `.01<D<2`, `ED/D<1e-10`, and
`s0>2^-350`. The maximum transported normalized residual is allowed to exceed
the historical rational-Householder cutoff; the new cubic rounding ledger
does not use that denominator or its cutoff.

## Finite inverse remainder

Write `B(s)` for the physical normalized Black price, `V=B'`,
`n=(b-B(s0))/(s0*V(s0))`, `A=H-T` and `Q=3H+T`. The cubic inverse update is

```text
delta3 = n - A*n²/2 + (2*A²+Q)*n³/6.
```

For the local inverse `G=B^-1`, direct differentiation gives

```text
B''/V   = A/s
B'''/V  = (A²-Q)/s²
B''''/V = (A³-3*A*Q+12*H)/s³
G''''   = -(6*A³+7*A*Q+12*H)/(s³*V⁴).
```

The checker independently verifies these polynomial identities. Taylor's
theorem is applied over the entire monotone price segment between `B(s0)` and
the target `b`. If `C4` bounds `|6A³+7AQ+12H|/24` on the corresponding root
segment and `K=Hmax/(1-r)³+Tmax*(1+r)`, its dimensionless remainder is at most

```text
C4 * |n|⁴ * exp(4*K*r) / (1-r)³.
```

The exponential is enclosed rationally by `exp(x)<=1/(1-x)` with an explicit
`0<=x<1/2` check. This is a finite fourth-derivative bound, rather than an
unbounded asymptotic remainder. Source-order geometry, residual, cubic
coefficients, Horner operations, and the final FMA are then charged explicitly.
The old root's lower attainable-accuracy denominator remains valid because
the target root and target price did not change.

The checker audits all 5,474 imported wing/small/core continuous seed leaves;
the common-finish bound covers 5,457 after the existing guard excludes 17 wing
leaves. It is conditional on the archived seed cover, Mills models, compensated
reserve lemmas, and floating primitive contracts. It does not regenerate those
certificates, establish a whole-domain solver theorem, prove the compiler or
vendor libm, or make a speed claim. Empirical replay and host timing are separate
evidence.
