# Deferred central expansion and adaptive degrees

The deferred central route keeps its normal-shape tables and all 98 stored
Black-conversion coefficients, evaluating 48 conversion coefficients through
shorter Horner chains. Final pair rounding is unchanged. Its
accepted `q.hi` cap increases from `.5` to `.625`. The outer envelope becomes
`a<=.45, b<.26`, and the conservative `b+a/2` preflight becomes `.25`.
It still accepts only `z.hi<=.5`. There is no new residual evaluation, final
polish or fallback.

This is a modest extension of the existing analytic expansion. It does not
implement the proposed conformal coordinate or claim a new wing-wide method.

The exact coordinates are

```text
chi = sinh(a/2)/(a/2), d = a*chi,
q = sqrt(2*pi)*(b+d/2), z = (d/q)^2.
```

The true Bachelier normal shape satisfies `0<f(z)<=1`, so its preliminary
volatility `v=q*f/chi` satisfies `v<=q`. The source-order ledger also explicitly
checks that the polynomial normal shape, with its computed interval, remains
inside the enlarged conversion rectangle. Including the computed `q` guard's
high/low error gives `q<.625000000001`; the new model rectangle is
`a²<=.203, v²<=.391`.

The [checker](../scripts/verify_central_expansion.py) recomputes the shortened
Black-conversion polynomial's complete residual and finite analytic tails on
that rectangle. It then recomputes native pair operations, both shortened
Horner calculations, coefficient/coordinate rounding, the sinhc tail and final
binary64 rounding. The resulting continuous helper bound is
**`rho_J<.570599`**.

## Lower normal-shape degrees

The four normal-shape intervals now evaluate degrees **18, 20, 20, 25**,
instead of 25 throughout. Constant-degree specializations remove 7, 5 and 5
two-lane Horner stages in the first three intervals. Both compensated final
stages and the derivative correction for the coordinate low part remain.
Coefficient storage is unchanged.

For the deleted suffix `T(x)=sum(c_i*x^i, i>d)`, the exact rational bound is

```text
|T(x_exact)| <= sum(|c_i|*R^i)
                 + (coordinate_low_bound + coordinate_error_bound)
                   * sum(i*|c_i|*S^(i-1)),
```

where `R` bounds the computed high coordinate and `S` also includes its low
part and error. Both cover the archived padded cells, including values just
outside `[-1,1]`. Thus the bound includes the omitted derivative times the
coordinate low part. The suffix bounds in units of `u=2^-53` are below
`0.003164`, `0.002370`, `0.002828` and zero. The checker separately recomputes
the shortened source-order arithmetic, derivative-coefficient rounding and
the domination required by the compensated final additions.

Combined with the shorter conversion below, the branch bounds are below
`0.567129`, `0.567769`, `0.566887` and `0.570599`. The maximum remains set by
the degree-25 normal interval. This preserves the conditional precision budget;
it does not assert that every result is more accurate. The degree-25 normal
arithmetic ledger is checked for exact equality with its archived counterpart.

## Shorter Black-conversion polynomial

Write the represented conversion polynomial as
`C(A,w)=sum(c_ij*A^i*w^j)`, where `A=a²`, `w=v²`, and the resulting total
volatility is `v*(1+w*C)`. The retained row lengths are now
**12, 11, 9, 7, 5, 3, 1, 0, 0, 0**, down from
16, 15, 14, 13, 11, 10, 8, 6, 4, 1. Only the lowest powers in each retained
row are evaluated; row six retains its constant and rows seven through nine
are omitted. The original coefficient arrays and represented bits remain.

For the omitted polynomial `DeltaC=C_full-C`, exact represented coefficients
give the continuous diagnostic

```text
|DeltaC(A,w)| <= sum_omitted |c_ij| * .203^i * .391^j,
.391 * |DeltaC| < 0.314710*u,  u=2^-53.
```

The final certificate recomputes the analytic residual for the **shortened
polynomial itself**. It does not reuse the old conversion error or insert the
diagnostic above as an unconditioned relative-output allowance. The archived
degree-24 bivariate interval arithmetic, exponential recurrence,
shifted-Chebyshev coefficient norm and Cauchy tail at radius 16 enclose
`E=|K(C)-exp(w/8)*T|`. On the entire root-localization segment of radius
`r=2^-28`, the lower derivative bound `m` satisfies `E<r*m`. The exact
residual-to-root and relative-volatility transfer used by the checker is

```text
relative_conversion_error
  <= WMAX * E / (m * (1 - WMAX * (Cmax+r))).
```

Here `Cmax` bounds the shortened polynomial on the whole rectangle. The
recomputed values are `E<1.719557e-17`, `m>0.989235`, and relative conversion
error below `6.913492e-18`. Native rounding is independently recomputed for
every retained row and the shortened outer Horner chain. Coordinate-high
errors use exact finite-polynomial derivative bounds for the new polynomial.

The conversion now evaluates 48 coefficients instead of 98. Its three paired
row evaluations and degree-six outer chain remove **22 two-lane Horner FMAs
and six scalar Horner FMAs** per accepted central call. This is an arithmetic
count; measured latency is recorded separately in the follow-up report.

The full-degree normal-shape model and its padded `z` cells are unchanged.
Their reviewed certificate and every imported proof module are bound by
SHA-256. The exact [source bridge](../scripts/central_expansion.py) restores
only the reviewed guards, normal-degree specialization and conversion
truncation to the old whole-file hash. Both degree reductions receive their
separate error bounds above; the bridge
does not assert bit equivalence. The complete
small source and coefficient representation are independently bound through
the seed source bridge, including the `Pair` constructors and caller. Central
coefficients must equal the original represented bits. A changed seed accuracy
certificate is not a premise of this helper proof.

The certificate remains conditional on the inherited exact normal-shape
model, the reviewed baseline arithmetic correspondence, RN-even binary64
basic operations, single-rounding FMA, gradual underflow and preservation of
source order and independent SIMD lanes. It is a helper proof, not a
whole-domain solver theorem or an end-to-end Lean result.

The [seam fixtures](../tests/experimental_central_expansion.rs) contain 1,026
unique exact binary64 inputs beside old/new `q`, `a`, preflight and `b` guards,
and throughout all four normal-shape intervals and their seams. Independent
Black roots and `eta_J` agree at
160 and 240 decimal digits. Both FMA feature policies pass these final-output
checks. Nineteen deliberate mutations, including changed source/proof degrees
and conversion row lengths,
a doubled `Pair::from` constructor, an enlarged `q` cap and a changed
conversion FMA, are rejected, alongside one positive binding test.

```sh
rtk proxy .venv/bin/python scripts/verify_central_expansion.py --reference-root /path/to/implied-black-volatility
rtk proxy .venv/bin/python scripts/test_central_expansion.py --reference-root /path/to/implied-black-volatility
rtk proxy .venv/bin/python scripts/generate_central_expansion_seams.py --output tests/experimental_central_expansion.rs
rtk cargo test --locked --offline --features experimental --test experimental_central_expansion
rtk cargo test --locked --offline --features experimental,fma --test experimental_central_expansion
```

The [initial 2026-11-06 report](theoretical-cost-reduction.md) records the
normal-degree reduction's accuracy and timing. The
[follow-up report](further-theory-improvements.md) records the conversion
truncation's validation and timing.
