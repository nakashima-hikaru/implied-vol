# Large-route direct log-residual correction

`scripts/verify_large_direct.py` checks the replacement of the large-route
final correction and the precise handoff at `2^-16`. The certificate is
conditional on the same primitive and real-analysis contracts as the bound
reference proof. It is not an end-to-end Lean proof or a vendor-libm guarantee.
Finite reference tests are separate evidence.

## Exact correction

At the current point write

```
x = s/2 - a/s,   y = s/2 + a/s,
q = log(p),     z = q_target - q,
lambda = p/(s*p'),
alpha = -x*y*lambda,
beta = (x*x + y*y - x*y)*lambda^2.
```

Here `p` is the normalized price on the lower route and its complement on
the upper route. Thus `lambda = D/s` below and `lambda = -D/s` above. The
same polynomial handles both signs:

```
delta = lambda * (z + (1-alpha)*z^2/2
                   + (1-3*alpha+2*alpha^2+beta)*z^3/6).
```

The inverse volatility is `s*(1+delta)` through order three. This follows
from differentiating with respect to `q`; the checker verifies the resulting
polynomial identities exactly. In particular, with
`gamma = 12*(a/s)^2*lambda^3`,

```
s''''/s = lambda*(1-7*alpha+12*alpha^2+6*beta
                  -6*alpha^3-7*alpha*beta-gamma).
```

The preserved brackets give `|lambda|<1` below and `<2/5` above,
`-4<alpha<1`, `0<=beta<16`, and `|gamma|<64`. These bounds follow from the
Mills integral and moment inequalities; the derivation is recorded beside
`symbolic_and_shape`. Both the current point and the internal target root
are inside the bracket, so the same bounds hold along the Taylor segment.
For `r = gate + 2^-38`, a conservative relative remainder is

```
(1213/24) * r^4 / (1-r)^2.
```

## Source rounding and target transfer

The `2^-38` residual bound is used only to localize the Taylor interval.
Final accuracy instead uses the compensated evaluator's much smaller bound

```
|lambda| * |residual_pair - true_residual|
  <= mu*.376u + (u/64)/(1-u/8) + 200000u^2,
u = 2^-53,  mu = 1 below or 2/5 above.
```

The checker recomputes the inherited source geometry and residual ledgers.
It then accounts for the actual new source operations, in their order:

- The rounded values of compensated `D` and `s`, and their quotient for
  `lambda`; the relative quotient error is less than `4u`.
- Every multiply, add, FMA and division used for `alpha`, `beta`, the cubic
  coefficients and the nested polynomial. Magnitude bounds can be sharpened
  analytically, but accumulated operation errors are never discarded.
- Conversion of the compensated residual to binary64, the compensated
  geometry error, the omitted product `s.lo*delta`, the FMA constructing the
  output low part, and the final rounded addition.
- Absolute allowances for zero and subnormal intermediates. All divisors
  used by the new correction are positive normal values; `s>1/8` bounds
  amplification of the final absolute error.

The resulting internal-root budgets are transferred to the original binary64
price by recalculating the existing lower-price, ordinary-cap, fixed-cap and
expm1 formulas with the new budgets. No old finalizer accuracy bound is used.
The target logarithm ledger is also rerun.

The exact rational checks establish these upper bounds (decimal values below
are rounded upward):

| Precise gate | Lower price | Ordinary/fixed cap | expm1 | Overall |
| --- | ---: | ---: | ---: | ---: |
| `2^-22` | 0.883813869 | 0.658213055 | 0.777170983 | 0.883813869 |
| `2^-17` | 0.884627735 | 0.659001481 | 0.777703704 | 0.884627735 |
| `2^-16` | 0.896240684 | 0.670587768 | 0.785532276 | 0.896240684 |

These are bounds on the reference proof's original-input `rho_J`, not
observed errors. The selected `2^-16` gate remains below `0.9`.

## Source and dependency binding

The checker binds the exact finalizer expression and accepted gate constants.
It compares unchanged geometry, compensated Mills, residual, target, integer
cap and iteration functions with the bound reference source. Explicit reverse
bridges cover only the compact Mills coefficient layout, the new signed
`lambda`, and the normal-price logarithm helper. All consumed constant bits
are checked, and the SIMD implementation is pinned by version and source hash.

The fast component retains its independently recomputed degree-16 certificate.
Both its numerical payload and its source functions/constants are bound;
changing a reported maximum alone is rejected. The old whole-file/finalizer
hash is not used as evidence for the new correction.

`verify_large_fast.py` certifies that component under explicit argument-bracket
and target-magnitude premises. It does not import this direct certificate.
This direct checker establishes the actual caller's target, bracket,
compensated evaluator and finish integration, then consumes the fast component
result. Thus the dependency is one-way. Removing the obsolete whole-file,
`solve` and `precise` bindings from the fast checker changes its declared scope;
no numerical formula or bound in its component calculation was changed.

The remaining imported premises are the Mills integral/continuous enclosure
lemmas, input/target bracket and finite-transfer lemmas, and the reference
proof's external Lefevre-Muller cap-separation theorem. Primitive contracts
require RN-even basic operations, one-rounding FMA and sqrt, gradual underflow,
no unsafe reassociation, reduced-log relative error at most `u`, fast-log
relative error at most `2u`, and reached expm1 relative error at most `2u`.

Run the checker with the reviewed reference checkout and fast certificate:

```sh
python3 scripts/verify_large_fast.py \
  --reference-root /path/to/implied-black-volatility \
  --out target/proof/large-fast-proof.json
python3 scripts/verify_large_direct.py \
  --reference-root /path/to/implied-black-volatility \
  --fast-certificate target/proof/large-fast-proof.json \
  --out target/proof/large-direct-proof.json
```

The source-binding negative controls deliberately changed a cubic coefficient,
the gate to `2^-15`, the lambda sign, a geometry expression, the normal-price
log helper, and the fast certificate's maximum. All six changes were rejected.
