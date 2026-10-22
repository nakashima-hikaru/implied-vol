# Large-route direct log-residual correction

`scripts/verify_large_direct.py` checks the replacement of the large-route
final correction and the precise handoff at `2^-14`. The certificate is
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
`-4<alpha<1`, `0<=beta<16`, and `|gamma|<64`. These broad bounds remain
checked. The current certificate tightens them on the actual accepted
Taylor segment, while retaining the same cubic expression.

## Accepted target and Taylor segment

The source lower-route decision is `ell.hi <= -RN(.6)`, and the exact
internal target is `ell.hi+ell.lo`. The source `add` returns a final `two`;
its represented low satisfies
`abs(lo) <= (u*abs(hi)+tau)/(1-u)`, where `tau=2^-1075`. Therefore

```
q_target <= -RN(.6)*(1-u/(1-u)) + tau/(1-u) < -.6+4u.
```

The recomputed precise high-component error bounds the difference from
`q(current)-q_target` by `2^-38`. With `r=gate+2^-38`, both endpoints
and every point on the inverse-Taylor segment thus satisfy
`q < -.6+4u+r < -.599`. The source's separate `2^-34` fast prefilter is
unchanged. Its bound includes rounding the DD target to binary64, and the
checker proves that its bypass also produces a precise high residual inside
each supported handoff gate.

For lower price `p=Phi(x)-phi(x)*M(sqrt(x*x+2a))`, a rational bound gives

```
p(a,.35) > 1582202077/2867200000 > .55,  a>10,
exp(-.599) < .55.
```

The first inequality uses `sqrt(.35^2+20)>112/25`, `phi(0)>.39`,
`M(y)<1/y`, and elementary integral/exponential Taylor bounds. The second
uses a positive Taylor lower bound for `exp(.599)`. These are continuous
inequalities, not sampled prices. Price increases with both `x` and `a`,
so the entire lower segment has `x<.35`.

For lower `x=-v<=0`, `lambda` is the average of the decreasing moment `I1`
on `[v,sqrt(v*v+2a)]`, hence `lambda<M(0)/sqrt(20)<1/3`.
For `0<=x<.35`, the Mills integral gives `M(-.35)<1288/751<7/4`;
this implies `lambda<2/5`, `abs(alpha)<1`, and `beta<5`.
The existing nonpositive-lower and upper inequalities complete the uniform
accepted-segment bounds

```
abs(lambda)<2/5,  abs(alpha)<1,  0<=beta<5,  abs(gamma)<8.
```

The last bound follows from `Q=3h^2+t^2>=3h^2`, so
`abs(gamma)<=4*beta*abs(lambda)`. Including the leading `lambda` in the
fourth derivative gives `abs(s''''/(24*s))<33/20`. A conservative relative
remainder is consequently

```
(33/20) * r^4 / (1-r)^2.
```

## Source rounding and target transfer

The `2^-38` residual bound is used only to localize the Taylor interval.
Final accuracy instead uses the compensated evaluator's much smaller bound

```
|lambda| * |residual_pair - true_residual|
  <= mu*.376u + (u/64)/(1-u/8) + 200000u^2,
u = 2^-53,  mu = 2/5 on both accepted routes.
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
The target logarithm ledger is also rerun. The lower original-price transfer
retains its older full `.376u` allowance; the new localization is not assumed
for the segment between the internal target root and the original-price root.

The exact rational checks establish these upper bounds (decimal values below
are rounded upward):

| Precise gate | Lower price | Ordinary/fixed cap | expm1 | Overall |
| --- | ---: | ---: | ---: | ---: |
| `2^-22` | 0.771012995 | 0.658212995 | 0.777170943 | 0.777170943 |
| `2^-17` | 0.771053541 | 0.658253541 | 0.777198338 | 0.777198338 |
| `2^-16` | 0.771447149 | 0.658647149 | 0.777464290 | 0.777464290 |
| `2^-14` | 0.874279206 | 0.761479206 | 0.846945409 | 0.874279206 |

These are bounds on the reference proof's original-input `rho_J`, not
observed errors. The selected `2^-14` gate remains below `0.9`.
Changing a gate changes when precise evaluation hands off to the cubic;
it does not itself establish a runtime speedup.

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

The component generator solves independent Gaussian-erfc roots at 160 and
240 decimal digits and requires identical stored reference bits. It includes
the `2^-14` handoff, both sides of its threshold, Mills seams, and the lower
`ell.hi=-RN(.6)` routing seam with both signs of `ell.lo`.

Eleven negative controls were rejected: a cubic coefficient change, the
unproved `2^-13` gate, lambda sign, geometry, target logarithm, lower routing
cutoff, dropped target low, fast prefilter, target-pair renormalization,
the fast certificate's numerical bound, and Python's disabled-assertion mode.
The finite component checks and these binding controls are separate from
the continuous conditional proof.
