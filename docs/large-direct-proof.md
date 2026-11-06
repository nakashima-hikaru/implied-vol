# Large-route direct log-residual correction

`scripts/verify_large_direct.py` checks the replacement of the large-route
final correction, degree-16 compensated Mills table evaluation, adaptive
upper predictor, third-reciprocal-derivative log iteration, and the precise
handoff at `2^-14`. The certificate is
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
  <= mu*.376u + .014590506u/(1-.087000966u) + 200000u^2,
u = 2^-53,  mu = 2/5 on both accepted routes.
```

The checker recomputes the inherited source geometry and residual arithmetic
ledgers, then replaces the inherited Mills model allowances with the new
continuous degree-16 bounds.
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
| `2^-22` | 0.770495746 | 0.657695746 | 0.776821450 | 0.776821450 |
| `2^-17` | 0.770536223 | 0.657736223 | 0.776848800 | 0.776848800 |
| `2^-16` | 0.770929761 | 0.658129761 | 0.777114704 | 0.777114704 |
| `2^-14` | 0.873761384 | 0.760961384 | 0.846595530 | 0.873761384 |

These are bounds on the reference proof's original-input `rho_J`, not
observed errors. The selected `2^-14` gate remains below `0.9`.
Changing a gate changes when precise evaluation hands off to the cubic;
it does not itself establish a runtime speedup.

## Shortened precise Mills evaluation

The table lanes now retain coefficients `0..=16`, with the same three
compensated final stages and the same correction for the argument low part.
An asymptotic lane retains degree 18. A table/table pair therefore removes
eight two-lane Horner stages; a mixed pair removes six. The table lane of a
mixed pair inserts exact zeros at degrees 18 and 17.

For each half-open cell `[k/2,(k+1)/2)`, the checker decodes every coefficient
as an exact binary64 fraction. It bounds the omitted tail by
`sum(j=17..24, abs(c[k][j])*r^j)`, with
`r=(1+u)/4+NU` rounded outward. The small enlargement includes local-centering
rounding. Together with `s=sqrt(x*x+2a)+x`, `a>10`, and the proved accepted
lower geometry `x<.35`, the continuous deleted-tail budget is

```
(abs(delta M1)+abs(delta M2))/(s*u) < .010460408.
```

This bounds the difference from the represented degree-24 polynomial.
The full new bound additionally includes the true-Mills Taylor enclosure,
stored coefficient defects, source Horner rounding, the argument low-part
correction, geometry defects and DD pair formation. Reusing the inherited
source-order checker with exactly zero-padded high coefficients yields

```
accepted |delta D|/(D*u) < .087000966,
accepted |delta D|/(s*u) < .014590506.
```

Leading zero Horner stages are exact for finite arguments. Their only effect
on the inherited bound is an extra harmless absolute-underflow allowance,
so the zero-padded bound covers the shortened source. The degree-18
asymptotic calculation is rerun unchanged. Every interval box is continuous;
no finite solver sample enters these inequalities.

Before acceptance, lower iteration still permits `x` up to 1. Its broader
model bound is `|delta D|/(D*u)<.630598273`. The checker uses that bound only
for the independently recomputed `2^-38` precise-high localization estimate,
then uses the tighter accepted bounds in the final error ledger. Exact cell
selection follows `floor(2*zhi)`; a lower endpoint at a seam cannot select the
preceding half-open cell.

The final `rho_J` bound is slightly tighter than the prior conservative
allowance because this ledger consumes the recomputed `.014590506u` bound
instead of the former broad `u/64` allowance. This is a proof-budget change,
not a claim that degree 16 is intrinsically more accurate than degree 24.

## Degree-13 fast Mills evaluation

The fast evaluator determines the bracket sign and whether to run the precise
finalizer. It does not need the finalizer's near-roundoff accuracy. Its table
polynomial now uses coefficients `0..=13`, matching the existing degree-13
asymptotic polynomial. Compared with degree 16, this removes three scalar
Horner stages from the upper predictor's table call and three two-lane stages
from each fast pair containing a table lane. Mixed pairs need no zero padding.
The compensated degree-16/18 evaluator is unchanged.

`verify_large_fast.py` independently recomputes all 34 continuous table cells
and 114 residual boxes. For each table cell the error budget includes
`sum(j=14..24, abs(c[j])*r^j)`, the original degree-24 continuous Taylor tail,
represented coefficient errors, all 13 source FMAs and centering error. The
largest omitted high tail increases to less than `1.562012e-13`; after the
appropriate cell/geometry conditioning, the worst overall fast log-residual
bound remains `7.752447e-13 < 2^-38`. Thus the existing sign, handoff and
forced-bisection inequalities retain the same bounds. The final `rho_J`
budget is unchanged because the precise evaluator and acceptance gates are
unchanged.

Degree 12 also passed an isolated error-budget probe, but it uses almost all
of the `2^-38` allowance and would require separate mixed-lane padding because
the asymptotic polynomial remains degree 13. Degree 11 failed that allowance.
Degree 13 therefore gives the selected simple reduction without spending more
of the worst-case residual budget. Unconditional order-three/four upper physical
predictors missed the precise gate on some saved inputs; the adaptive selection
below retains order five for the larger steps.

## Adaptive upper physical predictor

At the inverse-normal upper seed, write `y=sqrt(x*x+2*a)`, `s=x+y` and
`n=M(y)*y/s`. With an exact inverse-normal seed this is the Newton step for
the physical complement-price residual. For a general exact Newton ratio
`n=-f/f'` and `q_j=f^(j)*n^(j-1)/(j!*f')`,

```
f(x+n*t)/f(x) = 1-t-q2*t²-q3*t³-...,
f(x)/f(x+n*t) = c0+c1*t+c2*t²+c3*t³+...,
c0=c1=1,   c_k=c_(k-1)+sum(j=2..k, q_j*c_(k-j)).
```

Householder's `d`-th reciprocal derivative gives `n*c_(d-1)/c_d`. The checker
verifies the reciprocal polynomial identities through degree five exactly.
The source retains the original represented `g`, `w`, `q` and `c` recurrence
expressions, with their original arithmetic order and one-rounding FMAs.
It evaluates later expressions only when the selected order needs them:

| Computed measure | Reciprocal derivative `d` | Proposal |
| --- | ---: | --- |
| `abs(x*n) <= 1/8` | 3 | `n*c2/c3` |
| `1/8 < abs(x*n) <= 3/8` | 4 | `n*c3/c4` |
| otherwise | 5 | `n*c4/c5` |

The measure is the magnitude of the **linear term** in the Gaussian exponent
change, `-x*n-n²/2`. These thresholds are heuristic choices of computational
work, not proved one-step error thresholds. Individual ideal Householder
formulas with exact derivatives and Newton ratios have local orders 4, 5 and
6. Those are not convergence-rate claims for this single, rounded adaptive
predictor or for the complete solver.

The predictor is accepted only when finite and strictly inside the original
bracket. Its precise residual must still satisfy `2^-14`; otherwise the
guarded fast update and forced midpoint every fourth update run.
The continuous correctness and stopping argument therefore does not assume
that any selected order reaches the precise handoff in one step. Invalid
proposals or zero rational denominators remain covered by these same guards.
No universal one-step residual bound for inverse-normal plus this predictor is
claimed; finite one-step residual/count checks are performance evidence only.

`scripts/generate_large_predictor_seams.py` locates the two ideal thresholds,
then rounds the corresponding physical prices to binary64 and includes adjacent
quotes and offsets on both sides. It also covers representative order-three
and order-five benchmark inputs. All 95 final-output references independently
agree at 160 and 240 decimal digits using direct Gaussian-erfc inversion of
the original represented input. Runtime branch coverage is checked separately;
ideal threshold locations are not assumed to equal rounded source boundaries.

## Shortened logarithmic iteration

The guarded iteration now uses Householder's third reciprocal derivative,
which has fourth-order local convergence, in place of the fifth reciprocal
derivative. For the log-price residual `f(x)`, define

```
n = -f/f',  A = f''/f',  B = f'''/f'.
step = 3*(1/f)''/(1/f)'''
     = n*(1 + A*n/2)/(1 + A*n + B*n*n/6).
```

Writing `d=f'` and `H=p''/p'`, differentiation gives
`A=H-d` and `B=H'+A*(H-2*d)`. The source retains the old expressions for
`H`, `H'`, `A` and `B`, with a short fused rational evaluation. It eliminates
the `H''`, `H'''`, fourth- and fifth-residual-derivative expressions. The exact
rational checker verifies the reciprocal identity and agreement with the
inverse Taylor series

```
step = n - A*n²/2 + (A²/2 - B/6)*n³ + O(n⁴).
```

The local-order statement assumes a smooth residual and a simple root. It is
not used as a global convergence premise. A proposal is still accepted only
when finite and inside the current bracket; a rounded midpoint is forced
every fourth update. Residual signs used for bracket updates are certified
outside the fast filter, and the precise `2^-14` handoff is unchanged.
Consequently the existing stopping and final accuracy bounds apply even when
this shorter proposal needs additional iterations or its denominator is zero.
No uniform one-step handoff claim is made. A lower-order Halley probe missed
the first handoff on some saved inputs, so the third-derivative proposal is
the selected compromise; runtime observations are reported separately.

`scripts/generate_large_iteration_cases.py` adds direct Gaussian-erfc root
fixtures near the `a=10` and `ell=-.6` switches, the largest first residual
found in the saved-input probe, and lower tails down to the minimum normal
price. Each stored root high/low pair and `eta_J` agrees at 160 and 240 decimal
digits. These fixtures validate final output accuracy on finite exact inputs.

## Source and dependency binding

The checker binds the exact finalizer expression and accepted gate constants.
It compares unchanged geometry, compensated Mills tail stages, residual,
target, integer cap and safeguarded iteration functions with the bound reference
source. Explicit reverse bridges cover the compact Mills coefficient layout,
degree-16 table truncation, the new signed `lambda`, the normal-price logarithm
helper, and staging the physical predictor recurrence at orders three through five. All consumed constant bits
are checked, and the SIMD implementation is pinned by version and source hash.
The adaptive predictor bridge copies every retained expression from the bound
reference and requires the exact two guards and three return formulas. It also
binds the one-rounding `fma` helper. Only independent pure arithmetic is moved
between stages; no reached expression's rounding order is changed.
The new `log_step3` body is bound explicitly and its algebra is checked exactly;
only its call name is reversed when comparing the safeguarded caller with the
reference. Removing forced bisections or modifying the rational coefficients
must fail source validation.

The fast component uses its independently recomputed degree-13 certificate.
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

`scripts/test_large_reduction_binding.py` additionally runs a complete positive
control against an isolated source copy and rejects seven reduction defects:
table degree 17, missing table zero padding in mixed pairs, lowering the mixed
asymptotic degree to 16, removing a compensated stage, changing a physical
predictor coefficient, truncating it at order four, and dropping its strict
bracket guard. Each must fail at its intended source-order binding.

```sh
python3 scripts/test_large_reduction_binding.py \
  --reference-root /path/to/implied-black-volatility \
  --fast-certificate target/proof/large-fast-proof.json
```
