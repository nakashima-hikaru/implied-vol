# Numerical methods and accuracy

`calculate::<SpFn>()` selects Hybrid. Use `calculate_with::<S>()` on either
Black implied-volatility builder to select another solver. Features expose
optional types without changing the default. Pricing, Bachelier inversion,
and `calculate_explicit()` are separate from Black solver selection.

## Coordinates and boundaries

For the normalized Black problem, log-moneyness is `x=ln(F/K)`, normalized OTM
price is `b=price/sqrt(F*K)`, and total volatility is `s=sigma*sqrt(T)`.
The mathematical price cap is `exp(-abs(x)/2)`. Full builders validate prices,
normalize inputs, and annualize the result after inversion.

Hybrid, Jaeckel, and FlashIv classify against a rounded price cap. A price
equal to a rounded-down cap can therefore return infinity even when its exact
mathematical root is finite. Experimental compares against the mathematical
cap using compensated arithmetic and has its own ATM inversion.

Custom `SpecialFn` providers apply to Hybrid, Jaeckel, and FlashIv. Experimental
uses its own functions and explicit FMA. Provider accuracy and arithmetic policy
affect numerical results; the built-in kernels assume mathematical functions
with small errors at the arguments they reach.

## Solver methods

**Hybrid** uses a FlashIV iteration in selected lowest-price regions for
`abs(x)>=0.01`, and Let's Be Rational elsewhere. A normal cap and
`b<=0.0005*b_max` select the low-price route before constructing LBR's
interpolation nodes. That route uses stable Black expansions, one approximate
Householder step, two full-precision steps, and a third selected by the residual.
If the route declines, LBR handles the input.

**Jaeckel** uses Let's Be Rational interpolation and Householder corrections.
For `abs(x)<=2^-20`, Hybrid and Jaeckel use cancellation-free lower-node
evaluation and scaled lowest-price arithmetic. This protects microscopic
near-ATM inputs while retaining the Black pricing equation.

**FlashIv** follows [FlashIV Algorithm 1](https://arxiv.org/pdf/2605.29102v1).
The ordinary path uses a Li/asymptotic seed, one inexpensive Householder step,
two full-precision steps, and a third when the residual entering the second
full-precision step is at least `1e-4`. It evaluates the erfcx/log-price objective
directly. Microscopic near-ATM inputs use a terminal Bachelier/Mills branch;
prices at or above the `0.99` forward-normalized threshold use three Halley
steps on the complementary log-price objective.

FlashIv uses fixed iteration counts. Invalid arithmetic steps retain the
current iterate; a finite result does not certify convergence. Cancellation
in the erfcx difference and seed error in extreme tails can produce larger
errors than Jäckel's attainable-precision target. Its regression allowance
reflects this method's accuracy trade-off.

**Experimental** uses fitted seeds, compensated residuals, Mills evaluation,
and Householder corrections, with an internal LBR branch. Its common finish
uses seed coordinates `h=abs(x)/s0` and `t=s0/2`. It evaluates `D(h,t)` directly
for `0.001<t<=1/4`, and for `1/4<t<=1/2` when `h>=2`, selecting shorter
moment-polynomial rows at smaller `t`. It uses explicit FMA and SIMD lanes where
available. Its target is the normal-price core contract defined below.

The separate **explicit inverse-Gaussian formula**, selected with
`calculate_explicit::<SpFn>()`, checks relative quantile changes and uses a
rationalized reciprocal form for moderate parameters. Intermediate mean or
quantile overflow can prevent a representable volatility from being returned.

## Accuracy contracts

Jäckel's attainable-precision estimate accounts for the conditioning of the
inverse price problem. For an exact binary64 input `a=abs(x)` and price `b`,
let `s_star` be the exact mathematical root and `V_star` its normalized vega:

```text
V_star = exp(-a*a/(2*s_star*s_star) - s_star*s_star/8) / sqrt(2*pi)
eta_J  = 2^-52 * (1 + b/(s_star*V_star))
rho_J  = abs(s_hat/s_star - 1) / eta_J
```

`rho_J<1` means the result meets this conditioning-based target. It does not
mean correct rounding or uniformly 16 accurate decimal digits: sensitive
inputs permit a larger relative volatility error. Hybrid and Jaeckel target
attainable precision but exceed one estimate on some tested edge cases.
FlashIv has a different accuracy trade-off. The
[shared-input accuracy comparison](performance.md#accuracy) reports each
solver's finite-output error and boundary outputs separately.

Experimental's core target is `rho_J<1` for finite exact binary64 `a>=0` and
positive normal prices `2^-1022<=b<exp(-a/2)`, using the mathematical cap.
General subnormal prices, full-API normalization, and annualization lie outside
that core contract. Paper/adapter references for exact `(x,c)` inputs include
adapter rounding and are assessed separately from exact `(a,b)` core roots.

Validation uses independent high-precision roots, conditioned round trips,
boundary and seam inputs, and optional C++ comparisons. Near-ATM regression
fixtures include cutoff neighbors, minimum normal prices, and microscopic
volatility. These finite checks cover selected inputs; they do not establish
a universal accuracy guarantee for any solver or custom function provider.

## Conditional precision bounds

Experimental's AS1 component has a conditional bound `rho_AS1<0.916169`.
Its common finish paths have conditional bounds below `0.553` for wing,
`0.775` for small-rank, and `0.787` for finite-rank. These bounds concern
components under their seed/Mills certificate assumptions, not the complete
solver domain.

The certificates assume round-to-nearest-even, gradual underflow, preserved
operation order, correctly rounded explicit FMA/sqrt, relative errors at most
`u=2^-53` for exp/log and `2u` for log1p/expm1 at reached arguments, and the
required mathematical separation results. Rust's standard library does not
guarantee these platform math error bounds. End-to-end Lean verification is
unfinished; conditional component proofs and finite testing are distinct
evidence.

`scripts/verify_as1_atan.py` checks the represented AS1 polynomial and its
source-order rounding bound. `scripts/verify_direct_d.py --reference-root PATH`
checks the D0 polynomial, analytic moment tails, and source-order recurrence
rounding with the required certificates in an `implied-black-volatility` checkout.
These component checks bind the expected source and dependencies; they do not
establish vendor math-function accuracy.
[Current validation data](results.json) identifies the measured source and bounds.

## Algorithm references

- Peter Jäckel, [Let's Be Rational](https://www.jaeckel.org/LetsBeRational.pdf).
- Peter Jäckel, [Implied Normal Volatility](https://www.jaeckel.org/ImpliedNormalVolatility.pdf).
- [An Explicit Solution to Black-Scholes Implied Volatility](https://arxiv.org/abs/2604.24480),
  the inverse-Gaussian representation.
- [FlashIV](https://arxiv.org/abs/2605.29102v1), log-price Householder inversion.
- Li and Lee, [equation 42 and seed coefficients](https://mpra.ub.uni-muenchen.de/6867/1/MPRA_paper_6867.pdf).
