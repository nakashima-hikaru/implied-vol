# Numerical methods and accuracy

Black solver selection is explicit through `calculate_with::<S>()` on the full
and normalized implied-volatility builders. `calculate::<SpFn>()` always selects
Hybrid. Solver features expose additional types and do not change that default.
Pricing, Bachelier inversion, and the separate explicit Black formula are
independent of the selected Black inverse solver.

## Coordinates and boundaries

For the normalized Black problem, log-moneyness is `x=ln(F/K)`, normalized OTM
price is `b=price/sqrt(F*K)`, and total volatility is `s=sigma*sqrt(T)`.
The OTM upper limit is `exp(-abs(x)/2)`. Full builders perform price bounds and
normalization before inversion and annualize the result afterward.

Hybrid, Jaeckel, and FlashIv share the existing boundary handling. Experimental
retains its own ATM inversion and classification against the mathematical
normalized-price cap, including a price equal to a rounded-down exponential.
Custom `SpecialFn` providers apply to Hybrid, Jaeckel, and FlashIv. Experimental
uses its own functions and explicit FMA to preserve its numerical arithmetic.

## Solver methods

**Hybrid** uses an independently implemented FlashIV iteration in the qualified
lowest-price interpolation region, for `abs(x)>=0.01`. Other inputs use Let's Be
Rational. Its exact FlashIV steps use stable Black expansions in the asymptotic
and small-volatility regions. The route uses one approximate step, two exact
Householder steps, and a third when selected by the residual.

**Jaeckel** uses Let's Be Rational interpolation and its Householder corrections.
It does not switch to another inverse method.

**FlashIv** uses its own FlashIV seeds, scaled Mills ratios near ATM, a
complementary log-price objective near the upper bound, and safeguarded exact
Householder iterations. It is a numerical variant of FlashIV rather than a
bit-identical implementation of the fixed-count paper algorithm. It has no
Let's Be Rational inverse fallback or FlashIV+ final Newton alignment.
Numerical nonconvergence returns `None`.

**Experimental** is the author's own implementation. It uses fitted seeds,
compensated residuals, Mills evaluation, exact cap handling, and Householder
corrections. Its original dispatcher includes an internal LBR
branch. No additional output search or refinement wraps the selected solver.

## Accuracy contracts

The default hybrid targets Jäckel's attainable accuracy. Faster routes must pass
the existing accuracy checks without relaxed thresholds. Seeded round trips,
independent high-precision roots, and the optional C++ comparison exercise
regressions. These finite checks do not establish a universal error bound;
existing Jaeckel and forward-price evaluator limitations remain.

Experimental's core metric uses the exact mathematical root `s_star` and vega `V_star`:

```text
V_star = exp(-a*a/(2*s_star*s_star) - s_star*s_star/8) / sqrt(2*pi)
eta_J  = 2^-52 * (1 + b/(s_star*V_star)),  a=abs(x)
rho_J  = abs(s_hat/s_star - 1) / eta_J
```

Its accuracy target is `rho_J<1` for exact finite binary64 `a>=0` and positive
normal prices `2^-1022<=b<exp(-a/2)`, using the mathematical cap. This is a
conditioning-based bound rather than correct rounding or uniform 16-digit
relative accuracy. General subnormal prices, market normalization, and
annualization are outside the core contract.

The analytical bound `rho_J<0.999991` is conditional on round-to-nearest-even,
gradual underflow, preserved operation order, and correctly rounded explicit
FMA/sqrt. It assumes relative errors at most `u=2^-53` for exp/log and `2u` for
log1p/expm1 at reached arguments, plus an external
[mathematical separation theorem](https://perso.ens-lyon.fr/jean-michel.muller/TMDworstcases.pdf).
Rust's standard library does not guarantee those platform math error bounds.
Full end-to-end Lean verification remains unfinished; component proofs and
finite testing are separate evidence.

[Source correspondence](../performance/experimental-source-2026-10-03.json)
records unchanged numerical modules after normalizing module qualification and
renamed comment labels.
[Independent-reference validation](../performance/experimental-validation-2026-10-03.json)
replays 46,688 core references and 68,273 paper inputs through explicit
Experimental selection in a native release build with all features. The imported
implementation matched bit for bit, with zero
invalid outputs or `rho_J>1`, including 3,075 valid prices equal to a rounded-down
cap. Paper references target exact `(x,c)` inputs and include adapter rounding;
they are distinct from the exact `(a,b)` core contract.

## Algorithm references

- Peter Jäckel, [Let's Be Rational](https://www.jaeckel.org/LetsBeRational.pdf).
- Peter Jäckel, [Implied Normal Volatility](https://www.jaeckel.org/ImpliedNormalVolatility.pdf).
- [An Explicit Solution to Black-Scholes Implied Volatility](https://arxiv.org/abs/2604.24480),
  the separate inverse-Gaussian representation.
- [FlashIV](https://arxiv.org/abs/2605.29102v1), log-price Householder inversion.
- Li and Lee, [equation 42 and seed coefficients](https://mpra.ub.uni-muenchen.de/6867/1/MPRA_paper_6867.pdf).

The default build has no required runtime dependencies. Built-in inverse kernels
allocate no memory. The optional C++ comparison feature adds its own dependencies.
