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

For a normal cap, `abs(x)>=0.01` and `b<=0.0005*b_max` select that same
restricted path before constructing LBR's interpolation nodes. This guard
lies strictly inside both original lower boundaries. Writing `r=sqrt(abs(x))`,
the integral for `1-erfcx(r)` on `[0.5,1]` gives
`b_c/b_max > exp(-1)/(22*sqrt(pi)) > 0.00943` for `r>=0.1`.
For the first fitted `b_l` polynomial, its positive rational term and
`s_c^2*(0.0756099664-0.0967271929*s_c)` give bounds above `0.0005447`
on `[sqrt(0.02),0.5]` and above `0.0017304` on `[0.5,0.71]`.
In each remaining rational segment, `N-0.01D` has positive cubic-and-higher
coefficients and a positive, increasing quadratic remainder above the
segment's lower endpoint, hence `b_l/b_max>0.01`. A normal cap limits `s_c`
to less than 54. The smallest first-segment margin exceeds `4.47e-5` of
the cap, well beyond binary64 rounding even when the threshold is subnormal.
Subnormal caps retain the original classifier. If the early FlashIV attempt
declines, LBR runs once; the same attempt is not repeated. The iterations,
arithmetic policy, and default numerical outputs are preserved; calls to a
custom mathematical `SpecialFn` provider may occur in a different order.

**Jaeckel** uses Let's Be Rational interpolation and its Householder corrections.
It does not switch to another inverse method.

For `abs(x)<=2^-20`, both Hybrid and Jaeckel evaluate the lower tangent node
with a cancellation-free Taylor polynomial in `sqrt(abs(x))`. The first omitted
term is below `3e-30` at the cutoff. Lowest-branch interpolation uses `b/abs(x)`
and the similarly scaled lower map, while lower and middle corrections use
relative volatility changes and scaled derivative ratios. These are coordinate
changes to the Black solver, not a Bachelier approximation. The existing Black
price expansions and the two-correction limit are retained.

The separate explicit inverse-Gaussian implementation checks relative changes
in its standardized quantile. An absolute step tolerance would stop too early
as log-moneyness approaches zero. Its moderate-parameter mode is evaluated in
rationalized reciprocal form to avoid subtracting nearly equal terms. The
inverse-Gaussian method and iteration limits are unchanged; intermediate mean
or quantile overflow can still prevent a representable volatility from being
returned.

**FlashIv** follows [FlashIV Algorithm 1](https://arxiv.org/pdf/2605.29102v1).
The ordinary path uses the Li/asymptotic seed, one inexpensive Householder
step, two full-precision steps, and a third only when the residual entering
the second full-precision step is at least `1e-4`. It evaluates the paper's
erfcx/log-price objective directly. Microscopic near-ATM prices use a terminal
Bachelier/Mills branch; prices at or above the `0.99` forward-normalized
threshold use three Halley steps on the complementary log-price objective.
The ordinary path has no adaptive convergence loop or unconditional bracket
search. The solver has no Let's Be Rational inverse fallback or FlashIV+ final
Newton alignment. Invalid
arithmetic steps retain the current iterate as in the author's implementation;
a finite return value does not certify convergence. Invalid inputs return
`None`, and the microscopic guard may return its defensive zero limit.

This is an independent implementation of the paper's algorithm, using the
crate's sqrt-forward input convention, `SpecialFn` provider, and optional FMA
policy. The [author's May 2026 source archive](https://chasethedevil.github.io/post/thiophene-iv-rust-full.zip)
supplies the shared microscopic-branch formulas and rational seed coefficients.
That source predates the paper and uses two Householder steps near the upper
bound; Algorithm 1 specifies three Halley steps and takes precedence here.
The source calls its Bachelier seed LFK2026, whereas the FlashIV paper cites
LFK2016. These version differences and arithmetic choices preclude a claim
of bit-identical reference outputs.

FlashIv has the paper method's accuracy trade-off. It does not adopt the
default hybrid's Jäckel attainable-precision contract: cancellation in the
erfcx difference and an inadequate seed in extreme tails can cause larger
errors even when the returned result is finite. The paper's finite benchmark
results are not a whole-domain precision guarantee.

The regression tests retain the previous input population. Hybrid and Jaeckel
keep their existing attainable-precision thresholds; FlashIv uses a separate
empirical allowance for erfcx subtraction, log-price rounding, and the rational
Bachelier seed. For example, a round trip at `x=-1e-7` and total volatility
`1.25e-7`, just outside the microscopic guard, showed relative error around
`6e-10`. Passing the FlashIv checks therefore does not imply passing Jäckel's
precision target on those inputs.

**Experimental** is the author's own implementation. It uses fitted seeds,
compensated residuals, Mills evaluation, exact cap handling, and Householder
corrections. Its original dispatcher includes an internal LBR
branch. No additional output search or refinement wraps the selected solver.
Independent Horner chains share two SIMD lanes where available, while retaining
each chain's degree and FMA sequence and the scalar reduction order.

Conservative preflights skip attempts that fail the existing range checks.
For the deferred caller's `a<=0.36`, `b<a/16` implies `z>0.503`, and
`b+a/2>=0.2` implies `q>0.5013`, with margins exceeding rounding errors.
For validated `0.5<=a<=10`, `2^(-ceil(a)-1)` is an exact lower bound on half
the cap because `ln(2)>0.5`; a price at or below this bound cannot enter the upper
route. The large route specializes the initial logarithm for its already
validated positive normal price, preserving both compensated components.

The AS1 low-price route includes the next term of its inverse-logarithm seed.
In exact arithmetic, with `L=2*log(a/(b*sqrt(2*pi)))`, `l=log(L)`, and
`A=a*a/4`, the seed for `y=(a/s)^2` is

```text
y0 = L - 3*l + (9*l - 6 - A)/L
y  = y0 + (13.5*l*l - (45 + 3*A)*l + 39 + 5*A)/(L*L)
```

The additional polynomial uses explicit FMA and retains the single Householder
correction. Its changed geometry can send an input to the existing wing route;
the original guards remain in force. This improves precision margin at a small
AS1 timing cost, documented in [the performance comparison](performance.md#experimental-as1-seed-refinement).

AS1's quadrature evaluates `atan(sqrt(q))/sqrt(q)` with five explicit FMA
stages instead of the original seven. The coefficients were generated by
high-precision Remez exchange and rounded to binary64. An independent exact
rational checker encloses the represented polynomial with Bernstein bounds
and the alternating integral series; it does not assume that the numerical fit
is an exact minimax solution. The reached source argument satisfies
`0<=q<(1/120)*(1+16u)`, with `u=2^-53`, including rounding in the geometry
and quadrature.
The checker also bounds the wider interval up to `1/50` used by intermediate
proof rectangles. Run `python3 scripts/verify_as1_atan.py` to check the
polynomial and its five-FMA rounding bound against the recorded source implementation.
This helper check is separate from the final AS1 error certificate below.

The intermediate-`t` branch of `forward_D` selects
`min(mills_degree[cell]+1,18)` for its divided-difference polynomial. This
retains degree 18 for `h<1.5` and uses degrees 13 through 17 in the remaining
accepted cells. The last three compensated stages, the `t=0.001` and
`t=0.0625` branch boundaries, and the Householder correction are unchanged.
The degree reduction is justified by the error of the divided difference
itself, rather than by reusing the Mills-function degree without a new bound.

## Accuracy contracts

The default hybrid targets Jäckel's attainable accuracy. Faster routes must pass
the existing accuracy checks without relaxed thresholds. Seeded round trips,
independent high-precision roots, and the optional C++ comparison exercise
regressions. These finite checks do not establish a universal error bound;
existing Jaeckel and forward-price evaluator limitations remain.

Near-ATM regressions include 81 independent exact-input roots evaluated at two
precisions (at least 100 decimal digits), with cutoff neighbors, minimum normal
prices, microscopic volatility, and ordinary prices at tiny log-moneyness. The
fixtures on the stabilized path require one attainable-accuracy estimate; the
neighbor above its cutoff retains the existing four-estimate envelope. Both
explicit FMA policies are checked. This is finite regression evidence, not a
whole-domain precision proof.

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

The original implementation's analytical bound `rho_J<0.999991` is conditional on round-to-nearest-even,
gradual underflow, preserved operation order, and correctly rounded explicit
FMA/sqrt. It assumes relative errors at most `u=2^-53` for exp/log and `2u` for
log1p/expm1 at reached arguments, plus an external
[mathematical separation theorem](https://perso.ens-lyon.fr/jean-michel.muller/TMDworstcases.pdf).
Rust's standard library does not guarantee those platform math error bounds.
Full end-to-end Lean verification remains unfinished; component proofs and
finite testing are separate evidence.

Independent-reference validation during integration covered 46,688 core
references and 68,273 paper inputs through explicit Experimental selection in a
native release build with all features. The imported
implementation matched bit for bit, with zero
invalid outputs or `rho_J>1`, including 3,075 valid prices equal to a rounded-down
cap. Paper references target exact `(x,c)` inputs and include adapter rounding;
they are distinct from the exact `(a,b)` core contract.

Polynomial optimization revalidation expanded this population with 2,083
independent tiny/ATM/seam/cap roots. Both FMA feature policies passed all 117,044
inputs with bit-identical baseline outputs, zero invalid outputs, and zero
`rho_J>=1`. This population includes 3,121 valid core-coordinate prices equal
to a rounded-down cap and contains only positive normal adapted prices.

Dispatch/logarithm optimization passed the same 117,044 independent inputs
under both FMA policies with unchanged bits and the same strict `rho_J<1`
checks. A separate unfiltered 660-input seam replay also preserved all output
bits and classifications, including eight subnormal-price inputs and five
existing invalid cases. These seam inputs have no independent root references;
their parity does not extend the subnormal accuracy contract.

The AS1 next-term revision passed all 117,044 archived references under their
respective core or paper/adapter metrics with both FMA feature policies. One
archived output changed and improved. An additional 489 exact core-coordinate
inputs covered AS1 geometry, quadrature rules, tiny scaling, and dispatch seams:
91 errors decreased, 12 increased, and 386 outputs were unchanged. All remained
strictly below `rho_J=1`; the largest observed ratio fell from `0.742010` to
`0.493706`. Sixteen independent two-precision roots are checked in as regressions.

Conditional interval revalidation of the revised AS1 helper covered the original
606 rectangles on `L in [70,1424]`, `a^2 in [0,100]`: 554 accepted-region leaves
and 52 geometry exclusions, with no incomplete cells. Its bound is
`rho_AS1<0.915824`, using the original primitive assumptions and unchanged Mills
component certificates. The source-order seed error remains within the existing
preprocessing allowance. A separate routing bound preserves the wing/tiny
connection for declined AS1 attempts. These results concern AS1 and its
fallthrough; they are not a new whole-Experimental or end-to-end Lean certificate.
The [AS1 validation record](experimental-as1-2026-10-03.json) separates finite
measurements, conditional bounds, and the two input-coordinate contracts.

The subsequent degree-five AS1 revision reuses the same seed and geometry;
all 554 accepted leaves and 52 geometry exclusions of the 606-rectangle cover
were recomposed with the new polynomial and its native rounding errors.
The resulting conditional bound is `rho_AS1<0.916169`. The small increase
from the previous bound spends precision margin to remove two FMA stages per
quadrature atom while retaining `rho_J<1` as the acceptance criterion.
Old/new output equality is not required. Benchmark checksums must be stable
within each executable; numerical accuracy is checked against independent
roots. The component bounds, numerical results, and timing provenance are in
the [degree-reduction record](experimental-rho-2026-10-04.json).

For the reduced divided-difference polynomial, exact stored-coefficient tail
bounds and source-order rounding bounds add less than `u/32` to the absolute
error budget for `D`. Transport through the compensated residual and finite
Householder correction adds less than `1/32` to `rho_J`. The inherited 5,474
continuous seed leaves then give conditional bounds below `0.571` for wing,
`0.719` for small-rank, and `0.729` for finite-rank. These are component bounds
under the existing primitive and seed-certificate assumptions, not a new
whole-Experimental or end-to-end Lean proof.

`scripts/verify_dpoly_degree.py --reference-root PATH` checks the added error
budget and binds the inherited certificates in an explicitly supplied
`implied-black-volatility` checkout. Both verifiers reject changed source or
certificate fingerprints and Python's assertion-disabling `-O` mode. The
divided-difference checker audits every seed leaf's prerequisites; it does not
regenerate the inherited seed and Mills certificates. Coefficient fitting is
optional and uses `scripts/generate_as1_atan.py` with `mpmath`.

The combined degree reductions passed 49,936 independent normal-price core
inputs and 68,273 normal-price paper/adapter inputs under both FMA policies,
with finite positive outputs and strict ratios below one in each respective
metric. The largest observed core ratio was `0.604894`; the paper/adapter ratio
was `0.773289`, or `0.869034` after annualization and retotalization. Outputs
happened to match the baseline on this finite population. A further 82 inputs
outside the normal-price interior retained their existing outputs and
classifications, including non-finite results; this does not expand the core
contract. Twenty independent two-precision exact-input roots covering the AS1
quadrature rules, tiny scaling, and each reduced divided-difference degree
with adjacent prices are retained in `tests/experimental_polynomials.rs`.

## Algorithm references

- Peter Jäckel, [Let's Be Rational](https://www.jaeckel.org/LetsBeRational.pdf).
- Peter Jäckel, [Implied Normal Volatility](https://www.jaeckel.org/ImpliedNormalVolatility.pdf).
- [An Explicit Solution to Black-Scholes Implied Volatility](https://arxiv.org/abs/2604.24480),
  the separate inverse-Gaussian representation.
- [FlashIV](https://arxiv.org/abs/2605.29102v1), log-price Householder inversion.
- Li and Lee, [equation 42 and seed coefficients](https://mpra.ub.uni-muenchen.de/6867/1/MPRA_paper_6867.pdf).

The default build has no required runtime dependencies. Built-in inverse kernels
allocate no memory. The optional C++ comparison feature adds its own dependencies.
