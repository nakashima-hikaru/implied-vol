# Weighted direct-D feasibility probe

The production direct-D moment recurrence remains unchanged. The new
[probe](../scripts/probe_weighted_direct_d.py) makes the proposal's
`|Dhat-D|/(1+D)` weight part of approximation generation by fitting

```text
g = D/(1+D), Dhat = ghat/(1-ghat).
```

The exact identity

```text
|Dhat-D|/(1+D) = |ghat-g|/|1-ghat|
```

connects that transformed fit to the requested objective. The probe compares
it with ordinary direct-D interpolation on three `(h,t²)` cells, four degree
pairs and an independent 17-by-17 uniform evaluation grid. Coefficients are
rounded once to binary64 and evaluated in high precision. An exact rational
coefficient 1-norm also gives a continuous positive denominator bound for each
represented transformed fit. This denominator result does not certify its
approximation error.

At 90 decimal digits, representative sampled errors are as follows. Here
`u=2^-53`, and all error figures are `|Dhat-D|/((1+D)*u)`.

| Cell | h/t² degrees | Ordinary D | Transformed g | Coefficient bytes per fit |
|---|---:|---:|---:|---:|
| h in [0,.5], t² in [0,.25] | 12/10 | .4851 | .8355 | 1,144 |
| h in [2,2.5], t² in [0,1] | 10/8 | 1.8615 | .1760 | 792 |
| h in [7.5,8], t² in [0,1] | 8/6 | .1190 | .0800 | 504 |

The transformed objective helps some cells and worsens the highest-degree
near-ATM example. These finite results do not establish weighted minimax
optimality, continuous approximation remainder, binary64 source-order error,
final `rho_J`, full dispatch coverage or improved runtime. The coefficient
payloads and additional rational division also need comparison with the
existing single fitted D0 row and short moment recurrence. Consequently this
candidate is retained as a reproducible research diagnostic, without replacing
the qualified runtime implementation.

```sh
rtk proxy .venv/bin/python scripts/probe_weighted_direct_d.py --output /tmp/weighted-direct-d-probe.json
```

The report includes every represented coefficient and both favorable and
unfavorable fits. The approximation generator uses tensor Chebyshev
interpolation; it does not claim a weighted minimax solution.
