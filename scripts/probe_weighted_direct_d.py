#!/usr/bin/env python3
"""Finite, ideal-arithmetic feasibility probe for a weighted direct-D fit.

Fits g=D/(1+D) on (h,t^2) cells, then reconstructs D=g/(1-g).
The exact identity |Dhat-D|/(1+D)=|ghat-g|/|1-ghat| makes the final-IV
weight part of the fitted function. A coefficient l1 bound proves denominator
positivity for the represented fit, but sampled approximation errors are NOT
a continuous accuracy certificate or binary64/FMA/final-solver validation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import struct

import mpmath as mp

CELLS = [('near_atm', '0', '.5', '0', '.25'),
         ('middle', '2', '2.5', '0', '1'),
         ('wing', '7.5', '8', '0', '1')]
DEGREES = [(6, 4), (8, 6), (10, 8), (12, 10)]


def mills(x):
    return mp.sqrt(mp.pi/2)*mp.exp(x*x/2)*mp.erfc(x/mp.sqrt(2))


def direct_d(h, T):
    if T == 0:
        return 1-h*mills(h)
    t = mp.sqrt(T)
    return (mills(h-t)-mills(h+t))/(2*t)


def chebyshev_values(degree, x):
    values = [mp.mpf(1)]
    if degree:
        values.append(x)
    for _ in range(2, degree+1):
        values.append(2*x*values[-1]-values[-2])
    return values


def fit(cell, dx, dy, transformed):
    _, hl, hh, tl, th = cell
    tx = [mp.pi*(i+mp.mpf('.5'))/(dx+1) for i in range(dx+1)]
    ty = [mp.pi*(j+mp.mpf('.5'))/(dy+1) for j in range(dy+1)]
    values = []
    for a in tx:
        row = []
        h = (hl+hh)/2+(hh-hl)*mp.cos(a)/2
        for b in ty:
            T = (tl+th)/2+(th-tl)*mp.cos(b)/2
            D = direct_d(h, T)
            row.append(D/(1+D) if transformed else D)
        values.append(row)
    cx = [[mp.cos(i*a) for a in tx] for i in range(dx+1)]
    cy = [[mp.cos(j*b) for b in ty] for j in range(dy+1)]
    # Two independent DCT passes avoid an unnecessary four-index reduction.
    left = [[sum(values[k][l]*cx[i][k] for k in range(dx+1))
             * (1 if i == 0 else 2)/(dx+1) for l in range(dy+1)] for i in range(dx+1)]
    return [[mp.mpf(float(sum(left[i][l]*cy[j][l] for l in range(dy+1))
                         * (1 if j == 0 else 2)/(dy+1)))
             for j in range(dy+1)] for i in range(dx+1)]


def evaluate(coefficients, x, y):
    xs = chebyshev_values(len(coefficients)-1, x)
    ys = chebyshev_values(len(coefficients[0])-1, y)
    return sum(c*xs[i]*ys[j] for i, row in enumerate(coefficients) for j, c in enumerate(row))


def probe(dps, grid):
    mp.mp.dps = dps
    unit = mp.mpf(2)**-53
    records = []
    for raw in CELLS:
        cell = (raw[0], *(mp.mpf(x) for x in raw[1:]))
        label, hl, hh, tl, th = cell
        references = []
        for i in range(grid):
            x = 2*mp.mpf(i)/(grid-1)-1
            for j in range(grid):
                y = 2*mp.mpf(j)/(grid-1)-1
                h, T = (hl+hh)/2+(hh-hl)*x/2, (tl+th)/2+(th-tl)*y/2
                references.append((x, y, direct_d(h, T)))
        for dx, dy in DEGREES:
            for transformed in (False, True):
                coefficients = fit(cell, dx, dy, transformed)
                norm = sum(abs(F(float(c))) for row in coefficients for c in row)
                denominator = 1-norm if transformed else F(1)
                weighted, relative = mp.mpf(0), mp.mpf(0)
                for x, y, true in references:
                    value = evaluate(coefficients, x, y)
                    estimate = value/(1-value) if transformed else value
                    error = abs(estimate-true)
                    weighted = max(weighted, error/(1+true)/unit)
                    relative = max(relative, error/true/unit)
                records.append({
                    'cell': label, 'h_interval': [str(hl), str(hh)], 't_squared_interval': [str(tl), str(th)],
                    'fit': 'D/(1+D), reconstructed by a positive-denominator rational' if transformed else 'direct D polynomial',
                    'degree_h': dx, 'degree_t_squared': dy,
                    'coefficients': (dx+1)*(dy+1), 'coefficient_bytes': 8*(dx+1)*(dy+1),
                    'sampled_max_weighted_error_in_u': float(weighted),
                    'sampled_max_relative_D_error_in_u': float(relative),
                    'continuous_exact_denominator_lower_bound': str(denominator) if denominator > 0 else None,
                    'coefficient_bits': [[f'0x{struct.unpack("<Q", struct.pack("<d", float(c)))[0]:016x}' for c in row] for row in coefficients],
                })
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps', type=int, default=90)
    parser.add_argument('--grid', type=int, default=17)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 80 or args.grid < 3:
        parser.error('use --dps>=80 and --grid>=3')
    records = probe(args.dps, args.grid)
    result = {'scope': 'Finite ideal-arithmetic research diagnostic; no runtime replacement adopted',
              'fit_objective': 'Fit D/(1+D); reconstruct D so the exact weighted-error identity applies',
              'precision_digits': args.dps, 'held_out_grid': [args.grid, args.grid], 'records': records,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'limits': 'Chebyshev interpolation, not a weighted minimax optimum. Denominator positivity alone does not certify approximation error. Model remainder, binary64/FMA source rounding, final rho_J, whole-domain coverage, and paired timing remain unqualified.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(args.output), 'fits': len(records),
                      'best_sampled_weighted_error_in_u_by_cell': {
                          label: min(row['sampled_max_weighted_error_in_u'] for row in records if row['cell'] == label)
                          for label, *_ in CELLS}}, indent=2))


if __name__ == '__main__':
    main()
