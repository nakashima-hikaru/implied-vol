#!/usr/bin/env python3
"""Fit the D(h,0) coefficient rows for the direct D recurrence; requires mpmath.

The output is coefficient data, not a precision certificate. The stdlib-only
verify_direct_d.py checks the represented production table independently of
these fitting diagnostics and of the fitter's numerical moment evaluations.
"""
import argparse
import json
import struct
from pathlib import Path

import mpmath as mp


def d0(h):
    """I_1(h), evaluated at high precision before binary64 rounding."""
    mills = mp.sqrt(mp.pi / 2) * mp.exp(h * h / 2) * mp.erfc(h / mp.sqrt(2))
    return 1 - h * mills


def add(a, b):
    result = [mp.mpf(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        result[i] += value
    for i, value in enumerate(b):
        result[i] += value
    return result


def fit(dps):
    mp.mp.dps = dps
    nodes = 24
    target = mp.mpf(2) ** -62 / 8
    basis = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(4)]]
    for _ in range(2, nodes):
        basis.append(add([mp.mpf(0)] + [8 * v for v in basis[-1]], [-v for v in basis[-2]]))
    theta = [mp.pi * (i + mp.mpf('.5')) / nodes for i in range(nodes)]
    cells = []
    for cell in range(17):
        center = mp.mpf(2 * cell + 1) / 4
        values = [d0(center + mp.cos(q) / 4) for q in theta]
        cheb = [sum(v * mp.cos(n * q) for v, q in zip(values, theta)) * 2 / nodes for n in range(nodes)]
        cheb[0] /= 2
        degree = next(n for n in range(nodes) if sum(abs(v) for v in cheb[n + 1:]) < target)
        assert degree < 16
        coefficients = [mp.mpf(0)] * (degree + 1)
        for n in range(degree + 1):
            for i, value in enumerate(basis[n]):
                coefficients[i] += value * cheb[n]
        highs = [float(value) for value in coefficients]
        lows = [float(coefficients[i] - mp.mpf(highs[i])) for i in range(3)]
        cells.append({
            'degree': degree,
            'high_bits': [bits(value) for value in highs] + ['0x0000000000000000'] * (16 - len(highs)),
            'low_bits': [bits(value) for value in lows],
            'chebyshev_tail_estimate': mp.nstr(sum(abs(v) for v in cheb[degree + 1:]), 20),
        })
    return {'dps': dps, 'cells': cells, 'scope': 'D0 fitting diagnostics only; run verify_direct_d.py for the conditional source error certificate.'}


def bits(value):
    return f"0x{struct.unpack('>Q', struct.pack('>d', value))[0]:016x}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps', type=int, default=90)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 80:
        parser.error('use at least 80 decimal digits for moment recurrence')
    result = fit(args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(f'Saved 17 D0 coefficient rows to {args.output}')


if __name__ == '__main__':
    main()
