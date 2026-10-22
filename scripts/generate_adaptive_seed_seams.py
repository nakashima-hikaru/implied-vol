#!/usr/bin/env python3
"""Independent Black inversions around small/core seed rectangle boundaries.

Root high/low and eta_J must agree at 160 and 240 decimal digits. An optional
binary-pair replay executable supplies starting guesses only; independent
high-precision residual checks validate every root. Without it, use bisection
to obtain starting guesses. This finite corpus is not a continuous proof.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess

import mpmath as mp

if not __debug__:
    raise RuntimeError('Reference generation requires assertions; do not use python -O.')


def bits(value):
    return struct.unpack('<Q', struct.pack('<d', float(value)))[0]


def adjacent(value):
    return math.nextafter(value, -math.inf), value, math.nextafter(value, math.inf)


def inputs():
    seen = set()
    centers = []
    with mp.workdps(240):
        for q in ('0.5', '1', '1.5', '2'):
            for w in ('0.25', '2', '3.5', '5', '6.75', '8.5'):
                b = mp.mpf(q) / (mp.sqrt(2*mp.pi) * mp.exp(mp.mpf(w)))
                a = 2*b*mp.expm1(mp.mpf(w))
                centers.append(('small', float(a), float(b)))
        for a in ('0.125', '0.25', '0.5', '1', '2', '4', '7', '10'):
            for z in ('0.35', '0.7', '1.2', '2', '3', '4.5', '6', '8', '10', '11.5', '13'):
                centers.append(('core', float(a), float(mp.exp(-mp.mpf(a)/2-mp.mpf(z)))))
        for label, a0, b0 in centers:
            for a in adjacent(a0):
                for b in adjacent(b0):
                    key = bits(a), bits(b)
                    if key in seen or not (0 < a <= 10 and 0 < mp.mpf(b) < mp.exp(-mp.mpf(a)/2)):
                        continue
                    seen.add(key)
                    yield label, a, b


def invert(a, b, initial, dps):
    with mp.workdps(dps):
        a, b = mp.mpf(a), mp.mpf(b)
        sqrt2, sqrt2pi = mp.sqrt(2), mp.sqrt(2*mp.pi)

        def evaluate(s):
            h, t = a/s, s/2
            price = (mp.exp(-a/2)*mp.erfc((h-t)/sqrt2)
                     - mp.exp(a/2)*mp.erfc((h+t)/sqrt2))/2
            vega = mp.exp(-(h*h+t*t)/2)/sqrt2pi
            return price, vega

        left, right = mp.mpf(0), mp.sqrt(2*a)+4
        while evaluate(right)[0] <= b:
            right *= 2
        if initial is not None and math.isfinite(initial) and 0 < initial < right:
            s = mp.mpf(initial)
        else:
            for _ in range(20):
                s = (left+right)/2
                if evaluate(s)[0] < b:
                    left = s
                else:
                    right = s
            s = (left+right)/2
        tolerance = mp.power(10, -dps+25)
        for _ in range(40):
            price, vega = evaluate(s)
            correction = (price-b)/vega
            if abs(correction/s) < tolerance:
                break
            if price < b:
                left = s
            else:
                right = s
            candidate = s-correction
            s = candidate if left < candidate < right else (left+right)/2
        else:
            raise ArithmeticError(('root did not converge', a, b, dps))
        price, vega = evaluate(s)
        assert abs((price-b)/(s*vega)) < tolerance
        eta = mp.mpf(2)**-52 * (1+b/(s*vega))
        hi = float(s)
        return bits(hi), bits(s-hi), bits(eta)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--initial-executable', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    cases = list(inputs())
    raw = b''.join(struct.pack('<dd', a, b) for _, a, b in cases)
    guesses = [None]*len(cases)
    if args.initial_executable:
        output = subprocess.check_output([str(args.initial_executable.resolve())], input=raw)
        assert len(output) == len(cases)*8
        guesses = [row[0] for row in struct.iter_unpack('<d', output)]
    references = []
    for index, ((label, a, b), guess) in enumerate(zip(cases, guesses)):
        low = invert(a, b, guess, 160)
        high = invert(a, b, guess, 240)
        assert low == high, ('precision disagreement', label, a, b, low, high)
        references.append((bits(a), bits(b), *high))
        if (index+1) % 100 == 0:
            print(f'Confirmed {index+1}/{len(cases)} roots', flush=True)
    lines = [
        '//! Exact binary64 neighbors of small/core seed rectangle boundaries.',
        '//! Independent Black root high/low and eta agree at 160/240 digits.',
        '//! Regenerate with scripts/generate_adaptive_seed_seams.py.',
        '#![cfg(feature = "experimental")]',
        'use implied_vol::solver::{BlackSolver, Experimental};',
        '',
        '// a, b, root high, root low, eta_J.',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    for row in references:
        lines.append('    [')
        for value in row:
            digits = f'{value:016x}'
            lines.append('        0x'+'_'.join(digits[i:i+4] for i in range(0, 16, 4))+',')
        lines.append('    ],')
    lines += [
        '];', '', '#[test]',
        'fn adaptive_seed_boundaries_retain_attainable_precision() {',
        '    for row in REFERENCES {',
        '        let [a, b, hi, lo, eta] = row.map(f64::from_bits);',
        '        let s = Experimental::implied_total_volatility(-a, b).unwrap();',
        '        let rho = (((s - hi) - lo).abs() / hi) / eta;',
        '        assert!(s.is_finite() && rho < 1.0, "a={a:e}, b={b:e}, rho_J={rho}");',
        '    }', '}', '',
    ]
    args.output.write_text('\n'.join(lines))
    report = {'cases': len(references), 'precision_digits': [160, 240],
              'input_sha256': hashlib.sha256(raw).hexdigest(),
              'fixture_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
              'scope': 'Physical binary64 price inputs; boundary neighborhoods and interiors, not a continuous proof',
              'initial_guess_executable': str(args.initial_executable) if args.initial_executable else None}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
