#!/usr/bin/env python3
"""Independent binary64 price fixtures for the shortened large log iteration.

Cover the a=10 and ell=-.6 switches, the observed slowest first handoff near
log relative price -2, and deep lower tails including the minimum normal price.
Each direct Gaussian-erfc root/eta must agree at 160 and 240 decimal digits.
The generator does not invoke or differentiate the candidate iteration.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

import mpmath as mp
from generate_adaptive_seed_seams import adjacent, bits, invert


def inputs():
    seen = set()
    minimum = float.fromhex('0x1p-1022')
    with mp.workdps(240):
        for a in [math.nextafter(10.0, math.inf), 10.125, 20.0, 100.0, 1000.0, 1400.0]:
            centers = [minimum]
            centers += [float(mp.exp(-mp.mpf(a)/2-mp.mpf(q)))
                        for q in ['0.6', '2', '2.05', '20', '100', '500']]
            if a == math.nextafter(10.0, math.inf):
                centers.append(float.fromhex('0x1.c7aaf6968fa7fp-11'))
            for center in centers:
                for b in adjacent(center):
                    key = bits(a), bits(b)
                    if key in seen or not (minimum <= b < mp.exp(-mp.mpf(a)/2)):
                        continue
                    seen.add(key)
                    yield a, b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    assert __debug__, 'Reference generation requires assertions'
    cases = list(inputs())
    references = []
    for a, b in cases:
        low, high = invert(a, b, None, 160), invert(a, b, None, 240)
        assert low == high, ('precision disagreement', a, b, low, high)
        references.append((bits(a), bits(b), *high))
    lines = [
        '//! Lower large-route handoff and tail fixtures for Householder d=3.',
        '//! Independent Gaussian-erfc roots and eta_J agree at 160/240 digits.',
        '//! Regenerate with scripts/generate_large_iteration_cases.py.',
        '#![cfg(feature = "experimental")]',
        'use implied_vol::solver::{BlackSolver, Experimental};',
        '', '// a, b, root high, root low, eta_J.',
        '#[rustfmt::skip]',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    lines += ['    ['+', '.join(f'0x{v:016x}' for v in row)+'],' for row in references]
    lines += [
        '];', '', '#[test]',
        'fn shortened_log_iteration_retains_attainable_precision() {',
        '    let mut maximum_rho = 0.0_f64;',
        '    for row in REFERENCES {',
        '        let [a, b, hi, lo, eta] = row.map(f64::from_bits);',
        '        let s = Experimental::implied_total_volatility(-a, b).unwrap();',
        '        let rho = (((s - hi) - lo).abs() / hi) / eta;',
        '        assert!(',
        '            s.is_finite() && s > 0.0 && rho < 1.0,',
        '            "a={a:e}, b={b:e}, rho_J={rho}"',
        '        );',
        '        maximum_rho = maximum_rho.max(rho);',
        '    }',
        '    eprintln!(',
        '        "large log iteration: {} cases, maximum rho_J={maximum_rho}",',
        '        REFERENCES.len()',
        '    );',
        '}', '',
    ]
    args.output.write_text('\n'.join(lines))
    raw = b''.join(struct.pack('<dd', a, b) for a, b in cases)
    report = {'cases': len(cases), 'precision_digits': [160, 240],
              'input_sha256': hashlib.sha256(raw).hexdigest(),
              'fixture_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
              'scope': 'Finite exact-input whole-solver accuracy; no one-step or whole-domain guarantee'}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
