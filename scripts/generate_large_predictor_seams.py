#!/usr/bin/env python3
"""Independent exact-input fixtures around adaptive upper predictor choices.

The ideal x*n thresholds locate physical prices. The saved inputs are those
prices rounded to binary64, neighboring quotes, and interior offsets on both
sides. References invert the original represented price using the direct
Gaussian-erfc formula, independently of the candidate predictor. Root high/low
and eta_J must agree at 160 and 240 decimal digits. Runtime branch counts are
separate empirical evidence: ideal threshold locations are not source bits.
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
        centers = []
        for a0 in [math.nextafter(10.0, math.inf), 10.5, 20.0, 40.0, 100.0, 1000.0, 1400.0]:
            a = mp.mpf(a0)
            for cutoff in [mp.mpf(1)/8, mp.mpf(3)/8]:
                # For larger a, the 3/8 threshold usually lies closer to the
                # cap than binary64 prices can represent. Retain the ordinary
                # physical seam inputs instead of silently rounding to the cap.
                if cutoff > mp.mpf(1)/8 and a0 > 40:
                    continue
                def measure(x):
                    y = mp.sqrt(x*x+2*a)
                    mills = mp.sqrt(mp.pi/2)*mp.exp(y*y/2)*mp.erfc(y/mp.sqrt(2))
                    return x*mills*y/(x+y)
                guess = mp.sqrt(2*a)*cutoff/mp.sqrt(1-2*cutoff)
                seam = mp.findroot(lambda x: measure(x)-cutoff, (guess, guess*mp.mpf('1.1')))
                assert seam > 0 and abs(measure(seam)-cutoff) < mp.mpf('1e-200')
                for side in [-1, 0, 1]:
                    x = seam*(1+side*mp.mpf(2)**-20)
                    price = mp.exp(-a/2)*(1-mp.erfc(x/mp.sqrt(2))/2)
                    centers.append((a0, float(price)))
        # Same represented-input construction as the path benchmarks.
        centers.extend([(1000.0, math.exp(-500.0)*0.9),
                        (10.5, math.exp(-5.25)*(1.0-1e-12))])
        for a, center in centers:
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
        '//! Adaptive upper physical-predictor threshold and interior fixtures.',
        '//! Independent Gaussian-erfc roots and eta_J agree at 160/240 digits.',
        '//! Regenerate with scripts/generate_large_predictor_seams.py.',
        '#![cfg(feature = "experimental")]',
        'use implied_vol::solver::{BlackSolver, Experimental};',
        '', '// a, b, root high, root low, eta_J.',
        '#[rustfmt::skip]',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    lines += ['    ['+', '.join(f'0x{v:016x}' for v in row)+'],' for row in references]
    lines += [
        '];', '', '#[test]',
        'fn adaptive_physical_predictor_retains_attainable_precision() {',
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
        '        "adaptive physical predictor: {} cases, maximum rho_J={maximum_rho}",',
        '        REFERENCES.len()',
        '    );',
        '}', '',
    ]
    args.output.write_text('\n'.join(lines))
    raw = b''.join(struct.pack('<dd', a, b) for a, b in cases)
    report = {'cases': len(cases), 'precision_digits': [160, 240],
              'input_sha256': hashlib.sha256(raw).hexdigest(),
              'fixture_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
              'scope': 'Finite exact-input whole-solver accuracy near ideal thresholds and represented quote neighbors; no one-step or whole-domain guarantee'}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
