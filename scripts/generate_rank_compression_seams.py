#!/usr/bin/env python3
"""Independent final-output references inside compressed small rectangle 2.

All 625 physical binary64 inputs reach the rectangle with either the original
central guard or the proposed q<=.625 guard. Existing adaptive-seed fixtures
separately cover exact binary64 neighbors of rectangle seams. Root high/low
and eta_J must agree at 160 and 240 decimal digits. This finite replay does
not replace the continuous represented-coefficient perturbation proof.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct

import mpmath as mp
from generate_adaptive_seed_seams import bits, invert


def inputs():
    with mp.workdps(240):
        for i in range(25):
            w = 2 + 3*(mp.mpf(i)+mp.mpf('.5'))/25
            for j in range(25):
                q = mp.mpf('.65') + mp.mpf('.35')*(mp.mpf(j)+mp.mpf('.5'))/25
                b = q/(mp.sqrt(2*mp.pi)*mp.exp(w))
                a = 2*b*mp.expm1(w)
                yield float(a), float(b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    assert __debug__, 'Reference generation requires assertions'
    cases = list(inputs())
    references = []
    for a, b in cases:
        low, high = invert(a,b,None,160), invert(a,b,None,240)
        assert low == high, ('precision disagreement',a,b,low,high)
        references.append((bits(a),bits(b),*high))
    lines = [
        '//! Binary64 interior grid for compressed small rectangle 2.',
        '//! Independent Black roots and eta_J agree at 160/240 decimal digits.',
        '//! Regenerate with scripts/generate_rank_compression_seams.py.',
        '#![cfg(feature = "experimental")]',
        'use implied_vol::solver::{BlackSolver, Experimental};',
        '', '// a, b, root high, root low, eta_J.',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    for row in references:
        lines.append('    [')
        lines.extend(f'        0x{value:016x},' for value in row)
        lines.append('    ],')
    lines += [
        '];', '', '#[test]',
        'fn compressed_rank_two_retain_attainable_precision() {',
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
        '        "small rectangle 2: {} cases, maximum rho_J={maximum_rho}",',
        '        REFERENCES.len()',
        '    );',
        '}', '',
    ]
    args.output.write_text('\n'.join(lines))
    raw = b''.join(struct.pack('<dd',a,b) for a,b in cases)
    report = {'cases':len(cases), 'precision_digits':[160,240],
              'input_sha256':hashlib.sha256(raw).hexdigest(),
              'fixture_sha256':hashlib.sha256(args.output.read_bytes()).hexdigest(),
              'scope':'Physical binary64 interior grid for small rectangle 2, with q>.65; finite final-output validation'}
    if args.report:
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
