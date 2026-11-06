#!/usr/bin/env python3
"""Independent Black references beside central guards and shape-degree seams.

Roots and eta_J must agree at 160/240 digits for exact binary64 inputs. This
finite replay complements verify_central_expansion.py's continuous ledger.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import struct

import mpmath as mp

from generate_adaptive_seed_seams import adjacent, bits, invert


def inputs():
    centers = []
    with mp.workdps(240):
        for q in ('0.5', '0.625'):
            for z in ('.0001', '.01', '.16', '.30', '.40', '.5'):
                d = mp.mpf(q)*mp.sqrt(mp.mpf(z))
                a = 2*mp.asinh(d/2)
                b = mp.mpf(q)/mp.sqrt(2*mp.pi)-d/2
                centers.append(('q_' + q, float(a), float(b)))
        for a in (.36, .45):
            for b in (.025, .04, .10):
                centers.append(('a_' + str(a), a, b))
        for cut in (.2, .25):
            for a in (.1, .36):
                centers.append(('preflight_' + str(cut), a, cut-a/2))
        for a in (.1, .4):
            centers.append(('b_outer', a, .26))
        centers.extend([('new_interior', .4, .04), ('new_interior', .3, .075)])
        # Exercise every normal-shape row near both ends of its local
        # coordinate, and at its center, across the central q range. Convert
        # mathematical (q,z) to binary64 (a,b) before solving independently;
        # adjacent inputs also exercise rounded dispatch at the row seams.
        for q in ('.00001', '.01', '.1', '.35', '.5', '.625'):
            for z in ('.000001', '.08', '.159', '.16', '.161', '.23',
                      '.299', '.30', '.301', '.35', '.399', '.40',
                      '.401', '.45', '.499', '.5'):
                d = mp.mpf(q)*mp.sqrt(mp.mpf(z))
                a = 2*mp.asinh(d/2)
                b = mp.mpf(q)/mp.sqrt(2*mp.pi)-d/2
                centers.append(('shape_degree_' + z, float(a), float(b)))
    seen = set()
    for label, a0, b0 in centers:
        for a in adjacent(a0):
            for b in adjacent(b0):
                key = bits(a), bits(b)
                if key in seen:
                    continue
                seen.add(key)
                yield label, a, b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Reference generation requires assertions; do not use -O')
    cases = list(inputs())
    references = []
    for label, a, b in cases:
        low, high = invert(a, b, None, 160), invert(a, b, None, 240)
        assert low == high, ('independent precision disagreement', label, a, b)
        references.append((bits(a), bits(b), *high))
    lines = [
        '//! Exact binary64 neighbors of central guards and normal-shape degree seams.',
        '//! Independent Black roots/eta_J agree at 160/240 decimal digits.',
        '//! Regenerate with scripts/generate_central_expansion_seams.py.',
        '#![cfg(feature = "experimental")]',
        'use implied_vol::solver::{BlackSolver, Experimental};',
        '',
        '// a, b, root high, root low, eta_J.',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    for row in references:
        lines.append('    [')
        lines.extend(f'        0x{v:016x},' for v in row)
        lines.append('    ],')
    lines.extend([
        '];', '', '#[test]',
        'fn central_expansion_seams_retain_attainable_precision() {',
        '    for row in REFERENCES {',
        '        let [a, b, hi, lo, eta] = row.map(f64::from_bits);',
        '        let s = Experimental::implied_total_volatility(-a, b).unwrap();',
        '        let rho = (((s - hi) - lo).abs() / hi) / eta;',
        '        assert!(s.is_finite() && rho < 1.0, "a={a:e}, b={b:e}, rho_J={rho}");',
        '    }', '}', '',
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(lines))
    result = {'scope': 'Finite exact-input Black references, not a continuous proof',
              'cases': len(cases), 'groups': dict(Counter(label for label, _, _ in cases)),
              'precision_digits': [160, 240], 'output': str(args.output),
              'input_sha256': hashlib.sha256(b''.join(struct.pack('<dd', a, b) for _, a, b in cases)).hexdigest(),
              'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
