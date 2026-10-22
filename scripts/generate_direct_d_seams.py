#!/usr/bin/env python3
"""Regenerate final-root regression cases for the h-adaptive direct-D selector.

Requires mpmath. This independently inverts the Black price at 160 and 240
decimal digits using exact binary64 a and b inputs. It is a regression corpus,
not a continuous precision certificate; run verify_direct_d.py for that.
"""
import argparse
import math
import struct
from pathlib import Path

import mpmath as mp


def bits(value):
    return struct.unpack('>Q', struct.pack('>d', float(value)))[0]


def solve(a, b, initial, dps):
    with mp.workdps(dps):
        a, b, s = mp.mpf(a), mp.mpf(b), mp.mpf(initial)
        sqrt2 = mp.sqrt(2)
        sqrt2pi = mp.sqrt(2 * mp.pi)
        for _ in range(8):
            h, t = a / s, s / 2
            price = (mp.exp(-a / 2) * mp.erfc((h - t) / sqrt2)
                     - mp.exp(a / 2) * mp.erfc((h + t) / sqrt2)) / 2
            vega = mp.exp(-(h * h + t * t) / 2) / sqrt2pi
            s -= (price - b) / vega
        h, t = a / s, s / 2
        vega = mp.exp(-(h * h + t * t) / 2) / sqrt2pi
        eta = mp.mpf(2) ** -52 * (1 + b / (s * vega))
        hi = float(s)
        return (bits(hi), bits(s - hi), bits(eta))


def cases():
    # Each h cutoff in both applicable t regimes, plus t boundaries.
    geometry = [(5., t) for t in (.8, .9, 1.)]
    geometry += [(6., t) for t in (.55, .625, .75, .8)]
    geometry += [(7., t) for t in (.55, .625, .7)]
    # Interior inputs exercise all four reductions even when the seed h differs
    # from the true-root h by much more than one binary64 ulp.
    geometry += [(6.25, .6), (7.25, .6), (5.25, .875), (6.1, .8)]
    for h, t in geometry:
        s = 2 * t
        center = h * s
        for a in (math.nextafter(center, -math.inf), center,
                  math.nextafter(center, math.inf)):
            if not 0 < a <= 10:
                continue
            with mp.workdps(240):
                ah, st = mp.mpf(a) / mp.mpf(s), mp.mpf(s) / 2
                price = float((mp.exp(-mp.mpf(a) / 2) * mp.erfc((ah - st) / mp.sqrt(2))
                               - mp.exp(mp.mpf(a) / 2) * mp.erfc((ah + st) / mp.sqrt(2))) / 2)
            for b in (math.nextafter(price, -math.inf), price,
                      math.nextafter(price, math.inf)):
                low = solve(a, b, s, 160)
                high = solve(a, b, s, 240)
                assert low == high, (a, b, low, high)
                yield (bits(a), bits(b), *high)


def component_cases():
    for cutoff in (5., 6., 7.):
        for h in (math.nextafter(cutoff, -math.inf), cutoff,
                  math.nextafter(cutoff, math.inf)):
            for h_low in (-h * 2**-54, h * 2**-54):
                ts = [v for t in (.5, .75) for v in
                      (math.nextafter(t, -math.inf), t,
                       math.nextafter(t, math.inf))]
                ts += [math.nextafter(1., -math.inf), 1.]
                for t in ts:
                    references = []
                    for dps in (160, 240):
                        with mp.workdps(dps):
                            x, y = mp.mpf(h) + mp.mpf(h_low), mp.mpf(t)
                            def mills(z):
                                return mp.sqrt(mp.pi / 2) * mp.exp(z*z/2) * mp.erfc(z/mp.sqrt(2))
                            value = (mills(x-y) - mills(x+y)) / (2*y)
                            hi = float(value)
                            references.append((bits(hi), bits(value-hi)))
                    assert references[0] == references[1], (h, h_low, t)
                    yield (bits(h), bits(h_low), bits(t), *references[1])


def write_component(path):
    rows = list(component_cases())
    lines = [
        '//! Component references at the exact h.high selector boundaries.',
        '//! Independent Mills differences agree at 160 and 240 decimal digits.',
        '//! The h lows are nonzero; each t cutoff is tested on both sides.',
        '//! Regenerate with `scripts/generate_direct_d_seams.py --component-output`.',
        'use super::{Pair, evaluate};',
        '',
        '// h.high, h.low, t, D high, D low.',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    for row in rows:
        lines.append('    [')
        for v in row:
            digits = f'{v:016x}'
            lines.append('        0x' + '_'.join(digits[i:i+4] for i in range(0,16,4)) + ',')
        lines.append('    ],')
    lines += [
        '];',
        '',
        '#[test]',
        'fn selector_boundaries_retain_the_certified_absolute_budget() {',
        '    for row in REFERENCES {',
        '        let [h, h_low, t, hi, lo] = row.map(f64::from_bits);',
        '        let d = evaluate(Pair::new(h, h_low), Pair::from(t));',
        '        let error = ((d.hi - hi) + (d.lo - lo)).abs();',
        '        assert!(',
        '            d.hi.is_finite() && error < f64::EPSILON / 2.0,',
        '            "h={h:e}, h_low={h_low:e}, t={t:e}, error={error:e}"',
        '        );',
        '    }',
        '}',
        '',
    ]
    path.write_text('\n'.join(lines))
    print(f'Saved {len(rows)} independently checked component cases to {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--component-output', type=Path)
    args = parser.parse_args()
    lines = [
        '//! Exact binary64 prices around h=5,6,7 and inside all reduced-count regions.',
        '//! Root high/low and eta bits agree at 160 and 240 decimal digits.',
        '//! Regenerate with `scripts/generate_direct_d_seams.py`.',
        '#![cfg(feature = "experimental")]',
        'use implied_vol::solver::{BlackSolver, Experimental};',
        '',
        '// a, b, root high, root low, eta_J.',
        'const REFERENCES: &[[u64; 5]] = &[',
    ]
    rows = list(cases())
    for row in rows:
        lines.append('    [')
        for v in row:
            digits = f'{v:016x}'
            lines.append('        0x' + '_'.join(digits[i:i+4] for i in range(0,16,4)) + ',')
        lines.append('    ],')
    lines += [
        '];',
        '',
        '#[test]',
        'fn adaptive_direct_d_root_neighborhoods_retain_attainable_precision() {',
        '    for row in REFERENCES {',
        '        let [a, b, hi, lo, eta] = row.map(f64::from_bits);',
        '        let s = Experimental::implied_total_volatility(-a, b).unwrap();',
        '        let rho = (((s - hi) - lo).abs() / hi) / eta;',
        '        assert!(s.is_finite() && rho < 1.0, "a={a:e}, b={b:e}, rho_J={rho}");',
        '    }',
        '}',
        '',
    ]
    args.output.write_text('\n'.join(lines))
    print(f'Saved {len(rows)} independently checked final-root cases to {args.output}')
    if args.component_output:
        write_component(args.component_output)


if __name__ == '__main__':
    main()
