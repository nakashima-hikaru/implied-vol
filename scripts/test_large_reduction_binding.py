#!/usr/bin/env python3
"""Reject source defects in the fast/precise large-route certificates.

Requires the reviewed reference checkout and a freshly generated fast component
certificate. A complete positive control runs first against an isolated copy;
negative cases must then fail at the intended source binding, not merely fail.
These controls test checker integrity, not the mathematical inequalities.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

if not __debug__:
    raise RuntimeError('Assertions are required; do not use Python -O')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--fast-certificate', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    source = args.source_root.resolve()
    current = (source/'src/experimental/large.rs').read_text()
    controls = [
        ('table_degree_17', 'if i <= 16 { ca[i] }', 'if i <= 17 { ca[i] }', 'source bridge pattern'),
        ('mixed_table_missing_zero_padding', 'if i <= 16 { ca[i] } else { 0.0 }', 'ca[i]', 'source bridge pattern'),
        ('mixed_asymptotic_degree_16', 'if A && B { 16 } else { 18 }', 'if A || B { 16 } else { 18 }', 'source bridge pattern'),
        ('missing_compensated_stage', 'for i in (0..3).rev()', 'for i in (0..2).rev()', 'unmodeled compensated Mills arithmetic'),
        ('physical_wrong_coefficient', 'fma(xn, g3, -nn * g2)', 'fma(xn, g3, nn * g2)', 'unmodeled physical predictor'),
        ('physical_order_four', 'n * c4 / c5', 'n * c3 / c4', 'unmodeled physical predictor'),
        ('physical_first_threshold', 'if xn.abs() <= 0.125', 'if xn.abs() <= 0.25', 'unmodeled physical predictor'),
        ('physical_second_threshold', 'if xn.abs() <= 0.375', 'if xn.abs() <= 0.5', 'unmodeled physical predictor'),
        ('physical_first_denominator', 'return n * c2 / c3;', 'return n * c2 / c2;', 'unmodeled physical predictor'),
        ('physical_second_denominator', 'return n * c3 / c4;', 'return n * c3 / c3;', 'unmodeled physical predictor'),
        ('physical_unfused_helper', '    a.mul_add(b, c)\n', '    a * b + c\n', 'inherited source order'),
        ('log_third_derivative_coefficient', '(a3 / 6.0).mul_add(n, a2)', '(a3 / 3.0).mul_add(n, a2)', 'unmodeled log iteration'),
        ('log_iteration_missing_newton_factor', 'n * num / den', 'num / den', 'unmodeled log iteration'),
        ('missing_forced_bisection', 'steps % 4 == 3 ||', 'steps % 400 == 399 ||', 'unreviewed target/routing/iteration source'),
        ('unchecked_predictor', 'next.is_finite() && next > left && next < right', 'next.is_finite()', 'unreviewed target/routing/iteration source'),
        ('fast_scalar_wrong_degree', 'let mut p = c[13];', 'let mut p = c[12];', 'fast source binding'),
        ('fast_pair_wrong_degree', 'let mut p = coef(13);', 'let mut p = coef(12);', 'fast source binding'),
        ('fast_pair_wrong_lane', 'if A { ca[i] } else { ASYM[i] }', 'if A { cb[i] } else { ASYM[i] }', 'fast source binding'),
    ]
    report = {'positive_control_passed': False, 'negative_controls': []}
    with tempfile.TemporaryDirectory(prefix='iv-large-binding-') as directory:
        isolated = Path(directory)
        for name in ('Cargo.toml', 'src/experimental/large.rs',
                     'src/experimental/large_constants.rs',
                     'src/experimental/mills_constants.rs', 'src/experimental/lanes.rs'):
            dest = isolated/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source/name, dest)
        command = [sys.executable, str(source/'scripts/verify_large_direct.py'),
                   '--source-root', str(isolated),
                   '--reference-root', str(args.reference_root.resolve()),
                   '--fast-certificate', str(args.fast_certificate.resolve()),
                   '--out', str(isolated/'proof.json')]
        positive = subprocess.run(command, capture_output=True, text=True)
        assert positive.returncode == 0, positive.stderr
        report['positive_control_passed'] = True
        for name, before, after, reason in controls:
            assert current.count(before) == 1, (name, 'ambiguous mutation')
            (isolated/'src/experimental/large.rs').write_text(current.replace(before, after, 1))
            result = subprocess.run(command, capture_output=True, text=True)
            assert result.returncode != 0 and reason in result.stderr, (name, result.stdout, result.stderr)
            if name.startswith('fast_'):
                # Both the component checker and its consumer must reject the
                # changed source; a stale numerical certificate is insufficient.
                fast_command = [sys.executable, str(source/'scripts/verify_large_fast.py'),
                                '--source-root', str(isolated),
                                '--reference-root', str(args.reference_root.resolve()),
                                '--out', str(isolated/'fast-proof.json')]
                fast_result = subprocess.run(fast_command, capture_output=True, text=True)
                assert fast_result.returncode != 0 and 'source function binding' in fast_result.stderr, (name, fast_result.stdout, fast_result.stderr)
            report['negative_controls'].append({'name': name, 'rejected': True, 'expected_binding': reason})
            print(name+': rejected at '+reason)
    report['all_rejected'] = all(c['rejected'] for c in report['negative_controls'])
    path = args.out or source/'target/proof/large-reduction-negative-controls.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    print('Positive control passed; '+str(len(controls))+' source defects rejected -> '+str(path))


if __name__ == '__main__':
    main()
