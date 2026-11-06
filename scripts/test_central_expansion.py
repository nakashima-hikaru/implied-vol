#!/usr/bin/env python3
"""Negative controls for the central domain/degree precision certificate."""
from pathlib import Path
import argparse
import sys
import tempfile
import unittest
from unittest.mock import patch

from central_expansion import BASELINE_SHA256, restore_central_source
from verify_central_expansion import digest, source_binding

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = None
FILES = ['src/experimental/small.rs', 'src/experimental/small/central.rs',
         'src/experimental/small/constants.rs', 'src/experimental/lanes.rs']


class SourceMutationControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        for name in FILES:
            path = self.root/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((ROOT/name).read_bytes())
        # Only the original coefficient file is required by the exact
        # representation bridge; no numerical certificate is consumed here.
        reference = self.root/'src/small/constants.rs'
        reference.parent.mkdir(parents=True)
        reference.write_bytes((REFERENCE/'src/small/constants.rs').read_bytes())

    def tearDown(self):
        self.temporary.cleanup()

    def reject(self, name, old, new):
        path = self.root/name
        text = path.read_text()
        self.assertEqual(text.count(old), 1)
        path.write_text(text.replace(old, new))
        with self.assertRaises(AssertionError):
            source_binding(self.root, self.root)

    def test_unmodified_binding_passes(self):
        source_binding(self.root, self.root)
        central = self.root/'src/experimental/small/central.rs'
        restored = restore_central_source(central.read_text())
        central.write_text(restored)
        self.assertEqual(digest(central), BASELINE_SHA256)

    def test_larger_q_guard_rejected(self):
        self.reject(FILES[1], 'q.hi <= 0.625', 'q.hi <= 0.75')

    def test_changed_preflight_rejected(self):
        self.reject(FILES[1], 'b + 0.5 * a >= 0.25', 'b + 0.5 * a >= 0.3')

    def test_changed_pair_arithmetic_rejected(self):
        self.reject(FILES[1], 'a.lo * b.hi)', 'a.lo * b.hi + a.lo * b.lo)')

    def test_changed_conversion_rounding_rejected(self):
        self.reject(FILES[1], 'y = fma(y, a2.hi, rows[i]);', 'y = y * a2.hi + rows[i];')

    def test_changed_conversion_retained_degree_rejected(self):
        self.reject(FILES[1], 'conversion_pair::<16, 15, 12, 11>',
                    'conversion_pair::<16, 15, 11, 10>')

    def test_changed_conversion_proof_lengths_rejected(self):
        with patch('verify_central_expansion.CONVERSION_LENGTHS', (16, 15, 14, 13, 11, 10, 1, 0, 0, 0)):
            with self.assertRaises(AssertionError):
                source_binding(self.root, self.root)

    def test_changed_conversion_final_row_proof_rejected(self):
        with patch('verify_central_expansion.CONVERSION_LENGTHS', (12, 11, 9, 7, 5, 3, 0, 0, 0, 0)):
            with self.assertRaises(AssertionError):
                source_binding(self.root, self.root)

    def test_changed_conversion_final_coefficient_rejected(self):
        self.reject(FILES[1], 'let mut y = conversion_6[7];',
                    'let mut y = conversion_6[6];')

    def test_changed_conversion_outer_order_rejected(self):
        self.reject(FILES[1], 'for i in (0..6).rev() {',
                    'for i in (0..5).rev() {')

    def test_changed_conversion_coefficient_offset_rejected(self):
        self.reject(FILES[1], 'let right_start = RIGHT - RIGHT_KEEP;',
                    'let right_start = RIGHT - RIGHT_KEEP - 1;')

    def test_changed_normal_degree_rejected(self):
        self.reject(FILES[1], '0 => normal_polynomial::<18>(k, t),',
                    '0 => normal_polynomial::<17>(k, t),')

    def test_changed_proof_degree_rejected(self):
        with patch('verify_central_expansion.NORMAL_DEGREES', (25, 25, 25, 25)):
            with self.assertRaises(AssertionError):
                source_binding(self.root, self.root)

    def test_changed_normal_branch_rejected(self):
        self.reject(FILES[1], '1 | 2 => normal_polynomial::<20>(k, t),',
                    '1 | 2 => normal_polynomial::<18>(k, t),')

    def test_changed_normal_horner_order_rejected(self):
        self.reject(FILES[1], 'for i in (2..DEGREE).rev() {',
                    'for i in (2..DEGREE) {')

    def test_missing_normal_low_correction_rejected(self):
        self.reject(FILES[1], 'Pair::new(p.hi, fma(dp, t.lo, p.lo))',
                    'Pair::new(p.hi, p.lo)')

    def test_changed_normal_compensation_rejected(self):
        self.reject(FILES[1], 'let residual = prod - (s - NC[k][i]);',
                    'let residual = 0.0;')

    def test_changed_outer_guard_rejected(self):
        self.reject(FILES[0], 'if a <= 0.45 && b < 0.26 {', 'if a <= 0.5 && b < 0.26 {')

    def test_changed_pair_constructor_rejected(self):
        self.reject(FILES[0], 'Self { hi, lo: 0.0 }', 'Self { hi: hi * 2.0, lo: 0.0 }')

    def test_changed_normal_coefficient_rejected(self):
        path = self.root/FILES[2]
        text = path.read_text()
        start = text.index('static NC:')
        pos = text.index('f64::from_bits(', start) + len('f64::from_bits(')
        text = text[:pos] + '0x3ff0000000000000' + text[text.index(')', pos):]
        path.write_text(text)
        with self.assertRaises(AssertionError):
            source_binding(self.root, self.root)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-root', type=Path, required=True)
    args = parser.parse_args()
    REFERENCE = args.reference_root.resolve()
    unittest.main(argv=[sys.argv[0]])
