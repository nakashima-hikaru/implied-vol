#!/usr/bin/env python3
"""Regenerate or check the exact degree-8/6 Chebyshev rank seeds."""
import argparse
import json
from pathlib import Path
import seed_reduction as seeds

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
parser.add_argument('--reference-root', type=Path, required=True)
parser.add_argument('--check', action='store_true')
parser.add_argument('--report', type=Path)
args = parser.parse_args()
if not __debug__:
    raise RuntimeError('Exact source checks require assertions; do not use -O')
source = args.source_root.resolve()
small_path = source/'src/experimental/small.rs'
constants_path = source/'src/experimental/small/constants.rs'
small, constants = small_path.read_text(), constants_path.read_text()
reference = (args.reference_root/'src/small/constants.rs').read_text()
assert seeds.sha256(reference) == seeds.REFERENCE_CONSTANTS_SHA256
if seeds.sha256(small) == seeds.BASELINE_SMALL_SHA256:
    reduced_small = seeds.reduce_small(small)
    reduced_constants, metadata = seeds.reduced_constants(constants)
else:
    small, constants, metadata = seeds.restore_and_validate(small, constants, reference)
    reduced_small, (reduced_constants, _) = seeds.reduce_small(small), seeds.reduced_constants(constants)
if args.check:
    assert small_path.read_text() == reduced_small
    assert constants_path.read_text() == reduced_constants
else:
    small_path.write_text(reduced_small)
    constants_path.write_text(reduced_constants)
report = {'method': 'Exact power-to-Chebyshev transform; truncate modes; exact inverse transform; one binary64 coefficient rounding',
          'x_degree': 8, 'y_degree': 6, 'rank_tables': len(metadata),
          'old_bytes': sum(v['old_bytes'] for v in metadata.values()),
          'new_bytes': sum(v['new_bytes'] for v in metadata.values()),
          'reference_constants_sha256': seeds.REFERENCE_CONSTANTS_SHA256, 'tables': metadata}
if args.report:
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({k: v for k, v in report.items() if k != 'tables'}, indent=2))
