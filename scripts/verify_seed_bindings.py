#!/usr/bin/env python3
"""Reject deliberate seed metadata, coefficient, shape and arithmetic defects.

Exercise both the shared direct-D bridge and the independent AS1 source
consumer in temporary source trees. No production source is modified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tempfile

import seed_reduction as seeds
import verify_as1_atan as as1


def change(source, old, new):
    assert source.count(old) >= 1, ('missing mutation target',old)
    return source.replace(old,new,1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--reference-root',type=Path,required=True)
    parser.add_argument('--report',type=Path)
    args = parser.parse_args()
    assert __debug__, 'Source verification requires assertions'
    root, reference_root = args.source_root.resolve(),args.reference_root.resolve()
    small = (root/'src/experimental/small.rs').read_text()
    constants = (root/'src/experimental/small/constants.rs').read_text()
    central = (root/'src/experimental/small/central.rs').read_text()
    reference = (reference_root/'src/small/constants.rs').read_text()
    match = next(m for m in seeds.DECLARATION.finditer(constants) if m.group(1)=='lr_small_2')
    value = re.search(r'f64::from_bits\((0x[0-9a-f]+)\)',match.group(3)).group(1)
    bit = f'0x{int(value,16)^1:016x}'
    coefficient = constants[:match.start(3)] + change(constants[match.start(3):],value,bit)
    cases = [
        ('rank metadata',small,change(constants,'(&lr_small_2, 2)','(&lr_small_2, 3)')),
        ('compressed coefficient bit',small,coefficient),
        ('compressed coefficient extent',small,change(constants,'static lr_small_2: [f64; 28]','static lr_small_2: [f64; 27]')),
        ('degree metadata',small,change(constants,'[[6, 5], [7, 5]','[[6, 6], [7, 5]')),
        ('rectangle selector',change(small,'SMALL_RANK_DEGREES[k]','SMALL_RANK_DEGREES[0]'),constants),
        ('compressed shape',change(small,'rank_seed_shape::<6, 5, 4, 28>','rank_seed_shape::<6, 5, 8, 56>'),constants),
        ('Horner axis degree',change(small,'if k < half { X_DEGREE } else { Y_DEGREE }','if k < half { Y_DEGREE } else { X_DEGREE }'),constants),
        ('FMA operation',change(small,'a.mul_add(b, c)','a * b + c'),constants),
        ('unqualified central guard',change(small,'if a <= 0.45 && b < 0.26 {','if a <= 0.5 && b < 0.26 {'),constants),
        ('unrelated constant',small,change(constants,'0x40076fcca4f7f770','0x40076fcca4f7f771')),
    ]
    records = []
    seeds.restore_and_validate(small,constants,reference)
    as1.verify_source(root/'src/experimental/small/central.rs',reference_root)
    with tempfile.TemporaryDirectory(prefix='rank-seed-negative-') as directory:
        source = Path(directory)/'experimental'
        local = source/'small'
        local.mkdir(parents=True)
        (local/'central.rs').write_text(central)
        for label, mutated_small, mutated_constants in cases:
            (source/'small.rs').write_text(mutated_small)
            (local/'constants.rs').write_text(mutated_constants)
            for consumer in ('direct-D shared bridge','AS1 source checker'):
                try:
                    if consumer=='direct-D shared bridge':
                        seeds.restore_and_validate(mutated_small,mutated_constants,reference)
                    else:
                        as1.verify_source(local/'central.rs',reference_root)
                except AssertionError:
                    records.append({'defect':label,'consumer':consumer,'rejected':True})
                else:
                    raise RuntimeError(f'Accepted {label}: {consumer}')
    report = {'controls':records,'scope':f'{len(cases)} deliberate defects, each rejected by the shared bridge and AS1 consumer',
              'source_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in [root/'src/experimental/small.rs',root/'src/experimental/small/constants.rs',
                                            root/'scripts/seed_reduction.py',root/'scripts/verify_seed_bindings.py']}}
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
