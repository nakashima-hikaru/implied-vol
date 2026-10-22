"""Exact Chebyshev reduction and a narrow source bridge for rank seeds.

This bridge changes an approximation. It does not claim bit equivalence.
Only the unchanged kernel source and archived seed cover are transported;
verify_seed_reduction separately bounds the changed approximation.
"""
from fractions import Fraction as F
import hashlib
import re
import struct

BASELINE_SMALL_SHA256 = '711bde18c91e67c18645581e61ce2b3a3fce95b79794281b783537a4aa19b7af'
BASELINE_CONSTANTS_SHA256 = '498b84a03308000fd71948ac554b4ae59cd914b9ff03e811aad735034c03e509'
REFERENCE_CONSTANTS_SHA256 = 'f9efd22b8d0c88f080855f61ab5778dedbf077c9788348f5ab4f199d283b0a8d'
DECLARATION = re.compile(r'pub\(super\) static (lr_(?:small|core)_\d+): \[f64; (\d+)\] = \[(.*?)\];', re.S)
ORIGINAL_DEGREES = '''let n = if DEGREE == 10 && LANES == 8 && k >= half {
            8
        } else {
            DEGREE
        };'''
REDUCED_DEGREES = '''let n = if DEGREE == 10 {
            if k < half { 8 } else { 6 }
        } else {
            DEGREE
        };'''
ORIGINAL_COMMENT = '// Preserve the original degree-eight y factors in wider degree-ten seeds.'
REDUCED_COMMENT = '// Small/core factors retain Chebyshev degrees eight in x and six in y.'
TEST_EDITS = [
    ('fn conservative_preflight_seams_preserve_original_results()', 'fn conservative_preflight_seams_retain_certified_outputs()'),
    ('// Frozen original experimental outputs at the new reject thresholds,\n        // including adjacent prices and the dyadic jump just above a=1.',
     '// Frozen outputs at reject thresholds, adjacent prices and the a=1 jump.\n        // The two reduced-seed changes were checked by independent Black roots\n        // at 100/160 decimal digits: rho < 0.053 and rho < 0.348 respectively.'),
    ('(0x3fb999999999999a, 0x3fc3333333333333, 0x3fdfbf7426673d99)', '(0x3fb999999999999a, 0x3fc3333333333333, 0x3fdfbf7426673d98)'),
    ('(0x4014000000000000, 0x3f8fffffffffffff, 0x400539daf0082b98)', '(0x4014000000000000, 0x3f8fffffffffffff, 0x400539daf0082b99)'),
]


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


def replace_once(source, old, new):
    assert source.count(old) == 1, ('seed source bridge pattern', old)
    return source.replace(old, new, 1)


def decode(value):
    return F(struct.unpack('<d', struct.pack('<Q', int(value, 16)))[0])


def bits(value):
    return f'0x{struct.unpack("<Q", struct.pack("<d", float(value)))[0]:016x}'


def represented_table(source, name):
    body = source.split('static ' + name + ':', 1)[1].split('=', 1)[1].split('];', 1)[0]
    return [decode(v) for v in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)', body)]


def polynomial_add(left, right):
    result = [F(0)] * max(len(left), len(right))
    for i, value in enumerate(left): result[i] += value
    for i, value in enumerate(right): result[i] += value
    return result


CHEBYSHEV = [[F(1)], [F(0), F(1)]]
for n in range(2, 11):
    CHEBYSHEV.append(polynomial_add([F(0)] + [2*v for v in CHEBYSHEV[-1]], [-v for v in CHEBYSHEV[-2]]))


def reduce_factor(values, degree):
    power = values[:]
    chebyshev = [F(0)] * len(power)
    for n in range(len(power)-1, -1, -1):
        chebyshev[n] = power[n] / CHEBYSHEV[n][n]
        for j, value in enumerate(CHEBYSHEV[n]): power[j] -= chebyshev[n] * value
    assert all(value == 0 for value in power)
    truncated = [F(0)] * len(values)
    for n in range(min(degree+1, len(values))):
        for j, value in enumerate(CHEBYSHEV[n]): truncated[j] += chebyshev[n] * value
    actual = [F(float(value)) for value in truncated]
    removed = sum(abs(v) for v in chebyshev[degree+1:])
    rounding = sum(abs(a-b) for a, b in zip(actual, truncated))
    return actual, removed, rounding


def reduced_constants(original):
    assert sha256(original) == BASELINE_CONSTANTS_SHA256
    source = original
    metadata = {}
    for match in list(DECLARATION.finditer(original))[::-1]:
        name, length, body = match.groups()
        values = [decode(v) for v in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)', body)]
        assert int(length) == len(values) and len(values) in (44, 88)
        lanes = 4 if len(values) == 44 else 8
        fields = []
        for lane in range(lanes):
            degree = 8 if lane < lanes//2 else 6
            new, removed, rounding = reduce_factor(values[lane::lanes], degree)
            fields.append({'lane': lane, 'removed': str(removed), 'rounding': str(rounding)})
            for i, value in enumerate(new): values[i*lanes+lane] = value
        values = values[:9*lanes]
        declaration = f'pub(super) static {name}: [f64; {len(values)}] = [\n'
        declaration += ''.join('    f64::from_bits(' + bits(value) + '),\n' for value in values) + '];'
        source = source[:match.start()] + declaration + source[match.end():]
        metadata[name] = {'name': name, 'lanes': lanes, 'fields': fields,
                          'old_bytes': int(length)*8, 'new_bytes': len(values)*8}
    assert len(metadata) == 26
    return source, metadata


def reduce_small(source):
    assert sha256(source) == BASELINE_SMALL_SHA256
    source = replace_once(source, ORIGINAL_DEGREES, REDUCED_DEGREES)
    source = replace_once(source, ORIGINAL_COMMENT, REDUCED_COMMENT)
    for lanes, old, new in [(4, 44, 36), (8, 88, 72)]:
        source = replace_once(source, f'rank_seed_shape::<10, {lanes}, {old}>', f'rank_seed_shape::<10, {lanes}, {new}>')
        source = replace_once(source, f'first_chunk::<{old}>()', f'first_chunk::<{new}>()')
    for old, new in TEST_EDITS:
        source = replace_once(source, old, new)
    return source


def restore_and_validate(small, constants, reference):
    """Return baseline source plus regenerated approximation metadata.

    No source is trusted merely because a newly recorded hash matches it.
    The narrow inverse must reproduce both reviewed baseline hashes, and
    exact coefficient regeneration must reproduce the actual candidate.
    """
    assert sha256(reference) == REFERENCE_CONSTANTS_SHA256
    references = {m.group(1): m.group(0) for m in DECLARATION.finditer(reference)}
    restored_constants = constants
    for match in list(DECLARATION.finditer(constants))[::-1]:
        restored_constants = restored_constants[:match.start()] + references[match.group(1)] + restored_constants[match.end():]
    assert sha256(restored_constants) == BASELINE_CONSTANTS_SHA256, 'non-seed constant change'
    expected, metadata = reduced_constants(restored_constants)
    assert constants == expected, 'candidate differs from exact Chebyshev reduction'
    restored_small = replace_once(small, REDUCED_DEGREES, ORIGINAL_DEGREES)
    restored_small = replace_once(restored_small, REDUCED_COMMENT, ORIGINAL_COMMENT)
    for lanes, old, new in [(4, 44, 36), (8, 88, 72)]:
        restored_small = replace_once(restored_small, f'rank_seed_shape::<10, {lanes}, {new}>(\n                c.first_chunk::<{new}>()', f'rank_seed_shape::<10, {lanes}, {old}>(\n                c.first_chunk::<{old}>()') if lanes == 4 else replace_once(restored_small, f'rank_seed_shape::<10, {lanes}, {new}>(\n            c.first_chunk::<{new}>()', f'rank_seed_shape::<10, {lanes}, {old}>(\n            c.first_chunk::<{old}>()')
    for old, new in TEST_EDITS:
        restored_small = replace_once(restored_small, new, old)
    assert sha256(restored_small) == BASELINE_SMALL_SHA256, 'non-degree small-kernel change'
    return restored_small, restored_constants, metadata
