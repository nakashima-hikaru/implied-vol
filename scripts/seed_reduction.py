"""Exact adaptive Chebyshev reduction and source binding for rank seeds.

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
# Fixed per-rectangle degrees qualified by the continuous inverse-finish ledger.
DEGREES = {
    'small': ((6, 5), (7, 5), (6, 5), (6, 5), (6, 5), (6, 5)),
    'core': ((6, 6), (6, 6), (6, 6), (7, 6), (6, 5),
             (6, 6), (6, 6), (6, 6), (6, 6), (6, 5),
             (6, 6), (6, 6), (6, 6), (6, 6), (6, 5),
             (6, 6), (6, 5), (6, 5), (6, 5), (6, 5)),
}
ORIGINAL_RUNTIME = """fn rank_seed<const DEGREE: usize>(t: (&[f64], usize), x: f64, y: f64) -> f64 {
    // The tables have four fixed shapes. Check the selected slice once so
    // Horner stages can use constant offsets without per-coefficient checks.
    const { assert!(DEGREE == 8 || DEGREE == 10) };
    let (c, rank) = t;
    if rank <= 2 {
        if DEGREE == 10 {
            rank_seed_shape::<10, 4, 44>(
                c.first_chunk::<44>().expect("rank coefficients"),
                rank,
                x,
                y,
            )
        } else {
            rank_seed_shape::<8, 4, 36>(
                c.first_chunk::<36>().expect("rank coefficients"),
                rank,
                x,
                y,
            )
        }
    } else if DEGREE == 10 {
        rank_seed_shape::<10, 8, 88>(
            c.first_chunk::<88>().expect("rank coefficients"),
            rank,
            x,
            y,
        )
    } else {
        rank_seed_shape::<8, 8, 72>(
            c.first_chunk::<72>().expect("rank coefficients"),
            rank,
            x,
            y,
        )
    }
}
#[inline]
fn rank_seed_shape<const DEGREE: usize, const LANES: usize, const LENGTH: usize>(
    c: &[f64; LENGTH],
    rank: usize,
    x: f64,
    y: f64,
) -> f64 {
    let half = LANES / 2;
    assert!(rank <= half, "rank factors");
    let mut values = [0.0; 8];
    for k in (0..LANES).step_by(2) {
        let z = F64x2::splat(if k < half { x } else { y });
        // Preserve the original degree-eight y factors in wider degree-ten seeds.
        let n = if DEGREE == 10 && LANES == 8 && k >= half {
            8
        } else {
            DEGREE
        };
        let mut v = F64x2::new(c[n * LANES + k], c[n * LANES + k + 1]);
        for j in (0..n).rev() {
            v = v.mul_add(z, F64x2::new(c[j * LANES + k], c[j * LANES + k + 1]));
        }
        values[k..k + 2].copy_from_slice(&v.to_array());
    }
    let mut v = values[0] * values[half];
    for k in 1..rank {
        v = fma(values[k], values[half + k], v);
    }
    v
}"""
ADAPTIVE_RUNTIME = """fn rank_seed<const DEGREE: usize>(t: (&[f64], usize), x: f64, y: f64) -> f64 {
    // The upper-gap tables retain their original degree-eight factors.
    const { assert!(DEGREE == 8) };
    let (c, rank) = t;
    if rank <= 2 {
        rank_seed_shape::<8, 8, 4, 36>(
            c.first_chunk::<36>().expect("rank coefficients"),
            rank,
            x,
            y,
        )
    } else {
        rank_seed_shape::<8, 8, 8, 72>(
            c.first_chunk::<72>().expect("rank coefficients"),
            rank,
            x,
            y,
        )
    }
}
#[inline]
fn adaptive_rank_seed(t: (&[f64], usize), degrees: [u8; 2], x: f64, y: f64) -> f64 {
    // Generated per-table degrees select fixed Horner chains and array sizes.
    let (c, rank) = t;
    match (rank <= 2, degrees) {
        (true, [6, 5]) => rank_seed_shape::<6, 5, 4, 28>(
            c.first_chunk::<28>().expect("rank coefficients"),
            rank,
            x,
            y,
        ),
        (false, [6, 5]) => rank_seed_shape::<6, 5, 8, 56>(
            c.first_chunk::<56>().expect("rank coefficients"),
            rank,
            x,
            y,
        ),
        (false, [6, 6]) => rank_seed_shape::<6, 6, 8, 56>(
            c.first_chunk::<56>().expect("rank coefficients"),
            rank,
            x,
            y,
        ),
        (false, [7, 5]) => rank_seed_shape::<7, 5, 8, 64>(
            c.first_chunk::<64>().expect("rank coefficients"),
            rank,
            x,
            y,
        ),
        (false, [7, 6]) => rank_seed_shape::<7, 6, 8, 64>(
            c.first_chunk::<64>().expect("rank coefficients"),
            rank,
            x,
            y,
        ),
        _ => unreachable!("rank coefficient degrees"),
    }
}
#[inline]
fn rank_seed_shape<
    const X_DEGREE: usize,
    const Y_DEGREE: usize,
    const LANES: usize,
    const LENGTH: usize,
>(
    c: &[f64; LENGTH],
    rank: usize,
    x: f64,
    y: f64,
) -> f64 {
    const {
        assert!(X_DEGREE >= Y_DEGREE);
        assert!(LENGTH == (X_DEGREE + 1) * LANES);
    }
    let half = LANES / 2;
    assert!(rank <= half, "rank factors");
    let mut values = [0.0; 8];
    for k in (0..LANES).step_by(2) {
        let z = F64x2::splat(if k < half { x } else { y });
        let n = if k < half { X_DEGREE } else { Y_DEGREE };
        let mut v = F64x2::new(c[n * LANES + k], c[n * LANES + k + 1]);
        for j in (0..n).rev() {
            v = v.mul_add(z, F64x2::new(c[j * LANES + k], c[j * LANES + k + 1]));
        }
        values[k..k + 2].copy_from_slice(&v.to_array());
    }
    let mut v = values[0] * values[half];
    for k in 1..rank {
        v = fma(values[k], values[half + k], v);
    }
    v
}"""
TEST_EDITS = [
    ('fn conservative_preflight_seams_preserve_original_results()', 'fn conservative_preflight_seams_retain_certified_outputs()'),
    ('// Frozen original experimental outputs at the new reject thresholds,\n        // including adjacent prices and the dyadic jump just above a=1.', '// Frozen outputs at reject thresholds, adjacent prices and the a=1 jump.\n        // Changed outputs were checked against independent Black roots at\n        // 160/240 decimal digits; their rho values are all below 0.363.'),
    ('(0x3fb999999999999a, 0x3fc3333333333333, 0x3fdfbf7426673d99)', '(0x3fb999999999999a, 0x3fc3333333333333, 0x3fdfbf7426673d98)'),
    ('(0x3fb999999999999a, 0x3fc3333333333334, 0x3fdfbf7426673d9a)', '(0x3fb999999999999a, 0x3fc3333333333334, 0x3fdfbf7426673d99)'),
    ('(0x3ff0000000000000, 0x3fd0000000000001, 0x3ffbc1183e525ee2)', '(0x3ff0000000000000, 0x3fd0000000000001, 0x3ffbc1183e525ee3)'),
    ('(0x3ff0000000000001, 0x3fbfffffffffffff, 0x3ff363b69fb0f1c0)', '(0x3ff0000000000001, 0x3fbfffffffffffff, 0x3ff363b69fb0f1bf)'),
    ('(0x4014000000000000, 0x3f90000000000000, 0x400539daf0082b99)', '(0x4014000000000000, 0x3f90000000000000, 0x400539daf0082b98)'),
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
        route, index = name.removeprefix('lr_').rsplit('_', 1)
        x_degree, y_degree = DEGREES[route][int(index)]
        fields = []
        for lane in range(lanes):
            degree = x_degree if lane < lanes//2 else y_degree
            new, removed, rounding = reduce_factor(values[lane::lanes], degree)
            fields.append({'lane': lane, 'removed': str(removed), 'rounding': str(rounding)})
            for i, value in enumerate(new): values[i*lanes+lane] = value
        values = values[:(x_degree+1)*lanes]
        declaration = f'pub(super) static {name}: [f64; {len(values)}] = [\n'
        declaration += ''.join('    f64::from_bits(' + bits(value) + '),\n' for value in values) + '];'
        source = source[:match.start()] + declaration + source[match.end():]
        metadata[name] = {'name': name, 'lanes': lanes, 'fields': fields,
                          'old_bytes': int(length)*8, 'new_bytes': len(values)*8,
                          'x_degree': x_degree, 'y_degree': y_degree}
    assert len(metadata) == 26
    return source + degree_declarations(), metadata


def degree_declarations():
    result = '\n'
    for route, degrees in DEGREES.items():
        prefix = f'pub(super) static {route.upper()}_RANK_DEGREES: [[u8; 2]; {len(degrees)}] ='
        if route == 'small':
            result += prefix + '\n    [' + ', '.join(f'[{x}, {y}]' for x, y in degrees) + '];\n'
        else:
            result += prefix + ' [\n'
            result += ''.join(f'    [{x}, {y}],\n' for x, y in degrees) + '];\n'
    return result


CALL_EDITS = [
    ('q * rank_seed::<10>(small_lr[k], (w - r[0]) * r[1], fma(q, q, -r[2]) * r[3]),',
     'q * adaptive_rank_seed(\n            small_lr[k],\n            SMALL_RANK_DEGREES[k],\n            (w - r[0]) * r[1],\n            fma(q, q, -r[2]) * r[3],\n        ),'),
    ('rank_seed::<10>(core_lr[k], (a - r[0]) * r[1], (z - r[2]) * r[3]),',
     'adaptive_rank_seed(\n            core_lr[k],\n            CORE_RANK_DEGREES[k],\n            (a - r[0]) * r[1],\n            (z - r[2]) * r[3],\n        ),'),
]


def reduce_small(source):
    assert sha256(source) == BASELINE_SMALL_SHA256
    source = replace_once(source, ORIGINAL_RUNTIME, ADAPTIVE_RUNTIME)
    for old, new in CALL_EDITS + TEST_EDITS:
        source = replace_once(source, old, new)
    return source


def restore_and_validate(small, constants, reference):
    """Regenerate every degree/coefficient and restore reviewed baseline hashes.

    Runtime dispatch, Horner degrees, array extents and table selection are
    bound by the exact inverse below. No newly recorded whole-file hash can
    admit a changed selector, coefficient or unrelated source operation.
    """
    assert sha256(reference) == REFERENCE_CONSTANTS_SHA256
    references = {m.group(1): m.group(0) for m in DECLARATION.finditer(reference)}
    degree_source = degree_declarations()
    assert constants.endswith(degree_source), 'adaptive degree metadata mismatch'
    restored_constants = constants.removesuffix(degree_source)
    for match in list(DECLARATION.finditer(restored_constants))[::-1]:
        restored_constants = restored_constants[:match.start()] + references[match.group(1)] + restored_constants[match.end():]
    assert sha256(restored_constants) == BASELINE_CONSTANTS_SHA256, 'non-seed constant change'
    expected, metadata = reduced_constants(restored_constants)
    assert constants == expected, 'candidate differs from exact adaptive Chebyshev reduction'
    restored_small = replace_once(small, ADAPTIVE_RUNTIME, ORIGINAL_RUNTIME)
    for old, new in CALL_EDITS + TEST_EDITS:
        restored_small = replace_once(restored_small, new, old)
    assert sha256(restored_small) == BASELINE_SMALL_SHA256, 'non-adaptive small-kernel change'
    return restored_small, restored_constants, metadata
