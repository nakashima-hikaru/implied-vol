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
# The two retained SVD factors of the degree-6/5 Chebyshev coefficient matrix
# for small rectangle 2, frozen once as binary64 Chebyshev coefficients. SVD
# is a search tool, not a proof premise: the checker bounds the exact difference
# between the represented original and generated rank-two polynomials.
# Each pair is (x factor, y factor); the exact power conversion below followed
# by one binary64 rounding defines the runtime coefficients reproducibly.
SMALL_2_RANK_TWO_CHEBYSHEV = (
    (
        ('0xbfe66f777266e96f', '0x3fc9a792d85df2f9', '0xbf9bb84e18878a6b',
         '0x3f65d704111dc031', '0xbf1b206d207761dc', '0xbefde5efa5b8085c',
         '0x3ee39cea670a5f9a'),
        ('0xbfe759a6039bfde8', '0xbf68697166c21bb9', '0xbef2f6264e59fef8',
         '0xbe83c457ab37298e', '0xbe17c8c76173c0eb', '0xbdaf2369f98e4bb5'),
    ),
    (
        ('0x3f83ff32948ee10a', '0x3fa080023860c07d', '0xbf8c76bff7f3c958',
         '0x3f6f9edf4ae78f8e', '0xbf4afe0cef5ae050', '0x3f224c9f4548f383',
         '0xbef2adc26bcb5f79'),
        ('0x3f239d2d4d75b0b8', '0xbfa2c289394ce8c1', '0xbf35d47b4ee3989b',
         '0xbecadaed345078e4', '0xbe619f597179b99e', '0xbdf85418787fd5ac'),
    ),
)
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


# Test-only accuracy cutover: retain exact inverse source binding while allowing
# platform-dependent output bits within the independent Black error budget.
PREFLIGHT_TEST_EDITS = [
    (
        'fn conservative_preflight_seams_retain_certified_outputs()',
        'fn conservative_preflight_seams_retain_attainable_precision()',
    ),
    (
        """        // Frozen outputs at reject thresholds, adjacent prices and the a=1 jump.
        // Changed outputs were checked against independent Black roots at
        // 160/240 decimal digits; their rho values are all below 0.363.""",
        """        // Reject thresholds, adjacent prices and the a=1 jump. The retained
        // historical outputs are reference anchors, not required output bits.
        // Independent Black roots and eta_J agree at 160/240 decimal digits;
        // see scripts/generate_adaptive_seed_seams.py::invert. Each correction
        // stores the high-precision root minus its historical anchor.""",
    ),
    (
        """            (0x4014000000000000, 0x3f90000000000001, 0x400539daf0082b99),
        ];
        for &(a, b, want) in cases {
            assert_eq!(
                solve(f64::from_bits(a), f64::from_bits(b)).to_bits(),
                want,
                "a={a:x}, b={b:x}"
            );
        }""",
        """            (0x4014000000000000, 0x3f90000000000001, 0x400539daf0082b99),
        ];
        // Root correction, eta_J (both encoded as binary64).
        let corrections: &[(u64, u64)] = &[
            (0xbc61b36975b0dc64, 0x3cb50ba378e54a32),
            (0xbc510a23f1b3675b, 0x3cb50ba378e54a32),
            (0x3c1528b07faea125, 0x3cb50ba378e54a32),
            (0x3c6818b23594cb0a, 0x3cbcc34e0871ec3c),
            (0x3c803b60b42fa9cc, 0x3cbcc34e0871ec3c),
            (0xbc663dac94177cac, 0x3cbcc34e0871ec3c),
            (0xbc6e44f7b265a44a, 0x3cb9f109ed6ec32f),
            (0x3c9eb4d7bb196026, 0x3cb9f109ed6ec32f),
            (0xbc9c503ae21a767b, 0x3cb9f109ed6ec330),
            (0x3ca27b92898d89c2, 0x3cb6fc8c09423134),
            (0xbc8433f32a1ad6fa, 0x3cb6fc8c09423134),
            (0x3c97c3c91aa19681, 0x3cb6fc8c09423134),
            (0x3cac4f81abe980a5, 0x3cb35c9b851805e8),
            (0x3cb29d74dee43753, 0x3cb35c9b851805e8),
            (0xbc91dc8c3cf36ab4, 0x3cb35c9b851805e8),
        ];
        assert_eq!(cases.len(), corrections.len());
        for (&(a, b, anchor), &(correction, eta)) in cases.iter().zip(corrections) {
            let s = solve(f64::from_bits(a), f64::from_bits(b));
            let hi = f64::from_bits(anchor);
            let lo = f64::from_bits(correction);
            let rho = (((s - hi) - lo).abs() / hi) / f64::from_bits(eta);
            assert!(
                s.is_finite() && s > 0.0 && rho < 1.0,
                "a={a:x}, b={b:x}, s={s:e}, rho_J={rho}"
            );
        }""",
    ),
]
TEST_EDITS += PREFLIGHT_TEST_EDITS

# Test-only migration for the central normal-degree reduction. Restore these
# exact edits rather than ignoring the test module in the source certificate.
ARCHIVED_CENTRAL_TEST = """    #[test]
    fn central_archived_input_keeps_attainable_precision() {
        // Exact binary64 inputs, with an independent Black root and eta_J
        // agreeing at 160/240 digits via generate_adaptive_seed_seams.py::invert.
        let a = f64::from_bits(0x3fb999999999999a);
        let b = f64::from_bits(0x3fa999999999999a);
        let hi = f64::from_bits(0x3fcd6496956a4ac3);
        let lo = f64::from_bits(0x3c6eb51a9e3bf1b6);
        let eta = f64::from_bits(0x3cb9aa2c0dd715c2);
        let s = solve(a, b);
        let rho = (((s - hi) - lo).abs() / hi) / eta;
        assert!(s.is_finite() && s > 0.0 && rho < 1.0, "rho_J={rho}");
    }
"""
TEST_EDITS += [
    (
        '        // Expected bit patterns come from the unchanged C++ experimental with its AVX2 lane order.',
        '        // Unchanged branches retain the archived C++/AVX2 regression anchors.\n'
        '        // The optimized central input is checked against an independent root below.',
    ),
    (
        '            (0x16687e92154ef7ac, 0x124212dd4de70913, 0x16364f16a8097a31),\n'
        '            (0x3fb999999999999a, 0x3fa999999999999a, 0x3fcd6496956a4ac3),\n'
        '            (0x3fb999999999999a, 0x39b4484bfeebc2a0, 0x3f82e778a810e85c),',
        '            (0x16687e92154ef7ac, 0x124212dd4de70913, 0x16364f16a8097a31),\n'
        '            (0x3fb999999999999a, 0x39b4484bfeebc2a0, 0x3f82e778a810e85c),',
    ),
    ('    #[test]\n    fn invalid_domain()',
     ARCHIVED_CENTRAL_TEST + '    #[test]\n    fn invalid_domain()'),
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


def power_to_chebyshev(values):
    power = values[:]
    chebyshev = [F(0)] * len(power)
    for n in range(len(power)-1, -1, -1):
        chebyshev[n] = power[n] / CHEBYSHEV[n][n]
        for j, value in enumerate(CHEBYSHEV[n]): power[j] -= chebyshev[n] * value
    assert all(value == 0 for value in power)
    return chebyshev


def chebyshev_to_power(values, length):
    result = [F(0)] * length
    for n, coefficient in enumerate(values):
        for j, value in enumerate(CHEBYSHEV[n]): result[j] += coefficient * value
    return result


def reduce_factor(values, degree):
    chebyshev = power_to_chebyshev(values)
    truncated = [F(0)] * len(values)
    for n in range(min(degree+1, len(values))):
        for j, value in enumerate(CHEBYSHEV[n]): truncated[j] += chebyshev[n] * value
    actual = [F(float(value)) for value in truncated]
    removed = sum(abs(v) for v in chebyshev[degree+1:])
    rounding = sum(abs(a-b) for a, b in zip(actual, truncated))
    return actual, removed, rounding


def rank_polynomial(values, lanes, rank):
    """Exact bivariate Chebyshev coefficients of represented power factors."""
    factors = [power_to_chebyshev(values[lane::lanes]) for lane in range(lanes)]
    n = len(factors[0])
    return [[sum(factors[k][i]*factors[lanes//2+k][j] for k in range(rank))
             for j in range(n)] for i in range(n)]


def rank_two_small_2():
    factors = []
    for side in range(2):
        for pair in SMALL_2_RANK_TWO_CHEBYSHEV:
            exact = chebyshev_to_power([decode(v) for v in pair[side]], 7)
            factors.append([F(float(v)) for v in exact])
    return [factors[lane][i] for i in range(7) for lane in range(4)]


def check_rank_two_search(original):
    """Reproduce the exploratory SVD at two precisions, outside proof trust."""
    import mpmath as mp
    values = represented_table(original, 'lr_small_2')
    for lane in range(8):
        factor, _, _ = reduce_factor(values[lane::8], 6 if lane<4 else 5)
        for i,value in enumerate(factor): values[8*i+lane] = value
    matrix = rank_polynomial(values, 8, 3)
    for precision in (90, 130):
        with mp.workdps(precision):
            exact = lambda value: mp.mpf(value.numerator)/value.denominator
            left,singular,right = mp.svd(mp.matrix([[exact(matrix[i][j]) for j in range(6)] for i in range(7)]))
            for k,pair in enumerate(SMALL_2_RANK_TWO_CHEBYSHEV):
                for side,factor in enumerate(pair):
                    generated = tuple(bits(mp.sqrt(singular[k])*(left[i,k] if side==0 else right[k,i]))
                                      for i in range(7 if side==0 else 6))
                    assert generated == factor, ('SVD search provenance',precision,k,side)
    return {'search_precision_digits':[90,130], 'frozen_chebyshev_factor_bits_reproduced':True,
            'scope':'Exploratory SVD provenance only; the proof uses exact represented coefficient differences'}


def rank_difference(original, old_lanes, old_rank, candidate, lanes, rank):
    before = rank_polynomial(original, old_lanes, old_rank)
    after = rank_polynomial(candidate, lanes, rank)
    return sum(abs(value-(after[i][j] if i<len(after) and j<len(after) else 0))
               for i,row in enumerate(before) for j,value in enumerate(row))


def reduced_constants(original):
    assert sha256(original) == BASELINE_CONSTANTS_SHA256
    source = original
    metadata = {}
    for match in list(DECLARATION.finditer(original))[::-1]:
        name, length, body = match.groups()
        values = [decode(v) for v in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)', body)]
        original_values = values[:]
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
        compression = {}
        if name == 'lr_small_2':
            values = rank_two_small_2()
            compression = {'old_lanes': lanes, 'old_rank': 3, 'rank': 2,
                           'fields': [],
                           'chebyshev_difference': str(rank_difference(original_values, lanes, 3, values, 4, 2))}
            lanes = 4
        declaration = f'pub(super) static {name}: [f64; {len(values)}] = [\n'
        declaration += ''.join('    f64::from_bits(' + bits(value) + '),\n' for value in values) + '];'
        source = source[:match.start()] + declaration + source[match.end():]
        metadata[name] = {'name': name, 'lanes': lanes, 'fields': fields,
                          'old_bytes': int(length)*8, 'new_bytes': len(values)*8,
                          'x_degree': x_degree, 'y_degree': y_degree, **compression}
    assert len(metadata) == 26
    source = replace_once(source, '(&lr_small_2, 3)', '(&lr_small_2, 2)')
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

# Independently qualified by scripts/verify_central_expansion.py. Restore this
# guard exactly before importing any unchanged small-kernel/AS1 certificate.
CENTRAL_GUARD_EDITS = [
    ('if a <= 0.36 && b < 0.21 {', 'if a <= 0.45 && b < 0.26 {'),
]


def reduce_small(source):
    assert sha256(source) == BASELINE_SMALL_SHA256
    source = replace_once(source, ORIGINAL_RUNTIME, ADAPTIVE_RUNTIME)
    for old, new in CALL_EDITS + CENTRAL_GUARD_EDITS + TEST_EDITS:
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
    restored_constants = replace_once(restored_constants, '(&lr_small_2, 2)', '(&lr_small_2, 3)')
    assert sha256(restored_constants) == BASELINE_CONSTANTS_SHA256, 'non-seed constant change'
    expected, metadata = reduced_constants(restored_constants)
    assert constants == expected, 'candidate differs from exact adaptive Chebyshev reduction'
    restored_small = replace_once(small, ADAPTIVE_RUNTIME, ORIGINAL_RUNTIME)
    for old, new in reversed(CALL_EDITS + CENTRAL_GUARD_EDITS + TEST_EDITS):
        restored_small = replace_once(restored_small, new, old)
    assert sha256(restored_small) == BASELINE_SMALL_SHA256, 'non-adaptive small-kernel change'
    return restored_small, restored_constants, metadata
