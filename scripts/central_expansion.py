"""Exact source bridge for the wider, lower-degree central expansion.

All coefficients and unchanged arithmetic remain source-bound. This bridge
restores the reviewed guards, normal-polynomial degree specialization and
conversion-polynomial truncation;
verify_central_expansion.py proves those changes separately, not as bit
equivalence to the original kernel.
"""
import hashlib

BASELINE_SHA256 = '61e04489aa1849afdbabe74273685f0b2062c436035c77c0d4bcae648f5a642f'
NORMAL_DEGREES = (16, 18, 18, 24)
CONVERSION_LENGTHS = (12, 11, 9, 7, 5, 3, 1, 0, 0, 0)
EDITS = (
    (
        '    // b<a/16 gives z>0.503; b+a/2>=0.2 gives q>0.5013.\n'
        '    if b < a / 16.0 || b + 0.5 * a >= 0.2 {',
        '    // b<a/16 gives z>0.503; b+a/2>=0.25 gives q>0.6266.\n'
        '    // The wider q cap is covered by verify_central_expansion.py.\n'
        '    if b < a / 16.0 || b + 0.5 * a >= 0.25 {',
    ),
    ('    if !(q.hi >= 1e-12 && q.hi <= 0.5) {',
     '    if !(q.hi >= 1e-12 && q.hi <= 0.625) {'),
    (
        '    // Independent Horner chains share the argument; keep each chain\'s degree\n'
        '    // and every original rounding point while evaluating in two NEON lanes.\n'
        '    let mut value_derivative = F64x2::new(fma(0.0, t.hi, NC[k][25]), 25.0 * NC[k][25]);\n'
        '    let arg = F64x2::splat(t.hi);\n'
        '    for i in (2..25).rev() {',
        '    // The represented coefficient tails are bounded on the padded cells in\n'
        '    // verify_central_expansion.py, including the final solver error budget.\n'
        '    match k {\n'
        '        0 => normal_polynomial::<16>(k, t),\n'
        '        1 | 2 => normal_polynomial::<18>(k, t),\n'
        '        _ => normal_polynomial::<24>(k, t),\n'
        '    }\n'
        '}\n'
        '#[inline]\n'
        'fn normal_polynomial<const DEGREE: usize>(k: usize, t: Pair) -> Pair {\n'
        '    const { assert!(DEGREE >= 2 && DEGREE <= 25) };\n'
        '    // Retain both compensated final stages and the low-argument correction.\n'
        '    let mut value_derivative =\n'
        '        F64x2::new(fma(0.0, t.hi, NC[k][DEGREE]), DEGREE as f64 * NC[k][DEGREE]);\n'
        '    let arg = F64x2::splat(t.hi);\n'
        '    for i in (2..DEGREE).rev() {',
    ),
    (
        '''fn conversion_pair<const LEFT: usize, const RIGHT: usize>(
    a: &[f64; LEFT],
    b: &[f64; RIGHT],
    x: f64,
) -> [f64; 2] {
    // The rows have different degrees. Evaluate the longer row's leading
    // stages first, then pair the remaining stages without padding either row.
    const { assert!(LEFT >= RIGHT && RIGHT > 0) };
    let head = LEFT - RIGHT;
    let mut p = a[0];
    for c in &a[1..=head] {
        p = fma(p, x, *c);
    }
    let mut p = F64x2::new(p, b[0]);
    let arg = F64x2::splat(x);
    for i in 1..RIGHT {
        p = p.mul_add(arg, F64x2::new(a[head + i], b[i]));
    }
    p.to_array()
}''',
        '''fn conversion_pair<
    const LEFT: usize,
    const RIGHT: usize,
    const LEFT_KEEP: usize,
    const RIGHT_KEEP: usize,
>(
    a: &[f64; LEFT],
    b: &[f64; RIGHT],
    x: f64,
) -> [f64; 2] {
    // Drop only the source-certified highest powers of each polynomial.
    // Pair the retained Horner stages without padding either row.
    const {
        assert!(LEFT >= LEFT_KEEP && RIGHT >= RIGHT_KEEP);
        assert!(LEFT_KEEP >= RIGHT_KEEP && RIGHT_KEEP > 0);
    };
    let left_start = LEFT - LEFT_KEEP;
    let right_start = RIGHT - RIGHT_KEEP;
    let head = LEFT_KEEP - RIGHT_KEEP;
    let mut p = a[left_start];
    for c in &a[left_start + 1..=left_start + head] {
        p = fma(p, x, *c);
    }
    let mut p = F64x2::new(p, b[right_start]);
    let arg = F64x2::splat(x);
    for i in 1..RIGHT_KEEP {
        p = p.mul_add(
            arg,
            F64x2::new(a[left_start + head + i], b[right_start + i]),
        );
    }
    p.to_array()
}''',
    ),
    (
        '''    let [p0, p1] = conversion_pair(&conversion_0, &conversion_1, w.hi);
    let [p2, p3] = conversion_pair(&conversion_2, &conversion_3, w.hi);
    let [p4, p5] = conversion_pair(&conversion_4, &conversion_5, w.hi);
    let [p6, p7] = conversion_pair(&conversion_6, &conversion_7, w.hi);
    let [p8, p9] = conversion_pair(&conversion_8, &conversion_9, w.hi);
    let rows = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9];
    let mut y = rows[9];
    for i in (0..9).rev() {''',
        '''    let [p0, p1] = conversion_pair::<16, 15, 12, 11>(&conversion_0, &conversion_1, w.hi);
    let [p2, p3] = conversion_pair::<14, 13, 9, 7>(&conversion_2, &conversion_3, w.hi);
    let [p4, p5] = conversion_pair::<11, 10, 5, 3>(&conversion_4, &conversion_5, w.hi);
    let rows = [p0, p1, p2, p3, p4, p5];
    let mut y = conversion_6[7];
    for i in (0..6).rev() {''',
    ),
)


def restore_central_source(source):
    """Restore the exact reviewed edits, rejecting all other mutations."""
    if hashlib.sha256(source.encode()).hexdigest() == BASELINE_SHA256:
        return source
    for old, new in reversed(EDITS):
        assert source.count(new) == 1, ('central source bridge', new)
        source = source.replace(new, old, 1)
    assert hashlib.sha256(source.encode()).hexdigest() == BASELINE_SHA256, (
        'central arithmetic changed outside qualified edits')
    return source
