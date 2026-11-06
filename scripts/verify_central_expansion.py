#!/usr/bin/env python3
"""Continuous error budget for the wider, lower-degree central expansion.

The hash-bound archive supplies exact interval/polynomial arithmetic and the
original normal-shape certificate. Exact represented coefficient tails cover
the degree reductions, including padded cells and low-coordinate corrections.
The Black conversion model and native source-order rounding are recomputed.
No sampled error, fitted remainder, fallback, or output polish is used.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib
import json
import math
from pathlib import Path
import re
import sys

from central_expansion import (
    BASELINE_SHA256, CONVERSION_LENGTHS, NORMAL_DEGREES, restore_central_source,
)
import seed_reduction

sys.dont_write_bytecode = True
DEPENDENCIES = {
    'scripts/prove_nbt_conversion.py': '3d5473dd453a901393b8c6f778c72580a5f7b453da411a97a9c90c59e57872fa',
    'scripts/prove_nbt_native.py': '39ddc57a82629048510c1b72f2452b85ca547171dbf925b8aa242e1eff6d2153',
    'scripts/prove_nbt_model.py': 'ae94894216b8cd2acb8e4e98b82af6f0f14585af5eb6eab281e06cbc9bc4dc4a',
    'scripts/atm_interval_tools.py': 'dfa3a7533f7a96ac2679a2b8eca96d4dd2a2d5a59fce433781073e4a850eccb1',
    'scripts/prove_atm_small_rho.py': '3ab901a87c24b66f165777d441d126957333736fd4445eef3bf711543a4f1f22',
    'scripts/prove_atm_large_rho.py': '332cc94e65cbe33a94289a3863bcdb943f5af61c548d44265092e123fb2ebdca',
    'scripts/verify_global_routing_lemmas.py': '0d179a16fb3bd218f48773ca28053e03b2815e99c43cffcd19e9c8bbb9a012aa',
    'scripts/verify_lbr_domain_bounds.py': '37d6f0fb6e62c283e0c96179a475def4afb22c71281fe639e9dd2d171d8fcd9f',
    'src/small/constants.rs': 'f9efd22b8d0c88f080855f61ab5778dedbf077c9788348f5ab4f199d283b0a8d',
    'results/global_rho_proof/nbt_normal_model.json': '1fb719503ab81826dc8864beec4a90cd8df262b7807e96b9d0f4a6e3b6d677fa',
}
AMAX = F.from_float(.45)
AMAX2 = F(203, 1000)
QCAP = F(5, 8)
QMAX = F(625000000001, 10**12)
WMAX = F(391, 1000)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_table(source, name):
    block = source.split('static ' + name + ':', 1)[1].split('=', 1)[1].split('];', 1)[0]
    return re.findall(r'(-?)f64::from_bits\((0x[0-9a-f]+)\)', block)


def source_binding(source, reference):
    central_path = source / 'src/experimental/small/central.rs'
    central = central_path.read_text()
    assert digest(central_path) != BASELINE_SHA256, 'candidate guards not installed'
    restore_central_source(central)
    assert NORMAL_DEGREES[1] == NORMAL_DEGREES[2]
    for arm, degree in [('0', NORMAL_DEGREES[0]),
                        ('1 | 2', NORMAL_DEGREES[1]), ('_', NORMAL_DEGREES[3])]:
        assert central.count(f'{arm} => normal_polynomial::<{degree}>(k, t),') == 1, (
            'normal proof/source degree mismatch', arm, degree)
    for k, left, right in [(0, 16, 15), (2, 14, 13), (4, 11, 10)]:
        call = (f'conversion_pair::<{left}, {right}, {CONVERSION_LENGTHS[k]}, '
                f'{CONVERSION_LENGTHS[k+1]}>(&conversion_{k}, &conversion_{k+1}, w.hi)')
        assert central.count(call) == 1, ('conversion proof/source lengths mismatch', k)
    assert CONVERSION_LENGTHS[6:] == (1, 0, 0, 0), 'conversion final-row proof mismatch'
    assert central.count('let mut y = conversion_6[7];') == 1
    assert central.count('for i in (0..6).rev() {') == 1
    small = (source / 'src/experimental/small.rs').read_text()
    assert small.count('if a <= 0.45 && b < 0.26 {') == 1, 'central outer guard'
    assert 'fn value(self) -> f64 {\n        self.hi + self.lo\n    }' in small
    assert 'fn fma(a: f64, b: f64, c: f64) -> f64 {\n    a.mul_add(b, c)\n}' in small
    assert digest(source / 'src/experimental/lanes.rs') == '40502418e5f024c8a4a522b52e1a9b226e126906bf9aca14496fd56c336d6f32'
    old = (reference / 'src/small/constants.rs').read_text()
    new = (source / 'src/experimental/small/constants.rs').read_text()
    # This binds Pair constructors, the caller, and every unchanged helper to
    # the reviewed complete source. It consumes no seed error certificate.
    seed_reduction.restore_and_validate(small, new, old)
    for name in ['NC', 'NL'] + ['conversion_' + str(k) for k in range(10)]:
        assert exact_table(old, name) == exact_table(new, name), ('central coefficient', name)
    # Positive b below .26 has a positive Black root throughout a<=.45:
    # exp(-a/2)>=1-a/2>.26. The preflight .25 cutoff is conservative because
    # q=sqrt(2pi)(b+a*sinhc(a/2)/2)>=sqrt(2pi)(b+a/2).
    assert F.from_float(.26) < 1 - AMAX / 2
    assert AMAX * AMAX < AMAX2 and QMAX * QMAX < WMAX
    return {str(p.relative_to(source)): digest(p) for p in [
        central_path, source / 'src/experimental/small.rs',
        source / 'src/experimental/small/constants.rs', source / 'src/experimental/lanes.rs']}


def conversion_rows(m):
    """Read exact source-bound coefficients and retain the selected low powers."""
    source = (m.ROOT/'src/small/constants.rs').read_text()
    full = [[(-1 if sign else 1)*m.from_bits(int(bits, 16))
             for sign, bits in exact_table(source, 'conversion_'+str(k))]
            for k in range(10)]
    assert tuple(map(len, full)) == (16, 15, 14, 13, 11, 10, 8, 6, 4, 1)
    assert len(CONVERSION_LENGTHS) == 10
    assert all(0 <= n <= len(row) for n, row in zip(CONVERSION_LENGTHS, full))
    rows = [row[-n:] if n else [] for n, row in zip(CONVERSION_LENGTHS, full)]
    return full, rows


def conversion_bound(m):
    """Recompute K(C)-exp(w/8)T for the shortened polynomial itself."""
    m.AMAX, m.WMAX, m.DEG = AMAX2, WMAX, 24
    P, B, out, exp_b = m.P, m.B, m.out, m.exp_b
    nterms = 24
    A = P.make({(1, 0): AMAX2})
    w = P.make({(0, 1): WMAX})
    full, rows = conversion_rows(m)
    C = P.make({(k, j): value*AMAX2**k*WMAX**j
                for k, row in enumerate(rows) for j, value in enumerate(reversed(row))})
    deleted = sum(abs(value)*AMAX2**k*WMAX**j
                  for k, row in enumerate(full) for j, value in enumerate(reversed(row))
                  if j >= CONVERSION_LENGTHS[k])
    assert C.bound() < F(1, 20)
    wc, cp = w * C, C
    eta = [P()]
    for j in range(1, nterms + 1):
        term = A * cp * ((1 if j % 2 else -1) * F(j + 1, 2))
        if j == 1:
            term = term - w * w * C / 4
        if j == 2:
            term = term - w * w * w * C * C / 8
        eta.append(term)
        cp = cp * wc
    g = [P.as_p(1)]
    for n in range(1, nterms + 1):
        g.append(sum((j * eta[j] * g[n - j] for j in range(1, n + 1)), P()) / n)
    K = C * sum((g[n] / (n + 1) for n in range(nterms + 1)), P())
    R, cm = F(16), C.bound()
    d = WMAX * cm * R
    assert d < 1
    psi = (AMAX2 * cm * R * (2 + d) / (2 * (1 - d)**2)
           + WMAX**2 * cm * R / 4 + WMAX**3 * cm**2 * R**2 / 8)
    ktail = cm * exp_b(B(psi)).hi / (R**nterms * (R - 1))
    K.e = out(K.e + ktail, True)
    tc = {}
    for i in range(nterms + 1):
        for j in range(nterms + 1 - i):
            k = i + j + 1
            tc[j, i] = (F((-1)**i * m.factorial2(2*i + 1),
                          8**k * math.factorial(k) * m.factorial2(2*k + 1))
                        * AMAX2**j * WMAX**i)
    T = P.make(tc)
    mag = max(AMAX2, WMAX) / 8
    ttail = F(1, 8) * mag**(nterms + 1) / math.factorial(nterms + 1) / (1 - mag / (nterms + 2))
    T.e = out(T.e + ttail, True)
    term = ew = P.as_p(1)
    for n in range(1, nterms + 1):
        term = term * w / (8*n)
        ew = ew + term
    ewtail = (WMAX/8)**(nterms + 1) / math.factorial(nterms + 1) / (1 - WMAX / (8*(nterms + 2)))
    ew.e = out(ew.e + ewtail, True)
    error = out((K - ew*T).chebyshev_bound(), True)
    radius = F(1, 2**28)
    segment = cm + radius
    d = WMAX * segment
    psisegment = (AMAX2 * segment * (2 + d) / (2 * (1 - d)**2)
                  + WMAX**2 * segment / 4 + WMAX**3 * segment**2 / 8)
    derivative = exp_b(B(-psisegment)).lo
    assert error < radius * derivative
    relative = out(WMAX * error / derivative / (1 - WMAX * segment), True)
    return {'a_squared_max': str(AMAX2), 'v_squared_max': str(WMAX),
            'conversion_relative_error_bound': str(relative), 'display': float(relative),
            'retained_lengths': CONVERSION_LENGTHS,
            'evaluated_coefficients': sum(map(len, rows)),
            'full_coefficients': sum(map(len, full)),
            'deleted_polynomial_norm': str(deleted),
            'deleted_relative_polynomial_bound_in_u': float(WMAX*deleted*2**53),
            'residual_bound': str(error), 'root_derivative_min': str(derivative),
            'cauchy_radius': str(R), 'cauchy_tail': str(ktail), 'T_tail': str(ttail)}


def conversion_native_bound(m, w, A):
    """Recompute source-order rounding for the shortened Horner chains."""
    _, all_rows = conversion_rows(m)
    rows = all_rows[:7]
    assert all(rows) and not any(all_rows[7:])
    wh, ah = w.high(), A.high()
    values = [m.ordinary_horner(row, wh) for row in rows]
    h, e = values[-1]
    for row, er in values[-2::-1]:
        raw = h*ah+row
        e = e*ah.mag()+er+m.U*raw.mag()+m.NU
        h = raw.widen(m.U*raw.mag()+m.NU)
    # Exact finite-polynomial Lipschitz transfer covers using only coordinate
    # highs. Its derivatives belong to the new retained polynomial.
    wa, aa = max(wh.mag(), w.t.mag()), max(ah.mag(), A.t.mag())
    dw = da = F(0)
    for k, row in enumerate(rows):
        for j, c in enumerate(reversed(row)):
            if j:
                dw += abs(c)*j*wa**(j-1)*aa**k
            if k:
                da += abs(c)*k*aa**(k-1)*wa**j
    e += dw*(w.e+w.l)+da*(A.e+A.l)
    return m.D(h.widen(e), F(0), e)


def normal_shape_bound(m, k, z):
    """Reduced polynomial arithmetic plus its exact deleted-tail enclosure.

    The target of the native arithmetic ledger is the truncated represented
    polynomial. The deleted suffix is then transported to the original model
    at the exact transformed coordinate. Bounding the suffix on the computed
    high and its derivative over the low/error extension explicitly includes
    both the omitted derivative*low term and coordinate error.
    """
    B, D, U, NU = m.B, m.D, m.U, m.NU
    degree = NORMAL_DEGREES[k]
    midhi = [0x3fb47ae147ae147b, 0x3fcd70a3d70a3d71,
             0x3fd6666666666666, 0x3fdccccccccccccd][k]
    ml = m.MID[k] - m.from_bits(midhi)
    scalehi = [F(25, 2), m.from_bits(0x402c924924924925), F(20), F(20)][k]
    sl = m.SCALE[k] - scalehi
    x = m.mul(m.add(z, D(-m.MID[k], abs(ml))), D(m.SCALE[k], abs(sl)))
    hx = x.high()
    assert hx.mag() < F(100000001, 100000000)
    coefficients = [(m.NC[26*k+i], m.NL[26*k+i] if i < 2 else F(0))
                    for i in range(degree+1)]
    # The final compensated stages use FastTwoSum; prove its domination
    # with the shortened high recurrence rather than inheriting it.
    hp = B(coefficients[-1][0])
    for i in range(degree-1, 1, -1):
        p = hp*hx + coefficients[i][0]
        hp = p.widen(U*p.mag()+NU)
    for i in (1, 0):
        prod = (hp*hx).widen(U*(hp*hx).mag()+NU)
        assert abs(coefficients[i][0]) >= prod.mag()
        value = prod + coefficients[i][0]
        hp = value.widen(U*value.mag()+NU)
    _, low, error = m.native_polynomial(coefficients, hx, hx, F(0), F(0), F(0), 2)
    # Source multiplies each derivative coefficient by its integer index
    # before the FMA. Include that multiplication's rounding at every stage.
    derivative = [i*m.NC[26*k+i] for i in range(1, degree+1)]
    dp = B(derivative[-1]).widen(U*abs(derivative[-1])+NU)
    derr = U*abs(derivative[-1])+NU
    for c in derivative[-2::-1]:
        cerror = U*abs(c)+NU
        raw = dp*hx+c
        derr = derr*hx.mag()+cerror+U*(raw.mag()+cerror)+NU
        dp = raw.widen(cerror+U*(raw.mag()+cerror)+NU)
    lownew = (dp.mag()*x.l+low)*(1+U)+NU
    represented = m.COEFFICIENTS[k][:degree+1]
    p1 = [i*represented[i] for i in range(1, degree+1)]
    p2 = [i*(i-1)*represented[i] for i in range(2, degree+1)]
    span = hx.widen(x.l+x.e)
    first = m.polynomial(p1, span).mag()
    second = m.polynomial(p2, span).mag()
    error += derr*x.l + abs(m.NL[26*k+1])*x.l + second*x.l*x.l/2
    error += U*(dp.mag()*x.l+low)+NU + first*x.e
    if degree == 25:
        # The unchanged branch reproduces the archived native ledger exactly.
        assert (lownew, m.out(error, True)) == m.normal_shape(k, z)
    # Coefficients are exact binary64 rationals, already compared with the
    # archive by source_binding. The low coefficient table is used only at
    # i=0,1, hence every deleted coefficient is exactly its NC high.
    tail_high = sum(abs(m.COEFFICIENTS[k][i])*hx.mag()**i
                    for i in range(degree+1, 26))
    tail_derivative = sum(i*abs(m.COEFFICIENTS[k][i])*span.mag()**(i-1)
                          for i in range(degree+1, 26))
    tail_low = tail_derivative*x.l
    tail_coordinate = tail_derivative*x.e
    tail = m.out(tail_high+tail_low+tail_coordinate, True)
    # This is a whole-interval coefficient-norm bound, not a sampled error.
    assert tail < U/100
    return lownew, m.out(error+tail, True), {
        'degree': degree,
        'saved_simd_horner_stages': 25-degree,
        'transformed_high_magnitude': str(hx.mag()),
        'transformed_span_magnitude': str(span.mag()),
        'deleted_tail_high_bound': str(m.out(tail_high, True)),
        'deleted_derivative_times_low_bound': str(m.out(tail_low, True)),
        'deleted_coordinate_error_bound': str(m.out(tail_coordinate, True)),
        'deleted_tail_bound': str(tail),
        'deleted_tail_in_u': float(tail/U),
        'native_error_before_tail': str(m.out(error, True)),
    }


def native_bound(m, model, conversion):
    """Recompute native helper error, including the final binary64 rounding."""
    m.AMAX, m.QMAX = AMAX, QMAX
    B, D, U, NU = m.B, m.D, m.U, m.NU
    # The archived sinhc() used a .006 upper reserve for a<=.36. Recompute
    # that same polynomial/rounding/tail ledger; only its bound is widened.
    w = (AMAX / 2)**2
    true = [F(1, math.factorial(2*k + 1)) for k in range(9)]
    c = [(F(1), F(0)),
         (m.from_bits(0x3fc5555555555555), m.from_bits(0x3c65555555555555)),
         (m.from_bits(0x3f81111111111111), m.from_bits(0x3c01111111111111)),
         (m.from_bits(0x3f2a01a01a01a01a), m.from_bits(0x3b6a01a01a01a01a))]
    c += [(m.rn(true[k]), F(0)) for k in range(4, 9)]
    X, dh = B(0, w), U*w + NU
    _, low, error = m.native_polynomial(c, X, X.widen(dh), dh, NU, dh, 4)
    coefficient_error = sum(abs(a+b-t)*w**i for i, ((a, b), t) in enumerate(zip(c, true)))
    tail = w**9 / math.factorial(19) / (1 - w/(20*21))
    upper = sum(t*w**k for k, t in enumerate(true)) + tail
    assert upper < F(1009, 1000)
    chi = m.renorm(D(B(1, upper), low, error + coefficient_error + tail))
    assert 2 / m.SQRT2PI.lo < F(4, 5) and 1 / m.SQRT2PI.lo < F(2, 5)
    a, b = D(B(0, F(4, 5))), D(B(0, F(2, 5)))
    d = m.mul(a, chi).retarget(B(0, F(4, 5)))
    kh, kl = m.from_bits(0x40040d931ff62706), -m.from_bits(0x3caa6a0d6f814637)
    kerror = max(abs(kh+kl-m.SQRT2PI.lo), abs(kh+kl-m.SQRT2PI.hi))
    q = m.mul(D(m.SQRT2PI, abs(kl), kerror), m.add(b, d.scale(F(1, 2)))).retarget(B(1))
    qerr = q.e + q.l
    assert F.from_float(1e-12)/(1+qerr) > m.QMIN
    assert QCAP/(1-qerr) < QMAX and QMAX**2 < WMAX
    assert m.SQRT2PI.lo/(4*(1+U)) > QCAP * (1+qerr), 'conservative .25 preflight including addition rounding'
    # If b<a/16, then b/d<1/16 and z>1/(sqrt(2pi)*9/16)^2>.503.
    assert 1/(m.SQRT2PI.hi*F(9, 16))**2 > F(503, 1000)
    zbase = m.mul(m.div(d, q), m.div(d, q))
    assert zbase.e + zbase.l < m.ZPAD
    rows = []
    for k in range(4):
        cells = model['branches'][str(k)]['cells']
        z = zbase.retarget(B(F(cells[0]['lo']), F(cells[-1]['hi'])))
        low, error, normal_record = normal_shape_bound(m, k, z)
        shape_error = max(F(c['shape_error_bound']) for c in cells)
        shape = B(min(F(c['shape_min']) for c in cells), max(F(c['shape_max']) for c in cells))
        normal = D(shape, low, error + shape_error)
        v = m.div(m.mul(q, normal), chi)
        wscaled = m.mul(v, v)
        assert wscaled.t.hi*QMAX**2 <= WMAX, 'enlarged conversion model covers computed normal shape range'
        w = D(B(0, WMAX), wscaled.l*QMAX**2, wscaled.e*QMAX**2)
        A = m.mul(D(B(0, AMAX)), D(B(0, AMAX)))
        C = conversion_native_bound(m, w, A)
        result = m.mul(v, m.add(D(1), m.mul(w, C)))
        assert result.high().lo > result.l and result.t.lo > 0
        relative = result.e / result.t.lo
        model_error = F(conversion['conversion_relative_error_bound'])
        pre = relative + model_error*(1+relative)
        rho = (U+(1+U)*pre)/(2*U)
        assert result.t.lo*m.QMIN > F(1, 2**1022)
        assert result.h()*QMAX < 10 and rho < 1
        rows.append({'normal_branch': k, 'rho_bound': str(m.out(rho, True)),
                     'display': float(rho), 'native_relative_error': str(m.out(relative, True)),
                     'normal_polynomial': normal_record})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Proof checks require assertions; do not use -O')
    source, reference = args.source_root.resolve(), args.reference_root.resolve()
    for name, expected in DEPENDENCIES.items():
        assert digest(reference/name) == expected, ('central inherited binding', name)
    source_hashes = source_binding(source, reference)
    sys.path.insert(0, str(reference/'scripts'))
    conversion_module = importlib.import_module('prove_nbt_conversion')
    native_module = importlib.import_module('prove_nbt_native')
    assert conversion_module.ROOT.resolve() == native_module.ROOT.resolve() == reference
    model = json.loads((reference/'results/global_rho_proof/nbt_normal_model.json').read_text())
    assert model['status'] == 'model_bound_complete' and len(model['branches']) == 4
    for name, expected in model['source_sha256'].items():
        assert digest(reference/name) == expected, ('normal-shape certificate', name)
    conversion = conversion_bound(conversion_module)
    rows = native_bound(native_module, model, conversion)
    result = {
        'status': 'All exact rational assertions passed',
        'scope': 'Deferred central helper, a<=binary64(.45), b<binary64(.26), q.hi<=.625, z.hi<=.5',
        'method': 'Exact represented normal-polynomial tails on padded cells; recomputed polynomial/Cauchy remainder for shortened Black conversion and all native source-order errors',
        'conversion_model': conversion, 'branches': rows, 'max_rho': max(r['display'] for r in rows),
        'source_sha256': source_hashes, 'dependency_sha256': DEPENDENCIES,
        'script_sha256': {p.name: digest(p) for p in [Path(__file__), Path(__file__).with_name('central_expansion.py'), Path(__file__).with_name('seed_reduction.py')]},
        'conditional_inherited_premises': [
            'Hash-bound exact normal-shape certificate; its four model cell covers are unchanged, with deleted coefficient tails separately enclosed.',
            'Reviewed pair arithmetic and shortened normal/conversion Horner chains; fixed independent SIMD lanes preserve each retained rounding point.',
            'IEEE RN-even binary64 basic operations and single-rounding FMA, gradual underflow, source order preserved.'],
        'limits': 'Conditional continuous helper proof; no whole-domain solver theorem, Lean completion, vendor libm claim, or speed claim.'}
    out = args.out or source/'target/proof/central-expansion-proof.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'max_rho': result['max_rho'], 'output': str(out)}, indent=2))


if __name__ == '__main__':
    main()
