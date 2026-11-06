#!/usr/bin/env python3
"""Source-bound conditional certificate for the large direct log-price inverse.

The new finalizer maps the compensated log-price residual directly to relative
volatility with a cubic Taylor polynomial.  All numerical assertions below use
exact Fractions (or explicit outward rational rounding), not sampled roots.
The Mills, target/cap and routing lemmas are inherited from explicitly bound
analytic certificates.  This is not a Lean theorem or a vendor-libm guarantee.

Write outputs only below --out.  The reference checkout is read-only.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import re
import struct
import sys
import tomllib

sys.dont_write_bytecode = True
if not __debug__:
    raise RuntimeError("Assertions are required; do not use Python -O")

U, NU, TAU = F(1, 2**53), F(1, 2**700), F(1, 2**1075)
REFERENCE_FILES = {
    'src/large.rs': '15ecde382af035a18841805d583dd9db9a3aa4ea0b0b41c8d44d22b8ef8d18fa',
    'src/large_constants.rs': 'a114f8123411950b1f92f3f6d3f447206dd7b3143d0a87361df03ffdf36ad02f',
    'scripts/prove_global_mills.py': '2dce1c923d5ba36d7109ffdc3a17f4c726e2996c2a50d6ee49663dfb436e7547',
    'scripts/prove_global_large_source.py': '23e8d407539fb2f489176ccef3486f81ae8e664665882acde8dd73b638b885f4',
    'scripts/prove_large_finish.py': 'aa196629febe0647cce08debaee629688e2c54891422ef5db63d054bf534fb43',
    'scripts/prove_large_targets.py': '8b9e1c73935860ca8296dd8402031714185ec3df33a37bf5bb35f4e23f78e646',
    'scripts/prove_large_global.py': '8129da0dfbd5482ca61f9f2ea2ae95a2b6dbc02701aa153b5b9c7fc959d86049',
    'scripts/prove_global_domain_and_integer.py': 'ba62d9a216c61e6a91c5c0b8b6ca50c7bce8f769dc7330426ff0851b377f9bce',
    'scripts/prove_large_cap.py': '4a103ff91f5c19129e4e78cc303678ab2c43d06f6fdc551e4dbc2838a51e26dc',
    'reference/iv_large_a_cap_table_v50/scripts/exact.py': 'cd66c85457348fabdc238e7395f6c4e10bc954f23a5af7ab709c43c08ac03559',
    'notes/global-rho-large-source.md': '3fd1fbae6b2d4ee027fe7a154a4984e8ba9a888bc68649b150f520f7ac08a026',
    'notes/global-rho-large-targets.md': '3555290780bf11e73e7d59ce710a4b371ca918100af6aaa97f38cdeef4cbcc39',
    'notes/global-rho-large-review.md': 'b40d7d5cb1f0f2cb4ae268788d89db2c16a71a6c5d6066b141540b9117359e4d',
    'results/global_proof/mills_budget.json': 'f23d4ae8275bd608c66a1bf201d79b81afc991c838d05332dc3dc18c81963ae3',
    'results/global_proof/large_domain_integer.json': '282c526b98d57ede99e34a877a224962d54917037360ac040a3b5372c3c32b0b',
    'results/global_rho_proof/large_cap_bound.json': 'd066fe4d7ffbd8a69e52bf063463fe9c12e0272fc08742fed54c82c794304d86',
    'results/global_rho_proof/large_global_bound.json': '5efe00d255665c7b6075bd8195ebf53f701cdce05bd494a55ea4b583244333a4',
    'results/global_rho_proof/TMDworstcases.pdf': '718e25ca37ab8d826586fb15dd127c578fcae9ec60c4d5caee82120698d99c4f',
}
SOURCE_FILES = {
    'src/experimental/mills_constants.rs': '0fda9d58f198c692908879ead2020d62bc0b2859fb2abbab1913fa8f06ed3063',
    'src/experimental/large_constants.rs': '9f4f8dc26f8c1313e60e7e2817a24438a609d73ab4ba130a7962902c35be3efd',
    'src/experimental/lanes.rs': '40502418e5f024c8a4a522b52e1a9b226e126906bf9aca14496fd56c336d6f32',
}
FAST_FUNCTIONS = {
    'mills_fast': '78af8b69bcbd6a5ced95d45a03ac7a61a5dece11da386eb9c336a0d64165b8ca',
    'fast_pair_impl': 'd252df351d482d6f6ec10f0b006b5a691bbcbb3ce588573e146f95723f9c096a',
    'mills_fast_pair': '8e1deeb5de3e24aaee8e2ee2a53fbf27daa2c6bf07c534fb822fa1bb18c2d28f',
    'eval_fast': '8873a86ffc761a8a1372b33da7c0555301f848597dd53bc02fd7935d6eed210a',
}
SIMD_HASH = '65450215ad42e632e66e97a0f2f8229e0934441b50b9fcca483c2cad80f5af75'
# This reviewed fast checker recomputes the degree-13 model/source arithmetic.
# Its certificate is consumed only for functions in FAST_FUNCTIONS, whose
# current source and constants are compared independently.  No old finalizer
# or whole-file hash is accepted as evidence for the new finalizer.
FAST_VERIFIER_HASHES = {
    # Degree-13 table source/model recomputed with the same residual envelope.
    'a996f269f620916a88804701d1d495b04cad1265c7f9c6c7504a6cbedbf72028',
    # Same fast numerical calculation; updated precise-component source binding.
    '8f0624fb88456589d5ee3f89bd447aeb3371e4a18bcfbbbcefa055944e8b81d7',
    '93325b560759bf6e8f4af17ec513109615ecb0a472195cab70764c69d41f3a59',
    # Same numerical calculation, scoped to fast-component source and caller
    # premises; finalizer/routing bindings are owned by this direct checker.
    'd82e349a0d95a3c679234ad1d938a7b0906ab405e10467404aebf3af0b1e5732',
}
# Bind the reviewed exact numerical certificate as well as its metadata.
# This payload excludes only source integration/routing metadata, allowing a
# reviewed source bridge to regenerate the same fast-component certificate.
FAST_CERTIFICATE_FIELDS = ('maximum_residual_error', 'maximum_residual_error_in_u',
                           'table_errors', 'asymptotic_error', 'boxes', 'coefficient_identity')
FAST_CERTIFICATE_PAYLOAD_HASH = 'b15b040380e65c999e17a3f7e1046e576a8cd6847284e42207cddd2d892e68c8'

FINISH_SOURCE = """fn finish_direct(x: f64, pe: Precise) -> f64 {
    let yy = val(pe.y);
    let xy = x * yy;
    let alpha = -xy * pe.lambda;
    let beta = (x * x + yy * yy - xy) * (pe.lambda * pe.lambda);
    let c2 = 0.5 * (1.0 - alpha);
    let c3 = (alpha.mul_add(2.0 * alpha - 3.0, 1.0) + beta) / 6.0;
    let z = -val(pe.residual);
    let delta = pe.lambda * (z * z.mul_add(z.mul_add(c3, c2), 1.0));
    let low = pe.s.hi.mul_add(delta, pe.s.lo);
    pe.s.hi + low
}"""

NORMAL_LOG_SOURCE = """fn log_positive_normal_price(b: f64) -> DD {
    debug_assert!(b.is_normal() && b > 0.0);
    let bits = b.to_bits();
    let mut e = ((bits >> 52) & 2047) as i32 - 1022;
    let mut m = f64::from_bits((bits & 0x800fffffffffffff) | (1022u64 << 52));
    if m < f64::from_bits(0x3fe6a09e667f3bcd) {
        m *= 2.0;
        e -= 1;
    }
    let t = two(m.ln(), 0.0);
    add(scale(LN2, e as f64), t)
}"""

LOG_STEP3_SOURCE = """fn log_step3(a: f64, x: f64, e: &Eval) -> f64 {
    let n = -e.f / e.df;
    let d = e.df;
    let h = e.h;
    let iy = 1.0 / e.y;
    let iy2 = iy * iy;
    let h1 = -1.0 - x * iy * iy2 + (x * x - 2.0 * a) * iy2 * iy2;
    let a2 = h - d;
    let a3 = h1 + a2 * (h - 2.0 * d);
    let num = (0.5 * a2).mul_add(n, 1.0);
    let den = n.mul_add((a3 / 6.0).mul_add(n, a2), 1.0);
    n * num / den
}"""


def digest(data):
    return hashlib.sha256(data).hexdigest()


def function_body(text, name):
    match = re.search(r'\bfn ' + re.escape(name) + r'(?:\(|<)', text)
    assert match, ('missing function', name)
    start = match.start()
    i = text.index('{', start) + 1
    depth = 1
    while depth:
        depth += (text[i] == '{') - (text[i] == '}')
        i += 1
    return text[start:i]


def canonical(text):
    # Compared numerical functions contain no string literals.  Permit only
    # whitespace/comments to vary, never operators, constants or source order.
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
    text = re.sub(r'//[^\n]*', '', text)
    return re.sub(r'\s+', '', text)


def replace_once(text, old, new):
    assert text.count(old) == 1, ('source bridge pattern', old)
    return text.replace(old, new, 1)


def upward(q, precision=512):
    """Directed rational rounding to at least precision significant bits."""
    q = F(q)
    assert q >= 0
    if not q:
        return q
    shift = precision - (q.numerator.bit_length() - q.denominator.bit_length())
    if shift >= 0:
        n, d = q.numerator << shift, q.denominator
        return F((n + d - 1) // d, 1 << shift)
    n, d = q.numerator, q.denominator << -shift
    return F((n + d - 1) // d) * (1 << -shift)


def pack(q):
    q = upward(q)
    return {'exact_upper': str(q), 'preview': float(q)}


def qfield(value):
    return F(value['exact'])


@dataclass
class Ball:
    """|exact real expression| <= mag; |source-expression| <= err."""
    mag: F
    err: F = F(0)

    def add(self, other):
        m = self.mag + other.mag
        e = self.err + other.err
        return Ball(upward(m), upward(e + U * (m + e) + NU))

    def mul(self, other):
        m = self.mag * other.mag
        e = self.mag * other.err + other.mag * self.err + self.err * other.err
        return Ball(upward(m), upward(e + U * (m + e) + NU))

    def fma(self, other, addend):
        m = self.mag * other.mag + addend.mag
        e = self.mag * other.err + other.mag * self.err + self.err * other.err + addend.err
        return Ball(upward(m), upward(e + U * (m + e) + NU))

    def scale(self, value):
        # Exact power-of-two scaling, plus absolute underflow.  Multiplication
        # by 2 cannot underflow; the same harmless NU also covers scaling by .5.
        assert value in (F(1, 2), F(2))
        return Ball(upward(abs(value) * self.mag), upward(abs(value) * self.err + NU))

    def divide_constant(self, value):
        m, e = self.mag / abs(value), self.err / abs(value)
        return Ball(upward(m), upward(e + U * (m + e) + NU))

    def sharpen(self, value):
        # This replaces only the analytically proved exact magnitude.
        # Source errors already accumulated above are never discarded.
        return Ball(F(value), self.err)


class Poly:
    """Tiny exact polynomial ring in four named-at-call-site variables."""
    def __init__(self, terms=None):
        self.terms = {k: F(v) for k, v in (terms or {}).items() if v}

    @staticmethod
    def constant(value):
        return Poly({(0, 0, 0, 0): F(value)})

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = Poly.constant(other)
        d = self.terms.copy()
        for k, v in other.terms.items():
            d[k] = d.get(k, F(0)) + v
        return Poly(d)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -other if isinstance(other, Poly) else self + -F(other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Poly):
            other = Poly.constant(other)
        d = {}
        for a, av in self.terms.items():
            for b, bv in other.terms.items():
                key = tuple(x + y for x, y in zip(a, b))
                d[key] = d.get(key, F(0)) + av * bv
        return Poly(d)

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        return self * (F(1) / F(divisor))

    def __pow__(self, power):
        result = Poly.constant(1)
        for _ in range(power):
            result = result * self
        return result

    def partial(self, variable):
        d = {}
        for key, value in self.terms.items():
            if key[variable]:
                new = list(key)
                new[variable] -= 1
                d[tuple(new)] = value * key[variable]
        return Poly(d)

    def __eq__(self, other):
        return isinstance(other, Poly) and self.terms == other.terms


def symbolic_log_step3():
    """Check Householder d=3 and its local inverse series in exact arithmetic.

    For g=1/f, 3*g''/g''' is the update. With n=-f/f',
    A=f''/f' and B=f'''/f', it is n*(1+A*n/2)/(1+A*n+B*n²/6).
    This establishes the ideal formula, not its global convergence or a
    one-step residual bound. The source-level safeguard supplies termination.
    """
    n, a, b, c = [Poly({tuple(int(i == j) for i in range(4)): F(1)})
                  for j in range(4)]
    numerator = n*(1+a*n/2)
    denominator = 1+a*n+b*n**2/6
    # After substituting f=-n*f', f^4*g'''=-f'^3*(6+6An+Bn²),
    # while f^3*g''=f'^2*(2+An). Cross-multiply the reciprocal formula.
    assert numerator*(6+6*a*n+b*n**2) == 3*n*(2+a*n)*denominator
    inverse_cubic = n-a*n**2/2+(a**2/2-b/6)*n**3
    residual = inverse_cubic*denominator-numerator
    assert all(key[0] >= 4 for key in residual.terms)
    # Normalize f(x+step)/f'=-n+step+A*step²/2+B*step³/6+C*step⁴/24+...
    root_residual = (-n+inverse_cubic+a*inverse_cubic**2/2
                     +b*inverse_cubic**3/6+c*inverse_cubic**4/24)
    assert all(key[0] >= 4 for key in root_residual.terms)
    # If d=f' and H=p''/p', then d'=d*(H-d). Differentiation gives
    # A=H-d and B=H'+A*(H-2d), as implemented by a2 and a3.
    d, h, h1 = a, b, c
    a2 = h-d
    assert h1-d*a2+a2**2 == h1+a2*(h-2*d)
    return {'method': 'Householder reciprocal derivative d=3',
            'ideal_update': 'n*(1+A*n/2)/(1+A*n+B*n*n/6)',
            'local_order': 4,
            'identities_checked': ['3*(1/f)second/(1/f)third',
                                   'inverse Taylor agreement through n^3',
                                   'log third derivative ratio'],
            'global_one_step_bound_claimed': False}


def symbolic_physical_predictor():
    """Exact reciprocal-series basis for every adaptive predictor branch.

    With n=-f/f' and q_j=f^(j)*n^(j-1)/(j!*f'),
    f(x+n*t)/f(x)=1-t-q_2*t²-q_3*t³-... . The reciprocal
    coefficients c_j obey the source recurrence. Householder's d-th
    reciprocal derivative therefore gives n*c_(d-1)/c_d.
    """
    q2, q3, q4, q5 = [Poly({tuple(int(i == j) for i in range(4)): F(1)})
                      for j in range(4)]
    one = Poly.constant(1)
    q = [Poly(), one, q2, q3, q4, q5]
    c = [one, one]
    c.append(q2*c[0]+c[1])
    c.append(q3*c[0]+q2*c[1]+c[2])
    c.append(q4*c[0]+q3*c[1]+q2*c[2]+c[3])
    c.append(q5*c[0]+q4*c[1]+q3*c[2]+q2*c[3]+c[4])
    assert c[2] == 1+q2
    assert c[3] == 1+2*q2+q3
    assert c[4] == 1+3*q2+q2**2+2*q3+q4
    assert c[5] == 1+4*q2+3*q2**2+3*q3+2*q2*q3+2*q4+q5
    for degree in (3, 4, 5):
        # Coefficients through t^degree of the series times its reciprocal
        # must be exactly (1,0,...,0). No finite numerical samples enter.
        assert c[0] == one
        for k in range(1, degree+1):
            assert c[k]-sum((q[j]*c[k-j] for j in range(1, k+1)), Poly()) == Poly()
    return {'method': 'Adaptive Householder reciprocal derivatives d=3/4/5',
            'ideal_updates': ['n*c2/c3', 'n*c3/c4', 'n*c4/c5'],
            'selection': ['abs(x*n)<=1/8', '1/8<abs(x*n)<=3/8', 'otherwise'],
            'ideal_local_orders_given_exact_newton_ratio': [4, 5, 6],
            'reciprocal_polynomial_identities_checked_exactly': True,
            'global_one_step_bound_claimed': False}


def adaptive_physical_source(original):
    """Retain exact source expressions, staging only independent arithmetic.

    The original physical recurrence is file-hash bound. Each retained let
    expression is copied unchanged, so the operation/rounding order within
    it is preserved. The reordered expressions are pure scalar arithmetic;
    the only function call, mills_fast(y), still runs first and exactly once.
    """
    body = function_body(original, 'physical_step7')

    def statement(name):
        matches = re.findall(r'    let '+name+r' = .*?;\n', body, flags=re.S)
        assert len(matches) == 1, ('physical recurrence expression', name)
        return matches[0]

    stages = [
        ('n', 'iy', 'k', 'v', 'nn', 'xn', 'g0', 'g1', 'w0', 'g2',
         'v1', 'w1', 'v2', 'w2', 'q2', 'q3', 'c0', 'c1', 'c2', 'c3'),
        ('g3', 'v3', 'w3', 'q4', 'c4'),
        ('g4', 'v4', 'w4', 'q5', 'c5'),
    ]
    return ('fn physical_step_adaptive(x: f64, y: f64, s: f64) -> f64 {\n'
            + ''.join(map(statement, stages[0]))
            + '    if xn.abs() <= 0.125 {\n        return n * c2 / c3;\n    }\n'
            + ''.join(map(statement, stages[1]))
            + '    if xn.abs() <= 0.375 {\n        return n * c3 / c4;\n    }\n'
            + ''.join(map(statement, stages[2]))
            + '    n * c4 / c5\n}')


def symbolic_and_shape(ex):
    # Write h=a/s, t=s/2. Then x=t-h, y=t+h,
    # A=-xy=h^2-t^2 and Q=x^2+y^2-xy=3h^2+t^2.
    # For either p=normalized price or its complement, p''/p'=A/s;
    # lambda=p/(s*p') has its sign. Since d/dq=s*lambda*d/ds,
    # lambda'=lambda*(1-lambda-alpha),
    # A'=(A-Q)*lambda and Q'=(-6h^2+2t^2)*lambda.
    # Therefore alpha'=alpha-alpha^2-beta and
    # beta'=2beta-2alpha*beta-12h^2*lambda^3. The polynomial
    # checks below derive the cubic and remainder from these identities.
    variables = []
    for i in range(4):
        key = [0] * 4
        key[i] = 1
        variables.append(Poly({tuple(key): F(1)}))
    lam, alpha, beta, gamma = variables
    deriv = (lam * (1 - lam - alpha), alpha - alpha**2 - beta,
             2 * beta - 2 * alpha * beta - gamma, Poly())

    def dq(poly):
        return sum((poly.partial(i) * deriv[i] for i in range(4)), Poly())

    t2 = lam * lam + dq(lam)
    t3 = lam * t2 + dq(t2)
    t4 = lam * t3 + dq(t3)
    assert t2 == lam * (1 - alpha)
    assert t3 == lam * (1 - 3 * alpha + 2 * alpha**2 + beta)
    assert t4 == lam * (1 - 7 * alpha + 12 * alpha**2 + 6 * beta
                         - 6 * alpha**3 - 7 * alpha * beta - gamma)

    # M(z)=integral exp(-zt-t^2/2)dt and I1=-M' are inherited real-analysis
    # lemmas. I1(t)<=1 for t>=0 and I1(t)<1/t^2 for t>0. M(t)<1/t.
    # M(0)=sqrt(2pi)/2<4/3; M(-1)<exp(.5)(M(0)+1)<4.
    assert ex.K[1] / 2 < F(4, 3)
    assert ex.expq(F(1, 2))[1] < F(5, 3)
    assert F(5, 3) * (F(4, 3) + 1) < 4
    # All reachable y=sqrt(x^2+2a)>sqrt(20)>4. Lower x in[-40,1]:
    # x=-v<=0, s=y-v, D=integral_v^y I1(t)dt, so 0<lambda<1.
    # v>0 additionally D<s/(vy), hence 0<alpha=vy*lambda<1.
    # If v>=1: beta=(v^2+y^2+vy)lambda^2<3.
    # If 0<=v<1: D<4/3, beta<(16/9)(y^2+y+1)/(y-1)^2<=112/27.
    # The last inequality follows from (y-4)(4y-1)>=0 for y>=4.
    assert F(16, 9) * F(7, 3) == F(112, 27) < 5
    # Lower 0<=x<=1: D<M(-1)<4, s>4, |alpha|=xyD/s<4,
    # Q=x^2+y^2-xy<s^2, so beta<D^2<16 and lambda<1.
    # Upper x>=0: |lambda|=D/s<(4/3+1/4)/4<2/5;
    # M(x)<1/x and M(y)<1/y give 0<=alpha=xyD/s<1 (x=0 is exact).
    assert (F(4, 3) + F(1, 4)) / 4 < F(2, 5)
    assert (F(4, 3) + F(1, 4))**2 == F(361, 144) < 16
    # gamma=12h^2*lambda^3; Q>=3h^2 gives |gamma|<=4 beta |lambda|<64.
    fourth = (1 + 7*4 + 12*4**2 + 6*16 + 6*4**3 + 7*4*16 + 64) / F(24)
    assert fourth == F(1213, 24)
    return fourth


def accepted_region_shape(ex, gate_bits):
    """Localize the WHOLE inverse-Taylor segment using the bound caller.

    The source binding below fixes ell=add(log_positive_normal_price(b),dd(a/2)),
    upper=ell.hi>-.6, target=ell on the lower route, and the final two() in
    add(). Its exact represented target is ell.hi+ell.lo, not just ell.hi.
    The primitive TwoSum lemma gives |lo| <= (u*|hi|+TAU)/(1-u). For hi<=-c
    this implies target<=-c*(1-u/(1-u))+TAU/(1-u), including underflow.
    Neither the primitive error in log(b) nor the original physical root is
    substituted here: this segment ends at the exact internal DD target.
    Original-input transfer keeps its older, wider lower-price allowance.
    """
    radius = F(1, 2**gate_bits) + F(1, 2**38)
    switch = F(float(.6))
    lower_target_upper = -switch*(1-U/(1-U))+TAU/(1-U)
    assert lower_target_upper < -F(3, 5)+4*U
    q_upper = lower_target_upper+radius
    assert radius+4*U < F(1, 1000)
    assert q_upper < -F(599, 1000)
    # Both endpoints q(current), q(target), and hence every intermediate q,
    # are <=q_upper. The precise high error is <2^-38. The fast bypass is
    # separately proved below to enter this same precise high-component gate.
    # exp(q_upper)<.55 follows from a positive Taylor lower bound for exp(.599).
    from math import factorial
    assert sum(F(599, 1000)**j/F(factorial(j)) for j in range(6)) > F(20, 11)

    # At x=c=.35, y=sqrt(c²+2a)>112/25 for every a>10. With
    # Phi(c)=1/2+phi(0)*integral_0^c exp(-t²/2)dt, M(y)<1/y,
    # exp(-w)>=1-w and exp(-w)<=1-w+w²/2, the following is a
    # rational lower bound for p=Phi(c)-phi(c)*M(y), not a sampled price.
    c, y_lower = F(7, 20), F(112, 25)
    assert y_lower*y_lower < c*c+20
    assert ex.K[1] < F(100, 39)  # phi(0)>39/100
    integral_lower = c-c**3/6
    exp_upper = 1-c*c/2+c**4/8
    bracket = integral_lower-exp_upper/y_lower
    assert bracket > 0
    price_lower = F(1, 2)+F(39, 100)*bracket
    assert price_lower > F(11, 20)
    # Lower p increases with x and a: M decreases and dp/dx=s*phi(x)/y>0.
    # Thus all points on the accepted lower segment have x<.35.

    # For lower x=-v<=0, lambda is the interval average of I1=-M' on
    # [v,sqrt(v²+2a)]. I1 decreases because I2>0. Parametrizing that interval
    # by theta gives argument>theta*sqrt20, so lambda<M(0)/sqrt20<1/3.
    assert ex.K[1]/2 < F(4, 3)
    assert 20 > 4**2
    # For lower 0<=x<.35, M(-x)<=M(-.35)<7/4, by its integral identity.
    assert ex.K[1]/2 < F(63, 50)
    mills_negative_upper = (F(63, 50)+c)/(1-c*c/2)
    assert mills_negative_upper == F(1288, 751) < F(7, 4)
    assert F(35, 8)**2 < 20
    assert F(7, 4)/F(35, 8) == F(2, 5)
    # Therefore lambda=D/s<M(-.35)/sqrt20<2/5 and
    # |alpha|=x*y*D/s<x*D<1. Q=x²+y²-xy<=s² gives beta<=D²<5.
    assert c*F(7, 4) < 1 and F(7, 4)**2 < 5
    # On lower x<=0, the unchanged Mills inequalities above give alpha<1,
    # beta<3 when v>=1 and beta<112/27<5 when v<1. On upper x>=0,
    # |lambda|<2/5, alpha<1 and beta<=D²<(4/3+1/4)²<5.
    assert F(112, 27) < 5 and (F(4, 3)+F(1, 4))**2 < 5
    mu, alpha, beta = F(2, 5), F(1), F(5)
    gamma = 4*beta*mu  # Q>=3h² => |12h²*lambda³|<=4*beta*|lambda|.
    fourth = mu*(1+7*alpha+12*alpha**2+6*beta+6*alpha**3+7*alpha*beta+gamma)/24
    assert gamma == 8 and fourth == F(33, 20)
    return {'lower_target_upper': lower_target_upper, 'true_residual_radius': radius,
            'lower_segment_q_upper': q_upper, 'price_lower_at_x_0_35': price_lower,
            'mills_upper_at_minus_0_35': mills_negative_upper,
            'lambda_abs': mu, 'alpha_abs': alpha, 'beta': beta, 'gamma_abs': gamma,
            'fourth_derivative_bound_over_factorial': fourth}



def reduced_precise_budget(mp):
    """Degree-16 model and exact existing compensated source rounding ledger.

    Leading zero Horner stages are exact for finite qhi: zero times qhi plus
    zero stays zero; inserting c16 is exact. Thus a zero-padded degree-24
    evaluation has the same nonzero source operations as the shortened
    table lane. The inherited ledger additionally charges harmless NU terms
    for those omitted stages. Mixed lanes similarly zero-pad 18 and 17.
    All three compensated stages, low-argument correction, degree-18
    asymptotic evaluation and DD pair formation retain their source order.
    """
    ex, up, u, nu = mp.ex, mp.up, mp.U, mp.NU
    rows, asym = mp.coefficients()
    # Recheck the inherited analytic magnitude premises consumed by
    # finish_low and the asymptotic/source-order ledgers.
    assert max(abs(v) for row in rows for pair in row for v in pair) < 2**80
    assert max(map(abs, asym)) < 2**80
    assert ex.expq(F(51, 100))[1] < F(5, 3)
    assert F(5, 3)*(F(4, 3)+F(1001, 1000)) < 4
    assert 1+F(1001, 1000)*4 < 6
    assert (1+F(1001, 1000)**2)*4+F(1001, 1000) < 12
    assert 68*u/(1-u)+u/4+nu < F(1, 1000)
    shortened = [row[:17]+[(F(0), F(0))]*8 for row in rows]
    plain = [mp.table_bound(k, shortened[k+12], False) for k in range(-2, 32)]
    low = [mp.table_bound(k, shortened[k+12], True) for k in range(-2, 32)]
    ap, al = mp.asymptotic_bound(asym, False), mp.asymptotic_bound(asym, True)
    # Exact coefficient difference to the represented degree-24 polynomial;
    # local centering rounding can enlarge the mathematical radius slightly.
    radius = up((1+u)/4+nu)
    tails = [up(sum(abs(rows[k+12][j][0])*radius**j for j in range(17, 25)))
             for k in range(-2, 32)]

    def reachable(bounds, lo, hi, asym_bound):
        # floor(2*zhi) selects [k/2,(k+1)/2). A lower endpoint at
        # a seam cannot reach the preceding cell. This matters at x=0.
        choices = [bounds[k+2] for k in range(-2, 32)
                   if hi >= F(k, 2) and lo < F(k+1, 2)]
        if hi >= 16:
            choices.append(asym_bound)
        assert choices
        return max(choices)

    def boxes(accepted):
        result = []
        for upper, count in ((False, 82), (True, 32)):
            for i in range(count):
                left = F(i, 2)+(0 if upper else -40)
                right = left+F(1, 2)
                if accepted and not upper:
                    if left >= F(7, 20):
                        continue
                    right = min(right, F(7, 20))
                xa = min(abs(left), abs(right)) if left*right > 0 else F(0)
                xb = max(abs(left), abs(right))
                yl, yh = ex.sqrtq(xa*xa+20)[0], ex.sqrtq(xb*xb+2834)[1]
                sl = ex.sqrtq(left*left+20)[0]+left
                zl, zr = (left, right) if upper else (-right, -left)
                ef = reachable(plain, zl, zr, ap)
                ey = reachable(low, yl*(1-2*u)-nu, yh*(1+2*u)+nu, al)
                ed = up(ef+ey+65536*u*u+nu)
                if upper:
                    dl = ex.mills(right)[0]+yh/(yh*yh+1)
                elif right <= 0:
                    v = -left
                    dl = 20/((ex.sqrtq(v*v+20)[1]+v)*(v*v+23))
                else:
                    dl = F(3, 4)
                tail = (reachable(tails, zl, zr, F(0))+
                        reachable(tails, yl*(1-2*u)-nu, yh*(1+2*u)+nu, F(0)))
                result.append({'upper': upper, 'x': [str(left), str(right)],
                               'D_error_over_u': ex.pack(up(ed/u)),
                               'error_over_D_in_u': ex.pack(up(ed/dl/u)),
                               'error_over_s_in_u': ex.pack(up(ed/sl/u)),
                               'deleted_tail_over_s_in_u': ex.pack(up(tail/sl/u))})
        return result

    global_boxes, accepted_boxes = boxes(False), boxes(True)
    value = lambda record, key: F(record[key]['exact'])
    gd = max(value(b, 'error_over_D_in_u') for b in global_boxes)*u
    ad = max(value(b, 'error_over_D_in_u') for b in accepted_boxes)*u
    ass = max(value(b, 'error_over_s_in_u') for b in accepted_boxes)*u
    tail = max(value(b, 'deleted_tail_over_s_in_u') for b in accepted_boxes)
    assert gd < u and ad < u/8 and ass < u/32
    assert tail < F(11, 1000)
    return {'global_max_error_over_D': gd,
            'accepted_max_error_over_D': ad, 'accepted_max_error_over_s': ass,
            'accepted_deleted_tail_over_s_in_u': tail,
            'table_degree': 16, 'asymptotic_degree': 18,
            'compensated_stages': 3,
            'method': 'exact binary64 coefficient fractions, continuous Taylor/moment enclosure and source-order rounding; no input samples',
            'table_plain_errors': [ex.pack(v) for v in plain],
            'table_low_errors': [ex.pack(v) for v in low],
            'asymptotic_plain_error': ex.pack(ap), 'asymptotic_low_error': ex.pack(al),
            'global_boxes': global_boxes, 'accepted_boxes': accepted_boxes}


def finish_ledger(gate_bits, shape, mills):
    """Full new source rounding, stringent residual transfer, output rounding."""
    u, nu = U, NU
    gate = F(1, 2**gate_bits)
    mu, alpha_bound, beta_bound = shape['lambda_abs'], shape['alpha_abs'], shape['beta']
    fourth = shape['fourth_derivative_bound_over_factorial']
    # The high-component bound is used only for localization/Taylor remainder.
    # It is never substituted for the much smaller weighted residual bound.
    radius = shape['true_residual_radius']
    assert radius == gate + F(1, 2**38)
    zmag = (1 + u) * (gate + nu) / (1 - u) + nu
    assert zmag < radius
    true_s_error, true_y_error, true_d_error = 512*u*u, 64*u*u, mills['accepted_max_error_over_D']
    # D and s, their rounded values, and their quotient are positive normal.
    # The source sign multiplication is exact. Bound both ratio directions.
    num_upper = (1 + true_d_error) * (1 + u)
    num_lower = (1 - true_d_error) * (1 - u)
    den_upper = (1 + true_s_error) * (1 + u)
    den_lower = (1 - true_s_error) * (1 - u)
    lambda_relative = max(num_upper*(1+u)/den_lower-1,
                          1-num_lower*(1-u)/den_upper) + 128*nu
    assert lambda_relative < 4*u
    y_relative = (1 + true_y_error) * (1 + u) - 1 + 8*nu
    yy = Ball(F(68), 68*y_relative)
    x, one = Ball(F(40)), Ball(F(1))
    lam = Ball(mu, mu*lambda_relative)
    # Negation is exact; add() is also the absolute-error ledger for a
    # subtraction. Each source multiply/add/divide/FMA is modeled in order.
    xy = x.mul(yy)
    alpha = xy.mul(lam).sharpen(alpha_bound)
    q = x.mul(x).add(yy.mul(yy)).add(xy)
    beta = q.mul(lam.mul(lam)).sharpen(beta_bound)
    c2_mag = (1+alpha_bound)/2
    c3_mag = (1+3*alpha_bound+2*alpha_bound**2+beta_bound)/6
    c2 = one.add(alpha).scale(F(1, 2)).sharpen(c2_mag)
    two_alpha_minus_three = alpha.scale(F(2)).add(Ball(F(3)))
    c3 = alpha.fma(two_alpha_minus_three, one).add(beta).divide_constant(F(6)).sharpen(c3_mag)
    # z here denotes the actual binary64 -val(residual), so no error is
    # attached to this Ball. Its perturbation from the true residual is
    # transferred separately below, with the true |lambda| weight retained.
    z = Ball(zmag)
    delta = lam.mul(z.mul(z.fma(z.fma(c3, c2), one)))
    delta_mag = mu*zmag*(1+c2_mag*zmag+c3_mag*zmag*zmag)
    # Ball operations outward-round their magnitude bounds. Replace only
    # that magnitude with this exact expanded polynomial bound; keep every
    # accumulated source-operation error, including outward-rounding reserve.
    delta = delta.sharpen(delta_mag)
    computed_delta_mag = delta_mag + delta.err
    # Exact-cubic derivative/|lambda| <= 1+2r+(11/2)r^2 here.
    lipschitz = 1+2*c2_mag*radius+3*c3_mag*radius*radius
    weighted_pair_error = mu*F(376, 1000)*u + mills['accepted_max_error_over_s']/(1-true_d_error) + 200000*u*u + 256*nu
    # val(residual) is one RN sum, including when its result is subnormal.
    weighted_val_error = mu * (u*(gate+nu)/(1-u) + nu)
    residual_to_delta = lipschitz * (weighted_pair_error + weighted_val_error)
    # True-root Taylor remainder. |d log s/dq|=|lambda|<1 gives both the
    # intermediate s/s0 and the final s0/s* factors <=1/(1-radius).
    truncation = fourth*radius**4/(1-radius)**2
    # Let S=sh+sl be the normalized source DD geometry. |S/s0-1|<=es,
    # |sl|/s0 <=u(1+es)+8nu and |sh|/s0 <=(1+u)(1+es)+8nu.
    # final low=FMA(sh,delta_hat,sl) omits sl*delta_hat and rounds once.
    es = true_s_error
    low_mag = u*(1+es)+8*nu
    high_mag = (1+u)*(1+es)+8*nu
    omitted_low_product = low_mag * computed_delta_mag
    low_fma_round = u*(high_mag*computed_delta_mag+low_mag)+8*nu
    geometry = es * (1+computed_delta_mag)
    pre = residual_to_delta + delta.err + geometry + omitted_low_product + low_fma_round
    # Last addition is rounded once. s* >1/8 so its absolute term <=8nu.
    raw = (1 + truncation + pre/(1-radius))*(1+u)-1+8*nu
    assert computed_delta_mag < F(1, 1000)
    assert raw < F(1, 10**14)
    return {
        'gate_bits': gate_bits, 'true_residual_radius': radius,
        'lambda_relative_error': lambda_relative,
        'alpha_source_error': alpha.err, 'beta_source_error': beta.err,
        'c2_source_error': c2.err, 'c3_source_error': c3.err,
        'delta_source_error': delta.err, 'delta_source_magnitude': computed_delta_mag,
        'weighted_pair_residual_error': weighted_pair_error,
        'weighted_val_residual_error': weighted_val_error,
        'residual_polynomial_lipschitz': lipschitz,
        'residual_to_delta': residual_to_delta,
        'truncation_relative_to_root': truncation,
        'geometry_output_error': geometry,
        'omitted_low_product': omitted_low_product,
        'output_low_fma_rounding': low_fma_round,
        'raw_relative_error': raw, 'raw_over_epsilon': raw/(2*u),
    }


def target_transfer(lower_raw, upper_raw, cap, ex):
    """Recompute original-input rho with new internal budgets, not old rho."""
    u, nu = U, NU
    el = F(376, 1000)*u
    lower = ((1+lower_raw)/(1-el)-1)/(2*u)
    cap_log = el+4096*u*u
    cap_internal = (1+upper_raw)/(1-F(2, 5)*cap_log)-1
    ec = qfield(cap['cap_absolute_error'])
    ea = 16*u*u+nu
    eg = ec+ea
    theta = (5*ec+10*ea)/(1-5*ec)
    gap_min = (1-u)*F(1, 2**80)-eg
    assert gap_min > 0
    rgap = (eg/gap_min+5*ec)/(1-5*ec)
    factor = 1/(1-rgap)
    assert factor < F(1001, 1000) and cap_internal < 1
    cap_final = (cap_internal+factor*theta)/(2*u)
    grid = F(1, 2**248)
    integer_ratio = ((1+2691*grid*2**118)*(1+2*u*u)
                     /((1-5*2691*grid)*(1-2*u*u))-1)
    assert integer_ratio < 8*u*u
    h = F(3, 5)+u
    low = u*h+nu
    dmin = F(1, 2**20)/(1+h)
    rl = (1+h)*u+nu/dmin
    gamma = (2*u*(1+low)+u*(1+2*u)*rl+nu/dmin+low*rl/(2*(1-low)))/(1-rl/(1-low))
    assert gamma < 2*u+32*u*u
    e = F(38, 100)*u
    b = (1+gamma)/(1-e)-1
    a = e/(1-e)*(1+b)
    rgap = a/F(1, 2**21)+b
    factor_exp = 1/(1-rgap)
    assert rgap < F(1, 2**32)
    assert ex.expq(F(601, 1000))[1] < F(11, 6)
    k = (1+upper_raw)*factor_exp

    def at_chi(chi):
        lam = min(F(2, 5), F(5, 6)*chi)
        return (upper_raw+k*(a*chi+b*lam))/(2*u*(1+chi))

    expm1 = max(at_chi(F(0)), at_chi(F(12, 25)), k*a/(2*u))
    return {'lower': lower, 'ordinary_cap': cap_final,
            'fixed_cap': cap_final, 'expm1': expm1,
            'maximum': max(lower, cap_final, expm1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--fast-certificate', type=Path, required=True)
    parser.add_argument('--simd-root', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    source, reference = args.source_root.resolve(), args.reference_root.resolve()
    for root, bindings in ((source, SOURCE_FILES), (reference, REFERENCE_FILES)):
        for name, wanted in bindings.items():
            assert digest((root/name).read_bytes()) == wanted, ('reviewed file binding', root, name)
    original = (reference/'src/large.rs').read_text()
    current = (source/'src/experimental/large.rs').read_text()
    unchanged = ('dd', 'two', 'neg', 'add', 'sub', 'mul', 'scale', 'divide', 'val',
                 'frexp', 'scalbn', 'logdd', 'raw_stage', 'reduced_cap',
                 'mills_precise_pair', 'geometry', 'residual', 'reduced_grid',
                 'grid_dd', 'fixed_scaled_cap', 'resolve_cap', 'fma')
    for name in unchanged:
        assert canonical(function_body(current, name)) == canonical(function_body(original, name)), ('inherited source order', name)
    # Every adaptive branch stops the same physical recurrence at d=3/4/5.
    # Accuracy is not a premise: the preserved bracket/bisection fallback and
    # guarded precise handoff certify completion regardless of the proposal.
    predictor = adaptive_physical_source(original)
    assert canonical(function_body(current, 'physical_step_adaptive')) == canonical(predictor), 'unmodeled physical predictor'
    physical_predictor = symbolic_physical_predictor()
    # Changing the proposal cannot weaken the bracket or precise handoff. Its
    # exact source/formula is still bound to make the reviewed algorithm explicit.
    assert canonical(function_body(current, 'log_step3')) == canonical(LOG_STEP3_SOURCE), 'unmodeled log iteration'
    log_iteration = symbolic_log_step3()
    # Retain the full integer implementation, not merely its caller hashes.
    def integer_source(text):
        return canonical(text[text.index('struct U256'):text.index('fn reduced_grid')])
    assert integer_source(current) == integer_source(original)
    assert canonical(function_body(current, 'finish_direct')) == canonical(FINISH_SOURCE), 'unmodeled direct finalizer source'
    precise = canonical(function_body(current, 'precise'))
    precise = replace_once(precise,
        'letlambda=(ifupper{-1.0}else{1.0})*val(d)/val(s);',
        'letderivative=(ifupper{-1.0}else{1.0})*val(s)/val(y)/val(d);')
    precise = replace_once(precise, 's,lambda,', 's,derivative,')
    assert precise == canonical(function_body(original, 'precise')), 'unmodeled precise evaluator change'

    pair = canonical(function_body(current, 'precise_pair_impl'))
    pair = replace_once(pair, 'ifA{ifi<=16{ca[i]}else{0.0}}', 'ifA{ca[i]}')
    pair = replace_once(pair, 'ifB{ifi<=16{cb[i]}else{0.0}}', 'ifB{cb[i]}')
    pair = replace_once(pair, 'letn=ifA&&B{16}else{18};', 'letn=ifA||B{24}else{18};')
    pair = replace_once(pair, '(z,q,iv,(k+12)asusize)', '(z,q,iv,&M_COEFF[(k+12)asusize])')
    pair = replace_once(pair,
        'let(za,qa,ia,ka)=prepare(a,A);let(zb,qb,ib,kb)=prepare(b,B);letca=&M_COEFF_HI[ka];letcb=&M_COEFF_HI[kb];letla=&M_COEFF_LO[ka];letlb=&M_COEFF_LO[kb];',
        'let(za,qa,ia,ca)=prepare(a,A);let(zb,qb,ib,cb)=prepare(b,B);')
    pair = pair.replace('ca[i]', 'ca[i].hi').replace('cb[i]', 'cb[i].hi')
    pair = replace_once(pair,
        'letcoef_lo=|i:usize|F64x2::new(ifA{la[i]}else{0.0},ifB{lb[i]}else{0.0});',
        'letcoef_lo=|i:usize|{F64x2::new(ifA{ca[i].lo}else{0.0},ifB{cb[i].lo}else{0.0},)};')
    assert pair == canonical(function_body(original, 'precise_pair_impl')), 'unmodeled compensated Mills arithmetic'
    # Bit identities justify the representation-only reverse bridge above.
    def bits(text, name):
        section = re.split(rf'(?:const|static) {name}:', text, maxsplit=1)[1].split('=', 1)[1]
        section = re.split(r'[}\]];', section, maxsplit=1)[0]
        return re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)', section)
    shared = (source/'src/experimental/mills_constants.rs').read_text()
    old_constants = (reference/'src/large_constants.rs').read_text()
    local_constants = (source/'src/experimental/large_constants.rs').read_text()
    for name in ('ASYM', 'CAP_TABLE', 'EXP_C', 'LN2', 'LOGK'):
        assert bits(local_constants, name) == bits(old_constants, name), ('inherited constant bits', name)
    old_bits = bits(old_constants, 'M_COEFF')
    assert bits(shared, 'M_COEFF_HI') == old_bits[::2]
    assert bits(shared, 'M_COEFF_LO') == [old_bits[2*(25*k+j)+1] for k in range(44) for j in range(3)]

    solve = canonical(function_body(current, 'solve'))
    matched_gate = re.search(r'pe.residual.hi.abs\(\)<=f64::from_bits\((0x[0-9a-f]+)\)', solve)
    assert matched_gate
    gate_bits = int(matched_gate.group(1), 16)
    gate_value = F(struct.unpack('<d', struct.pack('<Q', gate_bits))[0])
    choices = {F(1, 2**k): k for k in (22, 17, 16, 14)}
    assert gate_value in choices, 'unreviewed handoff gate'
    selected_gate = choices[gate_value]
    old_solve = canonical(function_body(original, 'solve'))
    prefix = solve.removesuffix('finish_direct(x,pe)}')
    assert prefix != solve, 'unreviewed finish call'
    prefix = replace_once(prefix, 'log_positive_normal_price(b)', 'logdd(dd(b))')
    prefix = replace_once(prefix, 'physical_step_adaptive(', 'physical_step7(')
    prefix = replace_once(prefix, 'log_step3(', 'log_step5(')
    prefix = replace_once(prefix, matched_gate.group(1), '0x3e90000000000000')
    assert prefix == old_solve.split('letyy=val(pe.y);', 1)[0], 'unreviewed target/routing/iteration source'
    # accepted_region_shape uses precisely this lower routing predicate and
    # normalized target construction. The whole-prefix comparison additionally
    # binds every branch between these operations and the finalizer call.
    assert 'letell=add(log_positive_normal_price(b),dd(a/2.0));letupper=ell.hi>-0.6;letmuttarget=ell;' in solve
    assert canonical(function_body(current, 'add')).endswith('two(t.hi,e)}')
    assert 'e.f.abs()<=f64::from_bits(0x3dd0000000000000)' in solve, 'fast prefilter changed'
    # Normal-positive b makes frexp's exceptional branches unreachable. Its
    # exact mantissa/exponent extraction and normalization equal this helper.
    assert canonical(function_body(current, 'log_positive_normal_price')) == canonical(NORMAL_LOG_SOURCE), 'normal-price logarithm helper changed'

    fast = json.loads(args.fast_certificate.read_text())
    fast_payload = json.dumps({k: fast[k] for k in FAST_CERTIFICATE_FIELDS}, sort_keys=True, separators=(',', ':')).encode()
    assert digest(fast_payload) == FAST_CERTIFICATE_PAYLOAD_HASH, 'unreviewed fast numerical certificate'
    assert fast['script_sha256'] in FAST_VERIFIER_HASHES, 'unreviewed fast checker certificate'
    for name, wanted in FAST_FUNCTIONS.items():
        assert digest(function_body(current, name).encode()) == wanted, ('fast source binding', name)
        assert fast['source_functions'][name] == wanted
    for name, wanted in SOURCE_FILES.items():
        assert fast['source_file_bindings'][name] == wanted
    fast_error = qfield(fast['maximum_residual_error'])
    assert fast_error < F(1, 2**38)

    manifest = tomllib.loads((source/'Cargo.toml').read_text())
    assert manifest['dependencies']['fearless_simd']['version'] == '=1.0.0'
    if args.simd_root:
        simd_root = args.simd_root.resolve()
    else:
        registry = Path(os.environ.get('CARGO_HOME', Path.home()/'.cargo'))/'registry/src'
        candidates = sorted(registry.glob('*/fearless_simd-1.0.0'))
        assert len(candidates) == 1, 'Pass --simd-root for the reviewed backend'
        simd_root = candidates[0]
    simd_hash = hashlib.sha256()
    for path in sorted((simd_root/'src').rglob('*.rs')):
        simd_hash.update(str(path.relative_to(simd_root)).encode()+b'\0')
        simd_hash.update(path.read_bytes()+b'\0')
    assert simd_hash.hexdigest() == SIMD_HASH

    sys.path.insert(0, str(reference/'scripts'))
    spec = importlib.util.spec_from_file_location('bound_large_source', reference/'scripts/prove_global_large_source.py')
    inherited = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inherited)
    assert inherited.ROOT.resolve() == reference
    geometry, residual, reciprocal = inherited.geometry(), inherited.residual(), inherited.reciprocal()
    ex = inherited.ex
    assert qfield(reciprocal['inverse_relative_error_in_u2']) < 4
    assert qfield(reciprocal['square_relative_error_in_u2']) < 32
    assert qfield(geometry['y_relative_error_in_u2']) < 64
    assert qfield(geometry['negative_s_relative_error_in_u2']) < 64
    assert qfield(geometry['positive_s_relative_error_in_u2']) < 64
    assert qfield(residual['residual_arithmetic_error_in_u2']) < 200000
    assert qfield(residual['precise_high_error']) < F(1, 2**38)
    target_spec = importlib.util.spec_from_file_location('bound_large_targets', reference/'scripts/prove_large_targets.py')
    inherited_targets = importlib.util.module_from_spec(target_spec)
    target_spec.loader.exec_module(inherited_targets)
    # This source-order target-log proof does not call its old finalizer.
    # Its result supplies the .376u target budget used by target_transfer.
    target_logarithm = inherited_targets.logdd()
    assert target_logarithm['logdd_error_u'] < F(376, 1000)
    assert target_logarithm['ell_error_u'] < F(376, 1000)
    # Preserve the same target guards and root brackets after the new
    # finalizer. These exact inequalities are the inherited routing proof's
    # numerical premises, independent of the final correction formula.
    assert ex.expq(-F(7, 10))[0] > F(46, 100)
    assert 22*ex.ln2[1] < 100
    assert ex.expq(-F(3, 2))[0] > F(1, 5)
    assert 115*ex.ln2[1] < 80
    mills_spec = importlib.util.spec_from_file_location('bound_compensated_mills', reference/'scripts/prove_global_mills.py')
    inherited_mills = importlib.util.module_from_spec(mills_spec)
    mills_spec.loader.exec_module(inherited_mills)
    assert inherited_mills.ROOT.resolve() == reference
    mills = reduced_precise_budget(inherited_mills)
    # High-component localization uses the wider whole-iteration bracket.
    # The source arithmetic term is unchanged; the D-model defect is new.
    global_d = mills['global_max_error_over_D']
    precise_high = 1600*U + global_d/(1-global_d) + F(376, 1000)*U + 200000*U*U
    assert precise_high < F(1, 2**38)
    assert mills['accepted_max_error_over_D'] < U/8
    assert mills['accepted_max_error_over_s'] < U/32
    cap = json.loads((reference/'results/global_rho_proof/large_cap_bound.json').read_text())
    domain = json.loads((reference/'results/global_proof/large_domain_integer.json').read_text())
    assert domain['k_integer_range'] == [6, 1021]
    assert qfield(domain['large_integer_cap']['final_error_in_grid_units']) < 2691
    assert qfield(cap['gap_absolute_error']) < F(1, 2**85)

    # Absolute terms: all newly divided quantities have safe positive normal
    # divisors. On the bracket D>2^-13, 1/8<s<128, 4<y<68 and
    # |lambda|>2^-20. Cancellation/zero/subnormal polynomial terms retain NU.
    # The sole normalization of an absolute output error amplifies by <8.
    assert F(20)/((ex.sqrtq(F(1620))[1]+40)*1623) > F(1, 2**13)
    assert F(20, 81) > F(1, 8)
    assert 2**32*TAU*2**120 < NU
    # Every fourth update is still a rounded midpoint; seed accuracy does
    # not enter this stopping argument. The new gate accepts at least all
    # old precise handoffs, and the fast handoff fits all supported gates.
    width = F(41, 2**49)+80*U
    assert 128*width+fast_error < F(1, 2**34)
    assert fast_error < F(1, 2**34)
    assert F(1, 2**34)+fast_error+F(1, 2**38) < gate_value
    assert len(range(3, 196, 4)) == 49 and 196 < 208

    global_fourth = symbolic_and_shape(ex)
    all_bounds = {}
    for gate in (22, 17, 16, 14):
        # A fast-filter acceptance also has precise.high inside the selected
        # gate, so both entrances to finish_direct share the localized segment.
        assert F(1, 2**34)+fast_error+F(1, 2**38) < F(1, 2**gate)
        shape = accepted_region_shape(ex, gate)
        lower = finish_ledger(gate, shape, mills)
        upper = finish_ledger(gate, shape, mills)
        transferred = target_transfer(lower['raw_relative_error'], upper['raw_relative_error'], cap, ex)
        assert transferred['maximum'] < F(9, 10), ('rho budget failure', gate)
        all_bounds[str(gate)] = {
            # Preserve exact signed values and each field's bound direction:
            # this group contains lower price bounds as well as upper errors.
            'accepted_region': {k: {'exact': str(v), 'preview': float(v)}
                                for k, v in shape.items()},
            'lower': {k: pack(v) if isinstance(v, F) else v for k, v in lower.items()},
            'upper': {k: pack(v) if isinstance(v, F) else v for k, v in upper.items()},
            'original_input_rho': {k: pack(v) for k, v in transferred.items()},
        }
    out = {
        'status': 'All exact rational assertions passed',
        'scope': 'Conditional large-route direct cubic source certificate; exact binary64 normal-price inputs with 10<a<1417 and 0<b<exp(-a/2)',
        'selected_gate_bits': selected_gate,
        'selected_bound': all_bounds[str(selected_gate)],
        'all_proved_gates': all_bounds,
        'algebra': {'lambda_lower': 'D/s', 'lambda_upper': '-D/s',
                    'alpha': '-x*y*lambda', 'beta': '(x*x+y*y-x*y)*lambda*lambda',
                    'relative_update': 'lambda*(z+(1-alpha)*z^2/2+(1-3*alpha+2*alpha^2+beta)*z^3/6)',
                    'global_bracket_fourth_derivative_bound_over_factorial': pack(global_fourth),
                    'accepted_region_fourth_derivative_bound_over_factorial': pack(F(33, 20)),
                    'identities_checked_symbolically': ['s2/s=lambda*(1-alpha)', 's3/s=lambda*(1-3*alpha+2*alpha^2+beta)', 's4/s=lambda*(1-7*alpha+12*alpha^2+6*beta-6*alpha^3-7*alpha*beta-gamma)']},
        'recomputed_reduced_precise_mills': {k: pack(v) if isinstance(v, F) else v for k, v in mills.items()},
        'recomputed_precise_high_error': pack(precise_high),
        'physical_predictor': {**physical_predictor, 'termination': 'Exact preserved bracket/gate guards and forced bisections; no proposal accuracy premise'},
        'log_iteration': {**log_iteration, 'termination': 'Preserved finite/in-bracket guards, reliable residual signs and forced bisections; no proposal accuracy premise'},
        'recomputed_inherited_source_reciprocal': reciprocal,
        'recomputed_inherited_source_geometry': geometry,
        'recomputed_inherited_source_residual_arithmetic': {k: residual[k] for k in ('LN2_constant_error_in_u2', 'LOGK_constant_error_in_u2', 'residual_arithmetic_error_in_u2')},
        'recomputed_inherited_target_logarithm': {k: pack(v) for k, v in target_logarithm.items()},
        'source_sha256': {'src/experimental/large.rs': digest(current.encode()), **SOURCE_FILES},
        'reference_bindings': REFERENCE_FILES,
        'fast_certificate_sha256': digest(args.fast_certificate.read_bytes()),
        'fast_residual_bound': pack(fast_error),
        'simd_source_tree_sha256': SIMD_HASH,
        'script_sha256': digest(Path(__file__).read_bytes()),
        'primitive_contracts': ['RN-even binary64 basic/FMA/sqrt; gradual underflow; no unsafe reassociation',
                                'reduced ln error <=u*abs(ln); fast ln error <=2u*abs(ln)',
                                'expm1 relative error <=2u at reached negative arguments'],
        'inherited_analytic_dependencies': ['Mills integral/moment identities and hash-bound continuous compensated-Mills certificate',
            'Positive decreasing I1=-M\' on nonnegative arguments; lower p increases with x and a; elementary exp Taylor inequalities',
            'Existing source-identical target/cap/routing lemmas place the internal root inside the preserved brackets',
            'Source-identical integer cap arithmetic and external Lefevre-Muller Property 1 separation theorem',
            'Target-to-original-price finite transfer via Black vega log concavity and correlated gap/condition-number bounds'],
        'limits': 'No end-to-end Lean or vendor-libm proof; source hashes and finite tests are not themselves mathematical guarantees. Degree-16 compensated-Mills, finalizer source rounding and target budgets are recomputed, not inferred from old finalizer rho. Adaptive predictor has no universal one-step handoff guarantee; safeguarded stopping is preserved.',
    }
    path = args.out or source/'target/proof/large-direct-proof.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': out['status'], 'selected_gate_bits': selected_gate,
                      'rho_bounds': {k: v['original_input_rho']['maximum']['preview'] for k, v in all_bounds.items()},
                      'out': str(path)}, indent=2))


if __name__ == '__main__':
    main()
