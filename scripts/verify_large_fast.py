#!/usr/bin/env python3
"""Verify shortened large-route fast Mills tables with exact rational bounds.

Recomputes continuous table/source-order rounding and 114 residual boxes.
The final compensated arithmetic is unchanged; only its coefficient layout changed. This is a conditional local source
certificate, not a whole-domain solver or platform-libm correctness theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
from math import factorial
import os
import re
import struct
import sys
import tomllib

sys.dont_write_bytecode = True
if not __debug__:
    raise RuntimeError("Assertions are required; do not use Python -O")

SOURCE_FUNCTIONS = {
    'mills_fast': 'a38617929a68dc6a799e5cb7c0c86bc43141849e9f69b673ba4b0cc6d914a119',
    'fast_pair_impl': 'c6242010d9d618829a704f102dc9ef878b0a0486aaa7f480c76b7734b8d3f7b5',
    'mills_fast_pair': '8e1deeb5de3e24aaee8e2ee2a53fbf27daa2c6bf07c534fb822fa1bb18c2d28f',
    'precise_pair_impl': '51e16a623f37ea4a84e34795097a6b27030ea999ad6000b1a32f9768882f07b6',
    'mills_precise_pair': '90981c1c1e96ce6b95ddf85878b5b532cbbc564edf22ec30019bb9fd7c2ab960',
    'geometry': 'e525b7b6cc2ac1af38c9a6be0979d6f1d132388b3225f6bc40315028b1a3694a',
    'residual': 'f56a99511350b50d415cd759e4d0dc229a037fedd3c8068d3556abd6943b6af0',
    'precise': 'a3835258201fe947c32acaadaffd7e021b98c2630da7ce1e04085210336029b8',
    'eval_fast': '8873a86ffc761a8a1372b33da7c0555301f848597dd53bc02fd7935d6eed210a',
    'solve': '2879acc1c5896c69c1ccc32f64a40787cfac97c8e997ca4b93d400eca8b9e29a',
}
SOURCE_FILES = {
    'src/experimental/mills_constants.rs': '0fda9d58f198c692908879ead2020d62bc0b2859fb2abbab1913fa8f06ed3063',
    'src/experimental/large.rs': 'fc8e8cfb49481800aff398fc7268ea7714c578413a644c4fd551d5dd5b42ff29',
    'src/experimental/large_constants.rs': '9f4f8dc26f8c1313e60e7e2817a24438a609d73ab4ba130a7962902c35be3efd',
    'src/experimental/lanes.rs': '40502418e5f024c8a4a522b52e1a9b226e126906bf9aca14496fd56c336d6f32',
}
REFERENCE_FILES = {
    'src/large_constants.rs': 'a114f8123411950b1f92f3f6d3f447206dd7b3143d0a87361df03ffdf36ad02f',
    'scripts/prove_global_mills.py': '2dce1c923d5ba36d7109ffdc3a17f4c726e2996c2a50d6ee49663dfb436e7547',
    'scripts/prove_global_fast.py': 'b29064dd9d707b9a7843e70e69b8ec1dc9c6415e349560ab1b355736f95fdef3',
    'reference/iv_large_a_cap_table_v50/scripts/exact.py': 'cd66c85457348fabdc238e7395f6c4e10bc954f23a5af7ab709c43c08ac03559',
    'notes/global-rho-large-source.md': '3fd1fbae6b2d4ee027fe7a154a4984e8ba9a888bc68649b150f520f7ac08a026',
    'notes/global-rho-large-audit.md': 'd228455451fb3ef0c7f5b7b3b646fb5df58085dbda444e511fe3d95da5d3206e',
}
SIMD_HASH = '65450215ad42e632e66e97a0f2f8229e0934441b50b9fcca483c2cad80f5af75'


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


def decode(bits):
    return F(struct.unpack('<d', struct.pack('<Q', int(bits, 16)))[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--simd-root', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    source, reference = args.source_root.resolve(), args.reference_root.resolve()
    text = (source / 'src/experimental/large.rs').read_text()
    for name, want in SOURCE_FUNCTIONS.items():
        assert digest(function_body(text, name).encode()) == want, ('source function binding', name)
    for root, bindings in ((source, SOURCE_FILES), (reference, REFERENCE_FILES)):
        for name, want in bindings.items():
            assert digest((root / name).read_bytes()) == want, ('file binding', root, name)
    manifest = tomllib.loads((source / 'Cargo.toml').read_text())
    assert manifest['dependencies']['fearless_simd']['version'] == '=1.0.0'
    if args.simd_root:
        simd_root = args.simd_root.resolve()
    else:
        registry = Path(os.environ.get('CARGO_HOME', Path.home() / '.cargo')) / 'registry/src'
        candidates = sorted(registry.glob('*/fearless_simd-1.0.0'))
        assert len(candidates) == 1, 'Pass --simd-root for reviewed fearless_simd sources'
        simd_root = candidates[0]
    simd_digest = hashlib.sha256()
    simd_files = sorted((simd_root / 'src').rglob('*.rs'))
    assert simd_files
    for file in simd_files:
        simd_digest.update(str(file.relative_to(simd_root)).encode() + b'\0')
        simd_digest.update(file.read_bytes() + b'\0')
    assert simd_digest.hexdigest() == SIMD_HASH

    sys.path.insert(0, str(reference / 'scripts'))
    spec = importlib.util.spec_from_file_location('bound_fast', reference / 'scripts/prove_global_fast.py')
    old_fast = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old_fast)
    mp, ex = old_fast.mp, old_fast.ex
    assert mp.ROOT.resolve() == reference
    u, nu, up = mp.U, mp.NU, mp.up
    constants = (source / 'src/experimental/large_constants.rs').read_text()
    shared = (source / 'src/experimental/mills_constants.rs').read_text()
    reference_constants = (reference / 'src/large_constants.rs').read_text()

    def bits(text, name):
        section = re.split(rf'(?:const|static) {name}:', text, maxsplit=1)[1].split('=', 1)[1].split('];', 1)[0]
        return re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)', section)

    high_bits, low_bits = bits(shared, 'M_COEFF_HI'), bits(shared, 'M_COEFF_LO')
    reference_bits = bits(reference_constants, 'M_COEFF')
    asym_bits = bits(constants, 'ASYM')
    assert len(high_bits) == 44*25 and len(low_bits) == 44*3
    assert len(reference_bits) == 44*25*2 and len(asym_bits) == 19
    # Verify represented bits, including signed zero: every high and every
    # consumed low is the original value. Only lows j>=3 were discarded.
    assert high_bits == reference_bits[::2], 'high coefficient identity'
    assert low_bits == [reference_bits[(r*25+j)*2+1] for r in range(44) for j in range(3)], 'used low coefficient identity'
    assert asym_bits == bits(reference_constants, 'ASYM'), 'asymptotic coefficient identity'
    highs, lows, asym = list(map(decode, high_bits)), list(map(decode, low_bits)), list(map(decode, asym_bits))
    rows = [[(highs[r*25+j], lows[r*3+j] if j < 3 else F(0)) for j in range(25)] for r in range(44)]
    assert max(abs(v) for row in rows for pair in row for v in pair) < 2**80
    assert max(map(abs, asym)) < 2**80

    def table(k, stored):
        radius = up((1 + u) / 4 + nu)
        truth = mp.true_taylor(F(2 * k + 1, 4))
        # Bound the actual degree-24 high polynomial first, then remove its
        # represented degrees 17..24. This avoids a loose degree-16 Taylor
        # derivative bound on negative cells while retaining a rigorous model.
        coefficient = sum(max(abs(stored[j][0] - v) for v in truth[j]) * radius**j for j in range(25))
        omitted = sum(abs(stored[j][0]) * radius**j for j in range(17, 25))
        zmin = F(k, 2) - 2 * u
        if zmin < 0:
            moment, i1 = ex.expq(F(4))[1] * 2**25 * factorial(12), F(6)
        else:
            moment, i1 = F(2**12 * factorial(12)), F(1)
            if zmin > 0:
                moment, i1 = min(moment, F(factorial(25)) / zmin**26), min(i1, 1 / zmin**2)
        taylor_tail = moment * radius**25 / factorial(25)
        model = up(coefficient + omitted + taylor_tail)
        # Source starts exactly at stored[16]. Each of the 16 following FMAs
        # is rounded once. h bounds the rounded intermediate magnitude.
        h, rounding = abs(stored[16][0]), F(0)
        for j in range(15, -1, -1):
            p = up(h * radius + abs(stored[j][0]))
            rounding = up(radius * rounding + u * p + nu)
            h = up((1 + u) * p + nu)
        centering = i1 * (u / 4 + nu)
        return up(model + rounding + centering), {
            'k': k, 'omitted_high_tail': ex.pack(up(omitted)),
            'coefficient_error': ex.pack(up(coefficient)),
            'degree24_taylor_tail': ex.pack(up(taylor_tail)),
            'degree16_horner_rounding': ex.pack(rounding),
            'centering_error': ex.pack(up(centering)),
            'total_error': ex.pack(up(model + rounding + centering)),
        }

    computed = [table(k, rows[k + 12]) for k in range(-2, 32)]
    bounds, table_records = [r[0] for r in computed], [r[1] for r in computed]
    # The asymptotic polynomial, inverse-square and division are unchanged.
    # For a mixed pair, degrees 16..14 are exact zeros and degree 13 inserts
    # ASYM[13] exactly. The table lane's final division by 1 is exact.
    asym_error = old_fast.asymptotic(asym)
    boxes = []
    for upper, count in ((False, 82), (True, 32)):
        for i in range(count):
            left = F(i, 2) + (0 if upper else -40)
            right = left + F(1, 2)
            xa = min(abs(left), abs(right)) if left * right > 0 else F(0)
            xb = max(abs(left), abs(right))
            yl, yh = ex.sqrtq(xa * xa + 20)[0], ex.sqrtq(xb * xb + 2834)[1]
            zl, zr = (left, right) if upper else (-right, -left)
            ef = mp.maximum_error(bounds, zl, zr, asym_error)
            ey = mp.maximum_error(bounds, yl * (1 - 2 * u) - nu, yh * (1 + 2 * u) + nu, asym_error)
            geometry = up(2 * u / (yl * (1 - 2 * u)**2) + nu)
            sum_bound = up(ex.mills(zl)[1] + 1 / yl)
            ed = up(ef + ey + geometry + u * (sum_bound + ef + ey + geometry) + nu)
            if upper:
                dl = ex.mills(right)[0] + yh / (yh * yh + 1)
            elif right <= 0:
                v = -left
                dl = 20 / ((ex.sqrtq(v * v + 20)[1] + v) * (v * v + 23))
            else:
                dl = F(3, 4)
            relative = up(ed / dl)
            assert relative < F(1, 100) and dl > F(1, 2**13)
            # Source order: (-.5*x)*x, subtract LOGK.hi, ln(Dhat), add ln,
            # subtract RN(target.hi+target.lo). The 706u reserve covers the
            # latter sum's distance to exact tau (rather than treating tau
            # as an exact input). The inherited ln contract is 2u*abs(ln).
            arithmetic = 800*u + u + u*(800*(1+u)+1) + 32*u + 818*u + 1524*u + 706*u + 16*nu
            residual = up(relative / (1-relative) + arithmetic)
            boxes.append({'upper': upper, 'x': [str(left), str(right)],
                          'D_lower': ex.pack(dl), 'D_error': ex.pack(ed),
                          'residual_error': ex.pack(residual)})
    maximum = max(F(b['residual_error']['exact']) for b in boxes)
    epsilon = F(1, 2**38)
    assert maximum < epsilon
    assert (1 + u/2) * (1 + u) < 1 + 2*u
    assert (1 - u)**2 > 1 - 2*u

    # The unchanged solve() only updates a sign after |fast|>2^-34. Thus
    # epsilon<2^-34 preserves its correct sign. A fast gate also implies
    # |precise.hi|<=2^-34+epsilon+epsilon<2^-22, the precise finish gate.
    # Import the existing target/root bracket and |F'|<=128 contracts. With
    # midpoint error <=40u, 49 forced bisections give the following width.
    fast_gate, precise_gate = F(1, 2**34), F(1, 2**22)
    width = F(41, 2**49) + 80*u
    assert epsilon < fast_gate
    assert fast_gate + maximum + epsilon < precise_gate
    assert 128*width + maximum < fast_gate
    assert len(range(3, 196, 4)) == 49 and 196 < 208

    out = {
        'status': 'All exact rational assertions passed',
        'scope': 'Large fast Mills degree 16, unchanged degree-13 asymptotic and compensated final kernel; 34 continuous table cells and 114 continuous residual boxes',
        'maximum_residual_error': ex.pack(maximum),
        'maximum_residual_error_in_u': ex.pack(up(maximum/u)),
        'fits_2m38': maximum < epsilon,
        'coefficient_identity': {'high_bits': 44*25, 'used_low_bits': 44*3, 'asymptotic_bits': 19, 'all_equal_to_bound_reference': True},
        'table_errors': table_records, 'asymptotic_error': ex.pack(asym_error), 'boxes': boxes,
        'stopping': {'fast_sign_margin': ex.pack(fast_gate-maximum),
                     'fast_handoff_precise_bound': ex.pack(fast_gate+maximum+epsilon),
                     'precise_gate': ex.pack(precise_gate),
                     'width_after_49_bisections': ex.pack(width),
                     'fast_after_49_bisections': ex.pack(128*width+maximum)},
        'primitive_contracts': ['RN-even binary64 basic/FMA/sqrt, gradual underflow, source order, one-rounding lane-wise FMA',
                                'abs(ln_impl(Dhat)-ln(Dhat))<=2u*abs(ln(Dhat)) on reached normal Dhat'],
        'conditional_imports': ['Mills integral, moment recurrence, CF/Gaussian series and degree-24 Taylor remainder lemmas from the bound reference; their rational enclosures are recomputed',
                                'Unchanged target formation produces normalized DD with abs(target.hi)<=705 and a root inside the initial bracket',
                                'Unchanged precise high residual error <=2^-38, |F derivative|<=128 on bracket, and precise final-correction contract'],
        'source_functions': SOURCE_FUNCTIONS, 'source_file_bindings': SOURCE_FILES,
        'reference_bindings': REFERENCE_FILES, 'simd_source_tree_sha256': SIMD_HASH,
        'actual_source_sha256': {'src/experimental/large.rs': digest(text.encode()), **{p: digest((source/p).read_bytes()) for p in SOURCE_FILES}},
        'script_sha256': digest(Path(__file__).read_bytes()),
        'limits': 'Conditional continuous source-component proof, not a whole-domain solver certificate, end-to-end Lean proof, vendor-libm guarantee or speed claim.',
    }
    path = args.out or source / 'target/proof/large-fast-proof.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'status': out['status'], 'maximum_residual_error': float(maximum),
                      'in_u': float(maximum/u), 'fits_2m38': True,
                      'maximum_omitted_high_tail': max(r['omitted_high_tail']['upper_float'] for r in table_records),
                      'table_cells': len(table_records), 'residual_boxes': len(boxes), 'out': str(path)}, indent=2))


if __name__ == '__main__':
    main()
