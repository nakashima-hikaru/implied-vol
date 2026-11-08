#!/usr/bin/env python3
"""Verify that the optimized Experimental common-finish source is narrowly scoped.

The inverse source bridge restores the *entire* reviewed small.rs before
optimization. Only exact, enumerated snippet replacements are permitted.
This verifies source provenance, NOT by itself a floating-point rounding
theorem or re-issuance of the inherited conditional precision certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path

PINNED_MAIN = "c6bb0c1313ac4dc3a23e18ac4752ac21fc87e46f"
PINNED_SMALL_GIT_BLOB = "6a6c5e55719c64de2134d8803482a892c02eb638"


def replace_once(source, old, new, name):
    occurrences = source.count(old)
    if occurrences != 1:
        raise ValueError(f"{name}: expected exactly one source occurrence, got {occurrences}")
    return source.replace(old, new, 1)


def restore_source(source):
    source = replace_once(
        source,
        "#[cfg(test)]\nmod finish_optimization_tests;\n",
        "",
        "test module",
    )
    scale_optimized = """#[inline]
fn scale(a: Pair, b: f64) -> Pair {
    // Specialized multiplication by a scalar: the second low component
    // is exactly zero. The two live FMA operations retain their original
    // order, as does the final Pair normalization.
    let hi = a.hi * b;
    let lo = fma(a.hi, b, -hi);
    let lo = fma(a.lo, b, lo);
    add(Pair::from(hi), Pair::from(lo))
}

// The exact binary64 square is already normalized to a high/low pair
// whenever its rounded high component is a normal float. Avoid a generic
// twofold multiplication and normalization on this very common path.
#[inline]
fn square_scalar(x: f64) -> Pair {
    let high = x * x;
    if high.is_normal() {
        Pair::new(high, fma(x, x, -high))
    } else {
        mul(Pair::from(x), Pair::from(x))
    }
}

// Scaling a normalized twofold number by 1/2 is exact if its components
// remain normal. The tiny-value branch retains the original computation.
#[inline]
fn scale_negative_half(a: Pair) -> Pair {
    if a.hi >= f64::from_bits(0x20b0_0000_0000_0000) {
        Pair::new(-0.5 * a.hi, -0.5 * a.lo)
    } else {
        scale(a, -0.5)
    }
}

"""
    scale_original = """#[inline]
fn scale(a: Pair, b: f64) -> Pair {
    mul(a, Pair::from(b))
}
"""
    source = replace_once(source, scale_optimized, scale_original, "scalar scaling and square")
    source = replace_once(
        source,
        """    div(
        dual_mills(sub(h, t), add(h, t), false),
        Pair::from(2.0 * t.hi),
    )""",
        """    div(dual_mills(sub(h, t), add(h, t), false), scale(t, 2.0))""",
        "exact doubling of positive t",
    )
    source = replace_once(source,
                          "    let tt = square_scalar(t.hi);",
                          "    let tt = mul(t, t);",
                          "normal square specialization")
    source = replace_once(source,
                          "    let exponent = scale_negative_half(add(hh, tt));",
                          "    let exponent = scale(add(hh, tt), -0.5);",
                          "exact normalized half scaling")
    source = replace_once(source,
        """    // Since the other high factor is exactly one and ev is finite and
    // positive, the generic twofold product reduces to this TwoSum.
    // Its coefficient rounding matches the old FMA with a zero addend.
    let exp_pair = add(Pair::from(ev), Pair::from(ev * exponent.lo));
    let v = mul(
        exp_pair,""",
        """    let v = mul(
        mul(Pair::from(ev), Pair::new(1.0, exponent.lo)),""",
        "unit-high exponential product")
    return source


def git_blob_sha(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1] /
                        "src/experimental/small.rs")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    optimized = args.source.read_text()
    restored = restore_source(optimized).encode()
    got = git_blob_sha(restored)
    if got != PINNED_SMALL_GIT_BLOB:
        raise ValueError(f"original source bridge hash mismatch: {got}")
    if args.baseline is not None and args.baseline.read_bytes() != restored:
        raise ValueError("restored source differs from pinned baseline file")
    report = {
        "result": "exact original source reconstructed",
        "pinned_main_commit": PINNED_MAIN,
        "pinned_small_blob": PINNED_SMALL_GIT_BLOB,
        "restored_blob": got,
        "optimized_sha256": hashlib.sha256(optimized.encode()).hexdigest(),
        "restored_sha256": hashlib.sha256(restored).hexdigest(),
        "limitations": [
            "The bridge does not prove arithmetic equivalence.",
            "Full continuous-source numerical certificate requalification remains necessary.",
        ],
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
