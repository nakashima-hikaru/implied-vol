#!/usr/bin/env python3
"""Reproducible independent Black roots for the volfi whole-book A/B domain.

This generates exact-input binary64 option prices from deliberately distributed
intrinsic coordinates, then *re-inverts those rounded prices* using mpmath
at 100 decimal digits. It does not call volfi or Experimental to produce roots.
Usage: python scripts/volfi_book_stress.py --output /tmp/inputs.tsv
Requires mpmath==1.3.0 in the test environment only.
"""
import argparse
import math
import random
import struct
from pathlib import Path

import mpmath as mp

mp.mp.dps = 100
SQRT2 = mp.sqrt(2)
SQRT2PI = mp.sqrt(2 * mp.pi)


def norm_cdf(z):
    return mp.erfc(-z / SQRT2) / 2


def black(h, s):
    d1 = -h / s + s / 2
    d2 = -h / s - s / 2
    return mp.exp(-h / 2) * norm_cdf(d1) - mp.exp(h / 2) * norm_cdf(d2)


def vega(h, s):
    return mp.exp(-h*h / (2*s*s) - s*s / 8) / SQRT2PI


def invert(h, b):
    lo = mp.mpf(0)
    hi = mp.mpf(2)
    while black(h, hi) <= b:
        hi *= 2
        if hi > 128:
            raise ValueError("unrepresentable upper bracket")
    # The root can be tiny for normal binary64 far-wing prices.
    # Bisection to a *relative* bracket is essential: unguarded Newton
    # starting from a near-cap volatility makes microscopic progress while
    # b is many orders of magnitude below the current Black price.
    for _ in range(128):
        if lo > 0 and hi / lo <= mp.mpf("1.25"):
            break
        midpoint = (lo + hi) / 2
        if black(h, midpoint) > b:
            hi = midpoint
        else:
            lo = midpoint
    else:
        raise ValueError("failed to establish a relative root bracket")
    s = (lo + hi) / 2
    for _ in range(180):
        p = black(h, s)
        err = p - b
        if abs(err) <= mp.mpf("1e-73") * (s * vega(h, s) + b):
            return s
        if err > 0:
            hi = s
        else:
            lo = s
        slope = vega(h, s)
        trial = s - err / slope if slope > 0 else mp.nan
        s = trial if (mp.isfinite(trial) and lo < trial < hi) else (lo + hi)/2
    raise ValueError(f"Black inversion failed at h={h} b={b}")


def bits(value):
    return f"{struct.unpack('>Q', struct.pack('>d', value))[0]:016x}"


def cases():
    # Regime A: 22 by 12 interior points, with both small and near-boundary h/a.
    a_values = [0.08 + (2*math.pi - 0.11) * (i + 0.5) / 22 for i in range(22)]
    h_fractions = (0.015, 0.06, 0.14, 0.27, 0.38, 0.50, 0.61,
                   0.70, 0.79, 0.88, 0.95, 0.995)
    for a in a_values:
        for frac in h_fractions:
            yield a * 0.35 * frac, a, "A"

    # Regime B: 26 by 12 points distributed logarithmically to a=680.
    for i in range(26):
        a = math.exp(math.log(2*math.pi + 0.035) +
                     (math.log(680) - math.log(2*math.pi + 0.035)) * (i + 0.5)/26)
        hcap = min(4.0, 0.16*a + 0.85)
        for frac in h_fractions:
            yield hcap*frac, a, "B"

    # Independent interior samples, not aliases of the reference seam dataset.
    rng = random.Random(20261108)
    for i in range(256):
        if i % 2:
            a = rng.uniform(0.1, 2*math.pi - 0.01)
            h = 0.35*a*rng.uniform(0.025, 0.995)
            yield h, a, "A_random"
        else:
            a = math.exp(rng.uniform(math.log(2*math.pi+0.04), math.log(680)))
            h = min(4.0, 0.16*a+0.85)*rng.uniform(0.025, 0.995)
            yield h, a, "B_random"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    skipped = 0
    with args.output.open("w") as out:
        out.write("# h_bits b_bits root_hi_bits root_lo_bits eta_bits input_region\n")
        for h_double, a_double, label in cases():
            h, a = mp.mpf(h_double), mp.mpf(a_double)
            exact_b = mp.expm1(h) / mp.expm1(a) * mp.exp(-h/2)
            b_double = float(exact_b)
            if not (math.isfinite(b_double) and b_double >= float.fromhex("0x1p-1022")):
                skipped += 1
                continue
            b = mp.mpf(b_double)
            root = invert(h, b)
            root_hi = float(root)
            root_lo = float(root - mp.mpf(root_hi))
            eta = float(mp.power(2, -52) * (1 + b / (root * vega(h, root))))
            if not all(map(math.isfinite, (root_hi, root_lo, eta))) or eta <= 0:
                raise ValueError("invalid reference root or conditioning")
            out.write(" ".join([bits(h_double), bits(b_double), bits(root_hi),
                                bits(root_lo), bits(eta), label]) + "\n")
            written += 1
    print(f"MPMATH_REFERENCE|generated={written}|skipped={skipped}|dps={mp.mp.dps}")
    if written < 750:
        raise ValueError("too few independent stress roots")


if __name__ == "__main__":
    main()
