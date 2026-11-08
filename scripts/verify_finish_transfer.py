#!/usr/bin/env python3
"""Check exact binary64 transfer obligations for optimized small::finish.

The mathematical argument is conditional on the inherited source-bound
s0>2^-350, 0<=h<=8, 0<t<5.5 and normalized DD Pair contracts.
This script verifies source bindings and exact rational domain inequalities,
and runs deterministic finite IEEE-754 checks with a negative control.
Finite checks are NOT a general SMT proof of every floating-point operation.
"""
import argparse
import ctypes
import ctypes.util
from fractions import Fraction as F
import json
import math
from pathlib import Path
import random
import struct

from verify_finish_scalar_bridge import (
    PINNED_SMALL_GIT_BLOB, git_blob_sha, restore_source,
)

ROOT = Path(__file__).resolve().parents[1]
SQ = F(1, 2**900)
SEED_LOWER = F(1, 2**350)
U = F(1, 2**53)
MIN_NORM = F(1, 2**1022)
MIN_SUB = F(1, 2**1074)


def bits(v):
    return struct.unpack(">Q", struct.pack(">d", v))[0]


if hasattr(math, "fma"):
    fma = math.fma
else:
    # Ubuntu 24.04 defaults to Python 3.12 without math.fma (added in
    # Python 3.13). Use the C99 correctly-rounded fused operation from libm.
    path = ctypes.util.find_library("m")
    if path is None:
        raise RuntimeError("need Python 3.13 math.fma or an IEEE-754 libm.fma")
    _libm = ctypes.CDLL(path)
    _libm.fma.argtypes = [ctypes.c_double] * 3
    _libm.fma.restype = ctypes.c_double
    fma = _libm.fma


def add(a,b):
    ah,al = a
    bh,bl = b
    sm = ah + bh
    z = sm - ah
    e = ((ah-(sm-z))+(bh-z)) + (al+bl)
    hi = sm+e
    return hi,(sm-hi)+e


def mul(a,b):
    ah,al = a
    bh,bl = b
    hi=ah*bh
    low=fma(ah,bh,-hi)
    low=fma(ah,bl,low)
    low=fma(al,bh,low)
    return add((hi,0.0),(low+al*bl,0.0))


def scale(p,b):
    ah,al=p
    hi=ah*b
    low=fma(ah,b,-hi)
    low=fma(al,b,low)
    return add((hi,0.0),(low,0.0))


def old_scale(p,b):
    return mul(p,(b,0.0))


def old_square(x):
    return mul((x,0.0),(x,0.0))


def new_square(x):
    hi=x*x
    if hi >= float(SQ) and math.isfinite(hi):
        return hi,fma(x,x,-hi)+0.0
    return old_square(x)


def old_half(p):
    return old_scale(p,-0.5)


def new_half(p):
    if p[0] >= 2**-500:
        return -0.5*p[0], -0.5*p[1]+0.0
    return scale(p,-0.5)


def old_exp_pair(ev,elo):
    return mul((ev,0.0),(1.0,elo))


def new_exp_pair(ev,elo):
    return add((ev,0.0),(ev*elo,0.0))


def same(a,b):
    return bits(a[0])==bits(b[0]) and bits(a[1])==bits(b[1])


def check_exact_envelopes():
    # Original finite seed cover lower bound; accepted t=s0/2.
    assert (SEED_LOWER/2)**2 > SQ
    # A binary64 x=M*2^e has at most 53 significant bits. If
    # RN(x*x)>=2^-900 and 0<x<5.5, x>=2^-450 and e>=-502;
    # thus the exact product and its RN residual live on a
    # 2^(2e) >= 2^-1004 grid, which cannot underflow. The
    # residual is <=0.5 ULP of the product and has <=53 bits.
    assert F(1,2**1004) > MIN_NORM
    assert 2*(-450-52) == -1004
    assert F(8)**2+F(11,2)**2 == F(377,4)
    exp_abs_max = F(377,8)
    # Elementary strict bound ln(2)>0.69 from positive terms
    # ln(2)=sum_{n>=1} 2^-n/n, no numerical libm is used.
    ln2_lower = sum((F(1,n*2**n) for n in range(1,9)),F())
    assert ln2_lower > F(69,100)
    assert exp_abs_max < 70*F(69,100)
    # Therefore ev=exp(-0.5*(h^2+t^2))>2^-70, normal.
    assert F(1,2**70)>MIN_NORM
    # Signed-zero normalization is explicitly restored by '+0.0'
    # in both specialized low components.
    return {"small_seed_lower":"2^-350",
            "min_accepted_t2":"strictly above 2^-702",
            "specialized_square_threshold":"2^-900",
            "min_exact_product_grid":"2^-1004",
            "dd_half_guard":"2^-500",
            "exp_lower_theorem":"exp(-377/8)>2^-70",
            "normal_min":"2^-1022",
            "inherited_conditions":["normalized Pair","round-to-nearest-even",
                "gradual underflow","FMA with one final rounding",
                "0<=h<=8","0<t<5.5","small/core/wing seed s0>2^-350"]}


def run_checks(samples):
    rng=random.Random(0x786165cc)
    checks={"scalar":0,"square":0,"half":0,"exponential":0}
    # Adversarial witness: a normal rounded square may still have
    # an underflowing residual. The old overly-broad shortcut fails.
    witness=struct.unpack(">d",int("200ab8bf40df0b0a",16).to_bytes(8,"big"))[0]
    old=old_square(witness)
    naive=(witness*witness,fma(witness,witness,-witness*witness)+0.0)
    assert not same(old,naive), "negative control unexpectedly passed"
    assert same(old,new_square(witness)), "guard failed on underflow witness"
    half_witness=(struct.unpack(">d",int("2505fec697039bba",16).to_bytes(8,"big"))[0],
                  math.ldexp(1.0,-1074))
    assert not same(old_half(half_witness),
                    (-0.5*half_witness[0],-0.5*half_witness[1]))
    assert same(old_half(half_witness),new_half(half_witness))
    for i in range(samples):
        exponent=rng.randrange(-550,6)
        hi=math.ldexp(1.0+rng.random(),exponent)
        # 1/16 of tests intentionally use subnormal low components;
        # otherwise make a normalized pair from the source's add().
        if i%16 ==0:
            tiny=math.ldexp(1.0,-1074)
            low=tiny if i%2 else -tiny
        else:
            low=(rng.random()-0.5)*hi*math.ldexp(1.0,-52)
        pair=add((hi,0.0),(low,0.0))
        scalar=(-0.5,2.0,0.125,0.7)[i%4]
        if not same(old_scale(pair,scalar),scale(pair,scalar)):
            raise ValueError(f"scalar counterexample at {i}: {pair!r}, {scalar}")
        checks["scalar"]+=1
        if not same(old_square(hi),new_square(hi)):
            raise ValueError(f"square counterexample at {i}: {hi.hex()}")
        checks["square"]+=1
        if not same(old_half(pair),new_half(pair)):
            raise ValueError(f"half counterexample at {i}: {pair!r}")
        checks["half"]+=1
        ev=math.exp(-47.125*rng.random())
        elo=(rng.random()-0.5)*1e-14
        if not same(old_exp_pair(ev,elo),new_exp_pair(ev,elo)):
            raise ValueError(f"exponential counterexample at {i}: {(ev,elo)!r}")
        checks["exponential"]+=1
    # Exact-grid neighbors in the guard binade and near underflow.
    for k in [-545,-530,-515,-512,-511,-505,-480,-455,-451,-450,-449,-448,-447,-1,0]:
        for j in range(512):
            x=math.ldexp(1.0+j*math.ldexp(1.0,-14),k)
            assert same(old_square(x),new_square(x)),(k,j)
            checks["square"]+=1
    return {"counts":checks,
            "adverse_old_shortcut_detected":True,
            "adverse_old_signed_zero_detected":True,
            "all_finite_bit_checks_passed":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=ROOT/"src/experimental/small.rs")
    parser.add_argument("--report",type=Path)
    parser.add_argument("--samples",type=int,default=80000)
    args=parser.parse_args()
    if not 0 < args.samples <= 1_000_000:
        parser.error("sample count must be in 1..1000000")
    source=args.source.read_text()
    got=git_blob_sha(restore_source(source).encode())
    assert got==PINNED_SMALL_GIT_BLOB, "source bridge no longer matches pinned main"
    envelopes=check_exact_envelopes()
    finite=run_checks(args.samples)
    report={
        "status":"conditional IEEE equivalence ledger and finite bit checks passed",
        "pinned_baseline_blob":got,
        "exact_envelopes":envelopes,
        "finite_checks":finite,
        "proof_scope":"theorem transfer under inherited normalized-DD and finite seed assumptions",
        "not_claimed":[
            "exhaustive binary64 enumeration",
            "SMT proof of all helper expressions",
            "freshly executed archived continuous proof (requires external reference archive)",
            "C++ compiler or vendor libm certification",
        ]
    }
    print(json.dumps(report,indent=2))
    if args.report is not None:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2)+"\n")


if __name__=="__main__":
    main()
