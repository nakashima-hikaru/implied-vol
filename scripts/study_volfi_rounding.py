#!/usr/bin/env python3
"""Isolate volfi whole-book chart approximation from binary64 evaluation error.

Use the exact binary64 reference inputs and independent 160/240-digit Black roots
already embedded in tests/experimental_adaptive_seeds.rs. The coefficient literals
are parsed from the Rust port, rounded to binary64, then evaluated in high-precision
decimal arithmetic (without the per-operation binary64 rounding of the Rust kernel).
This identifies approximation/represented-coefficient error, NOT a formal bound,
and does not replace source-order rounding or reference provenance checks.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import struct
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PI = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068")


def binary64(bits: str) -> float:
    return struct.unpack(">d", int(bits.replace("_", ""), 16).to_bytes(8, "big"))[0]


def represented(value: str) -> Decimal:
    return Decimal.from_float(float(value))


def parse_floats(source: str, name: str) -> list[Decimal]:
    pattern = rf"pub\(super\) const {name}: \[f64; \d+\] = \[(.*?)\];"
    match = re.search(pattern, source, re.DOTALL)
    if match is None:
        raise ValueError(f"missing float coefficients: {name}")
    values = [represented(s.strip()) for s in match.group(1).split(",") if s.strip()]
    if not values:
        raise ValueError(f"empty coefficient table {name}")
    return values


def parse_usizes(source: str, name: str) -> list[int]:
    pattern = rf"pub\(super\) const {name}: \[usize; \d+\] = \[(.*?)\];"
    match = re.search(pattern, source, re.DOTALL)
    if match is None:
        raise ValueError(f"missing index coefficients: {name}")
    return [int(t.strip()) for t in match.group(1).split(",") if t.strip()]


def parse_fixture_bits(source: str) -> list[tuple[str, ...]]:
    rows = re.findall(r"\[\s*(0x[\da-fA-F_]+)\s*,\s*(0x[\da-fA-F_]+)\s*,\s*(0x[\da-fA-F_]+)\s*,\s*(0x[\da-fA-F_]+)\s*,\s*(0x[\da-fA-F_]+)\s*,?\s*\]", source)
    if len(rows) < 900:
        raise ValueError(f"expected >=900 independent references, got {len(rows)}")
    return rows


def clenshaw(coeff: list[Decimal], x: Decimal) -> Decimal:
    e0 = Decimal(0)
    e1 = Decimal(0)
    x2 = x + x
    for c in reversed(coeff[1:]):
        e0, e1 = x2 * e0 + c - e1, e0
    return x * e0 + coeff[0] - e1


def horner(coeff: list[Decimal], x: Decimal) -> Decimal:
    p = coeff[0]
    for c in coeff[1:]:
        p = p * x + c
    return p


def chart(h: Decimal, b: Decimal, tables: dict, pi_exact: bool,
          coordinate_override: Decimal | None = None) -> tuple[str, Decimal] | None:
    if not (h >= 2**(-500) and h <= 4 and b > 0):
        return None
    c = b * (h / 2).exp()
    if not (0 < c < 1):
        return None
    a = ((1 + (h.exp() - 1) / c).ln() if coordinate_override is None
         else coordinate_override)
    two_pi = 2 * PI if pi_exact else represented("6.283185307179586")
    if a < two_pi:
        if h > represented("0.35") * a:
            return None
        t0 = a * clenshaw(
            tables["NW_G0"],
            a / PI - 1 if pi_exact else a * represented("0.3183098861837907") - 1,
        )
        t = h / t0
        x, y = t*t, t0*t0
        accum = Decimal(0)
        for row in range(15, 0, -1):
            off = tables["NWRowA_OFF"][row]
            deg = tables["NWRowA_DEG"][row]
            accum = accum * x + horner(tables["NW_P"][off:off+deg+1], y)
        return "A", t * x * accum + t
    if a > represented("700") or h > 4 or h > represented("0.16")*a+represented("0.85"):
        return None
    u = 1 / a.sqrt()
    f1 = clenshaw(tables["NW_F1"], u*represented("5.537929048854298") + represented("-1.2093140434512706"))
    pp = u*f1
    t, w = h*pp, pp*pp
    z = h*h * (1/(4*PI*PI) if pi_exact else represented("2.5330295910584443e-2"))
    sp = (1+z).sqrt()+1
    q = z/(sp*sp)
    accum = Decimal(0)
    for row in range(16, 0, -1):
        off = tables["NWRowB_OFF"][row]
        deg = tables["NWRowB_DEG"][row]
        accum = accum*q + horner(tables["NW_D"][off:off+deg], w)*w
    return "B", t*q*accum + t


def read_rust_outputs(path: Path | None) -> dict[tuple[int, int], Decimal]:
    if path is None:
        return {}
    results = {}
    for h, b, vol in re.findall(r"VOLFI\\|([0-9a-f]{16})\\|([0-9a-f]{16})\\|([0-9a-f]{16})", path.read_text()):
        key = int(h, 16), int(b, 16)
        if key in results:
            raise ValueError(f"duplicate runtime output {key}")
        results[key] = Decimal.from_float(binary64("0x" + vol))
    if len(results) != 495:
        raise ValueError(f"expected 495 runtime output records, got {len(results)}")
    return results


def study(rust_outputs_path: Path | None = None) -> dict:
    table_text = (ROOT / "src/volfi_book/tables.rs").read_text()
    tables = {name: parse_floats(table_text, name) for name in ("NW_G0", "NW_F1", "NW_P", "NW_D")}
    tables.update({name: parse_usizes(table_text, name) for name in (
        "NWRowA_OFF", "NWRowA_DEG", "NWRowB_OFF", "NWRowB_DEG"
    )})
    ref_text = (ROOT / "tests/experimental_adaptive_seeds.rs").read_text()
    references = parse_fixture_bits(ref_text)
    runtime_results = read_rust_outputs(rust_outputs_path)
    results = []
    with localcontext() as ctx:
        ctx.prec = 90
        for bits in references:
            h, b, high, low, eta = (Decimal.from_float(binary64(v)) for v in bits)
            approximate = chart(h, b, tables, pi_exact=False)
            if approximate is None:
                continue
            region, volatility = approximate
            exact_constant = chart(h, b, tables, pi_exact=True)
            if exact_constant is None:
                raise ValueError("unexpected coverage mismatch with exact pi")
            h64, b64 = binary64(bits[0]), binary64(bits[1])
            c64 = b64 * math.exp(h64 / 2.0)
            a64 = math.log1p(math.expm1(h64) / c64)
            rounded_coordinate = chart(
                h, b, tables, pi_exact=False,
                coordinate_override=Decimal.from_float(a64)
            )
            if rounded_coordinate is None:
                raise ValueError("rounded chart coordinate changed region acceptance")
            runtime = runtime_results.get((int(bits[0].replace("_",""), 16),
                                           int(bits[1].replace("_",""), 16)))
            root = high + low
            rho = abs(volatility/root - 1)/eta
            rho_exact_pi = abs(exact_constant[1]/root - 1)/eta
            coord_rounding_rho = abs(volatility - rounded_coordinate[1])/(root*eta)
            evaluator_rounding_rho = (abs(runtime-rounded_coordinate[1])/(root*eta)
                                      if runtime is not None else None)
            results.append({
                "h_bits":bits[0], "b_bits":bits[1], "h":float(h),
                "b":float(b), "region":region,
                "rho_exact_polynomial":float(rho),
                "rho_exact_pi_only":float(rho_exact_pi),
                "constant_rounding_rho":float(abs(volatility-exact_constant[1])/(root*eta)),
                "coordinate_rounding_rho":float(coord_rounding_rho),
                "rounded_coord_ideal_rho":float(abs(rounded_coordinate[1]/root-1)/eta),
                "evaluator_rounding_rho":(float(evaluator_rounding_rho)
                                          if evaluator_rounding_rho is not None else None),
                "runtime_rho":(float(abs(runtime/root-1)/eta)
                               if runtime is not None else None),
                "signed_root_error":float((volatility-root)/(root*eta)),
                "approximate_vol":str(volatility),
                "reference":str(root)
            })
    if len(results) != 495:
        raise ValueError(f"expected 495 covered references, got {len(results)}")
    largest=sorted(results,key=lambda x:x["rho_exact_polynomial"],reverse=True)
    return {
        "status":"high-precision evaluation of represented (not exact analytic) chart coefficients",
        "reference_count":len(references), "direct_count":len(results),
        "region_counts":{k:sum(x["region"]==k for x in results) for k in ("A","B")},
        "rho_gt_one_count":sum(x["rho_exact_polynomial"]>1.0 for x in results),
        "max_rho_exact_polynomial":largest[0]["rho_exact_polynomial"],
        "max_rho_exact_pi_only":max(x["rho_exact_pi_only"] for x in results),
        "runtime_available":bool(runtime_results),
        "max_coordinate_rounding_rho":max(x["coordinate_rounding_rho"] for x in results),
        "max_rounded_coord_ideal_rho":max(x["rounded_coord_ideal_rho"] for x in results),
        "max_evaluator_rounding_rho":(max(x["evaluator_rounding_rho"] for x in results)
                                       if runtime_results else None),
        "max_runtime_rho":(max(x["runtime_rho"] for x in results)
                           if runtime_results else None),
        "runtime_rho_gt_one":(sum(x["runtime_rho"]>1 for x in results)
                              if runtime_results else None),
        "worst_runtime":(sorted(results,key=lambda x:x["runtime_rho"],reverse=True)[:10]
                         if runtime_results else []),
        "worst_10":largest[:10],
        "all_results":results,
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--rust-outputs",type=Path)
    args=parser.parse_args()
    result=study(args.rust_outputs)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="all_results"},indent=2))
