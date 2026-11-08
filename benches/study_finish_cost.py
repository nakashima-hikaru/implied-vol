#!/usr/bin/env python3
"""Pinned-base same-host scalar benchmark for Experimental's common-finish change.

The old and candidate executables must be compiled with the same source
of benches/implied_black_vol_paths.rs and the same compiler flags.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import platform
import statistics
import subprocess
from pathlib import Path

CASES = (
    "rank_small2", "rank_small3", "rank_core4", "finite_seed", "wing_seed",
    "lower_middle", "upper_middle", "mixed_normalised", "mixed_full",
    "central_deferred", "central_shape1", "mixed_near_atm", "large_lower",
)


def order(block):
    return ("base", "candidate", "candidate", "base") if block % 2 == 0 else (
        "candidate", "base", "base", "candidate"
    )


def invoke(path, name, iterations, rounds):
    p = subprocess.run(
        [str(path), str(iterations), str(rounds), "--solver", "experimental",
         "--experimental-paths", "--case", name],
        capture_output=True, text=True, check=False,
    )
    if p.returncode:
        raise ValueError(f"{name}: failed executable {path}: {p.stderr or p.stdout}")
    csv_rows = list(csv.DictReader(io.StringIO(p.stdout)))
    if len(csv_rows) != rounds:
        raise ValueError(f"{name}: expected {rounds} rows, got {csv_rows!r}")
    results = []
    for j, row in enumerate(csv_rows):
        if set(row) != {"case", "round", "ns_per_call", "checksum"}:
            raise ValueError(f"bad fields: {row!r}")
        if row["case"] != name or int(row["round"]) != j:
            raise ValueError(f"unexpected output: {row!r}")
        ns = float(row["ns_per_call"])
        ck = float(row["checksum"])
        if not (math.isfinite(ns) and ns > 0 and math.isfinite(ck) and ck > 0):
            raise ValueError("invalid timing/checksum")
        results.append({"round": j, "ns": ns, "checksum": row["checksum"]})
    return results


def benchmark(exes, iterations, rounds, blocks, cases):
    raw, summary = [], []
    for name in cases:
        for block in range(blocks):
            for position, version in enumerate(order(block)):
                rows = invoke(exes[version], name, iterations, rounds)
                raw.extend({**r, "case": name, "block": block,
                            "position": position, "version": version} for r in rows)
        data = {v: [r for r in raw if r["case"] == name and r["version"] == v] for v in exes}
        for version, rows in data.items():
            if len(rows) != 2 * blocks * rounds:
                raise ValueError(f"{name}/{version}: missing timing samples")
            if len({r["checksum"] for r in rows}) != 1:
                raise ValueError(f"{name}/{version}: checksum drift")
        checksums_equal = data["base"][0]["checksum"] == data["candidate"][0]["checksum"]
        if not checksums_equal:
            raise ValueError(f"{name}: baseline and candidate output checksums DIFFER")
        base, cand = (statistics.median(r["ns"] for r in data[v]) for v in exes)
        blocks_pct = []
        for block in range(blocks):
            b, c = (statistics.median(r["ns"] for r in data[v] if r["block"] == block)
                    for v in exes)
            blocks_pct.append((c / b - 1) * 100)
        summary.append({"case": name, "baseline_ns": base, "candidate_ns": cand,
                        "change_percent": (cand / base - 1) * 100,
                        "blocks_percent": blocks_pct, "checksums_equal": checksums_equal})
    return raw, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True, type=Path)
    p.add_argument("--candidate", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--iterations", type=int, default=120000)
    p.add_argument("--rounds", type=int, default=2)
    p.add_argument("--blocks", type=int, default=4)
    p.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    args = p.parse_args()
    if min(args.iterations, args.rounds, args.blocks) <= 0:
        p.error("all repetitions must be positive")
    exes = {"base": args.base.resolve(), "candidate": args.candidate.resolve()}
    hashes = {v: hashlib.sha256(path.read_bytes()).hexdigest() for v, path in exes.items()}
    raw, summary = benchmark(exes, args.iterations, args.rounds, args.blocks, args.cases)
    if hashes != {v: hashlib.sha256(path.read_bytes()).hexdigest() for v, path in exes.items()}:
        raise ValueError("binary changed during measurement")
    report = {"method": "pinned base, same CPU/compiler/benchmark, rotating ABBA",
              "rustc": subprocess.check_output(["rustc", "-vV"], text=True),
              "platform": platform.platform(), "hashes": hashes,
              "iterations": args.iterations, "rounds": args.rounds,
              "blocks": args.blocks, "summary": summary, "samples": raw}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("| Case | Baseline ns | Candidate ns | Change | Four block changes |")
    print("|---|---:|---:|---:|---|")
    for r in summary:
        changes = ", ".join(f"{x:+.2f}%" for x in r["blocks_percent"])
        print(f"| {r['case']} | {r['baseline_ns']:.2f} | {r['candidate_ns']:.2f} | "
              f"{r['change_percent']:+.2f}% | {changes} |")


if __name__ == "__main__":
    main()
