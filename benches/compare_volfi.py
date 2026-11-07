#!/usr/bin/env python3
"""Same-executable Experimental vs Rust volfi-book ABBA/BAAB latency comparison."""
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

SOLVERS = ("experimental", "volfi-book", "volfi-sinh", "volfi-rsqrt", "volfi-compensated", "volfi-selective", "volfi-endpoint")
CASES = ("book_a", "book_b", "book_seam", "fallback", "mixed_book", "mixed_market", "mixed_full")
FIELDS = ("case", "round", "ns_per_call", "checksum", "covered", "total")


def rotation(block, reference="experimental", candidate="volfi-book"):
    return (reference, candidate, candidate, reference) if block % 2 == 0 else (
        candidate, reference, reference, candidate
    )


def parse_output(stdout, case, solver, rounds):
    reader = csv.reader(io.StringIO(stdout), strict=True)
    if tuple(next(reader, ())) != FIELDS:
        raise ValueError(f"{case}/{solver}: unexpected header")
    rows = []
    seen = set()
    for row in reader:
        if len(row) != len(FIELDS) or row[0] != case:
            raise ValueError(f"{case}/{solver}: malformed row {row!r}")
        try:
            round_index = int(row[1])
            ns = float(row[2])
            checksum = float(row[3])
            covered, total = int(row[4]), int(row[5])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{case}/{solver}: invalid number") from exc
        if (round_index in seen or not (0 <= round_index < rounds)
                or not (math.isfinite(ns) and ns > 0)
                or not (math.isfinite(checksum) and checksum > 0)
                or not (0 <= covered <= total)):
            raise ValueError(f"{case}/{solver}: invalid round, timing, checksum or coverage")
        seen.add(round_index)
        rows.append({"round": round_index, "ns": ns, "checksum": checksum,
                     "covered": covered, "total": total})
    if seen != set(range(rounds)):
        raise ValueError(f"{case}/{solver}: missing rounds")
    return rows


def compare(executable, cases, blocks, iterations, rounds, reference="experimental", candidate="volfi-book"):
    if reference not in SOLVERS or candidate not in SOLVERS or reference == candidate:
        raise ValueError("unknown or identical comparison solvers")
    pair = (reference, candidate)
    samples = []
    for case in cases:
        for block in range(blocks):
            for position, solver in enumerate(rotation(block, reference, candidate)):
                command = [str(executable), str(iterations), str(rounds),
                           "--case", case, "--solver", solver]
                result = subprocess.run(command, capture_output=True, text=True, check=False)
                if result.returncode:
                    raise RuntimeError(f"{case}/{block}/{solver}: {result.stderr or result.stdout}")
                for row in parse_output(result.stdout, case, solver, rounds):
                    samples.append({**row, "case": case, "solver": solver,
                                    "block": block, "position": position})
    summary = []
    for case in cases:
        group = {solver: [r for r in samples if r["case"] == case and r["solver"] == solver]
                 for solver in pair}
        for solver in pair:
            rows = group[solver]
            if len(rows) != 2 * blocks * rounds:
                raise ValueError(f"{case}/{solver}: missing samples")
            if len({(r["block"], r["position"], r["round"]) for r in rows}) != len(rows):
                raise ValueError(f"{case}/{solver}: duplicate rounds")
            if any(rotation(r["block"], reference, candidate)[r["position"]] != solver for r in rows):
                raise ValueError(f"{case}/{solver}: invalid ABBA order")
            if len({r["checksum"] for r in rows}) != 1:
                raise ValueError(f"{case}/{solver}: checksum drift")
            if len({(r["covered"], r["total"]) for r in rows}) != 1:
                raise ValueError(f"{case}/{solver}: case population drift")
        before, after = (statistics.median(r["ns"] for r in group[s]) for s in pair)
        control_coverage = {(r["covered"], r["total"]) for r in group[reference]}
        candidate_coverage = {(r["covered"], r["total"]) for r in group[candidate]}
        if control_coverage != candidate_coverage:
            raise ValueError(f"{case}: different input coverage between solvers")
        covered, total = control_coverage.pop()
        blocks_percent = []
        for block in range(blocks):
            a, b = [statistics.median(r["ns"] for r in group[s] if r["block"] == block)
                    for s in pair]
            blocks_percent.append(100 * (b / a - 1))
        summary.append({
            "case": case, "reference": reference, "candidate": candidate,
            "reference_ns": before, "candidate_ns": after,
            "candidate_change_percent": 100 * (after / before - 1),
            "relative_speedup": before / after,
            "block_changes_percent": blocks_percent,
            "direct_coverage_count": covered, "case_count": total,
            "direct_coverage_ratio": covered / total,
            "checksum_changed": group[reference][0]["checksum"] != group[candidate][0]["checksum"],
            "samples_per_solver": len(group[reference]),
        })
    return samples, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--blocks", type=int, default=4)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=150000)
    parser.add_argument("--features", required=True)
    parser.add_argument("--flags", required=True)
    parser.add_argument("--reference", choices=SOLVERS, default="experimental")
    parser.add_argument("--candidate", choices=SOLVERS, default="volfi-book")
    args = parser.parse_args()
    if min(args.blocks, args.rounds, args.iterations) <= 0:
        parser.error("nonpositive measurement parameters")
    path = args.executable.resolve()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    samples, summary = compare(path, args.cases, args.blocks, args.iterations,
                               args.rounds, args.reference, args.candidate)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError("executable modified during measurement")
    report = {
        "platform": platform.platform(),
        "rustc": subprocess.check_output(["rustc", "-vV"], text=True),
        "exe": str(path), "sha256": digest, "features": args.features,
        "reference": args.reference, "candidate": args.candidate,
        "flags": args.flags, "iterations": args.iterations,
        "rounds": args.rounds, "blocks": args.blocks,
        "method": "same executable; same prepared inputs; rotating ABBA/BAAB; warmup outside timing",
        "caveat": "Synthetic single-threaded scalar quote timings, not official C++ volfi AVX2/AVX-512 throughput. Chart precision requires independent validation.",
        "summary": summary, "samples": samples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"| Workload | {args.reference} ns | {args.candidate} ns | Change | Direct coverage |")
    print("|---|---:|---:|---:|---:|")
    for r in summary:
        print(f"| {r['case']} | {r['reference_ns']:.2f} | {r['candidate_ns']:.2f} | "
              f"{r['candidate_change_percent']:+.2f}% | "
              f"{r['direct_coverage_count']}/{r['case_count']} |")


if __name__ == "__main__":
    main()
