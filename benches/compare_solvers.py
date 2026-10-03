#!/usr/bin/env python3
"""Compare Hybrid, Jaeckel, and C++ in one immutable benchmark executable."""

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


SOLVERS = ("hybrid", "jaeckel", "cpp")
CASES = (
    "atm", "lowest", "lower_middle", "upper_middle", "highest",
    "near_atm", "near_atm_wider", "mixed_normalised", "mixed_full",
    "deep_otm_full", "near_atm_short_full", "legacy_otm_full",
)
FIELDS = ["case", "round", "ns_per_call", "checksum"]
CAVEAT = (
    "C++ retains its original -Ofast and -ffp-contract=fast build policy; "
    "Rust uses the caller-reported flags and explicit FMA feature policy. "
    "These arithmetic policies differ. Inputs are prepared before timing and "
    "shared across solvers; checksums may differ between solvers. "
    "Output preflight checks validity, not root accuracy; the solvers retain "
    "their different numerical contracts. "
    "Measurements cover warm, single-threaded synthetic workloads on this host."
)


class ComparisonError(ValueError):
    """The measurement does not satisfy the benchmark contract."""


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as executable:
        for chunk in iter(lambda: executable.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def measurement_order(block):
    offset = block % len(SOLVERS)
    order = SOLVERS[offset:] + SOLVERS[:offset]
    return tuple(reversed(order)) if block % 2 else order


def parse_samples(output, solver, rounds, block=0, position=0):
    if solver not in SOLVERS or rounds <= 0:
        raise ComparisonError("unknown solver or nonpositive round count")
    reader = csv.reader(io.StringIO(output), strict=True)
    samples = []
    seen = set()
    context = f"block {block}/{solver}"
    try:
        if next(reader, None) != FIELDS:
            raise ComparisonError(f"{context}: expected CSV header {FIELDS}")
        for line, row in enumerate(reader, 2):
            if len(row) != len(FIELDS):
                raise ComparisonError(f"{context}: invalid column count on line {line}")
            case, round_text, timing_text, checksum_text = row
            if case not in CASES:
                raise ComparisonError(f"{context}: unknown case {case!r}")
            try:
                round_number = int(round_text)
                timing = float(timing_text)
                checksum = float(checksum_text)
            except (ValueError, OverflowError) as error:
                raise ComparisonError(f"{context}: invalid number on line {line}") from error
            if round_number not in range(rounds):
                raise ComparisonError(f"{context}: unexpected round {round_number}")
            if not math.isfinite(timing) or timing <= 0.0:
                raise ComparisonError(f"{context}/{case}: timing must be finite and positive")
            if not math.isfinite(checksum):
                raise ComparisonError(f"{context}/{case}: checksum must be finite")
            key = (case, round_number)
            if key in seen:
                raise ComparisonError(f"{context}/{case}: duplicate round {round_number}")
            seen.add(key)
            samples.append({
                "solver": solver, "block": block, "position": position,
                "case": case, "round": round_number, "ns_per_call": timing,
                "checksum": checksum_text,
            })
    except csv.Error as error:
        raise ComparisonError(f"{context}: malformed CSV: {error}") from error
    expected = {(case, round_number) for case in CASES for round_number in range(rounds)}
    missing = sorted(expected - seen)
    if missing:
        raise ComparisonError(f"{context}: missing case/round rows: {missing}")
    return samples


def summarize(samples, blocks, rounds):
    grouped = {(solver, case): [] for solver in SOLVERS for case in CASES}
    seen = set()
    for sample in samples:
        solver, case = sample["solver"], sample["case"]
        if (solver, case) not in grouped:
            raise ComparisonError(f"unknown solver/case {solver}/{case}")
        block, round_number = sample["block"], sample["round"]
        if block not in range(blocks) or round_number not in range(rounds):
            raise ComparisonError(f"{solver}/{case}: unexpected block or round")
        key = (solver, case, block, round_number)
        if key in seen:
            raise ComparisonError(f"{solver}/{case}: duplicate block/round")
        seen.add(key)
        grouped[(solver, case)].append(sample)

    medians = {}
    for (solver, case), rows in grouped.items():
        if len(rows) != blocks * rounds:
            raise ComparisonError(f"{solver}/{case}: missing samples")
        checksums = {float(row["checksum"]) for row in rows}
        if len(checksums) != 1:
            raise ComparisonError(f"{solver}/{case}: checksum changed within solver")
        medians[(solver, case)] = statistics.median(row["ns_per_call"] for row in rows)

    return [{
        "case": case, "solver": solver,
        "median_ns_per_call": medians[(solver, case)],
        "time_change_vs_cpp_percent": 100.0 * (
            medians[(solver, case)] / medians[("cpp", case)] - 1.0),
        "samples": len(grouped[(solver, case)]),
        "checksum": grouped[(solver, case)][0]["checksum"],
    } for case in CASES for solver in SOLVERS]


def markdown_report(report):
    rows = {(row["solver"], row["case"]): row for row in report["summary"]}
    lines = [
        "# Hybrid and LBR comparison", "",
        f"{report['iterations']:,} calls per sample; {report['blocks']} blocks "
        f"of {report['rounds']} rounds, giving "
        f"{report['blocks'] * report['rounds']} samples per solver and case.", "",
        "Cells show median ns/call and time change versus C++. "
        "Negative percentages mean less time; positive percentages mean more time.", "",
        "| Case | Hybrid | Jaeckel | C++ |",
        "|---|---:|---:|---:|",
    ]
    for case in CASES:
        cells = [case]
        for solver in SOLVERS:
            row = rows[(solver, case)]
            value = f"{row['median_ns_per_call']:.2f} ns"
            if solver != "cpp":
                value += f" ({row['time_change_vs_cpp_percent']:+.2f}%)"
            else:
                value += " (reference)"
            cells.append(value)
        lines.append("| " + " | ".join(cells) + " |")
    build = report["reported_build"]
    lines.extend([
        "", CAVEAT, "",
        f"Caller-reported Rust flags: `{build['flags']}`. "
        f"Features: `{build['features']}`. Build provenance is supplied by the caller, "
        "not inferred from the executable.", "",
        f"Host: {report['machine']}. Executable SHA-256 before and after: "
        f"`{report['executable']['sha256_before']}`.", "",
    ])
    return "\n".join(lines)


def run_comparison(args):
    executable = args.executable.resolve()
    destinations = [args.output.resolve(), args.summary.resolve()]
    if len(set(destinations + [executable])) != 3:
        raise ComparisonError("executable, JSON output, and Markdown summary must be distinct files")
    initial_hash = sha256(executable)
    samples = []
    invocations = []
    for block in range(args.blocks):
        for position, solver in enumerate(measurement_order(block)):
            command = [str(executable), str(args.iterations), str(args.rounds), "--solver", solver]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            if result.returncode:
                raise ComparisonError(
                    f"block {block}/{solver}: benchmark exited with {result.returncode}: "
                    f"{result.stderr.strip() or result.stdout.strip()}")
            samples.extend(parse_samples(result.stdout, solver, args.rounds, block, position))
            invocations.append({"block": block, "position": position, "solver": solver,
                                "command": command})
    summary = summarize(samples, args.blocks, args.rounds)
    final_hash = sha256(executable)
    if final_hash != initial_hash:
        raise ComparisonError("executable changed during measurement")
    report = {
        "machine": platform.platform(), "architecture": platform.machine(),
        "reported_build": {"flags": args.flags, "features": args.features},
        "iterations": args.iterations, "rounds": args.rounds, "blocks": args.blocks,
        "executable": {"path": str(executable), "sha256_before": initial_hash,
                       "sha256_after": final_hash},
        "cases": list(CASES), "solvers": list(SOLVERS), "caveat": CAVEAT,
        "invocations": invocations, "summary": summary, "samples": samples,
    }
    markdown = markdown_report(report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    args.summary.write_text(markdown, encoding="utf-8")
    return markdown


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--flags", required=True,
                        help="caller-supplied Rust build flags, not inferred from the executable")
    parser.add_argument("--features", required=True,
                        help="caller-supplied build features, not inferred from the executable")
    parser.add_argument("--iterations", type=int, default=200_000)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--blocks", type=int, default=6)
    args = parser.parse_args()
    if min(args.iterations, args.rounds, args.blocks) <= 0:
        parser.error("iterations, rounds, and blocks must be positive")
    if not args.flags.strip() or not args.features.strip():
        parser.error("flags and features must include caller-supplied provenance")
    try:
        print(run_comparison(args), end="")
    except (ComparisonError, OSError) as error:
        parser.exit(1, f"comparison failed: {error}\n")


if __name__ == "__main__":
    main()
