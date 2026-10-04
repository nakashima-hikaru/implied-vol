#!/usr/bin/env python3
"""Compare immutable solver benchmarks in rotating ABBA order.

Checksums must remain stable within each executable. Algorithm changes may
change outputs between versions; independent root accuracy is a separate gate.
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

VERSIONS = ("baseline", "candidate")
SOLVERS = ("experimental", "hybrid")
CASES = (
    "atm", "lowest", "lower_middle", "upper_middle", "highest",
    "near_atm", "near_atm_wider", "mixed_normalised", "mixed_full",
    "deep_otm_full", "near_atm_short_full", "legacy_otm_full",
)
EXPERIMENTAL_CASES = (
    "central_deferred", "wing_seed", "finite_seed", "large_lower", "large_upper",
    "large_asymptotic", "large_near_cap", "microscopic", "mixed_large", "mixed_near_atm",
)
AS1_CASES = ("as1_high_a_guard", "as1_central", "as1_deep_tail", "as1_scaled", "as1_mixed")
FIELDS = ["case", "round", "ns_per_call", "checksum"]
CAVEAT = (
    "Inputs are prepared before timing from fixed and seeded populations. "
    "Each executable must report every requested case and round. "
    "Within-version checksum stability checks reproducibility, not root accuracy. "
    "Different baseline/candidate checksums are allowed; numerical changes require "
    "separate independent accuracy validation. Measurements cover warm, "
    "single-threaded synthetic workloads on this host."
)


class ComparisonError(ValueError):
    """The measurement does not satisfy the benchmark contract."""


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as executable:
        for chunk in iter(lambda: executable.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cases_for(solver, experimental_paths=True, as1_paths=False, dpoly_paths=False):
    if solver not in SOLVERS:
        raise ComparisonError(f"unknown solver {solver!r}")
    if solver == "experimental":
        return CASES + (EXPERIMENTAL_CASES if experimental_paths else ()) + (
            AS1_CASES if as1_paths else ()) + (("dpoly_cell", "dpoly_mixed") if dpoly_paths else ())
    return CASES


def version_order(block):
    return (("baseline", "candidate", "candidate", "baseline") if block % 2 == 0
            else ("candidate", "baseline", "baseline", "candidate"))


def parse_samples(output, version, solver, rounds, block=0, position=0,
                  experimental_paths=True, as1_paths=False, dpoly_paths=False):
    cases = cases_for(solver, experimental_paths, as1_paths, dpoly_paths)
    if version not in VERSIONS or rounds <= 0:
        raise ComparisonError("unknown version or nonpositive round count")
    context = f"block {block}/{version}/{solver}"
    reader = csv.reader(io.StringIO(output), strict=True)
    samples, seen = [], set()
    try:
        if next(reader, None) != FIELDS:
            raise ComparisonError(f"{context}: expected CSV header {FIELDS}")
        for line, row in enumerate(reader, 2):
            if len(row) != len(FIELDS):
                raise ComparisonError(f"{context}: invalid column count on line {line}")
            case, round_text, timing_text, checksum_text = row
            if case not in cases:
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
            if not math.isfinite(checksum) or checksum <= 0.0:
                raise ComparisonError(f"{context}/{case}: checksum must be finite and positive")
            key = case, round_number
            if key in seen:
                raise ComparisonError(f"{context}/{case}: duplicate round {round_number}")
            seen.add(key)
            samples.append({"version": version, "solver": solver, "block": block,
                            "position": position, "case": case, "round": round_number,
                            "ns": timing, "checksum": checksum_text})
    except csv.Error as error:
        raise ComparisonError(f"{context}: malformed CSV: {error}") from error
    expected = {(case, number) for case in cases for number in range(rounds)}
    missing = sorted(expected - seen)
    if missing:
        raise ComparisonError(f"{context}: missing case/round rows: {missing}")
    return samples


def summarize(samples, blocks, rounds, experimental_paths=True, as1_paths=False, dpoly_paths=False):
    if blocks <= 0 or rounds <= 0:
        raise ComparisonError("blocks and rounds must be positive")
    grouped = {(version, solver, case): [] for version in VERSIONS for solver in SOLVERS
               for case in cases_for(solver, experimental_paths, as1_paths, dpoly_paths)}
    seen = set()
    for row in samples:
        version, solver, case = row["version"], row["solver"], row["case"]
        if (version, solver, case) not in grouped:
            raise ComparisonError(f"unknown version/solver/case {version}/{solver}/{case}")
        block, position, number = row["block"], row["position"], row["round"]
        if block not in range(blocks) or position not in range(4) or number not in range(rounds):
            raise ComparisonError(f"{version}/{solver}/{case}: unexpected block, position, or round")
        if version_order(block)[position] != version:
            raise ComparisonError(f"{version}/{solver}/{case}: position violates ABBA order")
        key = version, solver, case, block, position, number
        if key in seen:
            raise ComparisonError(f"{version}/{solver}/{case}: duplicate invocation/round")
        seen.add(key)
        if not math.isfinite(row["ns"]) or row["ns"] <= 0.0:
            raise ComparisonError(f"{version}/{solver}/{case}: timing must be finite and positive")
        checksum = float(row["checksum"])
        if not math.isfinite(checksum) or checksum <= 0.0:
            raise ComparisonError(f"{version}/{solver}/{case}: checksum must be finite and positive")
        grouped[(version, solver, case)].append(row)
    checksums, medians = {}, {}
    for key, rows in grouped.items():
        version, solver, case = key
        if len(rows) != 2 * blocks * rounds:
            raise ComparisonError(f"{version}/{solver}/{case}: missing samples")
        values = {float(row["checksum"]) for row in rows}
        if len(values) != 1:
            raise ComparisonError(f"{version}/{solver}/{case}: checksum changed within version")
        checksums[key] = rows[0]["checksum"]
        medians[key] = statistics.median(row["ns"] for row in rows)
    summary = []
    for solver in SOLVERS:
        for case in cases_for(solver, experimental_paths, as1_paths, dpoly_paths):
            before, after = (medians[(version, solver, case)] for version in VERSIONS)
            pair = [checksums[(version, solver, case)] for version in VERSIONS]
            block_changes = []
            for block in range(blocks):
                a, b = (statistics.median(row["ns"] for row in grouped[(version, solver, case)]
                                         if row["block"] == block) for version in VERSIONS)
                block_changes.append(100.0 * (b / a - 1.0))
            summary.append({"solver": solver, "case": case, "baseline_ns": before,
                            "candidate_ns": after, "change_percent": 100.0 * (after / before - 1.0),
                            "paired_block_changes_percent": block_changes,
                            "samples_per_version": 2 * blocks * rounds,
                            "baseline_checksum": pair[0], "candidate_checksum": pair[1],
                            "checksum_changed_between_versions": float(pair[0]) != float(pair[1])})
    return summary


def run_comparison(args):
    executables = {"baseline": args.baseline.resolve(), "candidate": args.candidate.resolve()}
    if args.output.resolve() in executables.values():
        raise ComparisonError("output must be distinct from the benchmark executables")
    hashes = {name: sha256(path) for name, path in executables.items()}
    samples, invocations = [], []
    for block in range(args.blocks):
        solvers = SOLVERS if block % 2 == 0 else tuple(reversed(SOLVERS))
        for solver in solvers:
            for position, version in enumerate(version_order(block)):
                command = [str(executables[version]), str(args.iterations), str(args.rounds),
                           "--solver", solver]
                if solver == "experimental":
                    if args.experimental_paths:
                        command.append("--experimental-paths")
                    if args.as1_paths:
                        command.append("--as1-paths")
                    if args.dpoly_paths:
                        command.append("--dpoly-paths")
                result = subprocess.run(command, capture_output=True, text=True, check=False)
                if result.returncode:
                    raise ComparisonError(
                        f"block {block}/{version}/{solver}: benchmark exited with {result.returncode}: "
                        f"{result.stderr.strip() or result.stdout.strip()}")
                samples.extend(parse_samples(result.stdout, version, solver, args.rounds, block,
                                             position, args.experimental_paths, args.as1_paths, args.dpoly_paths))
                invocations.append({"version": version, "solver": solver, "block": block,
                                    "position": position, "command": command})
    summary = summarize(samples, args.blocks, args.rounds, args.experimental_paths, args.as1_paths, args.dpoly_paths)
    final_hashes = {name: sha256(path) for name, path in executables.items()}
    if final_hashes != hashes:
        raise ComparisonError("executable changed during measurement")
    report = {"machine": platform.platform(), "measurement_host_rustc": subprocess.check_output(
        ["rustc", "-vV"], text=True), "reported_build": {"flags": args.flags,
        "features": args.features}, "iterations": args.iterations, "rounds": args.rounds,
        "blocks": args.blocks, "executables": {name: {"path": str(path), "sha256_before": hashes[name],
        "sha256_after": final_hashes[name]} for name, path in executables.items()},
        "cases": {solver: list(cases_for(solver, args.experimental_paths, args.as1_paths, args.dpoly_paths)) for solver in SOLVERS},
        "population_flags": {"experimental_paths": args.experimental_paths, "as1_paths": args.as1_paths, "dpoly_paths": args.dpoly_paths},
        "caveat": CAVEAT, "invocations": invocations, "summary": summary, "samples": samples}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=200_000)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--blocks", type=int, default=4)
    parser.add_argument("--flags", default="-C target-cpu=native",
                        help="caller-supplied build flags; not inferred from binaries")
    parser.add_argument("--features", default="flashiv,experimental",
                        help="caller-supplied build features; not inferred from binaries")
    parser.add_argument("--experimental-paths", action=argparse.BooleanOptionalAction, default=True,
                        help="include Experimental's extended populations (enabled by default)")
    parser.add_argument("--as1-paths", action="store_true",
                        help="include Experimental's fixed and seeded AS1 populations")
    parser.add_argument("--dpoly-paths", action="store_true",
                        help="include Experimental's affected Mills polynomial populations")
    args = parser.parse_args()
    if min(args.iterations, args.rounds, args.blocks) <= 0:
        parser.error("iterations, rounds, and blocks must be positive")
    if not args.flags.strip() or not args.features.strip():
        parser.error("flags and features must include caller-supplied provenance")
    try:
        for row in run_comparison(args):
            print(f"{row['solver']:12s} {row['case']:20s} {row['baseline_ns']:8.2f} -> "
                  f"{row['candidate_ns']:8.2f} ns ({row['change_percent']:+.2f}%)")
    except (ComparisonError, OSError) as error:
        parser.exit(1, f"comparison failed: {error}\n")


if __name__ == "__main__":
    main()
