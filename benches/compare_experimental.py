#!/usr/bin/env python3
"""Compare separately built immutable solver benchmarks in rotating ABBA order."""

import argparse
import csv
import hashlib
import io
import json
import platform
import statistics
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=200_000)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--blocks", type=int, default=4)
    parser.add_argument("--flags", default="-C target-cpu=native",
                        help="build flags supplied by the caller; not inferred from binaries")
    parser.add_argument("--features", default="flashiv,experimental",
                        help="build features supplied by the caller; not inferred from binaries")
    args = parser.parse_args()
    if min(args.iterations, args.rounds, args.blocks) <= 0:
        parser.error("iterations, rounds, and blocks must be positive")
    executables = {"baseline": args.baseline.resolve(), "candidate": args.candidate.resolve()}
    hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest()
              for name, path in executables.items()}
    samples = []
    for block in range(args.blocks):
        # Reverse which executable gets the two central positions per block.
        order = ["baseline", "candidate", "candidate", "baseline"]
        if block % 2:
            order = ["candidate", "baseline", "baseline", "candidate"]
        solvers = ["experimental", "hybrid"] if block % 2 == 0 else ["hybrid", "experimental"]
        for solver in solvers:
            for version in order:
                command = [str(executables[version]), str(args.iterations), str(args.rounds),
                           "--solver", solver]
                if solver == "experimental":
                    command.append("--experimental-paths")
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                reader = csv.DictReader(io.StringIO(result.stdout))
                rows = list(reader)
                if not rows or reader.fieldnames != ["case", "round", "ns_per_call", "checksum"]:
                    raise RuntimeError(f"{version}/{solver}: missing benchmark CSV output")
                for row in rows:
                    samples.append({"version": version, "solver": solver, "block": block,
                                    "case": row["case"], "round": int(row["round"]),
                                    "ns": float(row["ns_per_call"]), "checksum": row["checksum"]})
    summary = []
    for solver in ["experimental", "hybrid"]:
        cases = dict.fromkeys(row["case"] for row in samples if row["solver"] == solver)
        for case in cases:
            pair = {version: [row for row in samples if row["solver"] == solver
                              and row["case"] == case and row["version"] == version]
                    for version in executables}
            checksums = {row["checksum"] for rows in pair.values() for row in rows}
            if len(checksums) != 1:
                raise RuntimeError(f"{solver}/{case}: checksums changed: {checksums}")
            before, after = (statistics.median(row["ns"] for row in pair[version])
                             for version in ["baseline", "candidate"])
            summary.append({"solver": solver, "case": case, "baseline_ns": before,
                            "candidate_ns": after, "change_percent": 100 * (after / before - 1),
                            "samples_per_version": len(pair["baseline"]),
                            "checksum": next(iter(checksums))})
            print(f"{solver:12s} {case:20s} {before:8.2f} -> {after:8.2f} ns "
                  f"({100 * (after / before - 1):+.2f}%)")
    for name, path in executables.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != hashes[name]:
            raise RuntimeError(f"{name}: executable changed during measurement")
    report = {"machine": platform.platform(), "measurement_host_rustc": subprocess.check_output(
        ["rustc", "-vV"], text=True), "reported_build": {"flags": args.flags,
        "features": args.features}, "iterations": args.iterations, "rounds": args.rounds,
        "blocks": args.blocks, "executables": {name: {"path": str(path), "sha256": hashes[name]}
        for name, path in executables.items()}, "summary": summary, "samples": samples}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
