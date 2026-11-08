#!/usr/bin/env python3
"""Re-run the *original* continuous Experimental certificate on an exact source bridge.

The original proof depends on an external, hash-bound reference archive named
'implied-black-volatility'. This script will NOT substitute synthesized data,
assume its result, or silently skip it if the archive is absent.

The existing direct-D checker is run against a temporary canonical source tree
obtained by reverting only the enumerated finish specializations; the real
optimized repository is not mutated. To apply that result to the optimized
binary64 implementation also requires the conditional IEEE transfer lemmas
in docs/experimental-finish-transfer.md. Replaying the old certificate alone
does not independently mechanize those lemmas.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from verify_finish_scalar_bridge import (
    PINNED_MAIN,
    PINNED_SMALL_GIT_BLOB,
    git_blob_sha,
    restore_source,
)

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root",type=Path,required=True,
                        help="local immutable archived implied-black-volatility proof tree")
    parser.add_argument("--simd-root",type=Path,
                        help="optional reviewed fearless_simd 1.0.0 crate source")
    parser.add_argument("--out",type=Path,default=ROOT/"target/proof/finish-source-transfer.json")
    parser.add_argument("--transfer-samples",type=int,default=80000,
                        help="deterministic f64 equivalence cases checked before the old proof")
    args=parser.parse_args()
    # A successful archived proof is necessary but not sufficient:
    # re-run this PR's source and FP transfer validation as part of the
    # same invocation, rather than leaving it to an external checklist.
    transfer_script=ROOT/"scripts/verify_finish_transfer.py"
    if not 0 < args.transfer_samples <= 1_000_000:
        parser.error("--transfer-samples must be in 1..1000000")
    with tempfile.TemporaryDirectory(prefix="implied-vol-fp-transfer-") as fp_tmp:
        transfer_output=Path(fp_tmp)/"transfer.json"
        transfer_cmd=[
            sys.executable,str(transfer_script),
            "--source",str(ROOT/"src/experimental/small.rs"),
            "--report",str(transfer_output),
            "--samples",str(args.transfer_samples),
        ]
        transfer=subprocess.run(transfer_cmd,capture_output=True,text=True,check=False)
        if transfer.returncode:
            raise RuntimeError(
                f"optimized arithmetic transfer failed, code={transfer.returncode}"
                f"\nSTDOUT:\n{transfer.stdout[-5000:]}\nSTDERR:\n{transfer.stderr[-5000:]}"
            )
        transfer_report=json.loads(transfer_output.read_text())
    if not transfer_report["finite_checks"]["all_finite_bit_checks_passed"]:
        raise ValueError("the optimized arithmetic transfer checks failed")
    archive=args.reference_root.resolve()
    required=[archive/"src/small.rs",
              archive/"results/global_proof/small_seed/coverage.json",
              archive/"results/global_proof/small_seed/finish.json"]
    missing=[str(x) for x in required if not x.is_file()]
    if missing:
        parser.error("reviewed external reference archive unavailable; exact missing paths: "
                     +", ".join(missing))
    candidate=(ROOT/"src/experimental/small.rs").read_text()
    canonical=restore_source(candidate)
    sha=git_blob_sha(canonical.encode())
    if sha!=PINNED_SMALL_GIT_BLOB:
        raise ValueError(f"optimized-to-original source bridge mismatch: {sha}")
    # Exclude external build state, Git metadata and generated artifacts.
    with tempfile.TemporaryDirectory(prefix="implied-vol-finish-proof-") as tmp:
        tmpdir=Path(tmp)/"canonical"
        shutil.copytree(ROOT,tmpdir,ignore=shutil.ignore_patterns(
            ".git","target","__pycache__", ".venv", ".pytest_cache",
            "benchmark-results", "finish-results"
        ))
        (tmpdir/"src/experimental/small.rs").write_text(canonical)
        checker=tmpdir/"scripts/verify_direct_d.py"
        output=Path(tmp)/"proof.json"
        cmd=[sys.executable,str(checker),"--source-root",str(tmpdir),
             "--reference-root",str(archive),"--out",str(output)]
        if args.simd_root is not None:
            cmd.extend(("--simd-root",str(args.simd_root.resolve())))
        process=subprocess.run(cmd,capture_output=True,text=True,check=False)
        if process.returncode:
            raise RuntimeError("archived continuous source proof failed; "
                               f"exit={process.returncode}\nSTDOUT:\n{process.stdout[-6000:]}\n"
                               f"STDERR:\n{process.stderr[-6000:]}")
        result=json.loads(output.read_text())
        bounds={name:float(v["candidate_display"])
                for name,v in result["route_bounds"].items()}
        if not bounds or not all(0 <= val < 1 for val in bounds.values()):
            raise ValueError(f"original conditional proof does not establish rho_J<1: {bounds}")
        report={
            "status":"original continuous proof replayed on exactly restored main source",
            "reviewed_main_commit":PINNED_MAIN,
            "canonical_small_blob":sha,
            "optimized_source_sha256":hashlib.sha256(candidate.encode()).hexdigest(),
            "reference_root":str(archive),
            "reference_bounds_rho_J":bounds,
            "proof_result_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
            "transfer_finite_counts":transfer_report["finite_checks"]["counts"],
            "transfer_source_blob":transfer_report["pinned_baseline_blob"],
            "transfer_envelopes":transfer_report["exact_envelopes"],
            "qualification":"conditional on the separately stated IEEE equivalence "
                            "lemmas and all original proof/environment assumptions",
            "not_claimed":["a newly mechanized complete IEEE proof",
                           "an unqualified whole-domain theorem",
                           "reproducibility without the external reference archive"],
        }
    dest=args.out.resolve()
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
