#!/usr/bin/env python3
"""Verify the strict abs(value) > threshold boundary with a round-trip double."""

import argparse
import subprocess
import tempfile
from pathlib import Path

from check_h0 import read_csr


def run(command):
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    if completed.returncode:
        raise SystemExit(completed.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abacus", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    args = parser.parse_args()

    checker = args.source_root / "tests" / "h0" / "check_h0.py"
    with tempfile.TemporaryDirectory(prefix="abacus-h0-threshold-") as directory:
        root = Path(directory)
        first = root / "unfiltered"
        common = ["python3", str(checker), "--abacus", str(args.abacus),
                  "--source-root", str(args.source_root), "--case", "si_soc",
                  "--h0-type", "core", "--skip-baseline"]
        run(common + ["--keep", str(first), "--extra-input",
                      "h0_sparse_threshold 0\nh0_precision 17"])
        _, _, blocks = read_csr(first / "OUT.h0" / "hrs1_nao.csr")
        maximum = max(abs(value) for values, _, _ in blocks.values() for value in values)
        if maximum == 0.0:
            raise AssertionError("baseline unexpectedly contains no nonzero values")

        second = root / "at-boundary"
        # 17 significant digits round-trip the exact binary64 matrix value.
        threshold = format(maximum, ".17g")
        run(common + ["--keep", str(second), "--expect-empty", "--extra-input",
                      "h0_sparse_threshold %s\nh0_precision 17" % threshold])
        print("PASS strict threshold boundary at %s" % threshold)


if __name__ == "__main__":
    main()
