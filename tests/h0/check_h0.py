#!/usr/bin/env python3
"""Run one H0 case and validate/compare ABACUS text CSR output."""

import argparse
import cmath
import gzip
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
COMPLEX = re.compile(r"\((%s),(%s)\)" % (FLOAT, FLOAT))


def _payload(lines, begin, end):
    return " ".join(line.strip() for line in lines[begin:end] if not line.lstrip().startswith("#"))


def _values(text):
    if "(" in text:
        return [complex(float(real), float(imag)) for real, imag in COMPLEX.findall(text)]
    return [complex(float(value), 0.0) for value in text.split()]


def read_csr(path):
    if path.suffix == ".gz":
        with gzip.open(str(path), "rt", encoding="utf-8") as stream:
            lines = stream.read().splitlines()
    else:
        lines = path.read_text(encoding="utf-8").splitlines()
    metadata = next(line.strip() for line in lines if "h0_type=" in line)
    nbasis = int(next(re.match(r"\s*(\d+)\s+# number of localized basis", line).group(1)
                      for line in lines
                      if re.match(r"\s*(\d+)\s+# number of localized basis", line)))
    nr = int(next(re.match(r"\s*(\d+)\s+# number of Bravais", line).group(1)
                  for line in lines
                  if re.match(r"\s*(\d+)\s+# number of Bravais", line)))
    start = next(i for i, line in enumerate(lines) if "CSR Format" in line)
    blocks = {}
    i = start + 1
    header = re.compile(r"\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(\d+)\s*$")
    while i < len(lines):
        match = header.match(lines[i])
        if not match:
            i += 1
            continue
        rvec = tuple(int(match.group(j)) for j in range(1, 4))
        nnz = int(match.group(4))
        value_mark = i + 1
        while "# CSR values" not in lines[value_mark]:
            value_mark += 1
        col_mark = value_mark + 1
        while "# CSR column indices" not in lines[col_mark]:
            col_mark += 1
        row_mark = col_mark + 1
        while "# CSR row pointers" not in lines[row_mark]:
            row_mark += 1
        next_header = row_mark + 1
        while next_header < len(lines) and not header.match(lines[next_header]):
            next_header += 1
        values = _values(_payload(lines, value_mark + 1, col_mark))
        cols_text = _payload(lines, col_mark + 1, row_mark)
        rows_text = _payload(lines, row_mark + 1, next_header)
        cols = [int(value) for value in cols_text.split()] if cols_text else []
        rows = [int(value) for value in rows_text.split()] if rows_text else []
        if rvec in blocks:
            raise AssertionError("duplicate R block: %s" % (rvec,))
        if len(values) != nnz or len(cols) != nnz:
            raise AssertionError("bad nnz arrays for R=%s" % (rvec,))
        if len(rows) != nbasis + 1 or rows[0] != 0 or rows[-1] != nnz:
            raise AssertionError("bad row pointer for R=%s" % (rvec,))
        if any(a > b for a, b in zip(rows, rows[1:])):
            raise AssertionError("non-monotonic row pointer for R=%s" % (rvec,))
        for row in range(nbasis):
            row_cols = cols[rows[row]:rows[row + 1]]
            if row_cols != sorted(set(row_cols)):
                raise AssertionError("columns are not sorted and unique for R=%s row=%d" % (rvec, row))
            if any(col < 0 or col >= nbasis for col in row_cols):
                raise AssertionError("column out of range for R=%s row=%d" % (rvec, row))
        blocks[rvec] = (values, cols, rows)
        i = next_header
    if len(blocks) != nr:
        raise AssertionError("header declares %d R blocks, parsed %d" % (nr, len(blocks)))
    return metadata, nbasis, blocks


def dense_entries(nbasis, blocks):
    entries = {}
    for rvec, (values, cols, rows) in blocks.items():
        for row in range(nbasis):
            for offset in range(rows[row], rows[row + 1]):
                entries[(rvec, row, cols[offset])] = values[offset]
    return entries


def check_hermitian(nbasis, blocks, tolerance):
    entries = dense_entries(nbasis, blocks)
    for (rvec, row, col), value in entries.items():
        reverse = ((-rvec[0], -rvec[1], -rvec[2]), col, row)
        other = entries.get(reverse, 0j)
        if abs(value - other.conjugate()) > tolerance:
            raise AssertionError("Hermitian relation failed at R=%s (%d,%d): %r vs %r"
                                 % (rvec, row, col, value, other.conjugate()))


def compare(actual, expected, tolerance):
    _, actual_nbasis, actual_blocks = read_csr(actual)
    _, expected_nbasis, expected_blocks = read_csr(expected)
    if actual_nbasis != expected_nbasis or actual_blocks.keys() != expected_blocks.keys():
        raise AssertionError("CSR dimensions or R support differ from baseline")
    for rvec in actual_blocks:
        av, ac, ar = actual_blocks[rvec]
        ev, ec, er = expected_blocks[rvec]
        if ac != ec or ar != er:
            raise AssertionError("CSR sparsity differs from baseline at R=%s" % (rvec,))
        for index, (a, e) in enumerate(zip(av, ev)):
            limit = tolerance * max(1.0, abs(e))
            if abs(a - e) > limit:
                raise AssertionError("value differs at R=%s offset=%d: %r vs %r" % (rvec, index, a, e))


def prepare_input(case_dir, h0_type, init_chg, source_root, extra_input,
                  set_inputs, remove_inputs):
    text = (case_dir / "INPUT").read_text(encoding="utf-8")
    text = re.sub(r"(?m)^h0_type\s+\S+", "h0_type %s" % h0_type, text)
    text += "\ninit_chg %s\n" % init_chg
    text += extra_input
    pp_orb = source_root / "tests" / "PP_ORB"
    text = re.sub(r"(?m)^pseudo_dir\s+\S+", "pseudo_dir %s" % pp_orb, text)
    text = re.sub(r"(?m)^orbital_dir\s+\S+", "orbital_dir %s" % pp_orb, text)
    for assignment in set_inputs:
        key, separator, value = assignment.partition("=")
        if not separator or not key or not value:
            raise ValueError("--set-input requires KEY=VALUE")
        pattern = r"(?m)^%s\s+.*$" % re.escape(key)
        replacement = "%s %s" % (key, value)
        if re.search(pattern, text):
            text = re.sub(pattern, replacement, text)
        else:
            text += "\n%s\n" % replacement
    for key in remove_inputs:
        text = re.sub(r"(?m)^%s\s+.*\n?" % re.escape(key), "", text)
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abacus", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--h0-type", required=True, choices=("core", "full"))
    parser.add_argument("--init-chg", default="atomic", choices=("atomic", "auto", "file"))
    parser.add_argument("--expected", type=Path,
                        help="optional expected file for the first output; the frozen suite baseline is used by default")
    parser.add_argument("--density-file", type=Path,
                        help="copy a chg.cube fixture and require that init_chg reads it")
    parser.add_argument("--skip-baseline", action="store_true")
    parser.add_argument("--expect-failure", type=str,
                        help="require a failed run whose output contains this text")
    parser.add_argument("--extra-input", default="")
    parser.add_argument("--set-input", action="append", default=[])
    parser.add_argument("--remove-input", action="append", default=[])
    parser.add_argument("--expect-empty", action="store_true",
                        help="require every retained R block to have zero stored entries")
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--keep", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()

    args.source_root = args.source_root.resolve()
    args.abacus = args.abacus.resolve()
    case_dir = args.source_root / "tests" / "h0" / args.case
    temporary = tempfile.TemporaryDirectory(prefix="abacus-h0-") if args.keep is None else None
    run_dir = Path(temporary.name) if temporary else args.keep
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case_dir / "STRU", run_dir / "STRU")
    extra_input = args.extra_input
    if args.density_file:
        density_dir = run_dir / "initial_density"
        density_dir.mkdir(exist_ok=True)
        shutil.copy2(args.density_file.resolve(), density_dir / "chg.cube")
        extra_input += "\nread_file_dir initial_density\n"
    (run_dir / "INPUT").write_text(prepare_input(case_dir,
                                                  args.h0_type,
                                                  args.init_chg,
                                                  args.source_root,
                                                  extra_input,
                                                  args.set_input,
                                                  args.remove_input),
                                         encoding="utf-8")
    env = dict(__import__("os").environ)
    env["OMP_NUM_THREADS"] = str(args.omp_threads)
    command = [str(args.abacus)]
    if args.mpi_ranks > 1:
        if args.mpiexec is None:
            raise ValueError("--mpiexec is required when --mpi-ranks is greater than one")
        command = [str(args.mpiexec), "-np", str(args.mpi_ranks)] + command
    completed = subprocess.run(command, cwd=str(run_dir), env=env,
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (run_dir / "console.log").write_text(completed.stdout, encoding="utf-8")
    if args.expect_failure:
        if completed.returncode == 0:
            raise AssertionError("ABACUS unexpectedly succeeded")
        if args.expect_failure not in completed.stdout:
            raise AssertionError("expected failure text not found: %s" % args.expect_failure)
        print("PASS expected failure: %s" % args.expect_failure)
        return
    if completed.returncode:
        raise SystemExit("ABACUS failed with exit code %d; see %s" % (completed.returncode,
                                                                      run_dir / "console.log"))
    log = (run_dir / "OUT.h0" / "running_get_h0.log").read_text(encoding="utf-8")
    combined_log = completed.stdout + "\n" + log
    if args.density_file and "Read electron density from file:" not in combined_log:
        raise AssertionError("init_chg did not read the supplied density file")
    if args.init_chg == "auto" and not args.density_file and "use atomic initialization instead" not in combined_log:
        raise AssertionError("init_chg=auto did not report its atomic fallback")
    if re.search(r"\b(?:SCF|diagonalization|iteration)\b", log, re.IGNORECASE):
        raise AssertionError("solver/iteration text found in one-shot H0 log")
    outputs = sorted((run_dir / "OUT.h0").glob("hrs*_nao.csr"))
    expected_count = 2 if (args.case == "fe2_collinear" and args.h0_type == "full") else 1
    if len(outputs) != expected_count:
        raise AssertionError("expected %d CSR outputs, found %d" % (expected_count, len(outputs)))
    for output in outputs:
        metadata, nbasis, blocks = read_csr(output)
        if "h0_type=%s" % args.h0_type not in metadata:
            raise AssertionError("incorrect H0 metadata")
        check_hermitian(nbasis, blocks, args.tolerance)
        if args.expect_empty:
            if not blocks or any(len(values) != 0 for values, _, _ in blocks.values()):
                raise AssertionError("thresholded output did not preserve only empty R blocks")
        if not args.skip_baseline:
            baseline = (args.source_root / "tests" / "h0" / "reference"
                        / ("%s_%s_%s.gz" % (args.case, args.h0_type, output.name)))
            if not baseline.is_file():
                raise AssertionError("missing frozen baseline: %s" % baseline)
            compare(output, baseline, args.tolerance)
    if args.expected:
        compare(outputs[0], args.expected, args.tolerance)
    print("PASS case=%s h0_type=%s outputs=%d" % (args.case, args.h0_type, len(outputs)))


if __name__ == "__main__":
    main()
