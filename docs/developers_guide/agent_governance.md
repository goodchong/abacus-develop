# ABACUS H0 Development Governance

This document defines the maintenance boundary for the H0-only executable.

## Scope

The program reads `INPUT`, `STRU`, pseudopotentials, and numerical atomic
orbitals, constructs a localized-orbital real-space Hamiltonian, writes text
CSR H(R), and exits.

Supported formulas are:

- `core`: `T + V_nl + V_loc`
- `full`: `T + V_nl + V_loc + V_H[rho0] + V_xc[rho0]`

The full path must preserve the legacy first-Hamiltonian initialization order:
initialize/read `rho0`, establish nonlinear core charge, normalize, construct
the initial local potential, project it into H(R), and exit before any solver or
density update.

## Architectural rules

- The control flow is one-way: validate input, initialize cell and basis,
  initialize grids and neighbors, construct H0, gather, write, exit.
- The core path must not allocate charge density, density matrices, or solver
  state.
- The full path may retain only the charge and potential facilities required to
  construct the initial Hartree and built-in LDA/GGA XC potentials.
- Keep only CPU code, with optional OpenMP and optional MPI plus ScaLAPACK.
- Preserve owning HContainer storage and union-add semantics across kinetic,
  nonlocal, and local-potential supports.
- Root rank writes after MPI accumulation. Preserve R blocks, including empty
  blocks, sorted columns, unique entries, and CSR row pointers.
- Do not reintroduce removed ABACUS workflows or dependencies as dormant code.

## Compatibility rules

- C++11 and LF are mandatory.
- Avoid new `GlobalV`, `GlobalC`, or `PARAM` dependencies. Any unavoidable
  retained global access must remain non-increasing and be documented.
- Keep headers narrow and source registration deterministic.
- Do not add default arguments to established interfaces.
- INPUT changes require generated YAML and Markdown documentation updates.

## Required tests

Changes affecting numerical construction or output must cover:

- `core/full` with physical `nspin=1,2,4`;
- exact comparison with frozen pre-crop first-Hamiltonian CSR references;
- `init_chg=atomic/file/auto`, including fallback and missing-file failure;
- Hermitian H(R), threshold boundary behavior, empty R blocks, sorted and
  duplicate-free columns, and valid row pointers;
- serial, OpenMP, and MPI accumulation behavior when those paths are affected;
- confirmation that no KPT file, SCF loop, diagonalization, or density update is
  used.

Reference files may be updated only for an intentional numerical-contract
change with an independently explained baseline.

## Verification commands

```bash
cmake -S . -B build-serial -DENABLE_MPI=OFF -DENABLE_OPENMP=ON -DBUILD_TESTING=ON
cmake --build build-serial -j2
ctest --test-dir build-serial --output-on-failure

cmake -S . -B build-mpi -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DBUILD_TESTING=ON
cmake --build build-mpi -j2
ctest --test-dir build-mpi --output-on-failure

git diff --check
python3 tools/03_code_analysis/agent_governance_check.py --staged
```

MPI configuration and runtime tests must be repeated outside a restricted
sandbox if OpenMPI reports socket or process-visibility errors.

## Review record

Every change should report the exact commands run, pass/fail counts, expected
environment limitations, INPUT/doc changes, and any deliberate numerical
reference changes. Do not claim verification that was not run freshly.

## Initial extraction rationale

The H0-only extraction moves a small set of initialization functions out of
the deleted solver tree while removing the overwhelming majority of legacy
global accesses. The governance checker reports 77 added and 6969 removed
`GlobalV`/`GlobalC`/`PARAM` references (net -6892). The retained accesses are
limited to legacy INPUT state, process/rank state, and the existing ABACUS log
streams needed by the surviving initialization and integration kernels. This
is a migration-neutral extraction at each moved call site and a substantial
repository-level dependency reduction; new H0 work must not increase it.

Headers flagged for review either own the referenced value type or expose a
template whose complete types are required at instantiation. The retained
`phi_operator.hpp` is narrowly scoped to that template implementation and is
included only by its declaring header.
