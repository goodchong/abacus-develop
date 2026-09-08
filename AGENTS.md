# ABACUS H0 Agent Instructions

This repository is the one-shot LCAO H0 matrix generator. Read
`docs/developers_guide/agent_governance.md` before changing code.

## Product boundary

- Keep `calculation=get_h0` as the only workflow.
- Keep `basis_type=lcao` as the only basis.
- `h0_type=core` means `T + V_nl + V_loc`.
- `h0_type=full` adds `V_H[rho0] + V_xc[rho0]` from the initial density.
- Never add diagonalization, occupations, wavefunctions, density matrices,
  mixing, SCF/NSCF, forces, stress, MD, relaxation, K points, H(k), or S(R).
- Do not add GPU, Libxc, EXX, ELPA, PEXSI, Python, or API dependencies.

## Coding rules

- Keep C++11 compatibility and LF line endings.
- Pass workflow dependencies explicitly where practical. Do not increase
  cross-layer control through `GlobalV`, `GlobalC`, or `PARAM`.
- Keep header dependencies minimal. Do not add default arguments to existing
  interfaces or implementation-only `.hpp` headers.
- New source files must be listed deterministically in CMake.
- H0 containers must own their storage; do not wrap borrowed buffers whose
  lifetime or zeroing is controlled elsewhere.

## Verification

- Set `OMP_NUM_THREADS=1` for reference runtime tests.
- Build serial and MPI variants from clean build directories.
- Run all CTest cases and the parameter-document generation check.
- Run MPI tests outside restricted sandboxes when sockets or process visibility
  matter.
- Update both `docs/parameters.yaml` and
  `docs/advanced/input_files/input-main.md` for INPUT behavior changes.
- Run `git diff --check` and the governance checker, and report exact commands
  and results.
