# ABACUS H0

This repository is a single-purpose ABACUS variant. It reads the usual
`INPUT`, `STRU`, norm-conserving pseudopotential, and numerical atomic-orbital
files; constructs the initial Hamiltonian in an LCAO basis; writes sparse
real-space `H(R)` in Ry; and exits. It never reads `KPT` and never enters a
solver, diagonalization, occupation, density-matrix, mixing, or SCF loop.

Two Hamiltonians are available:

- `h0_type core`: `T + V_nl + V_loc`
- `h0_type full`: `T + V_nl + V_loc + V_H[rho0] + V_xc[rho0]`

For `full`, `rho0` is initialized with the original ABACUS first-iteration
ordering, including nonlinear core charge and normalization. `init_chg` accepts
`atomic`, `file`, and `auto`; `auto` first tries the normal restart/cube files
and then falls back to atomic density.

## Build

Required dependencies are a C++11 compiler, FFTW3, BLAS, LAPACK, and pthreads.
OpenMP is optional. MPI builds additionally require MPI and ScaLAPACK.

```bash
# Serial
cmake -S . -B build-serial \
  -DENABLE_MPI=OFF -DENABLE_OPENMP=ON -DBUILD_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-serial -j
ctest --test-dir build-serial --output-on-failure

# MPI
cmake -S . -B build-mpi \
  -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DBUILD_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-mpi -j
ctest --test-dir build-mpi --output-on-failure
```

Only one executable is produced: `build-*/abacus`.

## Run

```text
INPUT_PARAMETERS
calculation         get_h0
h0_type             full
basis_type          lcao
suffix              h0
pseudo_dir          ./
orbital_dir         ./
nspin               1
ecutwfc             60
h0_sparse_threshold 1e-10
h0_precision        16
```

Place `STRU` and the files named by it beside `INPUT`, then run:

```bash
OMP_NUM_THREADS=1 /path/to/abacus
```

The result is written below `OUT.<suffix>`:

- `nspin=1`: `hrs1_nao.csr`
- `nspin=2, core`: `hrs1_nao.csr`
- `nspin=2, full`: `hrs1_nao.csr` and `hrs2_nao.csr`
- `nspin=4`: one complex spinor matrix, `hrs1_nao.csr`

See [docs/h0.md](docs/h0.md) for the exact contract and
[docs/advanced/input_files/input-main.md](docs/advanced/input_files/input-main.md)
for generated INPUT metadata.

## License and citation

The original ABACUS license is retained in [LICENSE](LICENSE). See
[CITATIONS.md](CITATIONS.md) for the LCAO and general ABACUS references.
