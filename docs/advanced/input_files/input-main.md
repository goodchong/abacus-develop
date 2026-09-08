# Full List of INPUT Keywords

<!-- This file is auto-generated from parameters.yaml -->
<!-- Do not edit manually - changes will be overwritten -->

<!-- Table of Contents -->
- [Full List of INPUT Keywords](#full-list-of-input-keywords)
  - [H0 workflow](#h0-workflow)
    - [calculation](#calculation)
    - [h0\_type](#h0_type)
    - [h0\_sparse\_threshold](#h0_sparse_threshold)
    - [h0\_precision](#h0_precision)
  - [Files](#files)
    - [suffix](#suffix)
    - [ntype](#ntype)
    - [stru\_file](#stru_file)
    - [pseudo\_dir](#pseudo_dir)
    - [orbital\_dir](#orbital_dir)
    - [read\_file\_dir](#read_file_dir)
  - [Basis](#basis)
    - [basis\_type](#basis_type)
    - [ecutwfc](#ecutwfc)
    - [ecutrho](#ecutrho)
    - [nx](#nx)
    - [ny](#ny)
    - [nz](#nz)
    - [nbspline](#nbspline)
    - [lcao\_ecut](#lcao_ecut)
    - [lcao\_dk](#lcao_dk)
    - [bx](#bx)
    - [by](#by)
    - [bz](#bz)
  - [Parallel](#parallel)
    - [nb2d](#nb2d)
  - [Density and spin](#density-and-spin)
    - [nspin](#nspin)
    - [noncolin](#noncolin)
    - [lspinorb](#lspinorb)
    - [soc\_lambda](#soc_lambda)
    - [init\_chg](#init_chg)
    - [dft\_functional](#dft_functional)
    - [nelec](#nelec)
    - [nupdown](#nupdown)
    - [pseudo\_rcut](#pseudo_rcut)
  - [Structure](#structure)
    - [min\_dist\_coef](#min_dist_coef)

## H0 workflow

### calculation

- **Type**: String
- **Description**: The H0-only executable accepts get_h0 and exits after writing H(R).
- **Default**: get_h0

### h0_type

- **Type**: String
- **Description**: Required. core constructs T + V_nl + V_loc. full additionally constructs V_H[rho0] + V_xc[rho0] from the initial density, without SCF.
- **Default**: Required

### h0_sparse_threshold

- **Type**: Real
- **Description**: Write a matrix entry only when abs(value) is strictly greater than this threshold.
- **Default**: 1e-10

### h0_precision

- **Type**: Integer
- **Description**: Number of significant digits written for matrix values; valid range is 1 through 17.
- **Default**: 16

[back to top](#full-list-of-input-keywords)

## Files

### suffix

- **Type**: String
- **Description**: Write results beneath OUT.&lt;suffix&gt;.
- **Default**: ABACUS

### ntype

- **Type**: Integer
- **Description**: Number of atomic species. Zero detects the value from STRU.
- **Default**: 0

### stru_file

- **Type**: String
- **Description**: Structure file containing species, pseudopotential and orbital filenames, lattice, and atomic positions.
- **Default**: STRU

### pseudo_dir

- **Type**: String
- **Description**: Directory prepended to pseudopotential filenames in STRU.
- **Default**: ./

### orbital_dir

- **Type**: String
- **Description**: Directory prepended to numerical atomic orbital filenames in STRU.
- **Default**: ./

### read_file_dir

- **Type**: String
- **Description**: Directory searched for ABACUS restart density and spin-resolved cube files when init_chg is file or auto.
- **Default**: OUT.&lt;suffix&gt;

[back to top](#full-list-of-input-keywords)

## Basis

### basis_type

- **Type**: String
- **Description**: Only lcao is supported.
- **Default**: lcao

### ecutwfc

- **Type**: Real
- **Description**: Plane-wave cutoff used for real-space integration grids and LCAO radial tables.
- **Default**: 100
- **Unit**: Ry

### ecutrho

- **Type**: Real
- **Description**: Charge/potential grid cutoff. The LCAO H0 workflow requires ecutrho = 4*ecutwfc unless nx, ny, and nz are supplied explicitly.
- **Default**: 4*ecutwfc
- **Unit**: Ry

### nx

- **Type**: Integer
- **Description**: Explicit FFT grid dimension; nx, ny, and nz must be supplied together.
- **Default**: 0

### ny

- **Type**: Integer
- **Description**: Explicit FFT grid dimension; nx, ny, and nz must be supplied together.
- **Default**: 0

### nz

- **Type**: Integer
- **Description**: Explicit FFT grid dimension; nx, ny, and nz must be supplied together.
- **Default**: 0

### nbspline

- **Type**: Integer
- **Description**: Non-negative values enable cardinal B-spline structure factors; -1 disables them.
- **Default**: -1

### lcao_ecut

- **Type**: Real
- **Description**: Upper reciprocal-space cutoff for LCAO two-center integral tables.
- **Default**: ecutwfc
- **Unit**: Ry

### lcao_dk

- **Type**: Real
- **Description**: Reciprocal-space spacing for two-center integral tables.
- **Default**: 0.01

### bx

- **Type**: Integer
- **Description**: Real-space integration block size; zero selects it automatically.
- **Default**: 0

### by

- **Type**: Integer
- **Description**: Real-space integration block size; zero selects it automatically.
- **Default**: 0

### bz

- **Type**: Integer
- **Description**: Real-space integration block size; zero selects it automatically.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Parallel

### nb2d

- **Type**: Integer
- **Description**: Block size used to distribute H(R) in MPI builds; zero selects it automatically.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Density and spin

### nspin

- **Type**: Integer
- **Description**: 1 is spin-degenerate, 2 is collinear spin, and 4 is a complex two-component spinor.
- **Default**: 1

### noncolin

- **Type**: Boolean
- **Description**: Enable non-collinear magnetism and select nspin=4.
- **Default**: false

### lspinorb

- **Type**: Boolean
- **Description**: Enable spin-orbit coupling, require a compatible fully relativistic pseudopotential, and select nspin=4.
- **Default**: false

### soc_lambda

- **Type**: Real
- **Description**: Scale between scalar-relativistic (0) and full pseudopotential SOC (1).
- **Default**: 1

### init_chg

- **Type**: String
- **Description**: atomic builds rho0 from isolated atoms; file requires existing restart/cube density; auto tries files first and falls back to atomic. Ignored by core.
- **Default**: auto

### dft_functional

- **Type**: String
- **Description**: Use the pseudopotential functional by default. Explicit built-in choices are LDA/PZ, PWLDA, PBE, PBESOL, REVPBE, WC, BLYP, BP, PW91, HCTH, and OLYP. Libxc, meta-GGA, hybrid, HF, and EXX forms are unavailable.
- **Default**: default

### nelec

- **Type**: Real
- **Description**: Override the electron count inferred from pseudopotentials; zero keeps the inferred value.
- **Default**: 0

### nupdown

- **Type**: Real
- **Description**: Constrain the collinear initial spin imbalance used to construct rho0.
- **Default**: 0

### pseudo_rcut

- **Type**: Real
- **Description**: Maximum pseudopotential radial integration range.
- **Default**: 15
- **Unit**: Bohr

[back to top](#full-list-of-input-keywords)

## Structure

### min_dist_coef

- **Type**: Real
- **Description**: Reject structures with distances below this fraction of the reference covalent distance.
- **Default**: 0.2

[back to top](#full-list-of-input-keywords)
