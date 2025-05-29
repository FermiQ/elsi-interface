# Overview

The `elsi_lapack.f90` file defines the `ELSI_LAPACK` Fortran module. This module provides an interface to standard LAPACK (Linear Algebra PACKage) routines, primarily for solving dense eigenvalue problems on a single node (i.e., for matrices that are not distributed across multiple MPI processes or are treated as local). It handles both standard and generalized eigenvalue problems. For generalized cases ($H \mathbf{x} = \lambda S \mathbf{x}$), it employs Cholesky factorization of the overlap matrix $S$ to transform the problem into a standard form ($H' \mathbf{y} = \lambda \mathbf{y}$). The module also incorporates routines for applying and reversing frozen core approximations. A notable aspect is its hybrid approach for the eigensolution: after LAPACK routines reduce the matrix to tridiagonal form, it calls an ELPA routine (`elsi_elpa_tridiag`) to solve the tridiagonal eigenproblem.

# Key Components

- **Module `ELSI_LAPACK`**: The main module encapsulating LAPACK-based functionalities.

- **Public Solver Interface**:
  - `elsi_solve_lapack` (generic interface for `_real` and `_cmplx` versions):
    Solves dense eigenvalue problems using LAPACK routines for matrices assumed to be stored locally.
    1.  **Overlap Singularity Check (Optional)**: If `ph%ill_check` is true and the problem is generalized, it calls `elsi_check_ovlp_sp` to check for ill-conditioning in the overlap matrix.
    2.  **Generalized to Standard Transformation**: If not a standard eigenvalue problem (`.not. ph%unit_ovlp`):
        - Calls `elsi_factor_ovlp_sp` to compute $S = U^H U$ and then $U^{-1}$ (if $S$ is well-conditioned).
        - Calls `elsi_reduce_evp_sp` to form $H' = (U^{-1})^H H U^{-1}$.
    3.  **Tridiagonalization**: Reduces the (transformed) Hamiltonian $H'$ to tridiagonal form using LAPACK's `dsytrd` (for real matrices) or `zhetrd` (for complex matrices).
    4.  **Tridiagonal Eigensolver**: Calls `elsi_elpa_tridiag` (from the `ELSI_ELPA` module) to compute eigenvalues and eigenvectors of the tridiagonal matrix.
    5.  **Back-transformation of Eigenvectors**: Transforms eigenvectors from the tridiagonal basis back to the full basis of $H'$ using LAPACK's `dormtr` (real) or `zunmtr` (complex).
    6.  **Back-transformation for Generalized Problem**: If originally a generalized problem, calls `elsi_back_ev_sp` to convert eigenvectors $Z'$ of $H'$ back to eigenvectors $Z$ of the original problem ($Z = U^{-1} Z'$).
    7.  **Frozen Core Adjustment**: If frozen core was applied, adjusts dimensions back to the full problem size.

- **Frozen Core Routines**:
  - `elsi_do_fc_lapack` (generic interface for `_real` and `_cmplx` versions): Applies the frozen core approximation by transforming the Hamiltonian and overlap matrices to the smaller valence subspace.
  - `elsi_undo_fc_lapack` (generic interface for `_real` and `_cmplx` versions): Transforms the eigenvectors obtained from the valence-only problem back to the full basis and reconstructs the eigenvalues of the frozen core states.

- **Private Helper Subroutines (with `_real` and `_cmplx` versions)**:
  - `elsi_factor_ovlp_sp`: Computes $S = U^H U$ using `dpotrf`/`zpotrf`, then $U^{-1}$ using `dtrtri`/`ztrtri`.
  - `elsi_reduce_evp_sp`: Computes $H' = (U^{-1})^H H U^{-1}$ using `dgemm`/`zgemm`.
  - `elsi_back_ev_sp`: Computes $Z = U^{-1} Z'$ using `dtrmm`/`ztrmm` (or `dgemm`/`zgemm` if ill-conditioned).
  - `elsi_check_ovlp_sp`: Checks overlap matrix singularity by computing its eigenvalues (using `dsytrd`/`zhetrd` then `elsi_elpa_tridiag`). If singular, $S$ is overwritten by its scaled eigenvectors.

# Important Variables/Constants

- **`ph` (type `elsi_param_t`) / `bh` (type `elsi_basic_t`)**: ELSI parameter and basic information handles, respectively. These contain control flags and problem dimensions.
- **Input/Output Matrices (`ham`, `ovlp`, `eval`, `evec`)**: These are treated as full, non-distributed matrices (e.g., dimension `ph%n_basis` x `ph%n_basis`).
- **`ph%n_basis`**: Global dimension of the matrices.
- **`ph%n_good`**: Number of basis functions considered well-conditioned.
- **`ph%n_states_solve`**: Number of eigenpairs to compute.
- **`ph%ill_ovlp` (logical)**: True if the overlap matrix is found to be ill-conditioned.
- **`ph%ill_check` (logical)**: If true, triggers the overlap singularity check.
- **`ph%unit_ovlp` (logical)**: True if the overlap matrix $S$ is identity (standard EVP).
- **`ph%fc_method`, `ph%fc_perm`**: Control parameters for frozen core calculations.
- **`perm` (integer array)**: Permutation vector for frozen core transformations.
- **LAPACK routine names**: `dpotrf`, `zpotrf`, `dtrtri`, `ztrtri`, `dgemm`, `zgemm`, `dtrmm`, `ztrmm`, `dsytrd`, `zhetrd`, `dormtr`, `zunmtr`, `dsyrk`, `zherk`. These are standard LAPACK routines for various linear algebra operations.

# Usage Examples

The routines in `ELSI_LAPACK` are primarily intended for internal use by the ELSI library, particularly when ELSI is configured to run in a serial mode or for problems small enough to be handled efficiently by direct LAPACK calls on a single node.

Conceptual workflow if `elsi_solve_lapack_real` is invoked by ELSI:
```fortran
! Assume 'elsi_handle_instance' is an initialized ELSI handle.
! Hamiltonian_local, Overlap_local, Eigenvalues_local, Eigenvectors_local are
! local Fortran arrays holding the full matrices/vectors.

! This call would be made internally by ELSI if LAPACK path is chosen.
call elsi_solve_lapack_real(elsi_handle_instance%ph, elsi_handle_instance%bh, &
                            Hamiltonian_local, Overlap_local, &
                            Eigenvalues_local, Eigenvectors_local)

! The Eigenvalues_local and Eigenvectors_local arrays are populated upon return.
```
The internal sequence involves transformations (if generalized EVP), reduction to tridiagonal form (LAPACK), solving the tridiagonal problem (`elsi_elpa_tridiag`), and back-transformations (LAPACK).

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Uses constants like `FC_BASIC`, `FC_PLUS_V` for frozen core methods.
- **`ELSI_DATATYPE`**: Relies on `elsi_param_t` and `elsi_basic_t` for parameters and control flags.
- **`ELSI_ELPA`**: Critically depends on `elsi_elpa_tridiag` for solving the eigenvalue problem of the tridiagonal matrix generated by LAPACK's reduction step. This represents a hybrid solution strategy.
- **`ELSI_MALLOC`**: Used for allocating temporary workspace arrays needed by LAPACK routines (e.g., `offd`, `tau`, `tmp`).
- **`ELSI_OUTPUT`**: For logging messages (`elsi_say`) and timing (`elsi_get_time`).
- **`ELSI_PRECISION`**: Defines numerical precision kinds (`r8`, `i4`).
- **`ELSI_SORT`**: Uses `elsi_heapsort` and `elsi_permute` for ordering eigenvalues/eigenvectors in the context of frozen core undo operations.
- **LAPACK (External Library)**: This is the core dependency. The module makes numerous direct calls to standard LAPACK routines (e.g., `dpotrf`, `zhetrd`, `dormtr`). A functional LAPACK library must be linked with ELSI.
- **BLAS (External Library)**: LAPACK routines themselves depend on an underlying BLAS (Basic Linear Algebra Subprograms) library for performing vector and matrix operations.
- **Scope**: These routines are designed for single-node, shared-memory execution as they use standard LAPACK calls that operate on matrices fully resident in memory, not distributed across MPI processes using ScaLAPACK.

The file `docs/elsi_lapack.md` has been successfully created.
The documentation includes:
-   An **Overview** of `elsi_lapack.f90`, explaining its role in providing LAPACK-based eigensolvers for non-distributed matrices, including its hybrid use of `elsi_elpa_tridiag`.
-   **Key Components**, listing the `ELSI_LAPACK` module, its main public interfaces (`elsi_solve_lapack`, `elsi_do_fc_lapack`, `elsi_undo_fc_lapack`), and the private helper subroutines for specific LAPACK operations and transformations.
-   **Important Variables/Constants**, highlighting key parameters from ELSI's data structures (`ph`, `bh`) and the names of standard LAPACK routines used.
-   **Usage Examples**, providing a conceptual Fortran example of how ELSI might internally invoke these LAPACK-based solvers.
-   **Dependencies and Interactions**, detailing the module's reliance on other ELSI modules (especially `ELSI_ELPA` for the tridiagonal solve) and the critical external dependencies on LAPACK and BLAS libraries.

All steps of the subtask have been completed.
