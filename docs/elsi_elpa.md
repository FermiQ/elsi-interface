# Overview

The `elsi_elpa.f90` file defines the `ELSI_ELPA` Fortran module, which serves as a comprehensive interface between the ELSI library and the ELPA (Eigenvalue SoLvers for Petaflop-Applications) library, specifically its ELPA-AEO (Algorithmic Extensions and Optimizations) version. This module orchestrates various operations using ELPA, including solving standard and generalized eigenvalue problems, Cholesky factorization and matrix inversion, transformation of generalized eigenproblems to standard form, back-transformation of eigenvectors, handling of frozen core approximations, and density matrix extrapolation. It supports both real and complex arithmetic, offers GPU acceleration options, and can utilize ELPA's autotuning features.

# Key Components

- **Module `ELSI_ELPA`**: The primary module housing all ELPA interfacing logic.

- **Initialization and Cleanup**:
  - `elsi_init_elpa(ph, bh)`: Initializes the ELPA library. This includes calling `elpa_init`, setting up MPI communicators (`ph%elpa_comm_row`, `ph%elpa_comm_col`) for ELPA's 2D processor grid, and creating and configuring ELPA solver instances (`ph%elpa_aux` for auxiliary operations, `ph%elpa_solve` for main eigenvalue solutions) via the internal `elsi_elpa_setup` routine.
  - `elsi_cleanup_elpa(ph)`: Finalizes ELPA operations by deallocating ELPA solver instances and freeing the associated MPI communicators.

- **Main Solver Interface**:
  - `elsi_solve_elpa` (generic interface for `_real` and `_cmplx` versions):
    Orchestrates the solution of eigenvalue problems.
    1.  Optionally checks the overlap matrix for ill-conditioning using `elsi_check_ovlp_elpa`.
    2.  For generalized problems (`.not. ph%unit_ovlp`):
        - Performs Cholesky factorization and inversion of the overlap matrix using `elsi_factor_ovlp_elpa` (if it's the first call and the matrix is not ill-conditioned).
        - Transforms the generalized problem to a standard one using `elsi_reduce_evp_elpa`.
    3.  Solves the standard eigenvalue problem using the internal `elsi_elpa_evec` routine (which calls the actual ELPA eigensolver).
    4.  For generalized problems, back-transforms the eigenvectors using `elsi_back_ev_elpa`.
    5.  Adjusts matrix dimensions and related parameters if a frozen core approximation was applied.

- **Core Operations (Interfaces with Real/Complex Specifics)**:
  - `elsi_factor_ovlp_elpa(ph, bh, ovlp)` (also `elsi_cholesky_inverse_inplace_elpa`): Computes the Cholesky factor $U$ of the overlap matrix `ovlp` (so $S = U^T U$ or $U U^H$) and then replaces `ovlp` with $U^{-1}$.
  - `elsi_reduce_evp_elpa(ph, bh, ham, ovlp, evec)`: Transforms a generalized eigenproblem $(H, S)$ into a standard form $H' = U^{-T} H U^{-1}$. The transformed $H'$ is stored in `ham`, and $U^{-1}$ (from `ovlp`) is used in the process, with `evec` potentially used as workspace.
  - `elsi_back_ev_elpa(ph, bh, ham, ovlp, evec)`: Transforms eigenvectors $Z'$ computed from the standard problem $H'$ back to the basis of the original generalized problem, $Z = U^{-1} Z'$.
  - `elsi_check_ovlp_elpa(ph, bh, ovlp, eval, evec)`: Computes all eigenvalues of the overlap matrix `ovlp` to check for singularity. If singular, `ovlp` is overwritten with scaled eigenvectors, and `ph%ill_ovlp` is set.
  - `elsi_update_dm_elpa(ph, bh, ovlp0, ovlp1, dm0, dm1)`: Extrapolates a density matrix `dm0` (consistent with `ovlp0`) to `dm1` (consistent with `ovlp1`) using the formula $D_1 = U_1^{-T} U_0 D_0 U_0^T U_1^{-1}$, where $U_0, U_1$ are Cholesky factors of `ovlp0`, `ovlp1`.
  - `elsi_do_fc_elpa(ph, bh, ham, ovlp, evec, perm, ham_v, ovlp_v, evec_v)`: Applies frozen core approximation, transforming full Hamiltonian and overlap matrices to their valence-only counterparts (`ham_v`, `ovlp_v`). Updates ELSI parameters to reflect the smaller valence problem size.
  - `elsi_undo_fc_elpa(ph, bh, ham, ovlp, evec, perm, eval_c, evec_v)`: Reconstructs the full eigenvectors and core state eigenvalues after solving the valence-only frozen core problem.

- **Low-Level Routines**:
  - `elsi_elpa_setup(ph, bh, is_aux)`: (Private) Configures ELPA solver handles (`elpa_t` objects) with parameters like matrix size, block size, MPI settings, solver type (1-stage/2-stage), and GPU options.
  - `elsi_elpa_evec(ph, bh, mat, eval, evec, sing_check)`: (Private, generic) The core routine that calls the ELPA `eigenvectors` method. It handles single-precision execution for initial iterations (if `ph%elpa_n_single` > 0), ELPA autotuning, and the singularity check logic.
  - `elsi_elpa_tridiag(ph, bh, diag, offd, evec)`: Solves a tridiagonal eigenvalue problem using ELPA, intended for serial contexts (uses `MPI_COMM_SELF`).

# Important Variables/Constants

- **`ph%elpa_started` (logical)**: Flag indicating if `elsi_init_elpa` has been successfully called.
- **`ph%elpa_first` (logical)**: Flag for operations done only on the first call per geometry (e.g., overlap factorization).
- **`ph%elpa_aux`, `ph%elpa_solve`, `ph%elpa_tune` (pointers to `elpa_t`, `elpa_autotune_t`)**: Opaque handles for ELPA library functionalities.
- **`ph%elpa_comm_row`, `ph%elpa_comm_col` (integer)**: MPI communicators for ELPA's 2D process decomposition.
- **`ph%elpa_solver` (integer)**: Specifies ELPA solver kernel (1-stage vs. 2-stage).
- **`ph%elpa_gpu` (integer)**: Enables (1) or disables (0) GPU acceleration.
- **`ph%elpa_n_single` (integer)**: Number of initial solver calls to run in single precision (often for performance with GPUs or 2-stage CPU solvers).
- **`ph%elpa_autotune` (integer)**: Controls the level of ELPA's autotuning feature.
- **`ph%ill_ovlp` (logical)**: Set to true if the overlap matrix is detected as ill-conditioned.
- **`ph%ill_check` (logical)**: If true, an explicit check for overlap matrix singularity is performed.
- **`ph%ill_tol` (real)**: Threshold for determining singularity.
- **`ph%n_good` (integer)**: Number of basis functions remaining after ill-conditioning treatment.
- **ELPA Library Constants**: e.g., `ELPA_2STAGE_REAL_GPU`, `ELPA_AUTOTUNE_DOMAIN_REAL` used in `elsi_elpa_setup`.
- **ScaLAPACK routines**: Various `p*` routines (e.g., `pdgemm`, `pzherk`) are used for distributed matrix operations.

# Usage Examples

The routines in `ELSI_ELPA` are primarily called internally by the ELSI infrastructure when ELPA is selected as the solver.
A conceptual sequence for solving a generalized eigenvalue problem $H Z = S Z E$:

```fortran
! Assume 'handle' is an initialized ELSI handle, and ELPA is the chosen solver.
! H_matrix, S_matrix are the Hamiltonian and Overlap matrices.
! Eigenvalues_vec, Eigenvectors_mat are for storing results.

! 1. Initialize ELPA (if not already done for this ELSI session)
call elsi_init_elpa(handle%ph, handle%bh)

! 2. ELSI calls the main solver routine
call elsi_solve_elpa(handle%ph, handle%bh, H_matrix, S_matrix, Eigenvalues_vec, Eigenvectors_mat)

!    Internal steps within elsi_solve_elpa (simplified):
!    a. Optional: call elsi_check_ovlp_elpa(handle%ph, handle%bh, S_matrix, ...)
!       to check S matrix singularity. S_matrix might be modified.
!    b. If problem is generalized (S is not identity):
!       If first call and S is well-conditioned:
!          call elsi_factor_ovlp_elpa(handle%ph, handle%bh, S_matrix) ! S_matrix becomes S_inv_chol
!       call elsi_reduce_evp_elpa(handle%ph, handle%bh, H_matrix, S_matrix, Eigenvectors_mat)
!       ! H_matrix becomes H_prime (standard form), Eigenvectors_mat used as temp workspace.
!    c. Call elsi_elpa_evec(handle%ph, handle%bh, H_matrix, Eigenvalues_vec, Eigenvectors_mat, .false.)
!       ! Solves H_prime Z' = Z' E. Eigenvalues_vec and Eigenvectors_mat (Z') are populated.
!    d. If problem was generalized:
!       call elsi_back_ev_elpa(handle%ph, handle%bh, H_matrix_dummy, S_matrix, Eigenvectors_mat)
!       ! Eigenvectors_mat (Z') transformed to original basis (Z).

! 3. Results are in Eigenvalues_vec and Eigenvectors_mat.

! 4. At the end of ELSI computations:
call elsi_cleanup_elpa(handle%ph)
```

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Uses constants like `LT_MAT`, `UT_MAT`, `UNSET`, `FC_BASIC`, `FC_PLUS_V`.
- **`ELSI_DATATYPE`**: Interacts extensively with `elsi_param_t` (`ph`) and `elsi_basic_t` (`bh`) for parameters and matrix/MPI information.
- **`ELSI_MALLOC`**: For dynamic memory allocation (`elsi_allocate`, `elsi_deallocate`).
- **`ELSI_MPI` (and underlying MPI library)**: Crucial for setting up ELPA's distributed environment (`MPI_Comm_split`) and for collective operations (`MPI_Allreduce`). ELPA itself is MPI-parallel.
- **`ELSI_OUTPUT`**: For logging messages (`elsi_say`) and timing sections of code (`elsi_get_time`).
- **`ELSI_PRECISION`**: Defines Fortran kinds for real and integer types (`r4`, `r8`, `i4`, `i8`).
- **`ELSI_SORT`**: Uses `elsi_heapsort` and `elsi_permute` in the context of frozen core calculations.
- **`ELSI_UTIL`**: For error checking (`elsi_check_err`), global index calculation (`elsi_get_gid`), and matrix formatting (`elsi_set_full_mat`).
- **`ELPA` (ELPA library module)**: This is the direct, low-level Fortran interface to the ELPA library. The `ELSI_ELPA` module makes calls to `elpa_init`, `elpa_allocate`, `elpa_deallocate`, `elpa_autotune_deallocate`, and uses ELPA-defined constants. Core ELPA functionalities (e.g., `cholesky`, `eigenvectors`) are accessed through the opaque `elpa_t` objects (`ph%elpa_aux`, `ph%elpa_solve`).
- **ScaLAPACK/PBLAS**: Many distributed linear algebra operations are performed using ScaLAPACK routines (e.g., `pdgemm`, `pzherk`, `pdtran`, `descinit`, `numroc`).
