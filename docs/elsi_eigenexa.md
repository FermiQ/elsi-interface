# Overview

The `elsi_eigenexa.f90` file implements the `ELSI_EIGENEXA` Fortran module, which serves as the ELSI library's interface to the EigenExa eigensolver. EigenExa is a library designed for solving large-scale dense eigenvalue problems. This ELSI module manages the lifecycle of the EigenExa solver (initialization, execution, finalization), handles data redistribution between ELSI's standard BLACS (Basic Linear Algebra Communication Subprograms) layout and EigenExa's internal data distribution, and orchestrates the solution process. For generalized eigenvalue problems (Ax = &lambda;Bx), this module leverages routines from ELSI's ELPA interface to transform the problem into a standard form (A'y = &lambda;y) before invoking EigenExa, and then back-transforms the resulting eigenvectors.

# Key Components

- **Module `ELSI_EIGENEXA`**:
  The main Fortran module that encapsulates all interfacing logic with the EigenExa library.

- **`elsi_init_eigenexa(ph, bh)`**:
  A public subroutine responsible for initializing the EigenExa library. It calls `eigen_init` from EigenExa, retrieves EigenExa's processor grid configuration and local matrix dimensions, and stores this information in the `ph` (ELSI parameters) data structure. Marks EigenExa as started by setting `ph%exa_started = .true.`.

- **`elsi_cleanup_eigenexa(ph)`**:
  A public subroutine that finalizes the EigenExa library by calling `eigen_free()`. It also resets EigenExa-related status flags in `ph`.

- **Interface `elsi_solve_eigenexa`**:
  A public generic interface for invoking the EigenExa solver. It resolves to one of the following type-specific subroutines:
  - **`elsi_solve_eigenexa_real(ph, bh, ham, ovlp, eval, evec)`**: Handles real-valued Hamiltonian (`ham`) and overlap (`ovlp`) matrices to compute real eigenvalues (`eval`) and real eigenvectors (`evec`).
  - **`elsi_solve_eigenexa_cmplx(ph, bh, ham, ovlp, eval, evec)`**: Handles complex-valued Hamiltonian and overlap matrices to compute real eigenvalues and complex eigenvectors.

  Both solver routines follow a similar workflow:
  1.  **Generalized to Standard Transformation**: If the problem is generalized (i.e., `ph%unit_ovlp` is false), ELPA routines (`elsi_factor_ovlp_elpa`, `elsi_reduce_evp_elpa`) are used to convert it to a standard eigenvalue problem. The transformed Hamiltonian is stored in `ham`, and transformation vectors (related to the inverse factors of `ovlp`) are temporarily stored in `evec`.
  2.  **Data Redistribution (ELSI to EigenExa)**: The Hamiltonian matrix (now in standard form, if applicable) is redistributed from ELSI's BLACS layout to EigenExa's specific layout using `elsi_blacs_to_eigenexa_h`. Workspace arrays (`ham_exa`, `evec_exa`) are allocated for EigenExa.
  3.  **EigenExa Solver Call**: The appropriate EigenExa routine is called:
      - For real problems: `eigen_s` or `eigen_sx` (depending on `ph%exa_method`).
      - For complex problems: `eigen_h`.
      The `mode` parameter ("A" for all eigenvectors, "N" for eigenvalues only) is determined by `ph%n_states`.
  4.  **Data Redistribution (EigenExa to ELSI)**: Computed eigenvectors are redistributed from EigenExa's layout back to ELSI's BLACS layout using `elsi_eigenexa_to_blacs_ev` into the `evec` array.
  5.  **Workspace Deallocation**: EigenExa work arrays are deallocated.
  6.  **Back Transformation**: If the problem was generalized, eigenvectors are back-transformed using `elsi_back_ev_elpa`.

# Important Variables/Constants

These variables, mostly members of the `ph` (ELSI parameters) and `bh` (ELSI basic info) data structures, are crucial for the EigenExa interface:

- **`ph%exa_started` (logical)**: True if `elsi_init_eigenexa` has been called successfully.
- **`ph%exa_first` (logical)**: True if it's the first call to the solver within a geometry step (used to control one-time factorization of the overlap matrix).
- **`ph%exa_n_prow`, `ph%exa_n_pcol` (integer)**: EigenExa's processor grid dimensions.
- **`ph%exa_my_prow`, `ph%exa_my_pcol` (integer)**: Current process's coordinates in EigenExa's grid.
- **`ph%exa_n_lrow`, `ph%exa_n_lcol` (integer)**: Local matrix dimensions for EigenExa on the current process.
- **`ph%n_basis` (integer)**: The global size of the matrices.
- **`ph%n_states` (integer)**: Number of eigenvalues/eigenvectors to compute.
- **`ph%exa_method` (integer)**: Selects the specific EigenExa solver kernel (e.g., standard vs. expert for real matrices).
- **`ph%exa_blk_fwd`, `ph%exa_blk_bkwd` (integer)**: Block sizes for EigenExa internal operations.
- **`ph%unit_ovlp` (logical)**: If true, the overlap matrix is identity (standard EVP); otherwise, it's a generalized EVP.
- **`bh%comm` (integer)**: MPI communicator passed to EigenExa for its initialization.
- **`UNSET` (integer constant)**: From `ELSI_CONSTANT`, used to check initialization status of certain variables like `bh%nnz_g`.

# Usage Examples

The EigenExa solver is invoked through the ELSI framework. A user would typically select EigenExa as the solver in ELSI's settings. The ELSI library then internally calls the routines from the `ELSI_EIGENEXA` module.

Conceptual internal workflow in ELSI:
```fortran
! Assuming 'elsi_handle' is an initialized ELSI handle structure
! and elsi_handle%ph%solver is set to EIGENEXA_SOLVER.

! ... (ELSI setup, matrix population in elsi_handle arrays) ...

! 1. Initialize EigenExa (typically done once per ELSI session or major setup)
call elsi_init_eigenexa(elsi_handle%ph, elsi_handle%bh)

! 2. Prepare Hamiltonian (H) and Overlap (S) matrices in ELSI's distributed format.
!    Let H_elsi and S_elsi be these matrices.
!    Let E_elsi (vector) and Z_elsi (matrix) be ELSI arrays for eigenvalues/vectors.

! 3. Call the solver via the generic interface
!    For a real problem:
call elsi_solve_eigenexa(elsi_handle%ph, elsi_handle%bh, H_elsi, S_elsi, E_elsi, Z_elsi)
!    The results (eigenvalues and eigenvectors) will be in E_elsi and Z_elsi.

! ... (Post-processing) ...

! 4. Clean up EigenExa resources (typically at the end of the ELSI session)
call elsi_cleanup_eigenexa(elsi_handle%ph)
```

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Uses the `UNSET` constant.
- **`ELSI_DATATYPE`**: Heavily relies on `elsi_param_t` (`ph`) and `elsi_basic_t` (`bh`) derived types for operational parameters, MPI context, matrix dimensions, and solver-specific settings.
- **`ELSI_ELPA`**: Depends on `elsi_factor_ovlp_elpa`, `elsi_reduce_evp_elpa`, and `elsi_back_ev_elpa` for handling generalized eigenvalue problems. This means ELPA functionality is a co-requisite when solving generalized problems with EigenExa through this interface.
- **`ELSI_MALLOC`**: Uses `elsi_allocate` and `elsi_deallocate` for managing memory for workspace arrays required by EigenExa.
- **`ELSI_MPI` (and underlying MPI library)**: Uses `MPI_Allreduce` for global operations (like summing non-zero elements). EigenExa itself is an MPI-parallel library.
- **`ELSI_OUTPUT`**: Utilizes `elsi_say` for logging progress and `elsi_get_time` for performance timing.
- **`ELSI_PRECISION`**: Defines numerical precision kinds (`r8`, `i4`, `i8`).
- **`ELSI_REDIST`**: Essential for data format compatibility, using `elsi_blacs_to_eigenexa_h` to convert Hamiltonian matrices from ELSI's BLACS distribution to EigenExa's layout, and `elsi_eigenexa_to_blacs_ev` to convert eigenvectors back.
- **`ELSI_UTIL`**: Uses `elsi_check_err` for robust error handling of MPI calls.
- **`EIGEN_LIBS_MOD`**: This Fortran module provides the direct low-level interface to the EigenExa C/C++ library. Key EigenExa routines called include `eigen_init`, `eigen_get_procs`, `eigen_get_id`, `eigen_get_matdims`, `eigen_s`, `eigen_sx`, `eigen_h`, and `eigen_free`.
