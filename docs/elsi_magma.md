# Overview

The `elsi_magma.f90` file defines the `ELSI_MAGMA` Fortran module, which serves as an interface to the MAGMA (Matrix Algebra on GPU and Multicore Architectures) library. This module enables the ELSI library to offload dense eigenvalue problem solutions to NVIDIA GPUs, potentially accelerating these computationally intensive steps. It provides routines for initializing and finalizing the MAGMA environment, and for solving both standard ($Ax = \lambda x$) and generalized ($Ax = \lambda Bx$) eigenvalue problems using MAGMA's GPU-accelerated functions. The interface supports both real and complex arithmetic and allows users to choose between MAGMA's one-stage and two-stage eigensolver algorithms.

# Key Components

- **Module `ELSI_MAGMA`**:
  The main Fortran module that encapsulates all interfacing logic with the MAGMA library.

- **`elsi_init_magma(ph)`**:
  A public subroutine responsible for initializing the MAGMA library via `magmaf_init()`. It also queries the number of available GPUs using `magmaf_num_gpus()` and stores this count in `ph%magma_n_gpus`. Sets the `ph%magma_started` flag to true.

- **`elsi_cleanup_magma(ph)`**:
  A public subroutine that finalizes the MAGMA library by calling `magmaf_finalize()`. It also resets the `ph%magma_started` flag.

- **Interface `elsi_solve_magma`**:
  A public generic interface for invoking the MAGMA eigensolver. It resolves to one of the following type-specific subroutines:
  - **`elsi_solve_magma_real(ph, bh, ham, ovlp, eval, evec)`**: Handles real-valued Hamiltonian (`ham`) and, if applicable, overlap (`ovlp`) matrices to compute real eigenvalues (`eval`) and real eigenvectors (`evec`).
  - **`elsi_solve_magma_cmplx(ph, bh, ham, ovlp, eval, evec)`**: Handles complex-valued matrices to compute real eigenvalues and complex eigenvectors.

  Both solver routines implement a two-call approach for MAGMA routines that require workspace:
  1.  **Workspace Query**: The appropriate MAGMA eigensolver is called with workspace size parameters set to -1. This call returns the optimal sizes for `work`, `iwork` (and `rwork` for complex types) arrays.
      - For standard problems (`ph%unit_ovlp = .true.`):
        - `magmaf_dsyevdx_m` (real, 1-stage) or `magmaf_dsyevdx_2stage_m` (real, 2-stage).
        - `magmaf_zheevdx_m` (complex, 1-stage) or `magmaf_zheevdx_2stage_m` (complex, 2-stage).
      - For generalized problems (`ph%unit_ovlp = .false.`):
        - `magmaf_dsygvdx_m` (real, 1-stage) or `magmaf_dsygvdx_2stage_m` (real, 2-stage).
        - `magmaf_zhegvdx_m` (complex, 1-stage) or `magmaf_zhegvdx_2stage_m` (complex, 2-stage).
      The choice between 1-stage and 2-stage algorithms is controlled by `ph%magma_solver`.
  2.  **Memory Allocation**: Workspace arrays are allocated with the sizes determined in the query step.
  3.  **Solver Execution**: The selected MAGMA eigensolver is called again, this time with the correctly sized workspace arrays, to perform the eigenvalue computation. MAGMA routines typically overwrite the input Hamiltonian matrix (`ham`) with the computed eigenvectors.
  4.  **Error Handling**: Checks for errors from MAGMA and ensures the required number of eigenvalues were found.
  5.  **Result Copying**: The eigenvectors, which are stored in `ham` by MAGMA, are copied to the output `evec` array.
  6.  **Workspace Deallocation**: The allocated workspace arrays are freed.

# Important Variables/Constants

- **`ph` (type `elsi_param_t`)**: ELSI parameter data structure.
  - `ph%magma_started` (logical): True if MAGMA has been initialized.
  - `ph%magma_n_gpus` (integer): Number of GPUs available, passed to MAGMA routines.
  - `ph%magma_solver` (integer): Determines MAGMA algorithm (1 for 1-stage, 2 for 2-stage).
  - `ph%unit_ovlp` (logical): If true, solves standard EVP; otherwise, solves generalized EVP.
  - `ph%n_basis` (integer): The global dimension of the input matrices.
  - `ph%n_states` (integer): The number of eigenvalue/eigenvector pairs to compute.
- **`bh` (type `elsi_basic_t`)**: ELSI basic information data structure.
  - `bh%n_lrow`, `bh%n_lcol`: Local matrix dimensions. For these MAGMA wrappers, these are expected to be equal to `ph%n_basis`, implying node-local matrices.
- **Input/Output Arrays**:
  - `ham` (real/complex, intent `inout`): Input Hamiltonian matrix; overwritten by MAGMA with eigenvectors.
  - `ovlp` (real/complex, intent `inout`): Input overlap matrix (for generalized problems).
  - `eval` (real, intent `out`): Output array for eigenvalues.
  - `evec` (real/complex, intent `out`): Output array where eigenvectors (from `ham`) are copied.
- **MAGMA Fortran API routine names**: `magmaf_init`, `magmaf_finalize`, `magmaf_num_gpus`, and various `magmaf_d*` (double real), `magmaf_z*` (double complex) eigensolvers for standard (`evdx`, `evdx_2stage`) and generalized (`gvdx`, `gvdx_2stage`) problems with multiple GPU support (`_m`).

# Usage Examples

The MAGMA solver is invoked through the ELSI framework when specified in ELSI's configuration. These routines are generally not called directly by the end-user but are managed by ELSI's solver dispatch mechanism. The matrices are assumed to be fully present on the host memory of the node where MAGMA is called, and MAGMA internally handles data transfers to and from the GPU(s).

Conceptual internal workflow in ELSI:
```fortran
! Assume 'elsi_handle_inst' is an initialized ELSI handle.
! elsi_handle_inst%ph%solver is set to MAGMA_SOLVER.
! elsi_handle_inst%ph%magma_solver is set to 1 (1-stage) or 2 (2-stage).

! ... (ELSI setup, Hamiltonian H_mat and Overlap S_mat are prepared as full matrices on the host) ...
! ... (Eigenvalue E_val and Eigenvector Z_vec arrays are allocated) ...

! 1. Initialize MAGMA (typically once per ELSI session using MAGMA)
call elsi_init_magma(elsi_handle_inst%ph)
! Now, elsi_handle_inst%ph%magma_n_gpus contains the number of GPUs.

! 2. Call the ELSI MAGMA solver interface (e.g., for a real, generalized problem)
call elsi_solve_magma(elsi_handle_inst%ph, elsi_handle_inst%bh, &
                      H_mat, S_mat, E_val, Z_vec)
! This internally calls elsi_solve_magma_real, which then calls the appropriate
! MAGMA routine (e.g., magmaf_dsygvdx_m or magmaf_dsygvdx_2stage_m).
! Results are stored in E_val and Z_vec.

! ... (Post-processing of results) ...

! 3. Clean up MAGMA resources (at the end of the ELSI session)
call elsi_cleanup_magma(elsi_handle_inst%ph)
```

# Dependencies and Interactions

- **`ELSI_DATATYPE`**: Utilizes `elsi_param_t` (`ph`) and `elsi_basic_t` (`bh`) for parameters, matrix dimensions, and MAGMA-specific settings.
- **`ELSI_MALLOC`**: Employs `elsi_allocate` and `elsi_deallocate` for managing MAGMA's workspace arrays.
- **`ELSI_MPI`**: While MAGMA can manage multi-GPU parallelism which may involve MPI, these specific ELSI wrappers do not make direct MPI calls for the solver part. The `elsi_stop` utility, used for error handling, is MPI-aware.
- **`ELSI_OUTPUT`**: Uses `elsi_say` for logging information and `elsi_get_time` for timing the solver execution. `elsi_stop` is used for error termination.
- **`ELSI_PRECISION`**: Defines numerical precision kinds (`r8` for double precision, `i4` for integers).
- **MAGMA Library (`MAGMA` module)**: This is the core external dependency. The module directly calls Fortran interface routines from the MAGMA library (e.g., `magmaf_init`, `magmaf_finalize`, `magmaf_dsyevdx_m`, `magmaf_zhegvdx_2stage_m`). A functional MAGMA installation, compiled with Fortran support and appropriate GPU drivers, must be linked with ELSI.
- **Matrix Locality**: The MAGMA routines called by this interface typically expect matrices to be fully available in the host memory of the node executing the call, rather than distributed in a ScaLAPACK fashion across multiple nodes. `bh%n_lrow` and `bh%n_lcol` are used as the full matrix dimensions (`ph%n_basis`).

The file `docs/elsi_magma.md` has been successfully created.
The documentation includes:
-   An **Overview** of `elsi_magma.f90`, explaining its role as an interface to the MAGMA library for GPU-accelerated dense eigenvalue solutions.
-   **Key Components**, listing the `ELSI_MAGMA` module, initialization (`elsi_init_magma`), cleanup (`elsi_cleanup_magma`), and the solver interface `elsi_solve_magma` with its real and complex implementations. The two-call (query and execute) MAGMA workflow is detailed.
-   **Important Variables/Constants**, highlighting key members from `elsi_param_t` and `elsi_basic_t` that control the MAGMA interface (e.g., `ph%magma_started`, `ph%magma_solver`, `ph%unit_ovlp`), and the MAGMA Fortran API routine names.
-   **Usage Examples**, providing a conceptual Fortran workflow of how ELSI would internally manage calls to the MAGMA interface routines.
-   **Dependencies and Interactions**, detailing the module's reliance on other ELSI modules (`ELSI_DATATYPE`, `ELSI_MALLOC`, `ELSI_OUTPUT`, `ELSI_PRECISION`) and the crucial dependency on the MAGMA library itself. The expectation of node-local matrices is also noted.

All steps of the subtask have been completed.
