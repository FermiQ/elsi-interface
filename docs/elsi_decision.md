# Overview

The `ELSI_DECISION` module in the ELSI library is responsible for automatically selecting an appropriate eigensolver or density matrix (DM) solver if the user has specified `AUTO_SOLVER`. This selection process relies on a set of heuristics that consider matrix characteristics (size, sparsity), physical system properties (energy gap, dimensionality), and the configuration of available parallel solvers.

# Key Components

- **`elsi_decide_ev(ph, bh)`**:
    - **Purpose**: Determines the eigensolver to be used when `ph%solver` is `AUTO_SOLVER`.
    - **Logic**: Currently, this subroutine defaults to selecting the `ELPA_SOLVER`.
    - **Arguments**:
        - `ph` (inout, type `elsi_param_t`): ELSI parameter bundle, `ph%solver` is updated if it was `AUTO_SOLVER`.
        - `bh` (in, type `elsi_basic_t`): Basic ELSI information.

- **`elsi_decide_dm_real(ph, bh, mat)`**:
    - **Purpose**: Determines the density matrix solver for real, dense input matrices (`mat`).
    - **Logic**: Calculates the sparsity of the input matrix `mat` by counting elements above `bh%def0`. The global non-zero count is obtained via `MPI_Allreduce`. The calculated sparsity is then broadcasted and passed to `elsi_decide_dm_smart`.
    - **Arguments**:
        - `ph` (inout, type `elsi_param_t`): ELSI parameters; `ph%solver` is updated.
        - `bh` (in, type `elsi_basic_t`): Basic ELSI information.
        - `mat` (in, real(r8)): The local part of the dense real matrix.

- **`elsi_decide_dm_cmplx(ph, bh, mat)`**:
    - **Purpose**: Determines the density matrix solver for complex, dense input matrices (`mat`).
    - **Logic**: Similar to `elsi_decide_dm_real`, but for complex matrices. Calculates sparsity and calls `elsi_decide_dm_smart`.
    - **Arguments**:
        - `ph` (inout, type `elsi_param_t`): ELSI parameters; `ph%solver` is updated.
        - `bh` (in, type `elsi_basic_t`): Basic ELSI information.
        - `mat` (in, complex(r8)): The local part of the dense complex matrix.

- **`elsi_decide_dm_sparse(ph, bh)`**:
    - **Purpose**: Determines the density matrix solver when the input matrix is already in a sparse format and its global non-zero count (`bh%nnz_g`) is known.
    - **Logic**: Calculates sparsity based on `bh%nnz_g` and `ph%n_basis`. The sparsity is broadcasted and passed to `elsi_decide_dm_smart`.
    - **Arguments**:
        - `ph` (inout, type `elsi_param_t`): ELSI parameters; `ph%solver` is updated.
        - `bh` (in, type `elsi_basic_t`): Basic ELSI information.

- **`elsi_decide_dm_smart(ph, bh, sparsity)`**:
    - **Purpose**: Contains the core decision-making logic for selecting a density matrix solver (PEXSI, NTPoly, or ELPA as a fallback).
    - **Logic**: It evaluates conditions based on `ph%n_basis` (matrix size), `sparsity`, `ph%energy_gap`, `ph%dimensionality`, PEXSI availability (`elsi_get_pexsi_enabled`), and MPI process counts relative to PEXSI configuration (`ph%pexsi_options%nPoints`, `ph%pexsi_np_per_pole`).
        - PEXSI is considered for large, sparse systems (e.g., `n_basis >= 20000`, `sparsity >= 0.95`) if enabled and MPI configuration is compatible, and typically for lower dimensionality systems.
        - NTPoly is considered for very large, very sparse systems (e.g., `n_basis >= 50000`, `sparsity >= 0.98`) if the energy gap is not too small.
        - If neither PEXSI nor NTPoly meets the criteria, or if `ph%solver` remains `AUTO_SOLVER` after checks, `ELPA_SOLVER` is chosen as the default.
    - **Arguments**:
        - `ph` (inout, type `elsi_param_t`): ELSI parameters; `ph%solver` is updated.
        - `bh` (in, type `elsi_basic_t`): Basic ELSI information.
        - `sparsity` (in, real(r8)): The calculated sparsity of the Hamiltonian matrix.

# Important Variables/Constants

- `ph%solver`: (Input/Output) If `AUTO_SOLVER`, it's updated to `ELPA_SOLVER`, `PEXSI_SOLVER`, or `NTPOLY_SOLVER`.
- `ph%n_basis`: (Input) Global dimension of the problem.
- `bh%def0`: (Input) Threshold to consider a matrix element non-zero when calculating sparsity for dense matrices.
- `bh%nnz_g`: (Input) Global number of non-zero elements, used by `elsi_decide_dm_sparse`.
- `ph%energy_gap`: (Input) Estimated fundamental energy gap of the material.
- `ph%dimensionality`: (Input) Dimensionality of the system (e.g., 1D, 2D, 3D).
- `AUTO_SOLVER`, `ELPA_SOLVER`, `PEXSI_SOLVER`, `NTPOLY_SOLVER`: Constants representing solver choices.

# Usage Examples

These subroutines are intended for internal use by the ELSI library. They are invoked at the beginning of a calculation if the user has requested automatic solver selection.

```fortran
! Inside an ELSI routine, before calling a specific solver
if (handle%ph%solver == AUTO_SOLVER) then
    ! For a density matrix calculation with a real, dense Hamiltonian:
    call elsi_decide_dm_real(handle%ph, handle%bh, handle%ham_real_den)
    ! Now handle%ph%solver will be set to ELPA_SOLVER, PEXSI_SOLVER, or NTPOLY_SOLVER
end if

! Subsequent logic will use the determined handle%ph%solver
select case (handle%ph%solver)
    case (ELPA_SOLVER)
        call elsi_solve_elpa_dm(...)
    case (PEXSI_SOLVER)
        call elsi_solve_pexsi_dm(...)
    ! ... and so on
end select
```

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Provides constants for solver types (`AUTO_SOLVER`, `ELPA_SOLVER`, etc.).
- **`ELSI_DATATYPE`**: Defines the `elsi_param_t` and `elsi_basic_t` derived types.
- **`ELSI_MPI`**: Used for MPI communication primitives like `MPI_Allreduce` and `MPI_Bcast` to gather global information for decision making.
- **`ELSI_OUTPUT`**: Used for logging the automatically selected solver via `elsi_say`.
- **`ELSI_PRECISION`**: Provides definitions for real and integer kinds.
- **`ELSI_UTIL`**: For utility functions like `elsi_check_err`.
- The routine `elsi_get_pexsi_enabled` is called to check if PEXSI is available. (Note: Its module of origin should be ensured via a `use` statement if not part of a larger `use ELSI`).
