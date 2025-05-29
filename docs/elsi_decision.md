# Overview

The `elsi_decision.f90` file contains the `ELSI_DECISION` Fortran module. This module is responsible for the automatic selection of an appropriate computational solver when the user has opted for `AUTO_SOLVER`. It includes logic to choose between different eigensolvers or density matrix solvers based on a set of heuristics. These heuristics consider factors like the problem size (number of basis functions), matrix sparsity, physical properties of the system (e.g., energy gap, dimensionality), and the configuration of available parallel resources and libraries (e.g., PEXSI availability and MPI process layout).

# Key Components

- **Module `ELSI_DECISION`**:
  The main module encapsulating the decision-making logic.

- **`elsi_decide_ev(ph, bh)`**:
  A public subroutine that determines the eigensolver if `ph%solver` is set to `AUTO_SOLVER`.
  - Currently, its default behavior is to select the `ELPA_SOLVER`.

- **Interface `elsi_decide_dm`**:
  A public generic interface for selecting the density matrix (DM) solver when `ph%solver` is `AUTO_SOLVER`. It resolves to one of the following specific subroutines based on the input matrix type:
  - **`elsi_decide_dm_real(ph, bh, mat)`**: For real, dense input matrices. It computes matrix sparsity and then delegates to `elsi_decide_dm_smart`.
  - **`elsi_decide_dm_cmplx(ph, bh, mat)`**: For complex, dense input matrices. Similar to the real version, it calculates sparsity before calling `elsi_decide_dm_smart`.
  - **`elsi_decide_dm_sparse(ph, bh)`**: For sparse input matrices where global sparsity information (`bh%nnz_g`) is already available. It directly calls `elsi_decide_dm_smart`.

- **`elsi_decide_dm_smart(ph, bh, sparsity)`**:
  A private subroutine that implements the core decision tree for DM solvers. It evaluates criteria to choose between `ELPA_SOLVER`, `PEXSI_SOLVER`, or `NTPOLY_SOLVER`. If no specialized solver is deemed suitable by the heuristics, it defaults to `ELPA_SOLVER`.

# Important Variables/Constants

The decision logic is influenced by several parameters, primarily from the `elsi_param_t` (`ph`) and `elsi_basic_t` (`bh`) derived types:

- **`ph%solver`**: (Input/Output) If `AUTO_SOLVER`, this module updates it to the ID of the chosen solver (e.g., `ELPA_SOLVER`, `PEXSI_SOLVER`).
- **`ph%n_basis`**: Total number of basis functions. Larger values might favor PEXSI or NTPoly.
- **`sparsity`**: Calculated matrix sparsity (fraction of zero elements). Higher sparsity is a key factor for selecting PEXSI or NTPoly. Thresholds like 0.95, 0.98, 0.99 are used in the logic.
- **`bh%nnz_g`**: Global number of non-zero elements (for sparse input).
- **`bh%def0`**: Threshold used to identify non-zero elements when calculating sparsity for dense matrices.
- **`ph%energy_gap`**: The system's energy gap. Very small gaps might make NTPoly less favorable.
- **`ph%dimensionality`**: Dimensionality of the system (1D, 2D, 3D). PEXSI might be preferred for lower-dimensional systems under certain conditions.
- **`bh%n_procs`**: Total number of MPI processes available, crucial for PEXSI's processor grid compatibility.
- **`ph%pexsi_options%nPoints`**: PEXSI parameter defining the number of points in its algorithm.
- **`ph%pexsi_np_per_pole`**: PEXSI parameter for processors per pole.
- **`pexsi_enabled`**: A flag (obtained via `elsi_get_pexsi_enabled`) indicating if PEXSI support is compiled in ELSI.
- **Constants from `ELSI_CONSTANT` module**: `AUTO_SOLVER`, `ELPA_SOLVER`, `PEXSI_SOLVER`, `NTPOLY_SOLVER`, `UNSET` are used extensively.

# Usage Examples

The subroutines in this module are intended for internal use by the ELSI library. They are invoked when the ELSI handle is configured with `AUTO_SOLVER`.

Conceptual workflow:
```fortran
! Assume 'handle' is an initialized elsi_handle
! User has set handle%ph%solver = AUTO_SOLVER

! For an eigenvalue problem:
! elsi_main_driver might call:
if (handle%ph%solver == AUTO_SOLVER .and. task_is_eigenvalue_problem) then
    call elsi_decide_ev(handle%ph, handle%bh)
    ! Now handle%ph%solver contains the ID of the chosen eigensolver (e.g., ELPA_SOLVER)
end if

! For a density matrix calculation with a real dense Hamiltonian 'h_matrix':
! elsi_main_driver might call:
if (handle%ph%solver == AUTO_SOLVER .and. task_is_density_matrix) then
    call elsi_decide_dm(handle%ph, handle%bh, h_matrix) ! Resolves to elsi_decide_dm_real
    ! Now handle%ph%solver contains the ID of the chosen DM solver (e.g., PEXSI_SOLVER)
end if
```
The `elsi_decide_dm_smart` routine uses heuristics like:
- If PEXSI is enabled AND `n_basis` > 20000 AND `sparsity` > 0.95 AND MPI processors are compatible with PEXSI layout, then PEXSI is a candidate.
- If `n_basis` > 50000 AND `sparsity` > 0.98 AND `energy_gap` > 0.5, then NTPoly is a candidate.
- Specific conditions (e.g., PEXSI for dimensionality < 3, NTPoly for very large and sparse systems) might lead to their selection.
- If no advanced solver (PEXSI, NTPoly) is selected, ELPA is the default.

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Relies heavily on constants like `AUTO_SOLVER`, specific solver IDs (`ELPA_SOLVER`, etc.), and `UNSET`.
- **`ELSI_DATATYPE`**: Functions receive and modify components of `elsi_param_t` (`ph`) and read from `elsi_basic_t` (`bh`).
- **`ELSI_MPI` (and MPI library)**: Uses MPI routines (`MPI_Allreduce`, `MPI_Bcast`) for collective operations like calculating global non-zero counts and broadcasting sparsity values. (Note: The code shows direct MPI calls, implying a dependency on a linked MPI library).
- **`ELSI_OUTPUT`**: Uses `elsi_say` to log the automatically selected solver.
- **`ELSI_PRECISION`**: Uses precision kinds `r8`, `i4`, `i8`.
- **`ELSI_UTIL`**: Depends on `elsi_check_err` for MPI error handling and potentially `elsi_get_pexsi_enabled` (which might reside in `ELSI_UTIL` or a more PEXSI-specific interface module) to determine if PEXSI is available.
- **Solver Modules (e.g., `elsi_elpa`, `elsi_pexsi`, `elsi_ntpoly`)**: The primary interaction is that the `ph%solver` variable, set by this decision module, dictates which of the actual solver execution routines in other ELSI modules will be invoked by the main ELSI control flow.

The file `docs/elsi_decision.md` has been successfully created/overwritten with the detailed documentation. This documentation was derived from analyzing `src/elsi_decision.f90` and using the `documentation_template.md`.

The document includes:
-   An **Overview** of `elsi_decision.f90` and its role in automatic solver selection.
-   **Key Components**, listing the `ELSI_DECISION` module, `elsi_decide_ev`, the `elsi_decide_dm` interface (and its specific procedures `elsi_decide_dm_real`, `elsi_decide_dm_cmplx`, `elsi_decide_dm_sparse`), and the core logic in `elsi_decide_dm_smart`.
-   **Important Variables/Constants**, detailing parameters from `elsi_param_t` and `elsi_basic_t` that influence the decision (e.g., `ph%solver`, `ph%n_basis`, `sparsity`, `ph%energy_gap`, PEXSI configuration) and constants from `ELSI_CONSTANT`.
-   **Usage Examples**, providing a conceptual Fortran workflow illustrating how these decision routines might be invoked internally by ELSI when `AUTO_SOLVER` is specified.
-   **Dependencies and Interactions**, noting dependencies on `ELSI_CONSTANT`, `ELSI_DATATYPE`, MPI, `ELSI_OUTPUT`, `ELSI_PRECISION`, and `ELSI_UTIL`. It also explains how the decisions made impact the subsequent calls to actual solver modules.

All steps of the subtask have been completed successfully.
