# Overview

The `elsi_datatype.f90` file defines the `ELSI_DATATYPE` Fortran module. This module's core purpose is to declare and specify several custom derived types. These types are crucial for encapsulating a wide array of data within the ELSI (Electronic Structure Infrastructure) library, including operational parameters, environmental settings (like MPI and BLACS configurations), solver-specific variables and handles, and the actual numerical data for matrices and vectors.

# Key Components

- **Module `ELSI_DATATYPE`**: The primary module containing all derived type definitions.

- **Derived Type `elsi_basic_t`**: This structure groups fundamental information regarding:
  - I/O: Print levels, output units, JSON logging status.
  - MPI Parallelism: Communicators, process IDs, process counts for both the ELSI-specific group and the global scope, and MPI readiness flags.
  - BLACS Grid: BLACS context, ScaLAPACK descriptors, block sizes, processor grid layout, local matrix dimensions, local non-zero counts, and BLACS readiness flag.
  - Common Sparse Matrix Info: Global and local non-zero counts, local column counts for sparse data, and a zero threshold.
  - Format-Specific Sparse Info: Non-zero counts, column counts, and readiness flags for PEXSI CSC, SIESTA CSC, and Generic COO formats.

- **Derived Type `elsi_param_t`**: A comprehensive structure holding parameters that control ELSI's solvers and functionalities. It includes:
  - General Info: Selected solver, matrix format, parallel execution mode, call counters.
  - Overlap Matrix (`S`) Properties: Flags for saving `S`, assuming `S` is an identity matrix, settings for ill-conditioning checks (flag, tolerance), and the eigenvalue range of `S`.
  - Physical System Parameters: Number of electrons, basis functions, spins, k-points, states to compute, local spin/k-point identifiers, k-point weights, spin degeneracy, band structure energy, energy gap, spectrum width, system dimensionality.
  - Chemical Potential (`mu`): Fermi level, entropy, broadening scheme and width, tolerance, and related parameters.
  - Frozen Core Treatment: Method, basis counts for core/valence, and permutation flags.
  - Matrix Redistribution Flags: Boolean flags that track if conversions between different distributed matrix layouts have been performed.
  - Solver-Specific Sections: Dedicated members for ELPA (e.g., `elpa_solver` choice, `elpa_aux` handle), libOMM (e.g., `omm_flavor`, `omm_tol`), PEXSI (e.g., `pexsi_plan` handle, `pexsi_options` structure), EigenExa, SLEPc-SIPs, NTPoly (e.g., `nt_options` structure, `nt_perm` permutation object), MAGMA, and BSEPACK. These often include handles from the underlying libraries and status flags like `_started` or `_first`.

- **Derived Type `elsi_handle`**: This is the main opaque handle for ELSI operations. It bundles:
  - `bh`: An instance of `elsi_basic_t`.
  - `ph`: An instance of `elsi_param_t`.
  - `jh`: A `fjson_handle` from the `FORTJSON` library, used for JSON input/output operations.
  - Allocatable Arrays: A collection of arrays for storing Hamiltonian, Overlap, Eigenvector, and Density matrices. These are provided for both real and complex data, and for dense (`_den`) and various sparse (`_sp`) formats (e.g., `ham_real_den`, `ovlp_cmplx_sp`, `evec_real`).
  - Sparse Matrix Data: Pointers and index arrays for different sparse storage schemes (e.g., `row_ind_sp1`, `col_ptr_sp1`).
  - NTPoly Matrix Types: `nt_ham`, `nt_ovlp`, `nt_dm` (instances of `Matrix_ps` type from NTPoly).
  - Frozen Core Arrays: Matrices and permutation vectors specific to frozen core calculations.
  - Auxiliary Arrays: Storage for copies of matrices, occupation numbers (`occ`), and other temporary or solver-specific data.
  - `handle_init`: A logical flag that tracks whether the ELSI handle has been properly initialized.

- **Derived Type `elsi_rw_handle`**: A specialized handle designed for matrix read and write operations. It contains:
  - `bh`: An instance of `elsi_basic_t`.
  - `rw_task`: Specifies the operation type (read or write).
  - `parallel_mode`, `matrix_format`, `n_electrons`, `n_basis`: Core parameters defining the context for the I/O task.
  - `header_user`: Array for user-defined values in the matrix file header.
  - `handle_init`: A logical flag indicating its initialization state.

# Important Variables/Constants (members of derived types)

The components of these derived types are variables that collectively define the state, configuration, and data for an ELSI calculation. For instance:
- `elsi_handle%ph%solver`: An integer (typically set using named constants from `ELSI_CONSTANT`) that specifies which solver (ELPA, PEXSI, etc.) will be used.
- `elsi_handle%bh%n_lrow`, `elsi_handle%bh%n_lcol`: Store the local row and column dimensions for distributed dense matrices on the current MPI process.
- `elsi_handle%ph%elpa_aux`: A Fortran pointer to an ELPA solver instance (type `elpa_t`).
- `elsi_handle%ham_real_den(:,:)`: An allocatable 2D array that holds the dense real Hamiltonian matrix.
- `elsi_handle%occ(:,:,:)`: A 3D allocatable array for storing occupation numbers, potentially indexed by state, k-point, and spin.

# Usage Examples

These derived types are primarily instantiated and managed by the ELSI library itself. User code, especially Fortran applications directly linking ELSI, would typically interact with an `elsi_handle` variable.

```fortran
! Import necessary ELSI modules
use ELSI_DATATYPE, only: elsi_handle
use ELSI_INIT, only: elsi_init     ! Procedure to initialize ELSI and get a handle
use ELSI_FINALIZE, only: elsi_finalize ! Procedure to clean up and release ELSI resources
use ELSI_SET, only: elsi_set_solver   ! Example procedure to set a parameter
use ELSI_CONSTANT, only: NTPOLY_SOLVER ! Example constant for solver selection

implicit none

type(elsi_handle) :: current_elsi_calculation  ! Declare an ELSI handle
integer :: error_code

! 1. Initialize the ELSI session and obtain a handle.
!    This populates internal structures within 'current_elsi_calculation'.
call elsi_init(current_elsi_calculation, NTPOLY_SOLVER, &
 & parallel_setting, matrix_storage_format, &
 & number_of_basis_functions, number_of_electrons, number_of_states, error_code)

! After this call, current_elsi_calculation%ph%solver would be NTPOLY_SOLVER.
! Default values for many other parameters in current_elsi_calculation%ph and current_elsi_calculation%bh
! would also be set.

! 2. Further configure the calculation by calling ELSI_SET routines.
!    These routines modify members within current_elsi_calculation%ph.
!    For example, to set a tolerance for the NTPoly solver:
!    call elsi_set_ntpoly_tolerance(current_elsi_calculation, 1.0d-7)

! 3. Pass the handle to ELSI solver routines. These routines will:
!    - Access parameters from current_elsi_calculation%ph and current_elsi_calculation%bh.
!    - Use pre-allocated arrays within current_elsi_calculation (e.g., for input matrices).
!    - Allocate and fill output arrays within current_elsi_calculation (e.g., eigenvalues, eigenvectors).
!    Example:
!    call elsi_solve_density_matrix(current_elsi_calculation, hamiltonian_matrix, overlap_matrix, density_matrix, total_energy)
!    This might allocate and populate current_elsi_calculation%dm_real_den.

! 4. Retrieve results using ELSI_GET routines or by directly accessing data if appropriate.

! 5. Finalize the ELSI session to release memory and clean up.
call elsi_finalize(current_elsi_calculation)
! This deallocates arrays within current_elsi_calculation (like ham_real_den, dm_real_den)
! and finalizes any associated library handles (ELPA, PEXSI, etc.).
```

# Dependencies and Interactions

- **`ISO_C_BINDING`**: Utilized for the `c_intptr_t` type, which is a component of `f_ppexsi_options` (from `F_PPEXSI_INTERFACE`) and thus part of `elsi_param_t`.
- **`ELSI_PRECISION`**: This module is crucial as it supplies the kind parameters (`r8` for double precision reals, `i4` for standard integers, `i8` for long integers) used for declaring nearly all numerical and integer members of the defined derived types.
- **External Libraries/Modules (Type Definitions)**:
  - `ELPA`: The `elpa_t` (ELPA solver instance) and `elpa_autotune_t` (autotuning handle) types are incorporated into `elsi_param_t`.
  - `F_PPEXSI_INTERFACE`: The `f_ppexsi_options` derived type is used within `elsi_param_t` to store PEXSI-specific configuration.
  - `FORTJSON`: The `fjson_handle` type is a component of `elsi_handle`, enabling JSON input/output capabilities.
  - `NTPOLY`: Several types from NTPoly (`Permutation_t`, `Matrix_ps`, `SolverParameters_t`, `ProcessGrid_t`) are used in `elsi_param_t` and `elsi_handle` for seamless integration with the NTPoly library.
- **Other ELSI Modules**:
  - The derived types defined in `ELSI_DATATYPE`, particularly `elsi_handle`, are the central data structures that are passed to and manipulated by the vast majority of other ELSI modules (e.g., initialization routines in `ELSI_INIT`, parameter setting routines in `ELSI_SET`, result retrieval routines in `ELSI_GET`, and various solver interface modules).
  - Integer members within these types that denote specific choices or modes (like `solver`, `matrix_format`, `parallel_mode`) are typically assigned values using the named constants defined in the `ELSI_CONSTANT` module.
  - Utility routines, likely found in modules such as `ELSI_UTIL`, are responsible for the proper initialization, allocation of internal allocatable members, and deallocation/cleanup of instances of these derived types, especially `elsi_handle`.
