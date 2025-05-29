# Overview

The `ELSI_DATATYPE` module is a cornerstone of the ELSI (Electronic Structure Infrastructure) library, defining the essential derived data types (analogous to structures in C/C++) used to manage and store data throughout ELSI operations. These types encapsulate a wide array of information, including MPI/BLACS configurations, physical system parameters, matrix data, solver-specific settings, and operational handles.

# Key Components (Derived Types)

This module's primary role is the definition of the following derived types:

- **`elsi_basic_t`**:
    - **Purpose**: Stores fundamental information about the computational environment and matrix distribution.
    - **Key Members**:
        - I/O settings: `print_info`, `print_unit`, `print_json`.
        - MPI details: `myid`, `n_procs`, `comm` (MPI communicator for the solver group), `comm_all` (global MPI communicator).
        - BLACS grid information: `blacs_ctxt`, `desc` (ScaLAPACK descriptor), `n_prow`, `n_pcol`, `my_prow`, `my_pcol`, `n_lrow`, `n_lcol` (local matrix dimensions).
        - Sparse matrix properties: `nnz_g` (global non-zeros), `nnz_l_sp` (local non-zeros for generic sparse format), `def0` (zero threshold), and format-specific non-zero counts and local column numbers for PEXSI CSC, SIESTA CSC, and Generic COO formats.

- **`elsi_param_t`**:
    - **Purpose**: Contains a comprehensive set of parameters that control the behavior of ELSI solvers and other functionalities.
    - **Key Members**:
        - General settings: `solver` (selected solver ID), `matrix_format`, `parallel_mode`.
        - Overlap matrix handling: `save_ovlp`, `unit_ovlp`, `ill_check` (ill-conditioning check toggle), `ill_tol` (ill-conditioning tolerance), `ovlp_ev_min`/`_max`.
        - Physical system details: `n_electrons`, `n_basis`, `n_spins`, `n_kpts`, `n_states` (number of states to find/use), `spin_degen`, `energy_gap`, `spectrum_width`.
        - Chemical potential (`mu`) calculation: `mu_scheme`, `mu_width`, `mu_tol`.
        - Frozen core options: `fc_method`, `n_basis_c` (core basis size), `n_basis_v` (valence basis size).
        - Matrix redistribution flags and solver-specific parameter subsections for:
            - ELPA: `elpa_solver` (1-stage/2-stage), `elpa_gpu`, `elpa_autotune`.
            - libOMM: `omm_flavor`, `omm_tol`.
            - PEXSI: `pexsi_np_per_pole`, `pexsi_mu_min`/`_max`, `pexsi_options` (native PEXSI options structure).
            - EigenExa: `exa_method`.
            - SLEPc-SIPs: `sips_n_slices`, `sips_interval`.
            - NTPoly: `nt_method`, `nt_tol`, `nt_options` (native NTPoly options structure).
            - MAGMA: `magma_solver`, `magma_n_gpus`.
            - BSEPACK: `bse_n_lrow`, `bse_n_lcol`.

- **`elsi_handle`**:
    - **Purpose**: The primary opaque handle used in ELSI's main interface. It aggregates basic information, parameters, and data arrays.
    - **Key Members**:
        - `bh`: An instance of `elsi_basic_t`.
        - `ph`: An instance of `elsi_param_t`.
        - `jh`: An instance of `fjson_handle` for JSON logging.
        - Allocatable arrays for dense matrices: `ham_real_den`, `ham_cmplx_den`, `ovlp_real_den`, `ovlp_cmplx_den`, `eval` (eigenvalues), `evec_real`, `evec_cmplx`, `dm_real_den`, `dm_cmplx_den`.
        - Allocatable arrays for sparse matrices (CSC values, COO values/indices): `ham_real_sp`, `ovlp_cmplx_sp`, `row_ind_sp1`, `col_ptr_sp1`, etc.
        - Native sparse matrix types for NTPoly: `nt_ham`, `nt_ovlp`, `nt_dm`.
        - Data for frozen core calculations: `perm_fc` (permutation vector), `ham_real_v` (valence Hamiltonian).
        - Auxiliary and temporary storage arrays used by various routines.
        - `handle_init`: A logical flag indicating if the handle has been initialized.

- **`elsi_rw_handle`**:
    - **Purpose**: A specialized handle for matrix read and write operations.
    - **Key Members**:
        - `bh`: An instance of `elsi_basic_t`.
        - `rw_task`: Specifies the operation (read or write).
        - `parallel_mode`, `matrix_format`, `n_electrons`, `n_basis`.
        - `header_user`: User-defined values in the matrix file header.
        - `handle_init`: A logical flag indicating if the read/write handle has been initialized.

# Important Variables/Constants

The members within these derived types are crucial for ELSI's operation. For example:
- `elsi_handle%bh%comm`: The MPI communicator for the current solver.
- `elsi_handle%ph%solver`: The integer code for the chosen eigensolver or density matrix method.
- `elsi_handle%ph%n_electrons`: The number of electrons in the system.
- `elsi_handle%ham_real_den`: Array storing the Hamiltonian matrix (real, dense case).

# Usage Examples

These data types are instantiated and populated by ELSI's setup and interface routines. Users interacting with ELSI via its Fortran API would typically work with a variable of type `elsi_handle`.

```fortran
! Declaration of an ELSI handle
type(elsi_handle) :: my_elsi_calculation

! Initialization (simplified representation of what elsi_init does)
call elsi_init(my_elsi_calculation, ELPA_SOLVER, MULTI_PROC, BLACS_DENSE, &
               n_basis_val, n_electrons_val, n_states_val)

! Setting a specific parameter after initialization
my_elsi_calculation%ph%mu_tol = 1.0d-8

! Accessing information from the handle
if (my_elsi_calculation%bh%myid == 0) then
  write(*,*) "ELSI calculation initialized for solver: ", my_elsi_calculation%ph%solver
end if

! Matrices are allocated and accessed via the handle, e.g.:
! my_elsi_calculation%ham_real_den(i,j) = ...
```

# Dependencies and Interactions

- **`ELSI_PRECISION`**: Provides kind parameters (`r8`, `i4`, `i8`) for defining the precision of numeric members.
- **External Library Modules**:
    - `ELPA`: For `elpa_t` (ELPA handle) and `elpa_autotune_t` (ELPA autotuning handle).
    - `F_PPEXSI_INTERFACE`: For `f_ppexsi_options` (PEXSI options structure).
    - `FORTJSON`: For `fjson_handle` (FortJSON logging handle).
    - `NTPOLY`: For `Permutation_t`, `Matrix_ps` (NTPoly sparse matrix type), `SolverParameters_t`, `ProcessGrid_t`.
- **Other ELSI Modules**: These data types are fundamental and are passed to and modified by most other modules within ELSI, particularly `ELSI_SETUP`, `ELSI_SET`, `ELSI_SOLVER`, and the various solver-specific interface modules.
