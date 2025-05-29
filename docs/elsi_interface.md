# Overview

The `elsi_interface.f90` file defines the `ELSI` Fortran module. This module serves as the **primary public Application Programming Interface (API)** for the entire ELSI (Electronic Structure Infrastructure) library. It acts as a high-level facade, consolidating and re-exporting a wide range of functionalities from various specialized `elsi_*` sub-modules. By `use`-ing this single `ELSI` module, Fortran programs gain access to all necessary routines and data types for initializing the library, setting parameters, performing calculations (like solving eigenvalue problems or computing density matrices), retrieving results, and managing data I/O, thus providing a unified and user-friendly entry point to the library's capabilities.

# Key Components

- **Module `ELSI`**:
  The sole public module defined in this file. It does not implement extensive logic itself but rather re-exports public entities from other, more specialized `elsi_*` modules.

The re-exported functionalities can be grouped as follows:

- **Core Data Types**:
  - `elsi_handle`: The fundamental opaque data structure (derived type) that encapsulates the entire state of an ELSI session, including parameters, MPI/BLACS setup, solver-specific data, and matrix storage. It is passed to nearly all ELSI routines.
  - `elsi_rw_handle`: A specialized handle for matrix read/write operations.

- **Setup and Finalization Subroutines**:
  - `elsi_init()`: Initializes an ELSI session, sets up basic parameters, and returns an `elsi_handle`.
  - `elsi_set_mpi()` / `elsi_set_mpi_global()`: Configures MPI communicators.
  - `elsi_set_spin()` / `elsi_set_kpoint()`: Sets spin and k-point information.
  - `elsi_set_blacs()`: Configures the BLACS context for distributed matrix operations.
  - `elsi_set_csc()` / `elsi_set_csc_blk()` / `elsi_set_coo()`: Configures sparse matrix storage formats.
  - `elsi_reinit()`: Reinitializes an existing `elsi_handle`.
  - `elsi_finalize()`: Terminates an ELSI session, cleans up resources, and deallocates ELSI-internal data.

- **Parameter Setting Subroutines (`elsi_set_*`)**:
  A comprehensive suite of routines allowing users to customize ELSI's behavior and solver-specific parameters. Examples include:
  - `elsi_set_input_file()`: Reads parameters from an external file.
  - General settings: `elsi_set_output()`, `elsi_set_zero_def()`, `elsi_set_illcond_check()`, `elsi_set_illcond_tol()`.
  - Solver tuning: `elsi_set_elpa_solver()`, `elsi_set_pexsi_n_pole()`, `elsi_set_ntpoly_method()`, `elsi_set_mu_broaden_scheme()`.
  - Physical parameters: `elsi_set_energy_gap()`, `elsi_set_spin_degeneracy()`.

- **Data Retrieval Subroutines (`elsi_get_*`)**:
  Functions to obtain information, results, and status flags from an `elsi_handle`. Examples include:
  - `elsi_get_version()`: Gets ELSI library version.
  - `elsi_get_mu()`: Gets the chemical potential.
  - `elsi_get_eval()`: Gets computed eigenvalues.
  - `elsi_get_evec_real()` / `elsi_get_evec_complex()`: Gets computed eigenvectors.
  - `elsi_get_edm_real()` / `elsi_get_edm_real_sparse()`: Gets the energy-density matrix.
  - `elsi_get_occ()`: Gets occupation numbers.

- **Main Solver Subroutines**:
  Routines to perform the core electronic structure calculations.
  - Eigenvalue problems: `elsi_ev_real()`, `elsi_ev_complex()`, `elsi_ev_real_sparse()`, `elsi_ev_complex_sparse()`.
  - Density matrix calculations: `elsi_dm_real()`, `elsi_dm_complex()`, `elsi_dm_real_sparse()`, `elsi_dm_complex_sparse()`.
  - Bethe-Salpeter Equation: `elsi_bse_real()`, `elsi_bse_complex()`.
  - Utility solvers: `elsi_inverse_cholesky_real()`.

- **Tool Subroutines**:
  Auxiliary routines for common tasks in electronic structure calculations.
  - `elsi_orthonormalize_ev_real()` / `elsi_orthonormalize_ev_real_sparse()`: Orthonormalize eigenvectors.
  - `elsi_extrapolate_dm_real()` / `elsi_extrapolate_dm_real_sparse()`: Extrapolate density matrices.
  - `elsi_compute_dm_real()`: Compute density matrix from eigenvalues/vectors.
  - `elsi_compute_mu_and_occ()`: Compute chemical potential and occupations.
  - `elsi_suggest_blacs_distribution()`: Suggests a BLACS grid layout.

- **Matrix I/O Subroutines**:
  For reading and writing matrices from/to files.
  - `elsi_init_rw()` / `elsi_finalize_rw()`: Initialize/finalize I/O handle.
  - `elsi_read_mat_dim()` / `elsi_read_mat_real()` / `elsi_read_mat_real_sparse()`.
  - `elsi_write_mat_real()` / `elsi_write_mat_real_sparse()`.

- **Deprecated Routines**:
  Some routines from previous versions are included for backward compatibility (e.g., `elsi_get_n_sing`).

# Important Variables/Constants

- **`elsi_handle` (derived type)**: This is the central data structure. An instance of `elsi_handle` (often named `eh` in examples) is passed as an argument to almost all ELSI API routines. It holds all configuration parameters, MPI/BLACS information, solver-specific data, and allocated memory for matrices and vectors.
- **`elsi_rw_handle` (derived type)**: A similar handle, but specifically tailored for matrix read/write operations.
- **Routine Arguments**:
  - Most `elsi_set_*` routines take `elsi_handle` and the value(s) to be set.
  - Solver routines (`elsi_ev_*`, `elsi_dm_*`) take `elsi_handle`, input matrices (Hamiltonian, Overlap), and output arguments for results (eigenvalues, eigenvectors, density matrix).
  - `elsi_get_*` routines take `elsi_handle` and output arguments to store the retrieved data.

# Usage Examples

A typical workflow for using the ELSI library through this Fortran interface involves the following conceptual steps:

1.  **Initialization**:
    ```fortran
    program use_elsi
      use ELSI, only: elsi_handle, elsi_init, elsi_finalize, ELPA_SOLVER ! ELPA_SOLVER from ELSI_CONSTANT via ELSI
      implicit none
      type(elsi_handle) :: eh
      integer :: err_stat
      integer :: n_basis = 100
      real(kind=selected_real_kind(15,300)) :: n_elec_real = 50.0
      integer :: n_states_calc = 50

      ! Define solver, parallel mode, matrix format (example values)
      call elsi_init(eh, solver_choice=ELPA_SOLVER, parallel_mode=0, matrix_format=0, &
       &              n_basis=n_basis, n_electron=n_elec_real, n_state=n_states_calc, &
       &              error_status=err_stat)
      if (err_stat /= 0) stop "ELSI Init failed!"
    ```

2.  **Set Options (Programmatically or via Input File)**:
    ```fortran
      use ELSI, only: elsi_set_input_file, elsi_set_illcond_tol
      ! Option A: Read from file
      call elsi_set_input_file(eh, "elsi_config.in")
      ! Option B: Set specific parameters
      call elsi_set_illcond_tol(eh, 1.0d-10)
    ```

3.  **Set Up Problem Details (MPI, BLACS, Matrices)**:
    ```fortran
      use ELSI, only: elsi_set_mpi, elsi_set_blacs
      real(kind=selected_real_kind(15,300)), allocatable :: hamiltonian(:,:), overlap(:,:)
      ! ... (Code to set up MPI communicator 'my_comm', BLACS context 'blacs_ctxt', block size 'blk_sz') ...
      ! call elsi_set_mpi(eh, my_comm)
      ! call elsi_set_blacs(eh, blacs_ctxt, blk_sz)
      ! ... (Allocate and populate Hamiltonian and Overlap matrices in distributed format) ...
    ```

4.  **Execute Solver**:
    ```fortran
      use ELSI, only: elsi_ev_real
      real(kind=selected_real_kind(15,300)), allocatable :: eigenvalues(:), eigenvectors(:,:)
      ! ... (Allocate result arrays based on dimensions in eh%ph and eh%bh) ...
      ! allocate(eigenvalues(eh%ph%n_states), eigenvectors(eh%bh%n_lrow, eh%bh%n_lcol))
      
      ! call elsi_ev_real(eh, hamiltonian, overlap, eigenvalues, eigenvectors)
    ```

5.  **Retrieve Additional Results (Optional)**:
    ```fortran
      use ELSI, only: elsi_get_mu
      real(kind=selected_real_kind(15,300)) :: fermi_energy
      ! call elsi_get_mu(eh, fermi_energy)
      ! print *, "Calculated Fermi energy: ", fermi_energy
    ```

6.  **Finalize ELSI Session**:
    ```fortran
      call elsi_finalize(eh)
      ! ... (Deallocate Hamiltonian, Overlap, Eigenvalues, Eigenvectors if allocated here) ...
    end program use_elsi
    ```

# Dependencies and Interactions

- The `ELSI` module serves as a high-level **facade** or **aggregator**. It does not implement most functionalities itself but rather `use`s other specialized `elsi_*` modules and re-exports their public entities.
- **Key `use`d modules whose functionalities are exposed via `ELSI`**:
  - `ELSI_DATATYPE`: Provides the definitions for `elsi_handle` and `elsi_rw_handle`.
  - `ELSI_SETUP`: Provides `elsi_init`, `elsi_finalize`, `elsi_reinit`, and initial configuration routines (MPI, BLACS, sparse formats).
  - `ELSI_SET`: Provides the extensive suite of `elsi_set_*` routines for parameter customization.
  - `ELSI_GET`: Provides the `elsi_get_*` routines for data retrieval.
  - `ELSI_INPUT`: Provides `elsi_set_input_file`.
  - `ELSI_SOLVER`: Provides the core computational routines (`elsi_ev_*`, `elsi_dm_*`, `elsi_bse_*`) and some `elsi_compute_*` utilities.
  - `ELSI_GEO`: Provides routines for eigenvector orthonormalization and density matrix extrapolation.
  - `ELSI_RW`: Provides routines for matrix read/write operations.
  - `ELSI_UTIL`: Provides `elsi_suggest_blacs_distribution`. (Other utilities from `ELSI_UTIL` are used internally by the aforementioned modules but not directly re-exported by `ELSI`).
- By design, a Fortran program using ELSI typically only needs to `use ELSI` to access the full public API of the library. This simplifies usage and decouples user code from the internal modular structure of ELSI.
