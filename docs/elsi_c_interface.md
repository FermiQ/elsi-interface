# Overview

This module provides a C language interface to the ELSI (Electronic Structure Infrastructure) library. It allows C or C++ programs to call ELSI's Fortran subroutines by using the `ISO_C_BINDING` feature for interoperability. The primary mechanism involves passing a C pointer (`type(c_ptr)`) as a handle (`h_c`) which corresponds to an internal Fortran handle (`elsi_handle` or `elsi_rw_handle`).

# Key Components

The module defines a large number of C wrapper functions, typically prefixed with `c_elsi_`. These can be broadly categorized as:

- **Initialization and Finalization:**
    - `c_elsi_init`: Initializes an ELSI handle.
    - `c_elsi_finalize`: Finalizes an ELSI handle and releases resources.
    - `c_elsi_reinit`: Reinitializes an ELSI handle.
- **Parameter Setting:**
    - Functions to set MPI communicators, BLACS context, matrix distribution and format (CSC, COO).
        - Example: `c_elsi_set_mpi`, `c_elsi_set_blacs`, `c_elsi_set_csc`
    - Functions to set problem-specific parameters like number of spins, k-points, electron count, basis size.
        - Example: `c_elsi_set_spin`, `c_elsi_set_kpoint`
    - Functions to configure solver-specific settings for various eigensolvers and density matrix solvers (ELPA, OMM, PEXSI, EigenExa, SIPS, NTPoly, MAGMA).
        - Example: `c_elsi_set_elpa_solver`, `c_elsi_set_pexsi_n_pole`, `c_elsi_set_ntpoly_method`
    - Functions to set general options like output verbosity, file names, tolerances.
        - Example: `c_elsi_set_input_file`, `c_elsi_set_output_log`, `c_elsi_set_mu_tol`
- **Computational Routines:**
    - Functions to perform eigenvalue calculations for real and complex, dense and sparse matrices.
        - Example: `c_elsi_ev_real`, `c_elsi_ev_complex_sparse`
    - Functions to calculate the density matrix for real and complex, dense and sparse matrices.
        - Example: `c_elsi_dm_real`, `c_elsi_dm_complex_sparse`
    - Functions to solve Bethe-Salpeter Equation (BSE) problems.
        - Example: `c_elsi_bse_real`, `c_elsi_bse_complex`
- **Data Retrieval:**
    - Functions to get version information, timestamps.
        - Example: `c_elsi_get_version`, `c_elsi_get_datestamp`
    - Functions to retrieve results like eigenvalues, eigenvectors, occupation numbers, chemical potential, energy, entropy.
        - Example: `c_elsi_get_eval`, `c_elsi_get_evec_real`, `c_elsi_get_mu`, `c_elsi_get_entropy`
    - Functions to get error status or specific solver metrics.
        - Example: `c_elsi_get_n_illcond`
- **Matrix Operations:**
    - Functions for orthonormalizing eigenvectors.
        - Example: `c_elsi_orthonormalize_ev_real`
    - Functions for extrapolating density matrices.
        - Example: `c_elsi_extrapolate_dm_complex`
- **Read/Write Functionality:**
    - `c_elsi_init_rw`: Initializes a handle for read/write operations.
    - `c_elsi_finalize_rw`: Finalizes the read/write handle.
    - Functions to read matrix dimensions and data from files for dense and sparse matrices.
        - Example: `c_elsi_read_mat_dim`, `c_elsi_read_mat_real_sparse`
    - Functions to write matrix data to files for dense and sparse matrices.
        - Example: `c_elsi_write_mat_complex`, `c_elsi_write_mat_real_sparse`
- **Utility Functions:**
    - `str_c2f`: Converts C-style strings (null-terminated) to Fortran-style strings. This is used internally when C functions pass string arguments to Fortran.

# Important Variables/Constants

- `h_c`: A `type(c_ptr)` which serves as an opaque handle for C programs to interact with the ELSI library. It is obtained from `c_elsi_init` (or `c_elsi_init_rw`) and passed to most other `c_elsi_...` functions.
- Various integer and real type arguments for specifying parameters, matrix dimensions, tolerances, and for passing data arrays (e.g., `ham_c`, `ovlp_c`, `eval_c`, `evec_c`). These are typically passed as C pointers (`type(c_ptr)`) for arrays/matrices, which are then converted to Fortran pointers within the interface.

# Usage Examples

A typical workflow from C would involve:

```c
// Conceptual C example
#include "elsi.h" // Assuming a C header file for ELSI C interface exists

elsi_handle_c h; // This would be a C-compatible type for the handle (e.g., void*)
int solver = ELSI_SOLVER_ELPA; // Example solver choice
int parallel_mode = ELSI_PARALLEL_MODE_MPI; // Example parallel mode
int matrix_format = ELSI_MATRIX_FORMAT_DENSE; // Example matrix format
int n_basis = 100;
double n_electron = 10.0;
int n_state = 10;

// Initialize ELSI
c_elsi_init(&h, solver, parallel_mode, matrix_format, n_basis, n_electron, n_state);

// Set MPI communicator (example, actual value depends on MPI setup)
// c_elsi_set_mpi(h, mpi_comm_c);

// ... set other parameters ...

// Allocate and prepare matrices (ham_c, ovlp_c, eval_c, evec_c)

// Perform eigenvalue calculation (real, dense example)
// c_elsi_ev_real(h, ham_c_ptr, ovlp_c_ptr, eval_c_ptr, evec_c_ptr);

// ... retrieve results ...
// c_elsi_get_eval(h, eval_c_ptr);

// Finalize ELSI
c_elsi_finalize(h);
```

# Dependencies and Interactions

- **ISO_C_BINDING:** This Fortran standard module is crucial for defining the C-compatible interfaces.
- **ELSI Fortran library:** This module acts as a wrapper around the core Fortran routines of the ELSI library. The actual computations and logic reside in the Fortran part of ELSI.
- **C/C++ calling code:** These C interface functions are intended to be called from user applications written in C or C++. The calling code needs to link against the compiled ELSI library.
- **MPI and BLACS (optional):** For parallel execution, the C calling code and ELSI must be properly configured with MPI and potentially ScaLAPACK/BLACS.
