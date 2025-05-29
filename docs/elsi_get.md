# Overview

The `elsi_get.f90` file defines the `ELSI_GET` Fortran module. This module provides a collection of public subroutines designed to retrieve various types of information from the ELSI (Electronic Structure Infrastructure) library after computations have been performed. Users can access ELSI version details, status indicators (like initialization state or number of ill-conditioned basis functions), calculated physical quantities (such as chemical potential, entropy, eigenvalues, eigenvectors, occupation numbers), and solver-specific results (e.g., bounds on the chemical potential from PEXSI). A key functionality of this module is also to compute and provide the energy-density matrix (EDM) in user-requested formats (dense or sparse), potentially involving internal conversions from different solver-native representations.

# Key Components

- **Module `ELSI_GET`**: The main module containing all data retrieval subroutines.

- **General Information Retrieval**:
  - `elsi_get_initialized(eh, initialized)`: Returns `1` if the ELSI handle `eh` is initialized, `0` otherwise.
  - `elsi_get_version(major, minor, patch)`: Retrieves the major, minor, and patch version numbers of the ELSI library.
  - `elsi_get_datestamp(datestamp)`: Retrieves the compilation datestamp of the ELSI library.

- **Status and Parameter Retrieval**:
  - `elsi_get_n_illcond(eh, n_illcond)`: Gets the number of basis functions identified as ill-conditioned and effectively removed. (Deprecated: `elsi_get_n_sing` is an alias).
  - `elsi_get_ovlp_ev_min(eh, ev_min)`: Returns the minimum eigenvalue of the overlap matrix (if computed, e.g., during an ill-conditioning check).
  - `elsi_get_ovlp_ev_max(eh, ev_max)`: Returns the maximum eigenvalue of the overlap matrix.
  - `elsi_get_pexsi_mu_min(eh, mu_min)`: Gets the lower bound of the chemical potential search range from PEXSI inertia counting.
  - `elsi_get_pexsi_mu_max(eh, mu_max)`: Gets the upper bound of the chemical potential search range from PEXSI.

- **Physical Quantity Retrieval**:
  - `elsi_get_mu(eh, mu)`: Retrieves the calculated chemical potential (Fermi level).
  - `elsi_get_entropy(eh, entropy)`: Retrieves the electronic entropy term ($T \times S$).
  - `elsi_get_eval(eh, eval)`: Retrieves the computed eigenvalues.
  - `elsi_get_evec_real(eh, evec)` / `elsi_get_evec_complex(eh, evec)`: Retrieves the computed eigenvectors in dense real or complex format.
  - `elsi_get_occ(eh, occ)`: Retrieves the computed electronic occupation numbers for the current spin and k-point.

- **Energy-Density Matrix (EDM) Retrieval**:
  These routines compute/retrieve the EDM, $E D = \sum_i \varepsilon_i f_i \psi_i \psi_i^\dagger$.
  - `elsi_get_edm_real(eh, edm)` / `elsi_get_edm_complex(eh, edm)`: Retrieve the EDM in dense real or complex format.
  - `elsi_get_edm_real_sparse(eh, edm)` / `elsi_get_edm_complex_sparse(eh, edm)`: Retrieve the EDM in the sparse format specified by `eh%ph%matrix_format`.
  *The EDM retrieval is sophisticated: it checks `eh%ph%edm_ready`. Depending on the solver that produced the density matrix (`eh%ph%solver`), it may construct the EDM from eigenvalues/vectors (ELPA, EigenExa), call solver-specific EDM computation routines (`elsi_compute_edm_omm`, `elsi_compute_edm_pexsi`, `elsi_compute_edm_ntpoly`, `elsi_build_dm_edm_sips`), and then use routines from `ELSI_REDIST` to convert the EDM to the user-requested output format (dense or specific sparse type) if necessary.*

# Important Variables/Constants

- **`eh` (type `elsi_handle`)**: The primary input to all routines, an ELSI handle containing the state and results of calculations.
- **Output Arguments**: Most routines have output arguments to return the requested data (e.g., `initialized`, `major`, `mu`, `eval`, `edm`).
- **Internal `elsi_handle` members accessed**:
  - `eh%handle_init`: For initialization checks.
  - Version/date info (indirectly via `elsi_version_info`).
  - `eh%ph%n_basis`, `eh%ph%n_good`: For ill-conditioned count.
  - `eh%ph%ovlp_ev_min`, `eh%ph%ovlp_ev_max`: For overlap eigenvalues.
  - `eh%ph%pexsi_options%muMin0`, `eh%ph%pexsi_options%muMax0`: For PEXSI $\mu$ range.
  - `eh%ph%mu`, `eh%ph%ts`: For chemical potential and entropy.
  - `eh%ph%edm_ready`, `eh%ph%eval_ready`, `eh%ph%evec_ready`, `eh%ph%occ_ready`: Status flags indicating data availability. These are often set to `.false.` after retrieval.
  - `eh%ph%solver`, `eh%ph%matrix_format`: Control logic within EDM retrieval for calling appropriate computation and redistribution routines.
  - `eh%eval(:)`, `eh%evec_real(:,:)`, `eh%evec_cmplx(:,:)`, `eh%occ(:,:,:)`: Arrays holding computed eigenvalues, eigenvectors, and occupations.
  - Various DM/EDM storage arrays: `eh%dm_real_den`, `eh%dm_cmplx_sp`, etc.
  - Solver-specific data: `eh%omm_c_real`, `eh%pexsi_ne_vec`, `eh%nt_ham`, `eh%nt_dm`.
  - Sparse matrix descriptors: `eh%row_ind_sp1`, `eh%col_ptr_sp1`, etc.
- **Solver and Format Constants (from `ELSI_CONSTANT`)**: Used internally, e.g., `ELPA_SOLVER`, `PEXSI_CSC`, `GET_EDM`.

# Usage Examples

After an ELSI calculation, these routines are used to extract results.

```fortran
use ELSI_GET
use ELSI_DATATYPE, only: elsi_handle
! ... other necessary ELSI modules ...

type(elsi_handle) :: my_elsi_h
! ... Assume my_elsi_h is initialized and calculations have been performed ...

integer :: maj_ver, min_ver, pat_ver
real(kind=r8) :: chemical_potential, electronic_entropy
real(kind=r8), allocatable :: eigenvalues_out(:)
real(kind=r8), allocatable :: edm_dense_out(:, :)

! Get ELSI version
call elsi_get_version(maj_ver, min_ver, pat_ver)
print *, "Using ELSI version: ", maj_ver, ".", min_ver, ".", pat_ver

! Get chemical potential and entropy
call elsi_get_mu(my_elsi_h, chemical_potential)
call elsi_get_entropy(my_elsi_h, electronic_entropy)

! Get eigenvalues
if (my_elsi_h%ph%eval_ready) then
  allocate(eigenvalues_out(my_elsi_h%ph%n_states))
  call elsi_get_eval(my_elsi_h, eigenvalues_out)
  ! Process eigenvalues_out
  deallocate(eigenvalues_out)
end if

! Get real dense energy-density matrix
if (my_elsi_h%ph%edm_ready) then
  allocate(edm_dense_out(my_elsi_h%bh%n_lrow, my_elsi_h%bh%n_lcol))
  call elsi_get_edm_real(my_elsi_h, edm_dense_out)
  ! Process edm_dense_out
  deallocate(edm_dense_out)
end if

! ... similar calls for elsi_get_evec_real, elsi_get_occ, etc. ...
```

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Uses various constants defining solver types, matrix formats, and task identifiers (e.g., `GET_EDM`).
- **`ELSI_DATATYPE`**: All routines operate on the `elsi_handle` derived type and its members.
- **`ELSI_MALLOC`**: `elsi_get_edm_real` and `elsi_get_edm_complex` (for dense output from diagonalization-based solvers) allocate a temporary `factor` array.
- **Solver-Specific Modules (`ELSI_NTPOLY`, `ELSI_OMM`, `ELSI_PEXSI`, `ELSI_SIPS`)**: The `elsi_get_edm_*` routines may call specialized computation functions from these modules (e.g., `elsi_compute_edm_pexsi`) to obtain the EDM if it was originally calculated by that solver or if it's the most direct way to compute it.
- **`ELSI_REDIST`**: This module is heavily utilized by the `elsi_get_edm_*` routines, especially when the requested output format (dense or a particular sparse type) differs from the internal storage format of the EDM after computation by a specific solver. Numerous `elsi_*_to_*_dm` conversion routines can be invoked.
- **`ELSI_UTIL`**:
  - `elsi_check_init`: Called by most routines to ensure the ELSI handle is valid.
  - `elsi_build_dm_edm`: Used by `elsi_get_edm_*` for constructing the EDM when results are from diagonalization-based solvers like ELPA or EigenExa.
  - `elsi_version_info`: (Likely a private routine elsewhere, e.g. `ELSI_MAIN` or `ELSI_UTIL`) provides raw version strings to `elsi_get_version` and `elsi_get_datestamp`.
  - `elsi_stop`: For error handling if data is not ready or a solver/format is unsupported.
- **Data Availability Flags**: Routines like `elsi_get_eval`, `elsi_get_evec_*`, `elsi_get_occ`, and the `elsi_get_edm_*` family rely on flags within `elsi_handle%ph` (e.g., `eval_ready`, `edm_ready`). These flags are typically set to `.false.` by the "get" routine after successful data retrieval, implying a one-time retrieval or that the data in the handle should be considered "consumed" by that get call.
