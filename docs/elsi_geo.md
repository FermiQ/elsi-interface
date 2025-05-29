# Overview

The `elsi_geo.f90` file defines the `ELSI_GEO` Fortran module. While "geo" might suggest geometry optimization, this module's primary functions are centered around operations that maintain the consistency and improve the efficiency of electronic structure calculations across successive steps, particularly when the atomic geometry changes or when reusing information from previous calculations. It provides routines for orthonormalizing eigenvectors against an overlap matrix and for extrapolating density matrices from a previous state to a new one, given a new overlap matrix. These procedures are vital for accelerating convergence in Self-Consistent Field (SCF) cycles and for providing good initial guesses in sequences of calculations.

# Key Components

- **Module `ELSI_GEO`**: The main container for the subroutines.

- **Eigenvector Orthonormalization**:
  - `elsi_orthonormalize_ev_real(eh, ovlp, evec)`: Orthonormalizes real eigenvectors (`evec`) with respect to a dense real overlap matrix (`ovlp`).
  - `elsi_orthonormalize_ev_complex(eh, ovlp, evec)`: Orthonormalizes complex eigenvectors with respect to a dense complex overlap matrix.
  - `elsi_orthonormalize_ev_real_sparse(eh, ovlp, evec)`: Orthonormalizes real eigenvectors with respect to a sparse real overlap matrix. The sparse overlap is internally converted to dense format before orthonormalization.
  - `elsi_orthonormalize_ev_complex_sparse(eh, ovlp, evec)`: Orthonormalizes complex eigenvectors with respect to a sparse complex overlap matrix, also via internal conversion to dense.
  *All orthonormalization routines utilize the `elsi_gram_schmidt` procedure from the `ELSI_UTIL` module.*

- **Density Matrix Extrapolation**:
  - `elsi_extrapolate_dm_real(eh, ovlp, dm)`: Extrapolates a previously stored real dense density matrix (from `eh%dm_real_copy`, consistent with `eh%ovlp_real_copy`) to a new real dense density matrix (`dm`) corresponding to a new real dense overlap matrix (`ovlp`). This typically uses `elsi_update_dm_elpa`.
  - `elsi_extrapolate_dm_complex(eh, ovlp, dm)`: Similar to the real version, but for complex matrices.
  - `elsi_extrapolate_dm_restart_real(eh, ovlp_old, ovlp_new, dm)`: Extrapolates a given real dense density matrix (`dm`, consistent with `ovlp_old`) to be consistent with `ovlp_new`. The result overwrites the input `dm`. This is useful when the old DM/overlap are not stored in the ELSI handle (e.g., read from restart files).
  - `elsi_extrapolate_dm_restart_complex(eh, ovlp_old, ovlp_new, dm)`: Complex counterpart to the restart extrapolation for dense matrices.
  - `elsi_extrapolate_dm_real_sparse(eh, ovlp, dm)`: Extrapolates a previously stored sparse real density matrix (internally stored in NTPoly format as `eh%nt_dm_copy`) to a new sparse real density matrix (`dm`) corresponding to a new sparse real overlap matrix (`ovlp`). This involves converting the input sparse overlap to NTPoly format, calling `elsi_update_dm_ntpoly`, and then converting the resulting NTPoly density matrix back to ELSI's sparse format.
  - `elsi_extrapolate_dm_complex_sparse(eh, ovlp, dm)`: Similar to the sparse real version, but for complex matrices.

# Important Variables/Constants

The subroutines in this module primarily operate on data passed as arguments or stored within the `elsi_handle` (`eh`):

- **`eh` (type `elsi_handle`)**: The main ELSI data structure.
  - `eh%bh%n_lrow`, `eh%bh%n_lcol`: Local dimensions for dense matrices.
  - `eh%bh%nnz_l_sp`: Number of local non-zero elements in input sparse matrices.
  - `eh%ovlp_real_den`, `eh%ovlp_cmplx_den`: Allocatable arrays within the handle used as temporary dense storage for overlap matrices, especially when converting from sparse formats for orthonormalization.
  - `eh%ovlp_real_copy`, `eh%dm_real_copy` (and complex versions): Store copies of the overlap and density matrix from a previous step, serving as the "old" state for dense density matrix extrapolation.
  - `eh%nt_ovlp`, `eh%nt_dm`: NTPoly sparse matrix objects (`Matrix_ps`) used during sparse density matrix extrapolation.
  - `eh%nt_ovlp_copy`, `eh%nt_dm_copy`: NTPoly sparse matrix objects storing the "old" overlap and density matrix for sparse extrapolation.
  - `eh%ph%matrix_format`: An integer (defined in `ELSI_CONSTANT`, e.g., `PEXSI_CSC`, `SIESTA_CSC`, `GENERIC_COO`) that dictates the specific sparse matrix format and thus the redistribution routines to be used.
  - `eh%row_ind_sp1`, `eh%col_ptr_sp1`, etc.: Arrays containing the index and pointer data for supported sparse matrix formats.

# Usage Examples

These routines are typically invoked by the host code between electronic structure calculation steps (e.g., SCF iterations or geometry optimization steps).

**Orthonormalizing Eigenvectors**:
```fortran
! Assuming 'my_elsi_handle' is an initialized elsi_handle.
! 'overlap_current_step' is the current overlap matrix S.
! 'eigenvectors_guess' contains eigenvectors from a previous calculation or a guess.

! For dense matrices:
call elsi_orthonormalize_ev_real(my_elsi_handle, overlap_current_step, eigenvectors_guess)
! Now, 'eigenvectors_guess' satisfies Z^T * S * Z = I.

! If 'overlap_current_step_sparse_vals' holds the non-zero values of a sparse S:
! call elsi_orthonormalize_ev_real_sparse(my_elsi_handle, overlap_current_step_sparse_vals, eigenvectors_guess)
```

**Extrapolating Density Matrix**:
```fortran
! Assuming 'my_elsi_handle' is an initialized elsi_handle.
! The previous S and DM are already stored in my_elsi_handle (e.g., ...%ovlp_real_copy).
! 'overlap_new_geometry' is the overlap matrix for the new atomic configuration.
! 'density_matrix_output' will store the extrapolated DM.

! For dense matrices:
call elsi_extrapolate_dm_real(my_elsi_handle, overlap_new_geometry, density_matrix_output)

! If 'overlap_new_geometry_sparse_vals' holds non-zero values of a sparse new S:
! call elsi_extrapolate_dm_real_sparse(my_elsi_handle, overlap_new_geometry_sparse_vals, density_matrix_output_sparse_vals)
```

# Dependencies and Interactions

- **`ELSI_CONSTANT`**: Uses constants like `PEXSI_CSC`, `SIESTA_CSC`, `GENERIC_COO` to identify sparse matrix formats.
- **`ELSI_DATATYPE`**: All routines take an `elsi_handle` as an argument and operate on its members.
- **`ELSI_ELPA`**: The `elsi_update_dm_elpa` subroutine is used for extrapolating dense density matrices.
- **`ELSI_MALLOC`**: Used to allocate temporary dense storage within the ELSI handle if it's not already allocated (e.g., `eh%ovlp_real_den`).
- **`ELSI_NTPOLY`**: The `elsi_update_dm_ntpoly` subroutine is used for extrapolating sparse density matrices.
- **`ELSI_PRECISION`**: Defines the kind `r8` for real variables.
- **`ELSI_REDIST`**: This module is critical for sparse operations. Routines like `elsi_sips_to_blacs_hs`, `elsi_generic_to_ntpoly_hs`, `elsi_ntpoly_to_siesta_dm`, etc., are used to convert matrices between different sparse formats (PEXSI_CSC, SIESTA_CSC, GENERIC_COO) and the dense (BLACS) or NTPoly formats required by the underlying computational kernels.
- **`ELSI_UTIL`**:
  - `elsi_check_init`: Called at the beginning of each public routine to ensure the ELSI handle is initialized.
  - `elsi_gram_schmidt`: Performs the actual Gram-Schmidt orthonormalization for all `elsi_orthonormalize_ev_*` routines.
  - `elsi_stop`: Used for error termination if an unsupported matrix format is encountered.
