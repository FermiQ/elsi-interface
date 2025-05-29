# Overview

This file, `elsi_bsepack.f90`, serves as an interface layer between the ELSI (Electronic Structure Infrastructure) library and the BSEPACK library. BSEPACK is used for solving dense generalized eigenvalue problems that arise from the Bethe-Salpeter Equation (BSE). This interface module provides routines to call BSEPACK for both real and complex matrices, handling workspace allocation and error checking.

# Key Components

- **Module `ELSI_BSEPACK`**:
  The primary Fortran module that encapsulates all BSEPACK interface functionality.

- **Interface `elsi_solve_bsepack`**:
  A public generic interface that directs calls to the appropriate real or complex solver routine based on the type of the input matrices.
  - **`elsi_solve_bsepack_real(ph, bh, mat_a, mat_b, eval, evec)`**:
    A subroutine that interfaces with the BSEPACK routine `pdbseig` to solve the generalized eigenvalue problem for real, distributed matrices `mat_a` and `mat_b`. It computes eigenvalues (`eval`) and eigenvectors (`evec`).
  - **`elsi_solve_bsepack_cmplx(ph, bh, mat_a, mat_b, eval, evec)`**:
    A subroutine that interfaces with the BSEPACK routine `pzbseig` to solve the generalized eigenvalue problem for complex, distributed matrices `mat_a` and `mat_b`. It computes real eigenvalues (`eval`) and complex eigenvectors (`evec`).

# Important Variables/Constants

- **`ph` (type `elsi_param_t`, intent `in`)**:
  An ELSI data structure holding various control parameters for the solver, such as the total number of basis functions (`ph%n_basis`), local and global dimensions for the eigenvector matrix (`ph%bse_n_lrow`, `ph%bse_n_lcol`), and its ScaLAPACK descriptor (`ph%bse_desc`).
- **`bh` (type `elsi_basic_t`, intent `in`)**:
  An ELSI data structure containing basic setup information, including MPI communicator details, local row/column counts for input matrices (`bh%n_lrow`, `bh%n_lcol`), and their ScaLAPACK descriptors (`bh%desc`).
- **`mat_a` (real/complex, intent `inout`)**:
  The distributed matrix A of the generalized eigenvalue problem A*x = lambda*B*x. It is overwritten by BSEPACK.
- **`mat_b` (real/complex, intent `in`)**:
  The distributed matrix B of the generalized eigenvalue problem A*x = lambda*B*x.
- **`eval` (real, intent `out`)**:
  Array where the computed eigenvalues are stored.
- **`evec` (real/complex, intent `out`)**:
  Array (distributed matrix) where the computed eigenvectors are stored.
- **`lwork` (integer)**:
  The size of the primary workspace array (`work`) for BSEPACK routines. Determined by a workspace query.
- **`liwork` (integer)**:
  The size of the integer workspace array (`iwork`) for BSEPACK routines. Determined by a workspace query.
- **`lrwork` (integer, complex version only)**:
  The size of the real workspace array (`rwork`) for the complex BSEPACK routine `pzbseig`. Determined by a workspace query.
- **`ierr` (integer)**:
  Error flag returned by the BSEPACK solver routines. Checked using `elsi_check_err`.
- **`caller` (character string, parameter)**:
  A constant string set to the name of the calling subroutine (e.g., "elsi_solve_bsepack_real"), used for logging and error reporting.

# Usage Examples

The subroutines in this module are typically called internally by the ELSI library when the chosen solver is BSEPACK. Direct user calls are less common but would conceptually follow this pattern:

```fortran
! Assume ph (elsi_param_t) and bh (elsi_basic_t) are initialized
! Assume mat_a and mat_b are distributed matrices set up according to bh%desc
! Assume eval and evec are allocated to appropriate dimensions

! For real matrices:
! call elsi_solve_bsepack_real(ph, bh, mat_a, mat_b, eval, evec)

! For complex matrices:
! call elsi_solve_bsepack_cmplx(ph, bh, mat_a, mat_b, eval, evec)

! After the call, eval will contain the eigenvalues and evec the eigenvectors.
```
The routines first perform a workspace query to BSEPACK (`pdbseig`/`pzbseig` with `lwork = -1`) to determine the required sizes for workspace arrays, allocate these arrays, and then call BSEPACK again to perform the actual computation.

# Dependencies and Interactions

- **ELSI Modules**:
  - `ELSI_DATATYPE`: Uses derived types `elsi_param_t` and `elsi_basic_t`.
  - `ELSI_MALLOC`: Uses `elsi_allocate` and `elsi_deallocate` for dynamic memory management of workspace arrays.
  - `ELSI_OUTPUT`: Uses `elsi_say` for logging messages and `elsi_get_time` for timing the solver execution.
  - `ELSI_PRECISION`: Uses `r8` (double precision real) and `i4` (integer kinds).
  - `ELSI_UTIL`: Uses `elsi_check_err` to check the error code returned by BSEPACK routines.

- **External Libraries**:
  - **BSEPACK**: This module is a direct interface to BSEPACK. It calls:
    - `pdbseig`: For solving real generalized eigenvalue problems.
    - `pzbseig`: For solving complex generalized eigenvalue problems.
  - **MPI (Message Passing Interface)**: While not explicitly calling MPI routines, the use of ScaLAPACK descriptors (`bh%desc`, `ph%bse_desc`) implies that the matrices are distributed and computations are performed in parallel using MPI.
  - **ScaLAPACK/PBLAS**: BSEPACK relies on ScaLAPACK and PBLAS for distributed linear algebra operations. The descriptors used are ScaLAPACK array descriptors.
