# Overview

This module provides an interface to the BSEPACK library, which is used to solve Bethe-Salpeter Equation (BSE) eigenproblems. It defines subroutines for both real and complex matrix types.

# Key Components

- `elsi_solve_bsepack_real`: Solves the BSE eigenproblem for real matrices using the `pdbseig` routine from BSEPACK.
- `elsi_solve_bsepack_cmplx`: Solves the BSE eigenproblem for complex matrices using the `pzbseig` routine from BSEPACK.

# Important Variables/Constants

- `ph`: Input parameter of type `elsi_param_t`, containing ELSI parameters.
- `bh`: Input parameter of type `elsi_basic_t`, containing basic ELSI information.
- `mat_a`, `mat_b`: Input matrices for the eigenproblem.
- `eval`: Output array for eigenvalues.
- `evec`: Output array for eigenvectors.

# Usage Examples

The subroutines are called with ELSI parameters, matrices, and output arrays for eigenvalues and eigenvectors.

```fortran
! Example (conceptual)
call elsi_solve_bsepack_real(ph, bh, mat_a, mat_b, eval, evec)
call elsi_solve_bsepack_cmplx(ph, bh, mat_a_cmplx, mat_b_cmplx, eval_cmplx, evec_cmplx)
```

# Dependencies and Interactions

- Depends on ELSI modules: `ELSI_DATATYPE`, `ELSI_MALLOC`, `ELSI_OUTPUT`, `ELSI_PRECISION`, `ELSI_UTIL`.
- Interacts with the BSEPACK library through the `pdbseig` and `pzbseig` routines.
