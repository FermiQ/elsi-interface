# Overview

The `ELSI_CONSTANT` module centralizes the definition of various integer and real constants that are used throughout the ELSI (Electronic Structure Infrastructure) library. These constants serve as standardized codes and parameters for selecting solvers, specifying matrix properties, controlling parallel execution, defining data types, and managing other operational aspects of ELSI.

# Key Components

This module does not define any functions, classes, or procedures. Its sole purpose is to declare named constants (`parameter` attributes).

# Important Variables/Constants

All constants defined in this module are significant for the configuration and operation of ELSI. They are grouped by their purpose:

- **Mathematical Constants:**
    - `SQRT_PI`: The square root of Pi.
    - `INVERT_SQRT_PI`: The inverse of the square root of Pi.

- **General Purpose Constants:**
    - `UNSET`: A default integer value indicating an unset or undefined parameter.
    - `N_SOLVERS`: Total number of available solver types.
    - `N_MATRIX_FORMATS`: Total number of supported matrix formats.
    - `N_PARALLEL_MODES`: Total number of parallelization modes.

- **Solver Types:** Integer codes to identify different solvers.
    - `AUTO_SOLVER`, `ELPA_SOLVER`, `OMM_SOLVER`, `PEXSI_SOLVER`, `EIGENEXA_SOLVER`, `SIPS_SOLVER`, `NTPOLY_SOLVER`, `MAGMA_SOLVER`, `BSEPACK_SOLVER`.

- **Data Types:** Integer codes for real or complex data.
    - `REAL_DATA`, `CMPLX_DATA`.

- **Matrix Formats:** Integer codes for different matrix storage schemes.
    - `BLACS_DENSE`: Dense matrix format using BLACS distribution.
    - `PEXSI_CSC`: Compressed Sparse Column (CSC) format, compatible with PEXSI.
    - `SIESTA_CSC`: CSC format, compatible with SIESTA.
    - `GENERIC_COO`: Generic Coordinate (COO) sparse matrix format.

- **Matrix Conversion Masks:** Integer codes for specifying which matrices (Hamiltonian H, Overlap S, or both HS) are involved in a conversion.
    - `MASK_HS`, `MASK_H`, `MASK_S`.

- **Triangular Matrix Types:** Integer codes for full, upper triangular, or lower triangular matrices.
    - `FULL_MAT`, `UT_MAT`, `LT_MAT`.

- **Parallelization Modes:** Integer codes for serial or parallel execution.
    - `SINGLE_PROC`, `MULTI_PROC`.

- **Broadening Schemes:** Integer codes for different electronic state broadening methods.
    - `GAUSSIAN`, `FERMI`, `METHFESSEL_PAXTON`, `CUBIC`, `COLD`.

- **Density Matrix Types:** Integer codes to specify the type of density matrix to compute or retrieve.
    - `GET_DM` (Density Matrix), `GET_EDM` (Energy Density Matrix), `GET_FDM` (Free Energy Density Matrix).

- **Frozen Core Methods:** Integer codes for different approaches to handling frozen core electrons.
    - `FC_BASIC`, `FC_PLUS_C`, `FC_PLUS_V`.

- **NTPoly Density Matrix Purification Methods:** Integer codes for purification algorithms in NTPoly.
    - `NTPOLY_PM`, `NTPOLY_TRS2`, `NTPOLY_TRS4`, `NTPOLY_HPCP`.

- **Density Matrix Extrapolation Methods:** Integer codes for density matrix extrapolation techniques.
    - `EXTRA_FACTOR`, `EXTRA_TRS2`.

- **Matrix Reading and Writing Parameters:**
    - `HEADER_SIZE`: Size of the header in matrix files.
    - `FILE_VERSION`: Version number for ELSI's matrix file format.
    - `READ_FILE`: Task identifier for reading a file.
    - `WRITE_FILE`: Task identifier for writing a file.

# Usage Examples

These constants are primarily used internally within ELSI to make decisions and configure operations. For instance, when a user specifies a solver via an input parameter, ELSI routines would compare this input against the defined solver constants:

```fortran
! In an ELSI internal routine
if (user_input_solver == ELPA_SOLVER) then
    call elsi_solve_elpa(...)
else if (user_input_solver == PEXSI_SOLVER) then
    call elsi_solve_pexsi(...)
end if

! When setting matrix properties
call elsi_set_matrix_format(handle, GENERIC_COO)
```

# Dependencies and Interactions

- **`ELSI_PRECISION`:** This module is used to define the precision (kind parameters `r8` for real, `i4` for integer) of the constants.
- **Other ELSI Modules:** The constants defined here are extensively used throughout the ELSI library by other modules to ensure consistent parameter values and to control conditional logic.
