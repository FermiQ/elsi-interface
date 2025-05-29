# Overview

The `elsi_input.f90` file defines the `ELSI_INPUT` Fortran module. This module is responsible for reading and processing runtime parameters for the ELSI (Electronic Structure Infrastructure) library from an external text file. It enables users to customize a wide array of ELSI settings and solver-specific parameters without needing to recompile the library. This includes choices for solvers, numerical tolerances, control flags for various algorithms, and detailed tuning parameters for underlying libraries like ELPA, PEXSI, OMM, NTPoly, etc.

# Key Components

- **Module `ELSI_INPUT`**:
  The primary module that encapsulates all functionality related to parsing the input file.

- **`elsi_set_input_file(eh, f_name)`**:
  This is the main public subroutine. It takes an ELSI handle (`eh` of type `elsi_handle`) and the input file name (`f_name` as a character string) as arguments.
  Its operation involves:
  1.  Opening the specified file.
  2.  Reading the file line by line.
  3.  Ignoring comment lines (starting with `#`) and blank lines.
  4.  Parsing each valid line to extract a keyword (case-insensitive) and its corresponding value(s).
  5.  Calling the appropriate `elsi_set_*` subroutine (from the `ELSI_SET` module) to apply the parsed parameter to the `elsi_handle` (`eh`).

- **`elsi_check_read(bh, ierr, kwd)`**:
  A private utility subroutine used to verify the success of internal Fortran `read` statements during the parsing process. If a read error (`ierr /= 0`) occurs, it calls `elsi_stop` to terminate the program with an informative message.

- **`elsi_str_to_int(val_str, val_i4)`**:
  A private utility subroutine that converts string values commonly used for boolean-like options (e.g., "true", ".false.", "yes", "no", "0", "1") from the input file into integer representations (typically `0` for false/no and `1` for true/yes).

# Important Variables/Constants (Keywords in Input File)

The `elsi_set_input_file` routine recognizes a comprehensive list of keywords. Each keyword corresponds to a specific configurable parameter within ELSI. The module reads the keyword and its value, then uses a corresponding `elsi_set_*` routine to update the ELSI handle. Examples of supported keywords include:

- **General Settings**:
  - `output`: Controls general verbosity.
  - `output_log`: Specifies if logging output should be written (e.g., to a JSON file).
  - `zero_def`: Threshold for considering a number as zero.
  - `sparsity_mask`: Mask for matrix operations.
  - `illcond_check`: Enables/disables ill-conditioning check for the overlap matrix.
  - `illcond_tol`: Tolerance for ill-conditioning.

- **Physical and Algorithmic Parameters**:
  - `energy_gap`: User-provided estimate for the HOMO-LUMO gap.
  - `spectrum_width`: User-provided estimate for the spectral width.
  - `dimensionality`: Dimensionality of the system (1D, 2D, 3D).
  - `extrapolation`: Method for density matrix extrapolation.

- **Solver-Specific Parameters (examples)**:
  - **ELPA**: `elpa_solver`, `elpa_n_single`, `elpa_gpu`, `elpa_autotune`.
  - **OMM**: `omm_flavor`, `omm_n_elpa`, `omm_tol`.
  - **PEXSI**: `pexsi_method`, `pexsi_n_mu`, `pexsi_n_pole`, `pexsi_np_per_pole`, `pexsi_np_symbo`, `pexsi_temp`, `pexsi_inertia_tol`.
  - **EigenExa**: `eigenexa_method`.
  - **SLEPc-SIPs**: `sips_n_elpa`, `sips_n_slice`.
  - **NTPoly**: `ntpoly_method`, `ntpoly_tol`, `ntpoly_filter`.
  - **MAGMA**: `magma_solver`.

- **Chemical Potential Calculation**:
  - `mu_broaden_scheme`: Broadening scheme for occupation numbers.
  - `mu_broaden_width`: Width for broadening.
  - `mu_tol`: Tolerance for chemical potential search.
  - `mu_mp_order`: Order for Methfessel-Paxton broadening.

- **Frozen Core**:
  - `n_frozen`: Number of frozen core states.
  - `frozen_method`: Method for frozen core calculations.

The input file parsing is line-oriented. Keywords are converted to lowercase internally, making them case-insensitive. Values are parsed as integers, reals, or strings, with strings often converted to integers for boolean flags by `elsi_str_to_int`.

# Usage Examples

An input file, conventionally named `elsi.in` or similar, would contain keyword-value pairs.

**Example `elsi.in` content**:
```
# General ELSI Settings
output              1       ! Verbosity level (0=quiet, 1=default, 2=debug)
illcond_check       true    ! Enable check for ill-conditioned overlap matrix
illcond_tol         1.0e-8  ! Tolerance for ill-conditioning

# ELPA Solver Settings
elpa_solver         2       ! Use ELPA 2-stage solver
elpa_gpu            true    ! Enable GPU acceleration for ELPA if available

# PEXSI Solver Settings
pexsi_n_pole        80      ! Number of poles for PEXSI calculation
pexsi_temp          300.0   ! Temperature for PEXSI (in Kelvin)
```

**Fortran code to read this input file**:
```fortran
program main
  use ELSI_DATATYPE, only: elsi_handle
  use ELSI_INIT, only: elsi_init ! Assuming elsi_init exists
  use ELSI_INPUT, only: elsi_set_input_file
  use ELSI_FINALIZE, only: elsi_finalize ! Assuming elsi_finalize exists

  implicit none

  type(elsi_handle) :: my_handle
  character(len=256) :: input_file_name
  integer :: error_status

  ! Initialize ELSI handle (simplified)
  ! Actual elsi_init might take solver, problem size, etc.
  call elsi_init(my_handle, default_solver, ..., error_status)
  if (error_status /= 0) stop "ELSI Init Failed"

  input_file_name = "elsi.in"
  call elsi_set_input_file(my_handle, input_file_name)

  ! At this point, parameters within my_handle%ph have been updated
  ! based on the content of "elsi.in".

  ! ... proceed with ELSI calculations ...

  call elsi_finalize(my_handle)

end program main
```

# Dependencies and Interactions

- **`ELSI_DATATYPE`**: The primary interaction is with the `elsi_handle` derived type. The `elsi_set_input_file` routine modifies the parameter component (`ph`) of the passed `elsi_handle`.
- **`ELSI_MPI`**: While no direct MPI calls are made for parallel communication in this module, the `elsi_stop` utility (invoked via `elsi_check_read`) is MPI-aware and ensures proper termination in parallel runs.
- **`ELSI_OUTPUT`**:
  - `elsi_get_unit`: Used to obtain a free Fortran I/O unit for opening and reading the input file.
  - `elsi_stop`: Called via `elsi_check_read` to halt execution if parsing errors occur or if the input file cannot be opened.
- **`ELSI_PRECISION`**: Defines the kinds `r8` (double precision real) and `i4` (integer) for variables that temporarily store the parsed values.
- **`ELSI_SET`**: This module is a crucial dependency. For each recognized keyword in the input file, `elsi_set_input_file` calls a corresponding subroutine from the `ELSI_SET` module (e.g., `elsi_set_output`, `elsi_set_elpa_gpu`, `elsi_set_pexsi_n_pole`) to actually apply the parameter value to the `elsi_handle`.
- **`ELSI_UTIL`**:
  - `elsi_check_init`: Ensures that the ELSI handle passed to `elsi_set_input_file` has been properly initialized.
- **Standard Fortran I/O**: The module relies on intrinsic Fortran features like `open`, `read` (list-directed for parsing), and `close` for file handling, and string functions like `trim` and `adjustl`.
- **No external parsing libraries** (like FortJSON) are used in this module; parsing is handled with basic Fortran string and I/O operations.
