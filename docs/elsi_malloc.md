# Overview

The `elsi_malloc.f90` file defines the `ELSI_MALLOC` Fortran module. This module provides a centralized set of wrapper routines for dynamic memory allocation and deallocation within the ELSI (Electronic Structure Infrastructure) library. Instead of direct calls to Fortran's `allocate` and `deallocate` statements, ELSI uses these routines to ensure consistent error handling, optional logging of memory operations, and automatic initialization of allocated arrays to zero. This approach enhances robustness and aids in debugging memory-related issues.

# Key Components

- **Module `ELSI_MALLOC`**:
  The main module that encapsulates all custom memory management routines.

- **Generic Interface `elsi_allocate`**:
  A public interface that provides a unified way to allocate arrays of various data types, ranks (1D, 2D, 3D), and precisions. It maps to specific internal subroutines based on the array's characteristics. Supported types include:
  - `integer(kind=i4)`: 1D, 2D, 3D. Also 1D with `integer(kind=i8)` dimension.
  - `integer(kind=i8)`: 1D, 2D. Also 1D with `integer(kind=i8)` dimension.
  - `real(kind=r4)` (single precision real): 1D, 2D.
  - `real(kind=r8)` (double precision real): 1D, 2D, 3D. Also 1D with `integer(kind=i8)` dimension.
  - `complex(kind=r4)` (single precision complex, via `complex8_2d`): 2D.
  - `complex(kind=r8)` (double precision complex, via `complex16_*`): 1D, 2D, 3D. Also 1D with `integer(kind=i8)` dimension.

  Each allocation subroutine performs the following:
  1.  Optionally prints a message indicating the amount of memory being allocated and for which array (controlled by `bh%print_info`).
  2.  Calls the standard Fortran `allocate` statement.
  3.  Checks the status of the allocation. If an error occurs (`ierr > 0`), it prints detailed error messages (including the requested memory size and the caller routine) and then stops the program execution via `elsi_stop`.
  4.  Initializes all elements of the newly allocated array to zero (0, 0.0, or (0.0, 0.0) as appropriate for the type).

- **Generic Interface `elsi_deallocate`**:
  A public interface that provides a unified way to deallocate arrays previously allocated through `elsi_allocate`. It maps to specific internal subroutines matching the type and rank of the array.
  Each deallocation subroutine:
  1.  Optionally prints a message indicating which array is being deallocated (controlled by `bh%print_info`).
  2.  Calls the standard Fortran `deallocate` statement.

# Important Variables/Constants

The key parameters for the allocation and deallocation routines are:

- **`bh` (type `elsi_basic_t`, intent `in`)**: An ELSI basic information handle. Its member `bh%print_info` (integer) controls the verbosity: if `bh%print_info > 2`, detailed messages about allocation/deallocation are printed.
- **`array` (various types, intent `inout`, `allocatable`)**: The Fortran array to be allocated or deallocated.
- **`dim_1`, `dim_2`, `dim_3` (integer `i4` or `i8`, intent `in`)**: The dimensions for the array to be allocated.
- **`label` (character string, intent `in`)**: A user-provided string to identify the array in log messages (e.g., "Hamiltonian Matrix", "Workspace Array").
- **`caller` (character string, intent `in`)**: A string indicating the name of the subroutine from which the allocation is requested. This is used in error messages if allocation fails.

There are no module-level constants specific to memory sizes or limits defined within `ELSI_MALLOC` itself; limits are system-dependent.

# Usage Examples

These routines are intended for internal use by other ELSI modules to ensure robust and traceable memory management.

**Allocating a 2D double precision real array**:
```fortran
module my_calculation_module
  use ELSI_MALLOC, only: elsi_allocate, elsi_deallocate
  use ELSI_DATATYPE, only: elsi_basic_t ! For bh, assuming it's available
  use ELSI_PRECISION, only: r8, i4     ! For kind parameters

  implicit none

  subroutine perform_calculation(bh_handle, num_rows, num_cols)
    type(elsi_basic_t), intent(in) :: bh_handle
    integer(kind=i4), intent(in) :: num_rows, num_cols
    real(kind=r8), allocatable :: work_matrix(:,:)

    ! Allocate work_matrix
    call elsi_allocate(bh_handle, work_matrix, num_rows, num_cols, &
     &                 label="Work_Matrix_XYZ", caller="perform_calculation")
    
    ! work_matrix is now allocated and all its elements are zero.
    ! ... perform operations using work_matrix ...

    ! Deallocate work_matrix
    call elsi_deallocate(bh_handle, work_matrix, label="Work_Matrix_XYZ")
  end subroutine perform_calculation

end module my_calculation_module
```

# Dependencies and Interactions

- **`ELSI_DATATYPE`**: Uses the `elsi_basic_t` derived type, specifically the `bh%print_info` member, to control logging verbosity.
- **`ELSI_MPI`**: While MPI communication is not directly performed by this module, the `elsi_stop` function (which is called upon allocation failure) is MPI-aware and ensures that if one MPI process fails to allocate memory, all processes in the MPI group are properly terminated.
- **`ELSI_OUTPUT`**:
  - `elsi_say`: Used to print informative messages about memory allocation/deallocation if verbose logging is enabled.
  - `elsi_stop`: Called to terminate the program with an error message if a memory allocation fails.
- **`ELSI_PRECISION`**: Uses kind parameters `i4`, `i8`, `r4`, `r8` to define the data types of arrays for which specific allocation/deallocation routines are provided.
- **Standard Fortran**: Relies on the intrinsic `allocate` and `deallocate` statements for the actual memory operations, and `size` for calculating memory usage for logging.
- **Other ELSI Modules**: The `elsi_allocate` and `elsi_deallocate` routines are fundamental utilities that are expected to be used by virtually all other ELSI modules that require dynamic allocation of arrays. This promotes a uniform approach to memory management and error reporting throughout the ELSI library.

The file `docs/elsi_malloc.md` has been successfully created.
The documentation includes:
-   An **Overview** of `elsi_malloc.f90`, explaining its role in providing wrapped memory allocation and deallocation routines with error handling and logging.
-   **Key Components**, listing the `ELSI_MALLOC` module and the generic interfaces `elsi_allocate` and `elsi_deallocate`, along with a description of the types and ranks of arrays they support. The actions performed by these wrappers (logging, allocation, error check, zero-initialization) are detailed.
-   **Important Variables/Constants**, describing the parameters passed to the routines, such as the `elsi_basic_t` handle (`bh`), the array itself, dimensions, label for logging, and caller name for error reporting.
-   **Usage Examples**, providing a conceptual Fortran example of how to use `elsi_allocate` and `elsi_deallocate` for a 2D array.
-   **Dependencies and Interactions**, outlining the module's reliance on other ELSI modules (`ELSI_DATATYPE`, `ELSI_OUTPUT`, `ELSI_PRECISION`, `ELSI_MPI` for `elsi_stop`), standard Fortran `allocate`/`deallocate`, and its role as a fundamental utility for other ELSI modules.

All steps of the subtask have been completed.
