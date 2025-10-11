module parameters
use iso_fortran_env, only: stdin => input_unit
use iso_fortran_env, only: stdout => output_unit
use iso_fortran_env, only: stderr => error_unit
use iso_fortran_env, only: int32, int64, real32, real64
!use iso_c_binding, only: c_int, c_long, c_float, c_double
implicit none

! Enable/disable debug tests
logical, parameter :: DEBUGGING = .true.

! Numerical kinds
integer, parameter :: ik = int32 ! 32-bit integer kind
!integer, parameter :: rk = real32 ! Single precision kind
integer, parameter :: rk = real64 ! Double precision kind

! Convergence tolerance
!real(rk), parameter :: CONV_TOL = 1E-6 ! Single precision
real(rk), parameter :: CONV_TOL = 1E-10 ! Double precision

! Mean squared distance tolerance
real(rk), parameter :: MSD_TOL = 1E-6

! Common character lengths
integer, parameter :: wl = 32 ! Word length
integer, parameter :: ll = 256 ! Line length

! Maximum coordination number
integer, parameter :: MAX_COORD = 10

! Displayed decimal places
integer, parameter :: DECIMAL_PLACES = 6

end module
