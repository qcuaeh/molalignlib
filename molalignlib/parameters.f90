module parameters
use iso_fortran_env, only: stdin => input_unit
use iso_fortran_env, only: stdout => output_unit
use iso_fortran_env, only: stderr => error_unit
use iso_fortran_env, only: int32, int64, real32, real64
!use iso_c_binding, only: c_int, c_long, c_float, c_double

! Default kinds
integer, parameter :: ik = int32 ! Selected integer kind
integer, parameter :: rk = real64 ! Selected real kind

! Convergence tolerance for numerical methods
real(rk), parameter :: eps_tol = max(100*epsilon(1.0_rk), 1.0e-10_rk)
real(rk), parameter :: bias_scale = 0.001

! Fixed character string lengths
integer, parameter :: wl = 32 ! Word length
integer, parameter :: ll = 256 ! Line length
integer, parameter :: max_coord = 10

end module
