module parameters
use iso_fortran_env, only: stdin => input_unit
use iso_fortran_env, only: stdout => output_unit
use iso_fortran_env, only: stderr => error_unit
use iso_fortran_env, only: int32, int64, real32, real64
!use iso_c_binding, only: c_int, c_long, c_float, c_double

! Default kinds
integer, parameter :: ik = int32 ! Selected integer kind
integer, parameter :: rk = real64 ! Selected real kind

! Fixed character string lengths
integer, parameter :: wl = 32 ! Word length
integer, parameter :: ll = 256 ! Line length

! Integer parameters
integer, parameter :: DEC_PREC = 6
integer, parameter :: MAX_COORD = 10

! Real parameters
real(rk), parameter :: BIAS_SF = 0.001
real(rk), parameter :: CONV_TOL = max(100*epsilon(1.0_rk), 1.0e-10_rk)

end module
