module kinds
use iso_fortran_env, only: sp=>real32, dp=>real64
!use iso_c_binding, only: sp=>c_float, dp=>c_double

integer, parameter :: wp = dp ! Working precision
integer, parameter :: wl = 32 ! Word length
integer, parameter :: ll = 256 ! Line length

end module
