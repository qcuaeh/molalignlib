module parameters

use iso_fortran_env, only: int64
use iso_fortran_env, only: sp => real32
use iso_fortran_env, only: dp => real64
use iso_fortran_env, only: input_unit
use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

implicit none

integer, parameter :: wp = dp
integer, parameter :: maxlblen = 16
integer, parameter :: maxstrlen = 256

end module
