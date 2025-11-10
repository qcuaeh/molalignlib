! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module parameters
use iso_fortran_env, only: stdout => output_unit
use iso_fortran_env, only: stderr => error_unit
use iso_fortran_env, only: int32, int64, real32, real64
!use iso_c_binding, only: c_int, c_long, c_float, c_double
implicit none

! Debug flags
logical, parameter :: DEBUG_TESTS = .false.

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

! Maximum number of atoms (currently not used)
integer, parameter :: MAX_ATOMS = 1000000

! Maximum coordination number
integer, parameter :: MAX_COORD = 10

! Displayed decimal places
integer, parameter :: DECIMAL_PLACES = 6

end module
