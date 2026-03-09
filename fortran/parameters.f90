! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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
use, intrinsic :: iso_fortran_env, only: int64, real64
use, intrinsic :: iso_fortran_env, only: stdout => output_unit
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use kinds

implicit none

! Do debug tests
logical(lk), parameter :: DEBUG_TESTS = .FALSE.

! Exit if random trials exceeds MAX_TRIALS
logical(lk), parameter :: MAX_TRIALS_EXIT = .TRUE.

! Prune unfeasible atom assignments
logical(lk), parameter :: PRUNE_ASSIGNMENTS = .TRUE.

! Convergence tolerance
!real(rk), parameter :: CONV_TOL = 1E-6 ! Single precision
real(rk), parameter :: CONV_TOL = 1E-10 ! Double precision

! Squared distance tolerance
real(rk), parameter :: SQDIST_TOL = 1E-6

! Common character lengths
integer(ik), parameter :: ll = 256 ! Line length

! Maximum coordination number
integer(ik), parameter :: MAX_COORDNUM = 10

! Default number of random trials to exit early
integer(ik), parameter :: MAX_TRIALS_DEFAULT = 10000

! Displayed decimal places
integer(ik), parameter :: DECIMAL_PLACES = 6

end module
