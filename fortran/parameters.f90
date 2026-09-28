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

! Search strategy for conformer alignment. If neither is set, the strategy
! is chosen automatically from the ratio of total to partial assignment
! combinations in the assignment tree: stochastic fixed-orientation search
! when the ratio is high, exhaustive orientation-independent search when it
! is low. Setting one of them forces that strategy regardless of the tree
! topology. They are mutually exclusive; if both are set, the exhaustive
! search takes precedence.
logical(lk), parameter :: FORCE_EXHAUSTIVE = .FALSE.
logical(lk), parameter :: FORCE_STOCHASTIC = .FALSE.

! Convergence tolerance
!real(rk), parameter :: CONV_TOL = 1E-6 ! Single precision
real(rk), parameter :: CONV_TOL = 1E-10 ! Double precision

! Squared distance tolerance
real(rk), parameter :: SQDIST_TOL = 1E-6

! Common character lengths
integer(ik), parameter :: ll = 256 ! Line length

! Maximum coordination number
integer(ik), parameter :: MAX_COORDNUM = 10

! Bond types. The bond state of an atom pair is an integer bond type:
! NO_BOND (0) when the atoms are not bonded, and a positive type when they
! are. Without bond types every bond is ANY_BOND; with bond types the
! types read from file are compacted to 1, 2, ... Two atom pairs match
! when their bond types are equal, so "no bond" behaves as a bond of type
! zero. A partition signature entry encodes (neighbor part index, bond
! type) as
!    part_idx*BOND_TYPE_RADIX + bondtype
! Bond types above MAX_BOND_TYPE are merged into MAX_BOND_TYPE, which only
! makes the comparison coarser, never wrong.
integer(ik), parameter :: NO_BOND = 0
integer(ik), parameter :: ANY_BOND = 1
integer(ik), parameter :: BOND_TYPE_RADIX = 256
integer(ik), parameter :: MAX_BOND_TYPE = BOND_TYPE_RADIX - 1

! Default number of random trials to exit early
integer(ik), parameter :: MAX_TRIALS_DEFAULT = 10000

! Displayed decimal places
integer(ik), parameter :: DECIMAL_PLACES = 6

end module
