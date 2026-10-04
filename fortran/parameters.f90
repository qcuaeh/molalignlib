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

! Enable internal consistency checks
logical(lk), parameter :: DEBUG_TESTS = .FALSE.

! Branch-and-bound pruning in the conformer assignment search
logical(lk), parameter :: PRUNE_ASSIGNMENTS = .TRUE.

! Search strategy for conformer alignment. By default it is chosen from the
! ratio of total to partial assignment combinations of the assignment tree:
! stochastic fixed-orientation search when the ratio is high, exhaustive
! orientation-independent search when it is low. Setting one of these
! forces that strategy; if both are set, the exhaustive search wins.
logical(lk), parameter :: FORCE_EXHAUSTIVE = .FALSE.
logical(lk), parameter :: FORCE_STOCHASTIC = .FALSE.

! Convergence tolerance of the Jacobi eigensolver (use 1E-6 in single
! precision)
real(rk), parameter :: CONV_TOL = 1E-10

! Slack added to the distance budget of the pruned assignment search
real(rk), parameter :: SQDIST_TOL = 1E-6

! Line length
integer(ik), parameter :: ll = 256

! Maximum coordination number
integer(ik), parameter :: MAX_COORDNUM = 10

! Bond types. The bond state of an atom pair is an integer bond type:
! NO_BOND (0) when the atoms are not bonded, and a positive type when they
! are. Without bond types every bond is GENERIC_BOND; with bond types the
! types read from file are compacted to 1, 2, ... Two atom pairs match
! when their bond types are equal, so "no bond" behaves as a bond of type
! zero. A partition signature entry encodes (neighbor part index, bond
! type) as
!    part_idx*BOND_TYPE_RADIX + bondtype
! Bond types above MAX_BOND_TYPE are merged into MAX_BOND_TYPE, which only
! makes the comparison coarser, never wrong.
integer(ik), parameter :: NO_BOND = 0
integer(ik), parameter :: GENERIC_BOND = 1
integer(ik), parameter :: BOND_TYPE_RADIX = 256
integer(ik), parameter :: MAX_BOND_TYPE = BOND_TYPE_RADIX - 1

! Default maximum number of random trials
integer(ik), parameter :: MAX_TRIALS_DEFAULT = 10000

! Default maximum number of molecular fragments (connected components of
! the bond graph of the included atoms) allowed by conformsd
integer(ik), parameter :: MAX_FRAGS_DEFAULT = 1

! Default mapping frequencies
integer(ik), parameter :: ATO_FREQ_DEFAULT = 10
integer(ik), parameter :: CONFO_FREQ_DEFAULT = 100
integer(ik), parameter :: ISO_FREQ_DEFAULT = 100

! Displayed decimal places
integer(ik), parameter :: DECIMAL_PLACES = 6

end module
