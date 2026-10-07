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

! Bond types.
! A bond type is an opaque label: a short case-insensitive string that has
! no intrinsic meaning. Bond types are compared literally, as strings
! (through their integer codes): two bonds match when their labels are
! the same string and differ otherwise, whatever the strings might be read
! to mean. 3/2, 1.5 and 6/4 are three different types, and no type is
! "closer" to another. The only exception is that directed types are
! compared undirected (dr/dl as dv, up/dn as 1, see
! molecule::comparable_bondtype), because an adjacency matrix is symmetric.
!
! A bond type string has one of three forms, each stored as an integer code
! in its own block:
!   d     a digit 1-9                            code d (1..9)
!   ab    a letter, then a letter or a digit     code BOND_BLOCK2 + 36*a + b
!   nsd   a digit, a separator, a digit          code BOND_BLOCK3 + 100*s + 10*n + d
! with a = 0..25 for a..z, b = 0..35 for 0..9, a..z, n, d = 0..9, and s the
! position (from 0) of the separator in BONDTYPE_SEPARATORS (/ : - . ,).
! Every string of these forms is a valid bond type, and no other string is
! (see molecule::bondtype_code and molecule::bondtype_str). The code 0
! (NO_BOND) is not the code of any string, so it is not a bond type: it is
! only used in adjacency matrices, for the atom pairs that are not bonded.
!
! UNDEFINED_BOND_TYPE (BOND_TYPE_RADIX) is not the code of any string either:
! it is the type of a bond whose type is undefined. Bonds perceived from
! geometry, MOL2 bonds of type un and MOL/SDF bonds of type 8 (any) have an
! undefined type, and so has every bond when bond types are not used (see
! below).
!
! Conventional labels. These carry no meaning for the comparison either;
! they are the strings the readers produce for the bond types of each file
! format (MOL2 types other than un are read as they are, MOL/SDF numbers
! other than 8 are mapped by molecule::sdf_bondtype), so that bond types
! read from different formats, or passed through the C interface, come out
! as the same strings and can be compared. Use them when labeling bonds by
! hand:
!   1..6     single, double, triple, quadruple, quintuple, sextuple
!   n/d      fractional bond order (1/2, 3/2, 5/2, 4/3, 5/4, ...)
!   n:m      multicenter bond with n centers and m electrons (3:2 is a
!            3c-2e bond, 3:4 a 3c-4e bond)
!   ar am    aromatic, amide
!   dv       dative, direction unspecified
!   dr dl    dative, electrons on the first / second atom
!   h1..h9   one metal-atom contact of an eta-1 .. eta-9 (haptic) bond
!   hp       haptic contact, hapticity unspecified
!   co hb    coordination, hydrogen bond
!   io       ionic
!   sd sa da query types: single or double, single or aromatic, double or
!            aromatic
!   up dn    single bond, stereo direction up / down
! Entries that are not bonds (MOL2 du and nc, MOL/SDF type 0) are dropped
! by the readers, as are bonds to dummy atoms, so every entry of a bond list
! is a bond.
!
! When bond types are not used, they are not validated, and every bond
! becomes UNDEFINED_BOND_TYPE whatever its code (0 included), so only
! connectivity is compared. When bond types are used, the types between
! compared atoms must be valid codes (1 to MAX_BOND_TYPE), and
! UNDEFINED_BOND_TYPE is an error, because an undefined type cannot be
! compared (see molecule::check_bondtypes and molecule::adjacency_from_bonds).
! Two atom pairs match when their bond types are equal, so "no bond"
! behaves as a bond of type NO_BOND. A partition signature entry encodes
! (neighbor part index, bond type) as
!    part_idx*BOND_TYPE_RADIX + bondtype
! which is unique because the bond types in adjacency lists are 1 to
! BOND_TYPE_RADIX (NO_BOND never appears in them).
character(*), parameter :: BONDTYPE_DIGITS = '0123456789'
character(*), parameter :: BONDTYPE_LETTERS = 'abcdefghijklmnopqrstuvwxyz'
character(*), parameter :: BONDTYPE_ALNUM = BONDTYPE_DIGITS//BONDTYPE_LETTERS
character(*), parameter :: BONDTYPE_SEPARATORS = '/:-.,'
integer(ik), parameter :: NO_BOND = 0
integer(ik), parameter :: BOND_BLOCK2 = 10
integer(ik), parameter :: BOND_BLOCK3 = BOND_BLOCK2 &
      + len(BONDTYPE_LETTERS)*len(BONDTYPE_ALNUM)
integer(ik), parameter :: MAX_BOND_TYPE = BOND_BLOCK3 &
      + 100*len(BONDTYPE_SEPARATORS) - 1
integer(ik), parameter :: BOND_TYPE_RADIX = MAX_BOND_TYPE + 1
integer(ik), parameter :: UNDEFINED_BOND_TYPE = BOND_TYPE_RADIX

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
