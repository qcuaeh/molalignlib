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

module error_codes
! Single source of truth for the error codes returned by atormsd_calculate,
! conformsd_calculate and isormsd_calculate (via their C-bound `error_code`
! output argument, see clib_atormsd.f90 / clib_conformsd.f90 /
! clib_isormsd.f90), and by the internal routines they call.
!
! `enum, bind(c)` guarantees these enumerators share their representation
! with a C int, so this module and error_codes.h describe the same
! wire values. There is no build step that keeps their *numbers* in sync
! automatically, though: this module is the canonical numbering, and
! error_codes.h is its C-facing mirror, maintained by hand. If a value
! changes here, update error_codes.h (and the Cython wrapper's use of
! it) to match.
!
! Every value has exactly one meaning, whichever function returns it. Not
! every function can return every code; the "Returned by" notes below say
! which can.
implicit none
private

public :: MOLALIGN_SUCCESS
public :: MOLALIGN_ERROR_NOT_ISOMERS
public :: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH
public :: MOLALIGN_ERROR_MISSING_BONDS
public :: MOLALIGN_ERROR_BOND_MISMATCH
public :: MOLALIGN_ERROR_NOT_CONFORMERS
public :: MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED
public :: MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER

enum, bind(c)
   ! No error. Returned by: all.
   enumerator :: MOLALIGN_SUCCESS                  = 0
   ! Different atom counts or compositions. Returned by: all.
   enumerator :: MOLALIGN_ERROR_NOT_ISOMERS        = 1
   ! Atom types differ in input order. Returned by: all (remap_flag=false only).
   enumerator :: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH = 2
   ! One or both molecules have no bonds. Returned by: conformsd, isormsd.
   enumerator :: MOLALIGN_ERROR_MISSING_BONDS      = 3
   ! Bonds differ in input order. Returned by: conformsd (remap_flag=false only).
   enumerator :: MOLALIGN_ERROR_BOND_MISMATCH      = 4
   ! Same composition but non-isomorphic bond graphs. Returned by: conformsd
   ! (remap_flag=true only). Also raised internally during isormsd's
   ! conformer refinement, where it is handled and never returned.
   enumerator :: MOLALIGN_ERROR_NOT_CONFORMERS     = 5
   ! No valid assignment under the pruning constraints (pruning tolerance
   ! might be too tight). Returned by: atormsd (remap_flag=true only).
   enumerator :: MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED  = 6
   ! An atomic number is outside the element tables in chemdata
   ! (0:num_elems, where 0 is the dummy atom). Returned by: atormsd,
   ! conformsd.
   enumerator :: MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER     = 7
end enum

end module error_codes
