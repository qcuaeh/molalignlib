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
! Error codes returned by atormsd_calculate and conformsd_calculate (their
! error_code argument) and by the routines they call. This module is the canonical numbering; error_codes.h mirrors it by
! hand for C, so update it (and the Cython wrapper) when a value changes.
! Each value has one meaning; the notes say which functions can return it.
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
   ! No error. Returned by: both.
   enumerator :: MOLALIGN_SUCCESS                  = 0
   ! Different atom counts or compositions. Returned by: both.
   enumerator :: MOLALIGN_ERROR_NOT_ISOMERS        = 1
   ! Atom types differ in input order. Returned by: both (remap_flag=false only).
   enumerator :: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH = 2
   ! One or both molecules have no bonds. Returned by: conformsd.
   enumerator :: MOLALIGN_ERROR_MISSING_BONDS      = 3
   ! Bonds differ in input order. Returned by: conformsd (remap_flag=false only).
   enumerator :: MOLALIGN_ERROR_BOND_MISMATCH      = 4
   ! Same composition but non-isomorphic bond graphs. Returned by: conformsd
   ! (remap_flag=true only).
   ! conformer refinement, where it is handled and never returned.
   enumerator :: MOLALIGN_ERROR_NOT_CONFORMERS     = 5
   ! No valid assignment under the pruning constraints (pruning tolerance
   ! might be too tight). Returned by: atormsd (remap_flag=true only).
   enumerator :: MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED  = 6
   ! An atomic number is outside the element tables in chemdata
   ! (0:n_elems, where 0 is the dummy atom). Returned by: atormsd,
   ! conformsd.
   enumerator :: MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER     = 7
end enum

end module error_codes
