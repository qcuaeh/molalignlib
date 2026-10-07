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
! Error codes returned by atormsd_calculate, conformsd_calculate and
! isormsd_calculate (their error_code argument) and by the routines they
! call. This module is the canonical numbering; molalign.h mirrors it by
! hand for C, so update it (and the Cython wrapper) when a value changes.
implicit none
private

public :: MOLALIGN_SUCCESS
public :: MOLALIGN_ERROR_NOT_ISOMERS
public :: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH
public :: MOLALIGN_ERROR_TOO_MANY_FRAGMENTS
public :: MOLALIGN_ERROR_BOND_MISMATCH
public :: MOLALIGN_ERROR_NOT_CONFORMERS
public :: MOLALIGN_ERROR_ASSIGNMENT_FAILED
public :: MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER
public :: MOLALIGN_ERROR_INVALID_BOUND
public :: MOLALIGN_ERROR_INVALID_BOND_TYPE
public :: MOLALIGN_ERROR_UNDEFINED_BOND_TYPE

enum, bind(c)
   ! No error.
   enumerator :: MOLALIGN_SUCCESS                        = 0
   ! Different atom counts or compositions.
   enumerator :: MOLALIGN_ERROR_NOT_ISOMERS              = 1
   ! Atom types differ in input order.
   enumerator :: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH       = 2
   ! More molecular fragments than max_fragments.
   enumerator :: MOLALIGN_ERROR_TOO_MANY_FRAGMENTS       = 3
   ! Bonds differ in input order.
   enumerator :: MOLALIGN_ERROR_BOND_MISMATCH            = 4
   ! Same composition but non-isomorphic bond graphs.
   enumerator :: MOLALIGN_ERROR_NOT_CONFORMERS           = 5
   ! No valid assignment under the pruning constraints.
   enumerator :: MOLALIGN_ERROR_ASSIGNMENT_FAILED        = 6
   ! Atomic number outside the element tables.
   enumerator :: MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER    = 7
   ! Count parameter less than 1.
   enumerator :: MOLALIGN_ERROR_INVALID_BOUND            = 8
   ! Bond type code outside 1..MAX_BOND_TYPE (not the code of any label).
   enumerator :: MOLALIGN_ERROR_INVALID_BOND_TYPE        = 9
   ! Bond of undefined type (UNDEFINED_BOND_TYPE) while bond types are used.
   enumerator :: MOLALIGN_ERROR_UNDEFINED_BOND_TYPE      = 10
end enum

end module error_codes
