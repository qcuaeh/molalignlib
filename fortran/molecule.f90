! MolAlignLib
! Copyright (C) 2025 José M. Vásquez, Carlos Z. Gómez

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

module molecule
! Atoms and bonds of a molecule, atom selection and coordinate transforms
use parameters
use error_codes
use str_utils
use flags
use chemdata
use adjacency
implicit none
private
public set_coords
public get_coords
public get_centroid
public include_all_atoms
public include_heavy_atoms
public pad_atoms
public complete_mapping
public bonds_from_atoms
public check_bondtypes
public adjacency_from_bonds
public bondtype_code
public bondtype_str
public is_valid_bondtype
public is_mol2_nonbond
public is_sdf_nonbond
public mol2_bondtype
public mol2_typestr
public sdf_bondtype
public sdf_bondnumber
public print_atoms
public print_bonds

type, public :: atom_t
   integer(ik) :: elnum
   integer(ik) :: group
   real(rk) :: coords(3)
end type

type, public :: bond_t
   integer(ik) :: atomidx1
   integer(ik) :: atomidx2
   integer(ik) :: bondtype
end type

! Conventional labels of MOL/SDF bond type numbers 1..10 (9 and 10 are
! V3000 only), see parameters. Number 8 (any) has no label: it is read as
! UNDEFINED_BOND_TYPE.
integer(ik), parameter :: SDF_UNDEFINED = 8
character(2), dimension(*), parameter :: SDF_BONDTYPES = &
      [character(2) :: '1', '2', '3', 'ar', 'sd', 'sa', 'da', '', 'co', 'hb']

interface get_coords
   module procedure get_coords_atoms
   module procedure get_coords_array
end interface

contains

subroutine include_all_atoms(atoms, atomset)
! Indices of all atoms. Pass only the real atoms (not the padding atoms).
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: atomidx

   allocate (atomset(size(atoms)))

   do atomidx = 1, size(atoms)
      atomset(atomidx) = atomidx
   end do
end subroutine

subroutine include_heavy_atoms(atoms, atomset)
! Indices of all heavy (non-hydrogen) atoms. Pass only the real atoms (not
! the padding atoms).
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: nel, atomidx

   allocate (atomset(count(atoms%elnum > 1)))

   nel = 0
   do atomidx = 1, size(atoms)
      if (atoms(atomidx)%elnum > 1) then
         nel = nel + 1
         atomset(nel) = atomidx
      end if
   end do
end subroutine

subroutine pad_atoms(atoms, n_padding)
! Append padding atoms (elnum = PADDING_ELNUM) until the molecule has
! n_padding atoms. The original atoms keep their indices, so bond tables and
! file line numbers stay valid. Padding atoms are recognised by index
! (i > original size) everywhere else, and their elnum never indexes the
! element tables.
   type(atom_t), dimension(:), allocatable, intent(inout) :: atoms
   integer(ik), intent(in) :: n_padding
   ! Local variables
   type(atom_t), dimension(:), allocatable :: padded
   integer(ik) :: n_real, i

   n_real = size(atoms)
   if (n_real >= n_padding) return

   allocate (padded(n_padding))
   padded(1:n_real) = atoms
   do i = n_real + 1, n_padding
      padded(i)%elnum = PADDING_ELNUM
      padded(i)%group = 0
      padded(i)%coords = 0.0_rk
   end do

   call move_alloc(padded, atoms)
end subroutine

subroutine neighbor_table(n_atoms, bonds, nbr_start, nbr_list)
! Compressed neighbour lists: the neighbours of atom i are
! nbr_list(nbr_start(i) : nbr_start(i+1)-1). Every bond counts, whatever
! its type.
   integer(ik), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(ik), dimension(:), allocatable, intent(out) :: nbr_start, nbr_list
   ! Local variables
   integer(ik), dimension(:), allocatable :: fill
   integer(ik) :: i, a1, a2

   allocate (nbr_start(n_atoms + 1))
   allocate (fill(n_atoms))
   fill = 0

   do i = 1, size(bonds)
      a1 = bonds(i)%atomidx1
      a2 = bonds(i)%atomidx2
      fill(a1) = fill(a1) + 1
      fill(a2) = fill(a2) + 1
   end do

   nbr_start(1) = 1
   do i = 1, n_atoms
      nbr_start(i + 1) = nbr_start(i) + fill(i)
   end do

   allocate (nbr_list(nbr_start(n_atoms + 1) - 1))
   fill = nbr_start(1:n_atoms)

   do i = 1, size(bonds)
      a1 = bonds(i)%atomidx1
      a2 = bonds(i)%atomidx2
      nbr_list(fill(a1)) = a2
      fill(a1) = fill(a1) + 1
      nbr_list(fill(a2)) = a1
      fill(a2) = fill(a2) + 1
   end do
end subroutine

subroutine complete_mapping(atomset1, atomset2, mapping1, atoms1, atoms2, &
      n_real1, n_real2, bonds1, bonds2, coords1, coords2, full_atomperm1)
! Build a full permutation of 1..size(atoms1) from the atom mapping of the
! included atoms. mapping1 is a permutation of 1..size(atomset1) in the
! compact numbering of the included atoms: included atom i of molecule 1
! (atom atomset1(i)) maps to included atom mapping1(i) of molecule 2 (atom
! atomset2(mapping1(i))). Both molecules must already be padded to the same
! size. Excluded atoms are paired in this order of preference:
!   1. Same element, bonded to the image of one of its included neighbours
!      (e.g. an H follows its heavy atom), closest first.
!   2. Same element, closest remaining real atom.
!   3. Whatever is left, in ascending index order. Because padding atoms are
!      appended at the end, real atoms are used up before padding atoms.
! atoms, bonds and coords are full (padded) arrays; coords1 and coords2 must
! be in the same (aligned) frame.
   integer(ik), dimension(:), intent(in) :: atomset1, atomset2
   integer(ik), dimension(:), intent(in) :: mapping1
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   integer(ik), intent(in) :: n_real1, n_real2
   type(bond_t), dimension(:), intent(in) :: bonds1, bonds2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), dimension(:), allocatable, intent(out) :: full_atomperm1
   ! Local variables
   logical(lk), dimension(:), allocatable :: in_set1, used2
   integer(ik), dimension(:), allocatable :: nbr1_start, nbr1_list, nbr2_start, nbr2_list
   integer(ik) :: n_atoms, n_incl, i, j, k, p, q, h, best
   real(rk) :: dist, best_dist

   n_atoms = size(atoms1)
   n_incl = size(atomset1)

   if (size(atoms2) /= n_atoms) then
      error stop 'complete_mapping: molecules must have the padded size'
   end if

   if (size(atomset2) /= n_incl .or. size(mapping1) /= n_incl) then
      error stop 'complete_mapping: atom sets and mapping must have the same size'
   end if

   allocate (full_atomperm1(n_atoms))
   allocate (in_set1(n_atoms), used2(n_atoms))
   full_atomperm1 = 0
   in_set1 = .FALSE.
   used2 = .FALSE.

   ! Translate the mapping of the included atoms to full indices, checking
   ! that it is a valid permutation of the included atoms
   do i = 1, n_incl
      j = mapping1(i)
      if (j < 1 .or. j > n_incl) then
         error stop 'complete_mapping: included atom mapped out of range'
      end if
      k = atomset2(j)
      if (used2(k)) then
         error stop 'complete_mapping: atom of molecule 2 assigned twice'
      end if
      full_atomperm1(atomset1(i)) = k
      in_set1(atomset1(i)) = .TRUE.
      used2(k) = .TRUE.
   end do

   ! Pass 1: follow bonds from included neighbours
   call neighbor_table(n_atoms, bonds1, nbr1_start, nbr1_list)
   call neighbor_table(n_atoms, bonds2, nbr2_start, nbr2_list)

   do i = 1, n_real1
      if (full_atomperm1(i) /= 0) cycle
      best = 0
      best_dist = huge(best_dist)
      do p = nbr1_start(i), nbr1_start(i + 1) - 1
         h = nbr1_list(p)
         if (.not. in_set1(h)) cycle
         j = full_atomperm1(h)
         do q = nbr2_start(j), nbr2_start(j + 1) - 1
            k = nbr2_list(q)
            if (k > n_real2) cycle
            if (used2(k)) cycle
            if (atoms2(k)%elnum /= atoms1(i)%elnum) cycle
            dist = sum((coords1(:, i) - coords2(:, k))**2)
            if (dist < best_dist) then
               best_dist = dist
               best = k
            end if
         end do
      end do
      if (best > 0) then
         full_atomperm1(i) = best
         used2(best) = .TRUE.
      end if
   end do

   ! Pass 2: closest remaining real atom of the same element
   do i = 1, n_real1
      if (full_atomperm1(i) /= 0) cycle
      best = 0
      best_dist = huge(best_dist)
      do k = 1, n_real2
         if (used2(k)) cycle
         if (atoms2(k)%elnum /= atoms1(i)%elnum) cycle
         dist = sum((coords1(:, i) - coords2(:, k))**2)
         if (dist < best_dist) then
            best_dist = dist
            best = k
         end if
      end do
      if (best > 0) then
         full_atomperm1(i) = best
         used2(best) = .TRUE.
      end if
   end do

   ! Pass 3: fill the remaining slots in ascending order. The counts of
   ! free slots and free targets are equal, so k never runs past n_atoms.
   k = 0
   do i = 1, n_atoms
      if (full_atomperm1(i) /= 0) cycle
      do
         k = k + 1
         if (.not. used2(k)) exit
      end do
      full_atomperm1(i) = k
      used2(k) = .TRUE.
   end do
end subroutine

subroutine set_coords(atoms, coords)
! Copy a 3 x n coordinate array into atoms
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(rk), dimension(:,:), intent(in) :: coords
   ! Local variables
   integer(ik) :: i

   do i = 1, size(atoms)
      atoms(i)%coords = coords(:, i)
   end do
end subroutine

subroutine bonds_from_atoms(atoms, bonds)
! Bonds perceived from geometry: two atoms are bonded when their distance is
! below the sum of their covalent radii plus bond_tol. All bonds are
! UNDEFINED_BOND_TYPE, since geometry gives no bond type. Pass only the real
! atoms (not the padding atoms).
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   integer(ik) :: i, j, n_atoms, n_bonds
   real(rk) :: atom_dist
   real(rk), dimension(:), allocatable :: atom_radii
   logical(lk), dimension(:,:), allocatable :: is_bonded

   n_atoms = size(atoms)
   allocate (is_bonded(n_atoms, n_atoms))
   is_bonded = .FALSE.

   atom_radii = covalent_radii(atoms%elnum)

   do i = 1, n_atoms
      do j = i + 1, n_atoms
         atom_dist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
         is_bonded(i, j) = atom_dist < atom_radii(i) + atom_radii(j) + bond_tol
      end do
   end do

   n_bonds = count(is_bonded)
   allocate (bonds(n_bonds))

   n_bonds = 0
   do i = 1, n_atoms
      do j = i + 1, n_atoms
         if (is_bonded(i, j)) then
            n_bonds = n_bonds + 1
            bonds(n_bonds)%atomidx1 = i
            bonds(n_bonds)%atomidx2 = j
            bonds(n_bonds)%bondtype = UNDEFINED_BOND_TYPE
         end if
      end do
   end do

   deallocate (is_bonded)
end subroutine

subroutine check_bondtypes(atomset, n_atoms, bonds, error_code)
! Check the bond types before adjacency_from_bonds. With usebondtype_flag
! the type of every bond between atoms of atomset must be a valid bond
! type (see is_valid_bondtype): error_code is
! MOLALIGN_ERROR_UNDEFINED_BOND_TYPE for a bond of undefined type
! (UNDEFINED_BOND_TYPE), which cannot be compared, and
! MOLALIGN_ERROR_INVALID_BOND_TYPE for any other invalid code.
! Bonds with an end outside the set are never used, so their types are not
! checked. Without usebondtype_flag bond types are not used, so they are not
! validated: any value is accepted (see adjacency_from_bonds). bonds use full atom indices
! (1..n_atoms).
   integer(ik), dimension(:), intent(in) :: atomset
   integer(ik), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(ik), intent(out) :: error_code
   ! Local variables
   logical(lk), dimension(:), allocatable :: in_set
   integer(ik) :: i, bondtype

   error_code = MOLALIGN_SUCCESS
   if (.not. usebondtype_flag) return

   allocate (in_set(n_atoms))
   in_set = .FALSE.
   in_set(atomset) = .TRUE.

   do i = 1, size(bonds)
      if (.not. (in_set(bonds(i)%atomidx1) .and. in_set(bonds(i)%atomidx2))) cycle
      bondtype = bonds(i)%bondtype
      if (bondtype == UNDEFINED_BOND_TYPE) then
         error_code = MOLALIGN_ERROR_UNDEFINED_BOND_TYPE
         return
      end if
      if (.not. is_valid_bondtype(bondtype)) then
         error_code = MOLALIGN_ERROR_INVALID_BOND_TYPE
         return
      end if
   end do
end subroutine

subroutine adjacency_from_bonds(atomset, n_atoms, bonds, adjcs, error_code)
! Adjacency lists of the atoms of atomset, in the compact numbering of the
! set (atom atomset(i) becomes atom i). bonds use full atom indices
! (1..n_atoms); bonds with an end outside the set are dropped. bonds are not
! modified. Every entry of bonds is a bond: NO_BOND only marks the unbonded
! pairs of the adjacency matrix.
! Without usebondtype_flag only connectivity counts: every bond becomes
! UNDEFINED_BOND_TYPE, whatever its type code (0 included). With it each
! neighbor carries the type of its bond as stored in bonds, which
! check_bondtypes must have found valid, except that directed types are
! made undirected (see comparable_bondtype). Bond types are otherwise
! compared literally, for equality of their labels, never interpreted (see
! parameters). A bond of undefined type (UNDEFINED_BOND_TYPE) within the set
! cannot be compared, so with usebondtype_flag error_code is then
! MOLALIGN_ERROR_UNDEFINED_BOND_TYPE and adjcs is left unallocated.
   integer(ik), dimension(:), intent(in) :: atomset
   integer(ik), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   integer(ik), intent(out) :: error_code
   ! Local variables
   integer(ik), dimension(:), allocatable :: newidx
   integer(ik), dimension(:,:), allocatable :: adjmat
   integer(ik) :: n_set, atomidx1, atomidx2, bondtype, i

   error_code = MOLALIGN_SUCCESS
   n_set = size(atomset)

   ! Map full indices to compact indices (0 = not in the set)
   allocate (newidx(n_atoms))
   newidx = 0
   do i = 1, n_set
      newidx(atomset(i)) = i
   end do

   ! Adjacency matrix of the bonds within the set (bond type of each bond,
   ! NO_BOND where there is none). A pair listed more than once keeps the
   ! type of its last occurrence.
   allocate (adjmat(n_set, n_set))
   adjmat = NO_BOND
   do i = 1, size(bonds)
      bondtype = bonds(i)%bondtype
      atomidx1 = newidx(bonds(i)%atomidx1)
      atomidx2 = newidx(bonds(i)%atomidx2)
      if (atomidx1 == 0 .or. atomidx2 == 0) cycle
      if (usebondtype_flag) then
         if (bondtype == UNDEFINED_BOND_TYPE) then
            error_code = MOLALIGN_ERROR_UNDEFINED_BOND_TYPE
            return
         end if
         if (.not. is_valid_bondtype(bondtype)) then
            error stop 'adjacency_from_bonds: invalid bond type (see check_bondtypes)'
         end if
         bondtype = comparable_bondtype(bondtype)
      else
         bondtype = UNDEFINED_BOND_TYPE
      end if
      adjmat(atomidx1, atomidx2) = bondtype
      adjmat(atomidx2, atomidx1) = bondtype
   end do

   call adjmat_to_adjcs(adjmat, adjcs)
end subroutine

function bondtype_code(typestr) result(code)
! Integer code of a bond type string, case insensitive, surrounding blanks
! ignored (see parameters):
!   d     a digit 1-9                         -> d
!   ab    a letter, then a letter or a digit  -> BOND_BLOCK2 + 36*a + b
!   nsd   a digit, a separator, a digit       -> BOND_BLOCK3 + 100*s + 10*n + d
! with a = 0..25 for a..z, b = 0..35 for 0..9, a..z, n, d = 0..9 and s the
! position (from 0) of the separator in BONDTYPE_SEPARATORS. The code
! identifies the string and nothing else: equal codes mean equal strings.
! -1 for a string of any other form (including blank and 0), which is not a
! bond type.
   character(*), intent(in) :: typestr
   integer(ik) :: code
   ! Local variables
   character(:), allocatable :: trimmed
   integer(ik) :: char1, char2, separator

   code = -1
   trimmed = lowercase(trim(adjustl(typestr)))

   select case (len(trimmed))
   case (1)
      ! The digit's value; 0 is not a bond type
      char1 = index(BONDTYPE_DIGITS, trimmed) - 1
      if (char1 >= 1) code = char1
   case (2)
      char1 = index(BONDTYPE_LETTERS, trimmed(1:1)) - 1
      char2 = index(BONDTYPE_ALNUM, trimmed(2:2)) - 1
      if (char1 >= 0 .and. char2 >= 0) then
         code = BOND_BLOCK2 + len(BONDTYPE_ALNUM)*char1 + char2
      end if
   case (3)
      char1 = index(BONDTYPE_DIGITS, trimmed(1:1)) - 1
      char2 = index(BONDTYPE_DIGITS, trimmed(3:3)) - 1
      separator = index(BONDTYPE_SEPARATORS, trimmed(2:2)) - 1
      if (separator >= 0 .and. char1 >= 0 .and. char2 >= 0) then
         code = BOND_BLOCK3 + 100*separator + 10*char1 + char2
      end if
   end select
end function

function bondtype_str(code) result(typestr)
! Bond type string of a valid code, 1..MAX_BOND_TYPE (the inverse of
! bondtype_code): one to three characters, lowercase
   integer(ik), intent(in) :: code
   character(:), allocatable :: typestr
   ! Local variables
   integer(ik) :: offset, i, j, k

   if (code >= 1 .and. code < BOND_BLOCK2) then
      typestr = BONDTYPE_DIGITS(code+1:code+1)
   else if (code >= BOND_BLOCK2 .and. code < BOND_BLOCK3) then
      offset = code - BOND_BLOCK2
      i = offset/len(BONDTYPE_ALNUM) + 1
      j = mod(offset, len(BONDTYPE_ALNUM)) + 1
      typestr = BONDTYPE_LETTERS(i:i)//BONDTYPE_ALNUM(j:j)
   else if (code >= BOND_BLOCK3 .and. code <= MAX_BOND_TYPE) then
      offset = code - BOND_BLOCK3
      k = offset/100 + 1
      i = mod(offset, 100)/10 + 1
      j = mod(offset, 10) + 1
      typestr = BONDTYPE_DIGITS(i:i)//BONDTYPE_SEPARATORS(k:k)//BONDTYPE_DIGITS(j:j)
   else
      error stop 'bondtype_str: bond type code out of range'
   end if
end function

function is_valid_bondtype(code) result(valid)
! Whether code is the code of a bond type string, i.e. 1..MAX_BOND_TYPE (see
! parameters). NO_BOND is not a bond type, and UNDEFINED_BOND_TYPE is not
! valid either: an undefined type cannot be compared.
   integer(ik), intent(in) :: code
   logical(lk) :: valid

   valid = code >= 1 .and. code <= MAX_BOND_TYPE
end function

function comparable_bondtype(code) result(comparable)
! Bond type of a valid code as compared between molecules: the code itself,
! so that types are compared literally, except for directed types. The
! adjacency matrix is symmetric, so a directed type would differ from its
! reverse (dr on bond a-b is dl on bond b-a): dative bonds dr and dl are
! compared as dv, and stereo single bonds up and dn as single bonds (1).
   integer(ik), intent(in) :: code
   integer(ik) :: comparable

   select case (bondtype_str(code))
   case ('dr', 'dl')
      comparable = bondtype_code('dv')
   case ('up', 'dn')
      comparable = bondtype_code('1')
   case default
      comparable = code
   end select
end function

function is_mol2_nonbond(typestr) result(nonbond)
! Whether a MOL2 bond type string marks an entry that is not a bond (du
! dummy, nc not connected), case insensitive. The readers drop such
! entries instead of storing them.
   character(*), intent(in) :: typestr
   logical(lk) :: nonbond

   select case (lowercase(trim(adjustl(typestr))))
   case ('du', 'nc')
      nonbond = .TRUE.
   case default
      nonbond = .FALSE.
   end select
end function

function is_sdf_nonbond(number) result(nonbond)
! Whether a MOL/SDF bond type number marks an entry that is not a bond (0,
! not a standard type). The readers drop such entries instead of storing
! them.
   integer(ik), intent(in) :: number
   logical(lk) :: nonbond

   nonbond = number == 0
end function

function mol2_bondtype(typestr) result(code)
! Bond type code of a MOL2 bond type string (not du or nc, see
! is_mol2_nonbond), case insensitive. The standard types 1, 2, 3, ar and am
! are already conventional labels (see parameters) and are kept as they
! are, and un (unknown) is UNDEFINED_BOND_TYPE; any other valid bond type
! string is accepted as an extension. -1 for a string that is not a bond
! type.
   character(*), intent(in) :: typestr
   integer(ik) :: code

   if (lowercase(trim(adjustl(typestr))) == 'un') then
      code = UNDEFINED_BOND_TYPE
   else
      code = bondtype_code(typestr)
   end if
end function

function mol2_typestr(code) result(typestr)
! MOL2 bond type string of a valid bond type code or UNDEFINED_BOND_TYPE.
! The standard MOL2 types (1, 2, 3, ar, am) are written as they are; an
! undefined type and the extensions are written as un (unknown), so that
! any MOL2 reader can read the file.
   integer(ik), intent(in) :: code
   character(2) :: typestr

   typestr = 'un'
   if (code == UNDEFINED_BOND_TYPE) return

   select case (bondtype_str(code))
   case ('1', '2', '3', 'ar', 'am')
      typestr = bondtype_str(code)
   end select
end function

function sdf_bondtype(number) result(code)
! Bond type code of a MOL/SDF bond type number 1..10 (not 0, see
! is_sdf_nonbond): UNDEFINED_BOND_TYPE for 8 (any), the code of its
! conventional label in SDF_BONDTYPES otherwise. -1 for anything else.
   integer(ik), intent(in) :: number
   integer(ik) :: code

   if (number == SDF_UNDEFINED) then
      code = UNDEFINED_BOND_TYPE
   else if (number >= 1 .and. number <= size(SDF_BONDTYPES)) then
      code = bondtype_code(SDF_BONDTYPES(number))
   else
      code = -1
   end if
end function

function sdf_bondnumber(code) result(number)
! MOL/SDF bond type number (1..10) of a valid bond type code or
! UNDEFINED_BOND_TYPE. An undefined type, and labels without a MOL/SDF
! number (anything not in SDF_BONDTYPES, e.g. quadruple and higher orders,
! fractional orders, amide, dative, haptic, ionic, multicenter and stereo
! bonds), are 8 (any).
   integer(ik), intent(in) :: code
   integer(ik) :: number
   ! Local variables
   integer(ik) :: i

   number = SDF_UNDEFINED
   if (code == UNDEFINED_BOND_TYPE) return
   do i = 1, size(SDF_BONDTYPES)
      if (bondtype_code(SDF_BONDTYPES(i)) == code) then
         number = i
         return
      end if
   end do
end function

function get_coords_atoms(atoms, weights, center, transmat) result(coords)
! Coordinates of all the atoms passed in, as a 3 x size(atoms) array. Pass
! the full atom array for all atoms, or atoms(atomset) for the included ones.
! See get_coords_array for the meaning of the other arguments.
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:), intent(in) :: weights
   real(rk), intent(in) :: center(3)
   real(rk), intent(in) :: transmat(3, 3)
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   real(rk), dimension(:,:), allocatable :: atomcoords
   integer(ik) :: i

   allocate (atomcoords(3, size(atoms)))
   do i = 1, size(atoms)
      atomcoords(:, i) = atoms(i)%coords
   end do

   coords = get_coords_array(atomcoords, weights, center, transmat)
end function

function get_coords_array(coords_in, weights, center, transmat) result(coords)
! Transformed copy of a 3 x n coordinate array. Point i becomes
!    sqrt(weights(i)) * (transmat . x_i - center)
! so the transformations are applied in this order:
!   transmat : linear transformation, e.g. IDENTITY_MATRIX or MIRROR_MATRIX
!   center   : subtracted after the transformation, so it must be given in
!              the transformed frame (see get_centroid); ORIGIN for none
!   weights  : per-atom weights, used as given (not normalised here). Pass
!              unit weights for plain coordinates, or weights normalised to
!              sum 1 for weighted coordinates.
   real(rk), dimension(:,:), intent(in) :: coords_in
   real(rk), dimension(:), intent(in) :: weights
   real(rk), intent(in) :: center(3)
   real(rk), intent(in) :: transmat(3, 3)
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   integer(ik) :: i

   if (size(weights) /= size(coords_in, dim=2)) then
      error stop 'get_coords: weights and coordinates must have the same size'
   end if

   allocate (coords(3, size(coords_in, dim=2)))

   do i = 1, size(coords_in, dim=2)
      coords(:, i) = sqrt(weights(i))*(matmul(transmat, coords_in(:, i)) - center)
   end do
end function

function get_centroid(coords, weights) result(centroid)
! Weighted center of the given coordinates (3 x n). Working on coordinates
! rather than atoms keeps the centroid in the same frame as the coordinates
! (e.g. mirrored when they are mirrored).
   real(rk), dimension(:,:), intent(in) :: coords
   real(rk), dimension(:), intent(in) :: weights
   ! Local variables
   real(rk) :: centroid(3)
   real(rk) :: total_weight, total_coords(3)
   integer(ik) :: i

   if (size(weights) /= size(coords, dim=2)) then
      error stop 'get_centroid: weights and coords must have the same size'
   end if

   total_weight = 0
   total_coords = 0
   do i = 1, size(coords, dim=2)
      total_weight = total_weight + weights(i)
      total_coords = total_coords + weights(i)*coords(:, i)
   end do
   centroid = total_coords/total_weight
end function

subroutine print_atoms(atoms)
! Debugging output
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   integer(ik) :: i
   character(:), allocatable :: fmtstr
   type(atom_t) :: atom

   write (stderr, '(A,2X,A,1X,A,4X,A,8X,A,8X,A,4X,A)') "idx", "sym", "type", &
         "X","Y","Z"

   do i = 1, size(atoms)
      atom = atoms(i)
      fmtstr = '(I3,3X,A2,1X,I3,3(1X,f8.4),2X)'
      write (stderr, fmtstr) i, element_symbol(atom%elnum), atom%group, atom%coords
   end do
end subroutine

subroutine print_bonds(bonds)
! Debugging output
   type(bond_t), dimension(:), intent(in) :: bonds
   ! Local variables
   integer(ik) :: i

   write (stderr, '(a)') "atomidx1 atomidx2"

   do i = 1, size(bonds)
      write (stderr, '(I3,2X,I3)') bonds(i)%atomidx1, bonds(i)%atomidx2
   end do
end subroutine

end module
