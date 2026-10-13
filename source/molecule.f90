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
public is_valid_bondtype
public print_atoms
public print_bonds

type, public :: atom_t
   integer(IK) :: elnum
   integer(IK) :: group
   real(RK) :: coords(3)
end type

type, public :: bond_t
   integer(IK) :: atomidx1
   integer(IK) :: atomidx2
   integer(IK) :: bondtype
end type

interface get_coords
   module procedure get_coords_atoms
   module procedure get_coords_array
end interface

contains

subroutine include_all_atoms(atoms, atomset)
! Indices of all atoms. Pass only the real atoms (not the padding atoms).
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(IK), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(IK) :: atomidx

   allocate (atomset(size(atoms)))

   do atomidx = 1, size(atoms)
      atomset(atomidx) = atomidx
   end do
end subroutine

subroutine include_heavy_atoms(atoms, atomset)
! Indices of all heavy (non-hydrogen) atoms. Pass only the real atoms (not
! the padding atoms).
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(IK), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(IK) :: nel, atomidx

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
   integer(IK), intent(in) :: n_padding
   ! Local variables
   type(atom_t), dimension(:), allocatable :: padded
   integer(IK) :: n_real, i

   n_real = size(atoms)
   if (n_real >= n_padding) return

   allocate (padded(n_padding))
   padded(1:n_real) = atoms
   do i = n_real + 1, n_padding
      padded(i)%elnum = PADDING_ELNUM
      padded(i)%group = 0
      padded(i)%coords = 0.0_RK
   end do

   call move_alloc(padded, atoms)
end subroutine

subroutine neighbor_table(n_atoms, bonds, nbr_start, nbr_list)
! Compressed neighbour lists: the neighbours of atom i are
! nbr_list(nbr_start(i) : nbr_start(i+1)-1). Every bond counts, whatever
! its type.
   integer(IK), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(IK), dimension(:), allocatable, intent(out) :: nbr_start, nbr_list
   ! Local variables
   integer(IK), dimension(:), allocatable :: fill
   integer(IK) :: i, a1, a2

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
   integer(IK), dimension(:), intent(in) :: atomset1, atomset2
   integer(IK), dimension(:), intent(in) :: mapping1
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   integer(IK), intent(in) :: n_real1, n_real2
   type(bond_t), dimension(:), intent(in) :: bonds1, bonds2
   real(RK), dimension(:,:), intent(in) :: coords1, coords2
   integer(IK), dimension(:), allocatable, intent(out) :: full_atomperm1
   ! Local variables
   logical(LK), dimension(:), allocatable :: in_set1, used2
   integer(IK), dimension(:), allocatable :: nbr1_start, nbr1_list, nbr2_start, nbr2_list
   integer(IK) :: n_atoms, n_incl, i, j, k, p, q, h, best
   real(RK) :: dist, best_dist

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
   real(RK), dimension(:,:), intent(in) :: coords
   ! Local variables
   integer(IK) :: i

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
   integer(IK) :: i, j, n_atoms, n_bonds
   real(RK) :: atom_dist
   real(RK), dimension(:), allocatable :: atom_radii
   logical(LK), dimension(:,:), allocatable :: is_bonded

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
   integer(IK), dimension(:), intent(in) :: atomset
   integer(IK), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(IK), intent(out) :: error_code
   ! Local variables
   logical(LK), dimension(:), allocatable :: in_set
   integer(IK) :: i, bondtype

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
! neighbor carries the type of its bond exactly as stored in bonds, which
! check_bondtypes must have found valid. Bond types are compared literally,
! for equality of their codes, never interpreted (see parameters). A bond
! of undefined type (UNDEFINED_BOND_TYPE) within the set
! cannot be compared, so with usebondtype_flag error_code is then
! MOLALIGN_ERROR_UNDEFINED_BOND_TYPE and adjcs is left unallocated.
   integer(IK), dimension(:), intent(in) :: atomset
   integer(IK), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   integer(IK), intent(out) :: error_code
   ! Local variables
   integer(IK), dimension(:), allocatable :: newidx
   integer(IK), dimension(:,:), allocatable :: adjmat
   integer(IK) :: n_set, atomidx1, atomidx2, bondtype, i

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
      else
         bondtype = UNDEFINED_BOND_TYPE
      end if
      adjmat(atomidx1, atomidx2) = bondtype
      adjmat(atomidx2, atomidx1) = bondtype
   end do

   call adjmat_to_adjcs(adjmat, adjcs)
end subroutine

function is_valid_bondtype(code) result(valid)
! Whether code is the code of a bond type string, i.e. 1..MAX_BOND_TYPE (see
! parameters). NO_BOND is not a bond type, and UNDEFINED_BOND_TYPE is not
! valid either: an undefined type cannot be compared.
   integer(IK), intent(in) :: code
   logical(LK) :: valid

   valid = code >= 1 .and. code <= MAX_BOND_TYPE
end function

function get_coords_atoms(atoms, weights, center, transmat) result(coords)
! Coordinates of all the atoms passed in, as a 3 x size(atoms) array. Pass
! the full atom array for all atoms, or atoms(atomset) for the included ones.
! See get_coords_array for the meaning of the other arguments.
   type(atom_t), dimension(:), intent(in) :: atoms
   real(RK), dimension(:), intent(in) :: weights
   real(RK), intent(in) :: center(3)
   real(RK), intent(in) :: transmat(3, 3)
   real(RK), dimension(:,:), allocatable :: coords
   ! Local variables
   real(RK), dimension(:,:), allocatable :: atomcoords
   integer(IK) :: i

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
   real(RK), dimension(:,:), intent(in) :: coords_in
   real(RK), dimension(:), intent(in) :: weights
   real(RK), intent(in) :: center(3)
   real(RK), intent(in) :: transmat(3, 3)
   real(RK), dimension(:,:), allocatable :: coords
   ! Local variables
   integer(IK) :: i

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
   real(RK), dimension(:,:), intent(in) :: coords
   real(RK), dimension(:), intent(in) :: weights
   ! Local variables
   real(RK) :: centroid(3)
   real(RK) :: total_weight, total_coords(3)
   integer(IK) :: i

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
   integer(IK) :: i
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
   integer(IK) :: i

   write (stderr, '(a)') "atomidx1 atomidx2"

   do i = 1, size(bonds)
      write (stderr, '(I3,2X,I3)') bonds(i)%atomidx1, bonds(i)%atomidx2
   end do
end subroutine

end module
