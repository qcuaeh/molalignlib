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
public extract_bonds
public bonds_from_atoms
public adjacency_from_bonds
public distinct_bondtypes
public mol2_bondtype
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

interface get_coords
   module procedure get_coords_atoms
   module procedure get_coords_array
end interface

contains

subroutine include_all_atoms(atoms, atomset)
! Indices of all real (non-dummy) atoms
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: nel, atomidx

   ! Dummy atoms (elnum = 0) are never included
   allocate (atomset(count(atoms%elnum > 0)))

   nel = 0
   do atomidx = 1, size(atoms)
      if (atoms(atomidx)%elnum > 0) then
         nel = nel + 1
         atomset(nel) = atomidx
      end if
   end do
end subroutine

subroutine include_heavy_atoms(atoms, atomset)
! Indices of all heavy atoms
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: nel, atomidx

   ! elnum > 1 excludes both hydrogens and dummy atoms (elnum = 0)
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
! Append dummy atoms (elnum = 0) until the molecule has n_padding atoms. The
! original atoms keep their indices, so bond tables and file line numbers
! stay valid. Padding atoms are recognised by index (i > original size)
! everywhere else.
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
      padded(i)%elnum = 0
      padded(i)%group = 0
      padded(i)%coords = 0.0_rk
   end do

   call move_alloc(padded, atoms)
end subroutine

subroutine neighbor_table(n_atoms, bonds, nbr_start, nbr_list)
! Compressed neighbour lists: the neighbours of atom i are
! nbr_list(nbr_start(i) : nbr_start(i+1)-1).
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

function extract_bonds(atomset, n_atoms, bonds) result(subbonds)
! Bonds between atoms of atomset, renumbered to the compact numbering of
! the set (atom atomset(i) becomes atom i). Bonds with an end outside the
! set are dropped.
   integer(ik), dimension(:), intent(in) :: atomset
   integer(ik), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(bond_t), dimension(:), allocatable :: subbonds
   ! Local variables
   integer(ik), dimension(:), allocatable :: newidx
   integer(ik) :: i, n_bonds, a1, a2

   ! Map full indices to compact indices (0 = not in the set)
   allocate (newidx(n_atoms))
   newidx = 0
   do i = 1, size(atomset)
      newidx(atomset(i)) = i
   end do

   n_bonds = 0
   do i = 1, size(bonds)
      if (newidx(bonds(i)%atomidx1) > 0 .and. newidx(bonds(i)%atomidx2) > 0) then
         n_bonds = n_bonds + 1
      end if
   end do

   allocate (subbonds(n_bonds))

   n_bonds = 0
   do i = 1, size(bonds)
      a1 = newidx(bonds(i)%atomidx1)
      a2 = newidx(bonds(i)%atomidx2)
      if (a1 > 0 .and. a2 > 0) then
         n_bonds = n_bonds + 1
         subbonds(n_bonds)%atomidx1 = a1
         subbonds(n_bonds)%atomidx2 = a2
         subbonds(n_bonds)%bondtype = bonds(i)%bondtype
      end if
   end do
end function

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
! below the sum of their covalent radii plus bond_tol. All bonds get type 1.
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
            bonds(n_bonds)%bondtype = 1
         end if
      end do
   end do

   deallocate (is_bonded)
end subroutine

subroutine adjacency_from_bonds(atoms, bonds, adjcs, bondtypes)
! Adjacency lists of all atoms in atoms. To restrict them to a set of
! atoms, pass the atoms of the set and the bonds from extract_bonds.
! Each neighbor carries the compacted type of its bond: the position of
! the bond's type in bondtypes (see distinct_bondtypes). Bond types are not
! interpreted, only compared, so both molecules must come from the same
! source (file format and parser). Bonds whose type is not in bondtypes are
! GENERIC_BOND, so passing an empty bondtypes ignores bond types and only
! connectivity is compared.
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   integer(ik), dimension(:), intent(in) :: bondtypes
   ! Local variables
   integer(ik), dimension(:,:), allocatable :: adjmat
   integer(ik) :: n_atoms, atomidx1, atomidx2, bondtype, i, j

   n_atoms = size(atoms)

   ! Adjacency matrix of the bonds (bond type of each bond, NO_BOND where
   ! there is none). A pair listed more than once keeps the type of its
   ! last occurrence.
   allocate (adjmat(n_atoms, n_atoms))
   adjmat = NO_BOND
   do i = 1, size(bonds)
      bondtype = GENERIC_BOND
      do j = 1, size(bondtypes)
         if (bondtypes(j) == bonds(i)%bondtype) then
            bondtype = min(j, MAX_BOND_TYPE)
            exit
         end if
      end do
      atomidx1 = bonds(i)%atomidx1
      atomidx2 = bonds(i)%atomidx2
      adjmat(atomidx1, atomidx2) = bondtype
      adjmat(atomidx2, atomidx1) = bondtype
   end do

   call adjmat_to_adjcs(adjmat, adjcs)
end subroutine

function distinct_bondtypes(bonds1, bonds2) result(bondtypes)
! Distinct bond types found in either molecule, in ascending order. Bond
! type bondtypes(k) is compacted to k in adjacency_from_bonds, so the bond
! types of both molecules are consistent when they share this array.
   type(bond_t), dimension(:), intent(in) :: bonds1, bonds2
   integer(ik), dimension(:), allocatable :: bondtypes
   ! Local variables
   integer(ik), dimension(:), allocatable :: alltypes
   integer(ik) :: n_types, i, j, value

   alltypes = [bonds1%bondtype, bonds2%bondtype]
   allocate (bondtypes(size(alltypes)))

   ! Insertion into a sorted list of unique values
   n_types = 0
   do i = 1, size(alltypes)
      value = alltypes(i)
      j = n_types
      do while (j > 0)
         if (bondtypes(j) <= value) exit
         j = j - 1
      end do
      if (j > 0) then
         if (bondtypes(j) == value) cycle
      end if
      bondtypes(j+2:n_types+1) = bondtypes(j+1:n_types)
      bondtypes(j+1) = value
      n_types = n_types + 1
   end do

   bondtypes = bondtypes(:n_types)
end function

function mol2_bondtype(typestr) result(bondtype)
! Integer bond type (MOL convention) of a MOL2 bond type string. Amide
! bonds are single bonds; dummy, unknown and not connected are 0.
   character(*), intent(in) :: typestr
   integer(ik) :: bondtype

   select case (trim(adjustl(typestr)))
   case ('1', 'am', 'AM', 'Am')
      bondtype = 1
   case ('2')
      bondtype = 2
   case ('3')
      bondtype = 3
   case ('ar', 'AR', 'Ar')
      bondtype = 4
   case default
      bondtype = 0
   end select
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
      write (stderr, fmtstr) i, atomic_symbols(atom%elnum), atom%group, atom%coords
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
