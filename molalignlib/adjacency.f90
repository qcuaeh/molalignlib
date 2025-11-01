! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

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

module adjacency
use parameters
use types_basic
use permutation
use euclidean
use sorting
implicit none
private

public adjcs_to_adjmat
public adjmat_to_adjcs
public adjacencydiff
public adjacencydelta
public compute_differing_bonds
public match_bonds2
public add_bonds
public delete_bonds
public bond_modifier_interface

type, public :: adjcs_t
   integer, allocatable :: cns(:)
   integer, allocatable :: lists(:,:)
end type

interface adjmat_to_adjcs
   module procedure adjmat_to_adjcs_all
   module procedure adjmat_to_adjcs_atomset
end interface

interface adjacencydiff
   module procedure adjacencydiff_perm
end interface

abstract interface
   subroutine bond_modifier_interface(adjcs1, adjcs2, atomperm1, moldiffs, adjcs1_mod, adjcs2_mod)
      import
      type(adjcs_t), intent(in) :: adjcs1, adjcs2
      integer, dimension(:), intent(in) :: atomperm1
      integer, dimension(:,:), intent(in) :: moldiffs
      type(adjcs_t), intent(out) :: adjcs1_mod, adjcs2_mod
   end subroutine
end interface

contains

function adjcs_to_adjmat(adjcs) result(adjmat)
! Convert adjacency lists to adjacency matrix
   type(adjcs_t), intent(in) :: adjcs
   logical, dimension(:,:), allocatable :: adjmat
   integer :: i, j, k, atoms_size

   atoms_size = size(adjcs%cns)
   allocate(adjmat(atoms_size, atoms_size))
   adjmat = .false.

   do i = 1, atoms_size
      do j = 1, adjcs%cns(i)
         k = adjcs%lists(j, i)
         adjmat(i, k) = .true.
      end do
   end do
end function

subroutine adjmat_to_adjcs_all(adjmat, adjcs)
   logical, dimension(:,:), intent(in) :: adjmat
   type(adjcs_t), intent(out) :: adjcs
   integer :: i, j, atoms_size, nadj

   atoms_size = size(adjmat, 1)
   allocate(adjcs%cns(atoms_size))
   allocate(adjcs%lists(MAX_COORD, atoms_size))
   adjcs%lists = 0

   do i = 1, atoms_size
      nadj = 0
      do j = 1, atoms_size
         if (adjmat(i, j)) then
            nadj = nadj + 1
            if (nadj > MAX_COORD) then
               write (stderr, '(A,1X,I0,1X,A,1X,A)') &
                     'Error: Coordination number of atom', i, &
                     'exceeds', MAX_COORD
               stop 1
            end if
            adjcs%lists(nadj, i) = j
         end if
      end do
      adjcs%cns(i) = nadj
   end do
end subroutine

subroutine adjmat_to_adjcs_atomset(atomset, adjmat, adjcs)
   integer, dimension(:), intent(in) :: atomset
   logical, dimension(:,:), intent(in) :: adjmat
   type(adjcs_t), intent(out) :: adjcs
   integer :: i, nadj, atomidx, atoms_size

   atoms_size = size(adjmat, 1)
   allocate(adjcs%cns(atoms_size))
   allocate(adjcs%lists(MAX_COORD, atoms_size))
   adjcs%lists = 0

   ! Populate adjacency lists for all atoms
   do atomidx = 1, atoms_size
      nadj = 0
      ! Only include bonds to atoms in atomset
      do i = 1, size(atomset)
         if (adjmat(atomidx, atomset(i))) then
            nadj = nadj + 1
            if (nadj > MAX_COORD) then
               write (stderr, '(A,1X,I0,1X,A,1X,A)') &
                     'Error: Coordination number of atom', atomidx, &
                     'exceeds', MAX_COORD
               stop 1
            end if
            adjcs%lists(nadj, atomidx) = atomset(i)
         end if
      end do
      adjcs%cns(atomidx) = nadj
   end do
end subroutine

function adjacencydiff_perm(atomset1, atomperm1, adjcs1, adjcs2) result(diff)
!------------------------------------------------------------------------------
! Calculate connectivity difference from adjacency lists.
! Returns the number of differing edges.
!------------------------------------------------------------------------------
   integer, dimension(:), intent(in) :: atomset1
   integer, dimension(:), intent(in) :: atomperm1
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   integer :: diff
   integer :: i, j, idx1, idx2, mapped_idx1, neighbor_idx1, mapped_neighbor_idx1
   integer :: common_edges, nadjs, total_edges1, total_edges2

   ! Count common edges and total edges, counting each edge only once
   ! by only considering edges where idx1 < neighbor_idx1 (upper triangle)
   common_edges = 0
   total_edges1 = 0

   do i = 1, size(atomset1)
      idx1 = atomset1(i)
      mapped_idx1 = atomperm1(idx1)
      nadjs = adjcs1%cns(idx1)

      do j = 1, nadjs
         neighbor_idx1 = adjcs1%lists(j, idx1)

         ! Only count edge if idx1 < neighbor_idx1 to avoid double-counting
         if (idx1 < neighbor_idx1) then
            total_edges1 = total_edges1 + 1

            mapped_neighbor_idx1 = atomperm1(neighbor_idx1)

            ! Check if edge (mapped_idx1, mapped_neighbor_idx1) exists in structure 2
            if (any(adjcs2%lists(1:adjcs2%cns(mapped_idx1), mapped_idx1) == mapped_neighbor_idx1)) then
               common_edges = common_edges + 1
            end if
         end if
      end do
   end do

   ! Calculate total edges in structure 2
   ! Use atomperm1(atomset1) to get the corresponding atoms in molecule 2
   total_edges2 = 0
   do i = 1, size(atomset1)
      idx2 = atomperm1(atomset1(i))
      nadjs = adjcs2%cns(idx2)

      do j = 1, nadjs
         neighbor_idx1 = adjcs2%lists(j, idx2)

         ! Only count edge if idx2 < neighbor to avoid double-counting
         if (idx2 < neighbor_idx1) then
            total_edges2 = total_edges2 + 1
         end if
      end do
   end do

   ! Edge difference = total edges in both - 2*common_edges
   diff = total_edges1 + total_edges2 - 2*common_edges
end function

function adjacencydelta(adjcs1, adjmat2, atomperm1, k, l) result(delta)
!------------------------------------------------------------------------------
! Efficiently compute the change in adjacency difference when swapping
! atoms k and l in the permutation. Uses adjacency lists for structure 1 and
! adjacency matrix for structure 2.
!------------------------------------------------------------------------------
   type(adjcs_t), intent(in) :: adjcs1
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(in) :: atomperm1
   integer, intent(in) :: k, l
   integer :: i, nkk, nkl, nll, nlk, delta, nadjs_k, nadjs_l

   nadjs_k = adjcs1%cns(k)
   nadjs_l = adjcs1%cns(l)

   nkk = 0
   nkl = 0

   do i = 1, nadjs_k
      if (adjcs1%lists(i, k) /= l) then
         if (adjmat2(atomperm1(k), atomperm1(adjcs1%lists(i, k)))) nkk = nkk + 1
         if (adjmat2(atomperm1(l), atomperm1(adjcs1%lists(i, k)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, nadjs_l
      if (adjcs1%lists(i, l) /= k) then
         if (adjmat2(atomperm1(l), atomperm1(adjcs1%lists(i, l)))) nll = nll + 1
         if (adjmat2(atomperm1(k), atomperm1(adjcs1%lists(i, l)))) nlk = nlk + 1
      end if
   end do

   ! The change in adjacency difference when swapping k and l:
   ! delta = (new_diff_kl + new_diff_lk) - (old_diff_kk + old_diff_ll)
   ! After simplification: delta = 2*(nkk + nll - nkl - nlk)
   delta = 2*(nkk + nll - nkl - nlk)
end function

subroutine compute_differing_bonds(atomset1, atomperm1, adjmat1, adjmat2, moldiffs)

   integer, dimension(:), intent(in) :: atomset1
   integer, dimension(:), intent(in) :: atomperm1
   logical, dimension(:,:), intent(in) :: adjmat1, adjmat2
   integer, dimension(:,:), allocatable, intent(out) :: moldiffs

   ! Local variables
   integer :: i, j, idx1, idx2, mapped_idx1, mapped_idx2
   integer :: atoms_size, max_edges, bond_count
   integer, dimension(:,:), allocatable :: temp_bonds
   integer :: atom1, atom2
   logical :: bond_in_mol1, bond_in_mol2

   atoms_size = size(atomset1)
   ! Maximum possible differing edges
   max_edges = atoms_size * (atoms_size - 1) / 2

   allocate(temp_bonds(2, max_edges))
   bond_count = 0

   ! Compare all pairs of atoms in atomset1
   do i = 1, size(atomset1)
      idx1 = atomset1(i)
      mapped_idx1 = atomperm1(idx1)

      do j = i + 1, size(atomset1)
         idx2 = atomset1(j)
         mapped_idx2 = atomperm1(idx2)

         ! Check bond status in both structures
         bond_in_mol1 = adjmat1(idx1, idx2)
         bond_in_mol2 = adjmat2(mapped_idx1, mapped_idx2)

         ! If bond status differs, it's a differing bond
         if (bond_in_mol1 .neqv. bond_in_mol2) then
            ! Store atom pair with lower index first
            atom1 = min(mapped_idx1, mapped_idx2)
            atom2 = max(mapped_idx1, mapped_idx2)

            bond_count = bond_count + 1
            temp_bonds(1, bond_count) = atom1
            temp_bonds(2, bond_count) = atom2
         end if
      end do
   end do

   ! Allocate final array with exact size
   allocate(moldiffs(2, bond_count))
   moldiffs = temp_bonds(:, 1:bond_count)

   ! Sort the bonds for efficient comparison
   if (bond_count > 0) then
      call sort_pairs(moldiffs)
   end if

   deallocate(temp_bonds)
end subroutine

subroutine match_bonds2(adjcs1, adjcs2, atomperm1, moldiffs, adjcs1_mod, adjcs2_mod)
   ! Modify mol2's bonds to match mol1's connectivity
   ! If bond exists in mol2: remove it (exists in mol2 but not mol1)
   ! If bond doesn't exist in mol2: add it (exists in mol1 but not mol2)
   ! Note: adjcs1 and atomperm1 are not used but present for interface compatibility
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   integer, dimension(:), intent(in) :: atomperm1
   integer, dimension(:,:), intent(in) :: moldiffs
   type(adjcs_t), intent(out) :: adjcs1_mod, adjcs2_mod
   logical, dimension(:,:), allocatable :: adjmat2
   integer :: i, atom1, atom2, atoms_size

   atoms_size = size(adjcs2%cns)

   ! Convert adjcs2 to matrix, modify it, and convert back
   adjmat2 = adjcs_to_adjmat(adjcs2)

   do i = 1, size(moldiffs, 2)
      atom1 = moldiffs(1, i)
      atom2 = moldiffs(2, i)

      ! Toggle the bond: if it exists, remove it; if it doesn't exist, add it
      adjmat2(atom1, atom2) = .not. adjmat2(atom1, atom2)
      adjmat2(atom2, atom1) = .not. adjmat2(atom2, atom1)
   end do

   ! adjcs1 remains unchanged - copy structure
   allocate(adjcs1_mod%cns(size(adjcs1%cns)))
   allocate(adjcs1_mod%lists(MAX_COORD, size(adjcs1%cns)))
   adjcs1_mod%cns = adjcs1%cns
   adjcs1_mod%lists = adjcs1%lists

   ! Convert modified adjmat2 back to adjcs
   call adjmat_to_adjcs(adjmat2, adjcs2_mod)
end subroutine

subroutine add_bonds(adjcs1, adjcs2, atomperm1, moldiffs, adjcs1_mod, adjcs2_mod)
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   integer, dimension(:), intent(in) :: atomperm1
   integer, dimension(:,:), intent(in) :: moldiffs
   type(adjcs_t), intent(out) :: adjcs1_mod, adjcs2_mod
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer :: i, atom1_mol2, atom2_mol2, atom1_mol1, atom2_mol1
   integer, dimension(:), allocatable :: inv_perm

   ! Convert to matrices
   adjmat1 = adjcs_to_adjmat(adjcs1)
   adjmat2 = adjcs_to_adjmat(adjcs2)

   ! Get inverse permutation to map molecule 2 indices back to molecule 1
   inv_perm = inverse_permutation(atomperm1)

   do i = 1, size(moldiffs, 2)
      atom1_mol2 = moldiffs(1, i)
      atom2_mol2 = moldiffs(2, i)

      ! Add bond to molecule 2
      adjmat2(atom1_mol2, atom2_mol2) = .true.
      adjmat2(atom2_mol2, atom1_mol2) = .true.

      ! Map to molecule 1 and add bond
      atom1_mol1 = inv_perm(atom1_mol2)
      atom2_mol1 = inv_perm(atom2_mol2)
      adjmat1(atom1_mol1, atom2_mol1) = .true.
      adjmat1(atom2_mol1, atom1_mol1) = .true.
   end do

   ! Convert back to adjacency lists
   call adjmat_to_adjcs(adjmat1, adjcs1_mod)
   call adjmat_to_adjcs(adjmat2, adjcs2_mod)
end subroutine

subroutine delete_bonds(adjcs1, adjcs2, atomperm1, moldiffs, adjcs1_mod, adjcs2_mod)
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   integer, dimension(:), intent(in) :: atomperm1
   integer, dimension(:,:), intent(in) :: moldiffs
   type(adjcs_t), intent(out) :: adjcs1_mod, adjcs2_mod
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer :: i, atom1_mol2, atom2_mol2, atom1_mol1, atom2_mol1
   integer, dimension(:), allocatable :: inv_perm

   ! Convert to matrices
   adjmat1 = adjcs_to_adjmat(adjcs1)
   adjmat2 = adjcs_to_adjmat(adjcs2)

   ! Get inverse permutation to map molecule 2 indices back to molecule 1
   inv_perm = inverse_permutation(atomperm1)

   do i = 1, size(moldiffs, 2)
      atom1_mol2 = moldiffs(1, i)
      atom2_mol2 = moldiffs(2, i)

      ! Remove bond from molecule 2
      adjmat2(atom1_mol2, atom2_mol2) = .false.
      adjmat2(atom2_mol2, atom1_mol2) = .false.

      ! Map to molecule 1 and remove bond
      atom1_mol1 = inv_perm(atom1_mol2)
      atom2_mol1 = inv_perm(atom2_mol2)
      adjmat1(atom1_mol1, atom2_mol1) = .false.
      adjmat1(atom2_mol1, atom1_mol1) = .false.
   end do

   ! Convert back to adjacency lists
   call adjmat_to_adjcs(adjmat1, adjcs1_mod)
   call adjmat_to_adjcs(adjmat2, adjcs2_mod)
end subroutine

end module
