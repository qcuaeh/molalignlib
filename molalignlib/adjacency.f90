module adjacency
use parameters
use derived_types
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
public match_bonds_to_mol1
public bonds_union
public bonds_intersection

interface adjacencydiff
   module procedure adjacencydiff_perm
end interface

type, public :: adjc_t
   integer, allocatable :: adjlist(:)
end type

contains

function adjcs_to_adjmat(adjcs) result(adjmat)
! Convert adjacency lists to adjacency matrix
   type(adjc_t), dimension(:), intent(in) :: adjcs
   logical, dimension(:,:), allocatable :: adjmat
   integer :: i, j, k

   allocate(adjmat(size(adjcs), size(adjcs)))
   adjmat = .false.

   do i = 1, size(adjcs)
      do j = 1, size(adjcs(i)%adjlist)
         k = adjcs(i)%adjlist(j)
         adjmat(i, k) = .true.
      end do
   end do
end function

subroutine adjmat_to_adjcs(adjmat, adjcs)
   logical, dimension(:,:), intent(in) :: adjmat
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   integer :: i, j, n, count
   integer, dimension(:), allocatable :: temp_list

   n = size(adjmat, 1)
   allocate(adjcs(n))
   allocate(temp_list(n))

   do i = 1, n
      count = 0
      do j = 1, n
         if (adjmat(i, j)) then
            count = count + 1
            temp_list(count) = j
         end if
      end do
      allocate(adjcs(i)%adjlist(count))
      adjcs(i)%adjlist = temp_list(1:count)
   end do

   deallocate(temp_list)
end subroutine

function adjacencydiff_perm(atomset1, atomperm1, adjcs1, adjcs2) result(diff)
!------------------------------------------------------------------------------
! Calculate connectivity difference from adjacency lists.
! Returns the number of differing edges.
!------------------------------------------------------------------------------
   integer, dimension(:), intent(in) :: atomset1
   integer, dimension(:), intent(in) :: atomperm1
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
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
      nadjs = size(adjcs1(idx1)%adjlist)

      do j = 1, nadjs
         neighbor_idx1 = adjcs1(idx1)%adjlist(j)

         ! Only count edge if idx1 < neighbor_idx1 to avoid double-counting
         if (idx1 < neighbor_idx1) then
            total_edges1 = total_edges1 + 1

            mapped_neighbor_idx1 = atomperm1(neighbor_idx1)

            ! Check if edge (mapped_idx1, mapped_neighbor_idx1) exists in structure 2
            if (any(adjcs2(mapped_idx1)%adjlist(:) == mapped_neighbor_idx1)) then
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
      nadjs = size(adjcs2(idx2)%adjlist)

      do j = 1, nadjs
         neighbor_idx1 = adjcs2(idx2)%adjlist(j)

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
   type(adjc_t), dimension(:), intent(in) :: adjcs1
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(in) :: atomperm1
   integer, intent(in) :: k, l
   integer :: i, nkk, nkl, nll, nlk, delta, nadjs_k, nadjs_l

   nadjs_k = size(adjcs1(k)%adjlist)
   nadjs_l = size(adjcs1(l)%adjlist)

   nkk = 0
   nkl = 0

   do i = 1, nadjs_k
      if (adjcs1(k)%adjlist(i) /= l) then
         if (adjmat2(atomperm1(k), atomperm1(adjcs1(k)%adjlist(i)))) nkk = nkk + 1
         if (adjmat2(atomperm1(l), atomperm1(adjcs1(k)%adjlist(i)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, nadjs_l
      if (adjcs1(l)%adjlist(i) /= k) then
         if (adjmat2(atomperm1(l), atomperm1(adjcs1(l)%adjlist(i)))) nll = nll + 1
         if (adjmat2(atomperm1(k), atomperm1(adjcs1(l)%adjlist(i)))) nlk = nlk + 1
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
   integer :: num_atoms, max_edges, bond_count
   integer, dimension(:,:), allocatable :: temp_bonds
   integer :: atom1, atom2
   logical :: bond_in_mol1, bond_in_mol2

   num_atoms = size(atomset1)
   ! Maximum possible differing edges
   max_edges = num_atoms * (num_atoms - 1) / 2

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

subroutine match_bonds_to_mol1(adjmat2, moldiffs)
   ! Modify mol2's bonds to match mol1's connectivity
   ! If bond exists in mol2: remove it (exists in mol2 but not mol1)
   ! If bond doesn't exist in mol2: add it (exists in mol1 but not mol2)
   logical, dimension(:,:), intent(inout) :: adjmat2
   integer, dimension(:,:), intent(in) :: moldiffs
   integer :: i, atom1, atom2

   do i = 1, size(moldiffs, 2)
      atom1 = moldiffs(1, i)
      atom2 = moldiffs(2, i)

      ! Toggle the bond: if it exists, remove it; if it doesn't exist, add it
      adjmat2(atom1, atom2) = .not. adjmat2(atom1, atom2)
      adjmat2(atom2, atom1) = .not. adjmat2(atom2, atom1)
   end do
end subroutine

subroutine bonds_union(adjmat1, adjmat2, atomperm1, moldiffs)
   logical, dimension(:,:), intent(inout) :: adjmat1, adjmat2
   integer, dimension(:), intent(in) :: atomperm1
   integer, dimension(:,:), intent(in) :: moldiffs
   integer :: i, atom1_mol2, atom2_mol2, atom1_mol1, atom2_mol1
   integer, dimension(:), allocatable :: inv_perm

   ! Create inverse permutation to map molecule 2 indices back to molecule 1
   allocate(inv_perm(size(atomperm1)))
   do i = 1, size(atomperm1)
      inv_perm(atomperm1(i)) = i
   end do

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

   deallocate(inv_perm)
end subroutine

subroutine bonds_intersection(adjmat1, adjmat2, atomperm1, moldiffs)
   logical, dimension(:,:), intent(inout) :: adjmat1, adjmat2
   integer, dimension(:), intent(in) :: atomperm1
   integer, dimension(:,:), intent(in) :: moldiffs
   integer :: i, atom1_mol2, atom2_mol2, atom1_mol1, atom2_mol1
   integer, dimension(:), allocatable :: inv_perm

   ! Create inverse permutation to map molecule 2 indices back to molecule 1
   allocate(inv_perm(size(atomperm1)))
   do i = 1, size(atomperm1)
      inv_perm(atomperm1(i)) = i
   end do

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

   deallocate(inv_perm)
end subroutine

end module
