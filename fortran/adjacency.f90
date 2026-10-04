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

module adjacency
! Bond graphs as adjacency lists with bond types, and their comparison under
! an atom mapping
use parameters
use common_types
use permutation
use euclidean
use sorting
implicit none
private

public adjacencydiff
public adjacencydelta
public adjcs_to_adjmat
public adjmat_to_adjcs
public match_bonds_to_mol2
public match_bonds_to_mol1
public intersect_bonds
public find_differing_bonds
public edge_code
public count_fragments

! Adjacency list of an atom
type, public :: adjc_t
   integer(ik) :: cn
   integer(ik) :: list(MAX_COORDNUM)
   ! Type of the bond to each neighbor in list (GENERIC_BOND when bond
   ! types are not used, a compacted bond type otherwise)
   integer(ik) :: bondtype(MAX_COORDNUM) = GENERIC_BOND
end type

! Edit modes of edit_mismatched_bonds
integer(ik), parameter :: MATCH_TO_MOL1 = 1
integer(ik), parameter :: MATCH_TO_MOL2 = 2
integer(ik), parameter :: INTERSECT = 3

contains

elemental function edge_code(part_idx, bondtype) result(code)
! Signature entry of a neighbor in part part_idx reached through a bond
! of type bondtype. assignment_conformer::update_hna_part inlines this
! expression, so keep both in sync.
   integer(ik), intent(in) :: part_idx, bondtype
   integer(ik) :: code
   code = part_idx*BOND_TYPE_RADIX + bondtype
end function

function adjacencydiff(mapping1, adjcs1, adjcs2) result(diff)
! Adjacency difference: the number of atom pairs whose bond types differ
! under mapping1, "no bond" being the bond type NO_BOND. A pair counts once
! whether it is bonded in only one molecule or bonded in both with
! different types. Without bond types this is the number of differing
! edges. mapping1 must be a full permutation of the atoms of adjcs1 onto
! the atoms of adjcs2.
!
! With A = pair bonded in 1, B = bonded in 2, C = bonded in both and
! S = bonded in both with the same type, a pair differs by A + B - C - S,
! so diff = total_edges1 + total_edges2 - common_edges - same_edges.
   integer(ik), dimension(:), intent(in) :: mapping1
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer(ik) :: diff
   integer(ik) :: j, p, idx1, idx2, mapped_idx1, neighbor_idx1, mapped_neighbor_idx1
   integer(ik) :: common_edges, same_edges, nadjs, total_edges1, total_edges2

   ! Edges of molecule 1 (each counted once, idx1 < neighbor_idx1) and those
   ! also present in molecule 2, with and without equal type
   common_edges = 0
   same_edges = 0
   total_edges1 = 0

   do idx1 = 1, size(mapping1)
      mapped_idx1 = mapping1(idx1)
      nadjs = adjcs1(idx1)%cn

      do j = 1, nadjs
         neighbor_idx1 = adjcs1(idx1)%list(j)

         if (idx1 < neighbor_idx1) then
            total_edges1 = total_edges1 + 1

            mapped_neighbor_idx1 = mapping1(neighbor_idx1)

            do p = 1, adjcs2(mapped_idx1)%cn
               if (adjcs2(mapped_idx1)%list(p) == mapped_neighbor_idx1) then
                  common_edges = common_edges + 1
                  if (adjcs2(mapped_idx1)%bondtype(p) == adjcs1(idx1)%bondtype(j)) then
                     same_edges = same_edges + 1
                  end if
                  exit
               end if
            end do
         end if
      end do
   end do

   ! Edges of molecule 2, each counted once
   total_edges2 = 0
   do idx2 = 1, size(mapping1)
      nadjs = adjcs2(idx2)%cn

      do j = 1, nadjs
         neighbor_idx1 = adjcs2(idx2)%list(j)

         if (idx2 < neighbor_idx1) then
            total_edges2 = total_edges2 + 1
         end if
      end do
   end do

   diff = total_edges1 + total_edges2 - common_edges - same_edges
end function

function adjacencydelta(adjcs1, adjmat2, mapping1, k, l) result(delta)
! Change of the adjacency difference (see adjacencydiff) when the images of
! atoms k and l are swapped in mapping1, in O(cn) operations. Uses
! adjacency lists for structure 1 and the adjacency matrix of structure 2
! (bond type of each bond, NO_BOND where there is no bond).
!
! Only pairs (k,n) and (l,n), n /= k,l, change. Per pair the difference is
! A + B - C - S (see adjacencydiff); the A and B terms are unchanged by the
! swap, so only bonds of structure 1 contribute, each by the weight
! w = C + S of its partner pair in structure 2 before minus after the swap.
! Without bond types w = 2*C, which gives the classical 2*(nkk + nll - nkl - nlk).
   type(adjc_t), dimension(:), intent(in) :: adjcs1
   integer(ik), dimension(:,:), intent(in) :: adjmat2
   integer(ik), dimension(:), intent(in) :: mapping1
   integer(ik), intent(in) :: k, l
   integer(ik) :: delta
   ! Local variables
   integer(ik) :: i, n, bondtype

   delta = 0

   do i = 1, adjcs1(k)%cn
      n = adjcs1(k)%list(i)
      if (n /= l) then
         bondtype = adjcs1(k)%bondtype(i)
         delta = delta + pair_weight(bondtype, adjmat2(mapping1(k), mapping1(n))) &
                       - pair_weight(bondtype, adjmat2(mapping1(l), mapping1(n)))
      end if
   end do

   do i = 1, adjcs1(l)%cn
      n = adjcs1(l)%list(i)
      if (n /= k) then
         bondtype = adjcs1(l)%bondtype(i)
         delta = delta + pair_weight(bondtype, adjmat2(mapping1(l), mapping1(n))) &
                       - pair_weight(bondtype, adjmat2(mapping1(k), mapping1(n)))
      end if
   end do
end function

pure function pair_weight(bondtype1, bondtype2) result(weight)
! C + S for a bond of structure 1 of type bondtype1 against the pair of
! structure 2 of type bondtype2: 1 if that pair is bonded, plus 1 if its
! type is also equal.
   integer(ik), intent(in) :: bondtype1, bondtype2
   integer(ik) :: weight

   weight = 0
   if (bondtype2 /= NO_BOND) then
      weight = 1
      if (bondtype2 == bondtype1) weight = 2
   end if
end function

subroutine find_differing_bonds(mapping1, adjmat1, adjmat2, moldiffs)
! Atom pairs (numbering of molecule 2) whose bond types differ under
! mapping1: bonded in only one molecule, or bonded in both with different
! types. Their number is the adjacency difference (see adjacencydiff).
   integer(ik), dimension(:), intent(in) :: mapping1
   integer(ik), dimension(:,:), intent(in) :: adjmat1, adjmat2
   integer(ik), dimension(:,:), allocatable, intent(out) :: moldiffs

   ! Local variables
   integer(ik) :: idx1, idx2, mapped_idx1, mapped_idx2
   integer(ik) :: n_atoms, n_bonds, max_edges
   integer(ik), dimension(:,:), allocatable :: temp_bonds
   integer(ik) :: atom1, atom2

   n_atoms = size(mapping1)
   max_edges = n_atoms * (n_atoms - 1) / 2

   allocate(temp_bonds(2, max_edges))
   n_bonds = 0

   do idx1 = 1, n_atoms
      mapped_idx1 = mapping1(idx1)

      do idx2 = idx1 + 1, n_atoms
         mapped_idx2 = mapping1(idx2)

         if (adjmat1(idx1, idx2) /= adjmat2(mapped_idx1, mapped_idx2)) then
            ! Lower index first
            atom1 = min(mapped_idx1, mapped_idx2)
            atom2 = max(mapped_idx1, mapped_idx2)

            n_bonds = n_bonds + 1
            temp_bonds(1, n_bonds) = atom1
            temp_bonds(2, n_bonds) = atom2
         end if
      end do
   end do

   allocate(moldiffs(2, n_bonds))
   moldiffs = temp_bonds(:, 1:n_bonds)

   ! Sorted, so that lists can be compared element by element
   if (n_bonds > 0) then
      call sort_pairs(moldiffs)
   end if

   deallocate(temp_bonds)
end subroutine

subroutine match_bonds_to_mol2(adjcs1, adjcs2, mapping1, adjcs1_mod, adjcs2_mod)
! Edit mol1 so that each of its bonds matches mol2 under mapping1: every
! mismatched atom pair takes mol2's state, i.e. mol2's bond with its type,
! or no bond. mol2 is unchanged.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer(ik), dimension(:), intent(in) :: mapping1
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs1_mod, adjcs2_mod

   call edit_mismatched_bonds(adjcs1, adjcs2, mapping1, MATCH_TO_MOL2, adjcs1_mod, adjcs2_mod)
end subroutine

subroutine match_bonds_to_mol1(adjcs1, adjcs2, mapping1, adjcs1_mod, adjcs2_mod)
! Edit mol2 so that each of its bonds matches mol1 under mapping1: every
! mismatched atom pair takes mol1's state, i.e. mol1's bond with its type,
! or no bond. mol1 is unchanged.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer(ik), dimension(:), intent(in) :: mapping1
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs1_mod, adjcs2_mod

   call edit_mismatched_bonds(adjcs1, adjcs2, mapping1, MATCH_TO_MOL1, adjcs1_mod, adjcs2_mod)
end subroutine

subroutine intersect_bonds(adjcs1, adjcs2, mapping1, adjcs1_mod, adjcs2_mod)
! Edit both molecules so that their bonds match under mapping1, keeping
! only the bonds both molecules have with the same type: every mismatched
! bond is deleted from both.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer(ik), dimension(:), intent(in) :: mapping1
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs1_mod, adjcs2_mod

   call edit_mismatched_bonds(adjcs1, adjcs2, mapping1, INTERSECT, adjcs1_mod, adjcs2_mod)
end subroutine

subroutine edit_mismatched_bonds(adjcs1, adjcs2, mapping1, mode, adjcs1_mod, adjcs2_mod)
! Edit the bonds of both molecules at every mismatched atom pair under
! mapping1, so that the edited molecules match under mapping1. A pair is
! mismatched when its bond types differ ("no bond" being NO_BOND): it is
! bonded in only one molecule, or bonded in both with different types.
! Without bond types mismatches reduce to bonds present in only one
! molecule.
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer(ik), dimension(:), intent(in) :: mapping1
   integer(ik), intent(in) :: mode
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs1_mod, adjcs2_mod
   ! Local variables
   integer(ik), dimension(:,:), allocatable :: adjmat1, adjmat2
   integer(ik) :: i, j, mi, mj, bondtype1, bondtype2, bondtype

   adjmat1 = adjcs_to_adjmat(adjcs1)
   adjmat2 = adjcs_to_adjmat(adjcs2)

   do i = 1, size(mapping1)
      mi = mapping1(i)
      do j = i + 1, size(mapping1)
         mj = mapping1(j)
         bondtype1 = adjmat1(i, j)
         bondtype2 = adjmat2(mi, mj)
         if (bondtype1 == bondtype2) cycle

         select case (mode)
         case (MATCH_TO_MOL2)
            bondtype = bondtype2
         case (MATCH_TO_MOL1)
            bondtype = bondtype1
         case default
            bondtype = NO_BOND
         end select

         adjmat1(i, j) = bondtype
         adjmat1(j, i) = bondtype
         adjmat2(mi, mj) = bondtype
         adjmat2(mj, mi) = bondtype
      end do
   end do

   call adjmat_to_adjcs(adjmat1, adjcs1_mod)
   call adjmat_to_adjcs(adjmat2, adjcs2_mod)
end subroutine

function adjcs_to_adjmat(adjcs) result(adjmat)
! Adjacency matrix: the type of each bond, NO_BOND where there is none
   type(adjc_t), dimension(:), intent(in) :: adjcs
   integer(ik), dimension(:,:), allocatable :: adjmat
   integer(ik) :: i, j

   allocate(adjmat(size(adjcs), size(adjcs)))
   adjmat = NO_BOND

   do i = 1, size(adjcs)
      do j = 1, adjcs(i)%cn
         adjmat(i, adjcs(i)%list(j)) = adjcs(i)%bondtype(j)
      end do
   end do
end function

subroutine adjmat_to_adjcs(adjmat, adjcs)
! Adjacency lists, with bond types, from a symmetric adjacency matrix
   integer(ik), dimension(:,:), intent(in) :: adjmat
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   integer(ik) :: i, j, n_atoms, nadj

   n_atoms = size(adjmat, 1)
   allocate(adjcs(n_atoms))

   do i = 1, n_atoms
      nadj = 0
      ! adjmat is symmetric, so read column i for contiguous access
      do j = 1, n_atoms
         if (adjmat(j, i) /= NO_BOND) then
            nadj = nadj + 1
            if (nadj > MAX_COORDNUM) then
               write (stderr, '(A,1X,I0,1X,A,1X,A)') &
                     'Coordination number of atom', i, &
                     'exceeds', MAX_COORDNUM
               stop
            end if
            adjcs(i)%list(nadj) = j
            adjcs(i)%bondtype(nadj) = adjmat(j, i)
         end if
      end do
      adjcs(i)%cn = nadj
   end do
end subroutine

function count_fragments(adjcs) result(n_frags)
! Number of molecular fragments: connected components of the bond graph
! given by adjcs. An atom without bonds is a fragment of its own, so a
! graph without bonds has as many fragments as atoms; an empty graph has
! none. Iterative depth-first search, so deep graphs cannot overflow the
! call stack.
   type(adjc_t), dimension(:), intent(in) :: adjcs
   integer(ik) :: n_frags
   ! Local variables
   logical(lk), dimension(:), allocatable :: visited
   integer(ik), dimension(:), allocatable :: stack
   integer(ik) :: n_atoms, n_stack, start, node, neighbor, j

   n_atoms = size(adjcs)
   allocate(visited(n_atoms))
   allocate(stack(n_atoms))
   visited = .false.
   n_frags = 0

   do start = 1, n_atoms
      if (visited(start)) cycle

      ! New fragment: mark every atom reachable from start. Atoms are
      ! marked when pushed, so each is pushed at most once and the stack
      ! never holds more than n_atoms entries.
      n_frags = n_frags + 1
      visited(start) = .true.
      n_stack = 1
      stack(1) = start

      do while (n_stack > 0)
         node = stack(n_stack)
         n_stack = n_stack - 1
         do j = 1, adjcs(node)%cn
            neighbor = adjcs(node)%list(j)
            if (.not. visited(neighbor)) then
               visited(neighbor) = .true.
               n_stack = n_stack + 1
               stack(n_stack) = neighbor
            end if
         end do
      end do
   end do
end function

end module
