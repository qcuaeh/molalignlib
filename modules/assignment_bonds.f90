! MolAlignLib
! Copyright (C) 2022 José M. Vásquez, Carlos Z. Gómez

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

module assignment_bonds
use parameters
use types_basic
use permutation
use euclidean
use chemistry
use adjacency
use random
implicit none
private
public allocate_memory_pool
public deallocate_memory_pool
public remap_mismatched_bonds

! Global pool depth
integer :: pool_depth

! Memory pool for recursive work arrays
! Organized as (size, depth_level) to eliminate allocations in recursion
integer, dimension(:,:), allocatable, target :: pool_perm1, pool_perm2, pool_track
logical, dimension(:,:), allocatable, target :: pool_tracked, pool_matched1, pool_matched2
integer, dimension(:,:), allocatable, target :: pool_matches, pool_mismatches1, pool_mismatches2

contains

subroutine allocate_memory_pool(n_atoms)
!------------------------------------------------------------------------------
! Initialize memory pools for bond assignment
! Call once before using remap_mismatched_bonds
! MUST be called before any calls to remap_mismatched_bonds
!------------------------------------------------------------------------------
   integer, intent(in) :: n_atoms
   
   ! Set pool depth equal to maximum possible recursion depth
   pool_depth = n_atoms
   
   ! Allocate memory pools once for all recursive calls
   allocate(pool_perm1(n_atoms, pool_depth))
   allocate(pool_perm2(n_atoms, pool_depth))
   allocate(pool_matches(n_atoms, pool_depth))
   allocate(pool_mismatches1(n_atoms, pool_depth))
   allocate(pool_mismatches2(n_atoms, pool_depth))
   allocate(pool_track(n_atoms, pool_depth))
   allocate(pool_tracked(n_atoms, pool_depth))
   allocate(pool_matched1(n_atoms, pool_depth))
   allocate(pool_matched2(n_atoms, pool_depth))
end subroutine

subroutine deallocate_memory_pool()
!------------------------------------------------------------------------------
! Deallocate memory pools
! Call once after all remap operations are complete
!------------------------------------------------------------------------------
   deallocate(pool_perm1, pool_perm2, pool_track)
   deallocate(pool_tracked, pool_matched1, pool_matched2)
   deallocate(pool_matches, pool_mismatches1, pool_mismatches2)
end subroutine

subroutine remap_mismatched_bonds(atomset1, atomtypes, adjcs1, adjcs2, adjmat2, &
                  coords1, coords2, perm1)
!------------------------------------------------------------------------------
! Find best correspondence between points of graphs
! Randomly selects starting atoms - fragment identification is implicit
! OPTIMIZED: Uses pre-allocated memory pools with branch memory reuse
! NOTE: allocate_memory_pool must be called before using this subroutine
!------------------------------------------------------------------------------

   integer, dimension(:), intent(in) :: atomset1
   type(partition_t), intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(inout) :: perm1

   ! Local variables
   integer :: ntrack, permdiff
   integer, dimension(:), allocatable :: track
   logical, dimension(:), allocatable :: tracked
   integer, dimension(:), allocatable :: perm2
   real(rk) :: permdist
   integer :: i

   ! Variables for random selection
   integer, dimension(:), allocatable :: untracked
   integer :: nuntrack, random_idx, start_atom

   ! Variables
   integer :: n_atoms

   n_atoms = size(adjcs1)

   allocate(perm2(n_atoms))
   allocate(track(n_atoms))
   allocate(tracked(n_atoms))
   allocate(untracked(n_atoms))

   ! Initialization
   ntrack = 0
   tracked(:) = .FALSE.
   perm2 = inverse_permutation(perm1)
   permdiff = adjacencydiff(atomset1, perm1, adjcs1, adjcs2)
!   permdist = sqdistsum(atomset1, perm1, coords1, coords2)

   ! Process all atoms by randomly selecting untracked ones
   ! Each random selection implicitly starts a new fragment
   do while (ntrack < n_atoms)
      ! Build list of untracked atoms
      nuntrack = 0
      do i = 1, n_atoms
         if (.not. tracked(i)) then
            nuntrack = nuntrack + 1
            untracked(nuntrack) = i
         end if
      end do

      ! Pick a random untracked atom
      random_idx = random_uniform_integer(1, nuntrack)
      start_atom = untracked(random_idx)

      ! Process fragment starting from this random atom
      ! Start at depth 1 of the memory pool
      call recur_remap_mismatched_bonds(atomtypes, adjcs1, adjcs2, adjmat2, &
            coords1, coords2, start_atom, perm1, perm2, permdiff, permdist, &
            ntrack, track, tracked, 1)
   end do

   if (DEBUG_TESTS) then
      if (adjacencydiff(atomset1, perm1, adjcs1, adjcs2) /= permdiff) then
         error stop 'permdiff is not equal to the final adjacency difference'
      end if
   end if

   ! Deallocate local arrays
   deallocate(untracked)
   deallocate(track, tracked, perm2)
end subroutine

subroutine match_neighbors(node, adjcs1, adjcs2, perm1, tracked, nmatch, matches, &
                nmismatch1, mismatches1, nmismatch2, mismatches2)
! Classify the atoms connected to node as matches or unmatched
   integer, intent(in) :: node
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer, dimension(:), intent(in) :: perm1
   logical, dimension(:), intent(in) :: tracked
   integer, intent(out) :: nmatch, nmismatch1, nmismatch2
   integer, dimension(:), intent(out) :: matches, mismatches1, mismatches2
   ! Local variables
   integer :: i, mapped_node

   mapped_node = perm1(node)
   nmatch = 0
   nmismatch1 = 0
   nmismatch2 = 0

   ! Classify neighbors of node in structure 1
   do i = 1, adjcs1(node)%cn
      if (any(adjcs2(mapped_node)%list(:adjcs2(mapped_node)%cn) &
            == perm1(adjcs1(node)%list(i)))) then
         nmatch = nmatch + 1
         matches(nmatch) = adjcs1(node)%list(i)
      else
         nmismatch1 = nmismatch1 + 1
         mismatches1(nmismatch1) = adjcs1(node)%list(i)
      end if
   end do

   ! Find neighbors in structure 2 that don't match
   do i = 1, adjcs2(mapped_node)%cn
      if (.not. any(perm1(matches(:nmatch)) == adjcs2(mapped_node)%list(i))) then
         nmismatch2 = nmismatch2 + 1
         mismatches2(nmismatch2) = adjcs2(mapped_node)%list(i)
      end if
   end do
end subroutine

recursive subroutine recur_remap_mismatched_bonds(atomtypes, adjcs1, adjcs2, &
               adjmat2, coords1, coords2, node, perm1, perm2, permdiff, permdist, &
               ntrack, track, tracked, depth)
! Backtracks structure to find assignments that minimize permdiff
! OPTIMIZED: Reuses branch memory at same recursion level to minimize depth usage
   integer, intent(in) :: node
   type(partition_t), intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(inout) :: perm1, perm2
   logical, dimension(:), intent(inout) :: tracked
   integer, intent(inout) :: permdiff, ntrack
   integer, dimension(:), intent(inout) :: track
   real(rk), intent(inout) :: permdist
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, intent(in) :: depth

   ! Local variables
   integer :: nmatch, nmismatch1, nmismatch2
   integer :: branch_permdiff, branch_ntrack
   real(rk) :: branch_permdist
   integer :: i, j

   ! Pointers to memory pool slices for this recursion level
   integer, dimension(:), pointer :: matches, mismatches1, mismatches2
   logical, dimension(:), pointer :: matched1, matched2

   ! Pointers for branch workspace (shared by all branches at this level)
   integer, dimension(:), pointer :: branch_perm1, branch_perm2, branch_track
   logical, dimension(:), pointer :: branch_tracked

   ! Check recursion depth limit
   if (depth > pool_depth) then
      error stop 'Maximum recursion depth exceeded'
   end if

   ! Point to this level's slices in the memory pool
   matches => pool_matches(:, depth)
   mismatches1 => pool_mismatches1(:, depth)
   mismatches2 => pool_mismatches2(:, depth)
   matched1 => pool_matched1(:, depth)
   matched2 => pool_matched2(:, depth)

   ! Reserve node as tracked
   ntrack = ntrack + 1
   track(ntrack) = node
   tracked(node) = .TRUE.

   ! Classify neighbor atoms as matches or mismatched for coords1/coords2
   call match_neighbors(node, adjcs1, adjcs2, perm1, tracked, nmatch, matches, &
               nmismatch1, mismatches1, nmismatch2, mismatches2)

!   print *, "node:", node, "depth:", depth
!   print *, "matches:", matches(:nmatch)
!   print *, "mismatches1:", mismatches1(:nmismatch1)
!   print *, "mismatches2:", mismatches2(:nmismatch2)

   ! Shuffle indices
   call shuffle(matches(:nmatch))
   call shuffle(mismatches1(:nmismatch1))
   call shuffle(mismatches2(:nmismatch2))

   ! Run over matched neighbors
   do i = 1, nmatch
      if (.not. tracked(matches(i))) then
         call recur_remap_mismatched_bonds(atomtypes, adjcs1, adjcs2, adjmat2, &
               coords1, coords2, matches(i), perm1, perm2, permdiff, permdist, &
               ntrack, track, tracked, depth+1)
      end if
   end do

   matched1(:nmismatch1) = .FALSE.
   matched2(:nmismatch2) = .FALSE.

   ! Point to shared branch workspace
   branch_perm1 => pool_perm1(:, depth)
   branch_perm2 => pool_perm2(:, depth)
   branch_track => pool_track(:, depth)
   branch_tracked => pool_tracked(:, depth)

   ! Run over mismatched neighbors - these create branches
   do i = 1, nmismatch1
      if (.not. tracked(mismatches1(i))) then
         do j = 1, nmismatch2
            if (.not. matched2(j)) then
               if (atomtypes%itemdir1(mismatches1(i)) &
                     == atomtypes%itemdir2(mismatches2(j))) then

                  ! Setup branch state
                  branch_ntrack = ntrack
                  branch_track(:) = track(:)
                  branch_tracked(:) = tracked(:)
                  branch_perm1(:) = perm1(:)
                  branch_perm2(:) = perm2(:)

                  ! Apply swap to perm1 branch
                  branch_perm1(mismatches1(i)) = mismatches2(j)
                  branch_perm1(perm2(mismatches2(j))) = perm1(mismatches1(i))

                  ! Apply swap to perm2 branch
                  branch_perm2(mismatches2(j)) = mismatches1(i)
                  branch_perm2(perm1(mismatches1(i))) = perm2(mismatches2(j))

                  ! Update adjd with swap
                  branch_permdiff = permdiff + adjacencydelta(adjcs1, adjmat2, &
                                perm1, mismatches1(i), perm2(mismatches2(j)))

                  ! Update permdist with swap
!                  branch_permdist = permdist + ( &
!                     - sum((coords2(:,perm1(mismatches1(i))) - coords1(:,mismatches1(i)))**2) &
!                     - sum((coords2(:,mismatches2(j)) - coords1(:,perm2(mismatches2(j))))**2) &
!                     + sum((coords2(:,mismatches2(j)) - coords1(:,mismatches1(i)))**2) &
!                     + sum((coords2(:,perm1(mismatches1(i))) - coords1(:,perm2(mismatches2(j))))**2))

                  ! Backtrack swapped index
                  call recur_remap_mismatched_bonds(atomtypes, adjcs1, adjcs2, adjmat2, coords1, &
                        coords2, mismatches1(i), branch_perm1, branch_perm2, branch_permdiff, &
                        branch_permdist, branch_ntrack, branch_track, branch_tracked, depth+1)

                  if (branch_permdiff < permdiff) then
                     ntrack = branch_ntrack
                     track(:) = branch_track(:)
                     tracked(:) = branch_tracked(:)
                     perm1(:) = branch_perm1(:)
                     perm2(:) = branch_perm2(:)
                     permdiff = branch_permdiff
!                     permdist = branch_permdist
                     matched1(i) = .TRUE.
                     matched2(j) = .TRUE.
                     exit   ! exits inner do loop
                  end if
               end if
            end if
         end do
      end if
   end do

   ! Run over non-matched neighbors
   do i = 1, nmismatch1
      if (.not. matched1(i)) then
         if (.not. tracked(mismatches1(i))) then
            call recur_remap_mismatched_bonds(atomtypes, adjcs1, adjcs2, adjmat2, &
                  coords1, coords2, mismatches1(i), perm1, perm2, permdiff, permdist, &
                  ntrack, track, tracked, depth+1)
         end if
      end if
   end do
end subroutine

end module
