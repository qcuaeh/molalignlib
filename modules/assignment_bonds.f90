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
public remap_mismatched_bonds

! Memory pool for recursive work arrays
! Organized as (size, depth_level) to eliminate allocations in recursion
integer, dimension(:,:), allocatable, target :: pool_matches, pool_mismatches1, pool_mismatches2
integer, dimension(:,:), allocatable, target :: pool_mapping, pool_unmapping, pool_track
logical, dimension(:,:), allocatable, target :: pool_tracked, pool_matched1, pool_matched2

contains

subroutine remap_mismatched_bonds(atomset1, atomtypes, adjcs1, adjcs2, adjmat2, &
                  coords1, coords2, atomperm1)
!------------------------------------------------------------------------------
! Find best correspondence between points of graphs
! Randomly selects starting atoms - fragment identification is implicit
! OPTIMIZED: Uses pre-allocated memory pools to eliminate recursive allocations
!------------------------------------------------------------------------------

   integer, dimension(:), intent(in) :: atomset1
   type(partition_t), intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(inout) :: atomperm1

   ! Local variables
   integer :: ntrack, permdiff
   integer, dimension(:), allocatable :: track
   logical, dimension(:), allocatable :: tracked
   integer, dimension(:), allocatable :: atomperm2
   real(rk) :: permdist
   integer :: i

   ! Variables for random selection
   integer, dimension(:), allocatable :: untracked_atoms
   integer :: num_untracked, random_idx, start_atom

   ! Memory pool configuration
   integer :: natoms, max_depth

   natoms = size(atomperm1)

   ! Estimate maximum recursion depth
   ! Conservative estimate: assume worst case is full depth tree
   ! In practice, depth is typically much smaller (connected components, branching)
   ! A reasonable heuristic: natoms (worst case) but could be tuned lower
   max_depth = natoms

   ! Allocate memory pools once for all recursive calls
   ! Each "column" is used by one recursion level
   allocate(pool_matches(natoms, max_depth))
   allocate(pool_mismatches1(natoms, max_depth))
   allocate(pool_mismatches2(natoms, max_depth))
   allocate(pool_mapping(natoms, max_depth))
   allocate(pool_unmapping(natoms, max_depth))
   allocate(pool_track(natoms, max_depth))
   allocate(pool_tracked(natoms, max_depth))
   allocate(pool_matched1(natoms, max_depth))
   allocate(pool_matched2(natoms, max_depth))

   allocate(track(natoms))
   allocate(tracked(natoms))
   allocate(atomperm2(natoms))
   allocate(untracked_atoms(natoms))

   ! Initialization
   ntrack = 0
   tracked(:) = .false.
   atomperm2 = inverse_permutation(atomperm1)
   permdiff = adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2)
!   permdist = sqdistsum(atomset1, atomperm1, coords1, coords2)

   ! Process all atoms by randomly selecting untracked ones
   ! Each random selection implicitly starts a new fragment
   do while (ntrack < natoms)
      ! Build list of untracked atoms
      num_untracked = 0
      do i = 1, natoms
         if (.not. tracked(i)) then
            num_untracked = num_untracked + 1
            untracked_atoms(num_untracked) = i
         end if
      end do

      ! Pick a random untracked atom
      random_idx = random_uniform_integer(1, num_untracked)
      start_atom = untracked_atoms(random_idx)

      ! Process fragment starting from this random atom
      ! Start at depth 1 of the memory pool
      call recurse_remap_mismatched_bonds(start_atom, atomtypes, adjcs1, adjcs2, adjmat2, &
            atomperm1, atomperm2, tracked, permdiff, permdist, ntrack, track, coords1, &
            coords2, 1, max_depth)
   end do

   if (DEBUG_TESTS) then
      if (adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2) /= permdiff) then
         error stop 'permdiff is not equal to the final adjacency difference'
      end if
   end if

   ! Deallocate memory pools
   deallocate(pool_matches, pool_mismatches1, pool_mismatches2)
   deallocate(pool_mapping, pool_unmapping, pool_track)
   deallocate(pool_tracked, pool_matched1, pool_matched2)

   ! Deallocate other arrays
   deallocate(untracked_atoms)
end subroutine

subroutine match_neighbors(node, adjcs1, adjcs2, atomperm1, tracked, nmatch, matches, &
                nmismatch1, mismatches1, nmismatch2, mismatches2)
! Classify the atoms connected to node as matches or unmatched
   integer, intent(in) :: node
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer, dimension(:), intent(in) :: atomperm1
   logical, dimension(:), intent(in) :: tracked
   integer, intent(out) :: nmatch, nmismatch1, nmismatch2
   integer, dimension(:), intent(out) :: matches, mismatches1, mismatches2
   ! Local variables
   integer :: i, mapped_node

   mapped_node = atomperm1(node)
   nmatch = 0
   nmismatch1 = 0
   nmismatch2 = 0

   ! Classify neighbors of node in structure 1
   do i = 1, adjcs1(node)%cn
      if (any(adjcs2(mapped_node)%list(:adjcs2(mapped_node)%cn) &
            == atomperm1(adjcs1(node)%list(i)))) then
         nmatch = nmatch + 1
         matches(nmatch) = adjcs1(node)%list(i)
      else
         nmismatch1 = nmismatch1 + 1
         mismatches1(nmismatch1) = adjcs1(node)%list(i)
      end if
   end do

   ! Find neighbors in structure 2 that don't match
   do i = 1, adjcs2(mapped_node)%cn
      if (.not. any(atomperm1(matches(:nmatch)) == adjcs2(mapped_node)%list(i))) then
         nmismatch2 = nmismatch2 + 1
         mismatches2(nmismatch2) = adjcs2(mapped_node)%list(i)
      end if
   end do
end subroutine

recursive subroutine recurse_remap_mismatched_bonds(node, atomtypes, adjcs1, adjcs2, &
      adjmat2, atomperm1, atomperm2, tracked, permdiff, permdist, ntrack, track, &
      coords1, coords2, depth, max_depth)
! Backtracks structure to find assignments that minimize permdiff
! OPTIMIZED: Uses memory pool slices based on recursion depth to avoid allocations
   integer, intent(in) :: node
   type(partition_t), intent(in) :: atomtypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(inout) :: atomperm1, atomperm2
   logical, dimension(:), intent(inout) :: tracked
   integer, intent(inout) :: permdiff, ntrack
   integer, dimension(:), intent(inout) :: track
   real(rk), intent(inout) :: permdist
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, intent(in) :: depth, max_depth

   ! Local variables
   integer :: nmatch, nmismatch1, nmismatch2
   integer :: moldiff_branch, ntrack_branch
   real(rk) :: moldist_branch
   integer :: i, j
   integer :: natoms

   ! Pointers to memory pool slices for this recursion level
   integer, dimension(:), pointer :: matches, mismatches1, mismatches2
   logical, dimension(:), pointer :: matched1, matched2

   ! For branching, we need the next depth level
   integer, dimension(:), pointer :: mapping_branch, unmapping_branch, track_branch
   logical, dimension(:), pointer :: tracked_branch
   integer :: next_depth

   natoms = size(atomperm1)

   ! Check recursion depth limit
   if (depth > max_depth) then
      error stop 'Maximum recursion depth exceeded in remap_mismatched_bonds'
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
   tracked(node) = .true.

   ! Classify neighbor atoms as matches or mismatched for coords1/coords2
   call match_neighbors(node, adjcs1, adjcs2, atomperm1, tracked, nmatch, matches, &
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
         call recurse_remap_mismatched_bonds(matches(i), atomtypes, adjcs1, adjcs2, adjmat2, &
               atomperm1, atomperm2, tracked, permdiff, permdist, ntrack, track, coords1, &
               coords2, depth + 1, max_depth)
      end if
   end do

   matched1(:nmismatch1) = .false.
   matched2(:nmismatch2) = .false.

   ! Run over mismatched neighbors - these create branches
   do i = 1, nmismatch1
      if (.not. tracked(mismatches1(i))) then
         do j = 1, nmismatch2
            if (.not. matched2(j)) then
               if (atomtypes%itemdir1(mismatches1(i)) &
                     == atomtypes%itemdir2(mismatches2(j))) then

                  ! Use next depth level for branch arrays
                  next_depth = depth + 1
                  if (next_depth > max_depth) then
                     error stop 'Maximum recursion depth exceeded in branching'
                  end if

                  mapping_branch => pool_mapping(:, next_depth)
                  unmapping_branch => pool_unmapping(:, next_depth)
                  track_branch => pool_track(:, next_depth)
                  tracked_branch => pool_tracked(:, next_depth)

                  ntrack_branch = ntrack
                  track_branch(:) = track(:)
                  tracked_branch(:) = tracked(:)
                  mapping_branch(:) = atomperm1(:)
                  unmapping_branch(:) = atomperm2(:)

                  ! Apply swap to atomperm1 branch
                  mapping_branch(mismatches1(i)) = mismatches2(j)
                  mapping_branch(atomperm2(mismatches2(j))) = atomperm1(mismatches1(i))

                  ! Apply swap to atomperm2 branch
                  unmapping_branch(mismatches2(j)) = mismatches1(i)
                  unmapping_branch(atomperm1(mismatches1(i))) = atomperm2(mismatches2(j))

                  ! Update adjd with swap
                  moldiff_branch = permdiff + adjacencydelta(adjcs1, adjmat2, &
                                atomperm1, mismatches1(i), atomperm2(mismatches2(j)))

                  ! Update ssd with swap
!                  moldist_branch = permdist + ( &
!                     - sum((coords2(:,atomperm1(mismatches1(i))) - coords1(:,mismatches1(i)))**2) &
!                     - sum((coords2(:,mismatches2(j)) - coords1(:,atomperm2(mismatches2(j))))**2) &
!                     + sum((coords2(:,mismatches2(j)) - coords1(:,mismatches1(i)))**2) &
!                     + sum((coords2(:,atomperm1(mismatches1(i))) - coords1(:,atomperm2(mismatches2(j))))**2))

                  ! Backtrack swapped index
                  ! Branch recursion uses next_depth + 1 since we've already used next_depth for branch state
                  call recurse_remap_mismatched_bonds(mismatches1(i), atomtypes, adjcs1, adjcs2, adjmat2, &
                        mapping_branch, unmapping_branch, tracked_branch, moldiff_branch, moldist_branch, &
                        ntrack_branch, track_branch, coords1, coords2, next_depth + 1, max_depth)

                  if (moldiff_branch < permdiff) then
                     ntrack = ntrack_branch
                     track(:) = track_branch(:)
                     tracked(:) = tracked_branch(:)
                     atomperm1(:) = mapping_branch(:)
                     atomperm2(:) = unmapping_branch(:)
                     permdiff = moldiff_branch
!                     permdist = moldist_branch
                     matched1(i) = .true.
                     matched2(j) = .true.
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
            call recurse_remap_mismatched_bonds(mismatches1(i), atomtypes, adjcs1, adjcs2, &
                  adjmat2, atomperm1, atomperm2, tracked, permdiff, permdist, ntrack, track, &
                  coords1, coords2, depth + 1, max_depth)
         end if
      end if
   end do

   ! No deallocation needed - memory pool is managed by caller

end subroutine

end module
