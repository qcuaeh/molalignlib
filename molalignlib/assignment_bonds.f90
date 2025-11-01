module assignment_bonds
use parameters
use derived_types
use permutation
use euclidean
use chemistry
use adjacency
use random
implicit none
private

public minimize_adjdiff

! Module-level variables used by minimize_adjdiff and its helper procedures
! These are read-only after initialization and shared across recursive calls
integer, dimension(:), allocatable :: blkidx1, blkidx2
integer, dimension(:), allocatable :: eqvidx1, eqvidx2

contains

subroutine minimize_adjdiff(atomset1, atomtypes, scnatypes, adjcs1, adjcs2, adjmat2, &
                  coords1, coords2, atomperm1)
!------------------------------------------------------------------------------
! Find best correspondence between points of graphs
! Randomly selects starting atoms - fragment identification is implicit
!------------------------------------------------------------------------------

   integer, dimension(:), intent(in) :: atomset1
   type(partition_t), intent(in) :: atomtypes, scnatypes
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
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

   allocate(track(size(atomperm1)))
   allocate(tracked(size(atomperm1)))
   allocate(atomperm2(size(atomperm1)))
   allocate(untracked_atoms(size(atomperm1)))

   ! Set atoms block indices
   blkidx1 = atomtypes%itemdir1
   blkidx2 = atomtypes%itemdir2

   ! Set atoms equivalence indices
   eqvidx1 = scnatypes%itemdir1
   eqvidx2 = scnatypes%itemdir2

   ! Initialization
   ntrack = 0
   tracked(:) = .false.
   atomperm2 = inverse_permutation(atomperm1)
   permdiff = adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2)
!   permdist = sqdistsum(atomset1, atomperm1, coords1, coords2)

   ! Process all atoms by randomly selecting untracked ones
   ! Each random selection implicitly starts a new fragment
   do while (ntrack < size(atomperm1))
      ! Build list of untracked atoms
      num_untracked = 0
      do i = 1, size(atomperm1)
         if (.not. tracked(i)) then
            num_untracked = num_untracked + 1
            untracked_atoms(num_untracked) = i
         end if
      end do

      ! Pick a random untracked atom
      random_idx = random_uniform_integer(1, num_untracked)
      start_atom = untracked_atoms(random_idx)

      ! Process fragment starting from this random atom
      ! Recursion naturally explores the entire connected component
      call recurse_minimize_adjdiff(start_atom, adjcs1, adjcs2, adjmat2, atomperm1, &
              atomperm2, tracked, permdiff, permdist, ntrack, track, coords1, coords2)
   end do

   if (DEBUG_TESTS) then
      if (adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2) /= permdiff) then
         error stop 'permdiff is not equal to the final adjacency difference'
      end if
   end if

   ! Deallocate arrays
   deallocate(untracked_atoms)
   deallocate(blkidx1, blkidx2, eqvidx1, eqvidx2)
end subroutine

subroutine match_neighbors(node, adjcs1, adjcs2, atomperm1, tracked, nmatch, matches, &
                nmismatch1, mismatches1, nmismatch2, mismatches2)
! Classify the atoms connected to node as matches or unmatched
   integer, intent(in) :: node
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
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
   do i = 1, adjcs1%cns(node)
      if (any(adjcs2%lists(1:adjcs2%cns(mapped_node), mapped_node) == atomperm1(adjcs1%lists(i, node)))) then
         nmatch = nmatch + 1
         matches(nmatch) = adjcs1%lists(i, node)
      else
         nmismatch1 = nmismatch1 + 1
         mismatches1(nmismatch1) = adjcs1%lists(i, node)
      end if
   end do

   ! Find neighbors in structure 2 that don't match
   do i = 1, adjcs2%cns(mapped_node)
      if (.not. any(atomperm1(matches(:nmatch)) == adjcs2%lists(i, mapped_node))) then
         nmismatch2 = nmismatch2 + 1
         mismatches2(nmismatch2) = adjcs2%lists(i, mapped_node)
      end if
   end do
end subroutine

recursive subroutine recurse_minimize_adjdiff(node, adjcs1, adjcs2, adjmat2, atomperm1, &
                atomperm2, tracked, permdiff, permdist, ntrack, track, coords1, coords2)
! Backtracks structure to find assignments that minimize permdiff
   integer, intent(in) :: node
   type(adjcs_t), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(inout) :: atomperm1, atomperm2
   logical, dimension(:), intent(inout) :: tracked
   integer, intent(inout) :: permdiff, ntrack
   integer, dimension(:), intent(inout) :: track
   real(rk), intent(inout) :: permdist
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   ! Local variables
   integer :: nmatch, nmismatch1, nmismatch2
   integer, dimension(:), allocatable :: matches, mismatches1, mismatches2
   integer, dimension(:), allocatable :: mapping_branch, unmapping_branch, track_branch
   integer :: moldiff_branch, ntrack_branch
   logical, dimension(:), allocatable :: tracked_branch, matched1, matched2
   real(rk) :: moldist_branch
   integer :: i, j

   allocate(matches(size(atomperm1)))
   allocate(mismatches1(size(atomperm1)))
   allocate(mismatches2(size(atomperm1)))
   allocate(mapping_branch(size(atomperm1)))
   allocate(unmapping_branch(size(atomperm1)))
   allocate(track_branch(size(atomperm1)))
   allocate(tracked_branch(size(atomperm1)))
   allocate(matched1(size(atomperm1)))
   allocate(matched2(size(atomperm1)))

   ! Reserve node as tracked
   ntrack = ntrack + 1
   track(ntrack) = node
   tracked(node) = .true.

   ! Classify neighbor atoms as matches or mismatched for coords1/coords2
   call match_neighbors(node, adjcs1, adjcs2, atomperm1, tracked, nmatch, matches, &
               nmismatch1, mismatches1, nmismatch2, mismatches2)

!   print *, "node:", node
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
         call recurse_minimize_adjdiff(matches(i), adjcs1, adjcs2, adjmat2, atomperm1, &
                 atomperm2, tracked, permdiff, permdist, ntrack, track, coords1, coords2)
      end if
   end do

   matched1(:nmismatch1) = .false.
   matched2(:nmismatch2) = .false.

   ! Run over mismatched neighbors
   do i = 1, nmismatch1
      if (.not. tracked(mismatches1(i))) then
         do j = 1, nmismatch2
            if (.not. matched2(j)) then
               if (blkidx1(mismatches1(i)) == blkidx2(mismatches2(j))) then

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
!                     - sum((coords2(:, atomperm1(mismatches1(i))) - coords1(:, mismatches1(i)))**2) &
!                     - sum((coords2(:, mismatches2(j)) - coords1(:, atomperm2(mismatches2(j))))**2) &
!                     + sum((coords2(:, mismatches2(j)) - coords1(:, mismatches1(i)))**2) &
!                     + sum((coords2(:, atomperm1(mismatches1(i))) - coords1(:, atomperm2(mismatches2(j))))**2))

                  ! Backtrack swapped index
                  call recurse_minimize_adjdiff(mismatches1(i), adjcs1, adjcs2, adjmat2, mapping_branch, &
                     unmapping_branch, tracked_branch, moldiff_branch, moldist_branch, ntrack_branch, &
                     track_branch, coords1, coords2)

                  if ( &
                     moldiff_branch < permdiff &
                     .and. ( &
                        eqvidx1(mismatches1(i)) == eqvidx1(atomperm2(mismatches2(j))) &
                        .and. eqvidx2(atomperm1(mismatches1(i))) == eqvidx2(mismatches2(j)) &
                     ) &
                  ) then
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
            call recurse_minimize_adjdiff(mismatches1(i), adjcs1, adjcs2, adjmat2, atomperm1, &
                    atomperm2, tracked, permdiff, permdist, ntrack, track, coords1, coords2)
         end if
      end if
   end do
end subroutine

end module
