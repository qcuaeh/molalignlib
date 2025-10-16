module adjacency
use parameters
use derived_types
use random
use permutation
use euclidean
use chemistry
use molecule
implicit none
private

public adjcs_to_adjmat
public adjacencydiff
public adjacencydelta
public minimize_adjdiff

interface adjacencydiff
   module procedure adjacencydiff_perm
end interface

! Module-level variables used by minimize_adjdiff and its helper procedures
! These are read-only after initialization and shared across recursive calls
integer, dimension(:), allocatable :: blkidx1, blkidx2
integer, dimension(:), allocatable :: eqvidx1, eqvidx2

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

function adjacencydiff_perm(atomset1, atomperm, adjcs1, adjcs2) result(diff)
!------------------------------------------------------------------------------
! Calculate connectivity difference from adjacency lists.
! Returns the number of differing edges.
!------------------------------------------------------------------------------
   integer, dimension(:), intent(in) :: atomset1
   integer, dimension(:), intent(in) :: atomperm
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
      mapped_idx1 = atomperm(idx1)
      nadjs = size(adjcs1(idx1)%adjlist)
      
      do j = 1, nadjs
         neighbor_idx1 = adjcs1(idx1)%adjlist(j)
         
         ! Only count edge if idx1 < neighbor_idx1 to avoid double-counting
         if (idx1 < neighbor_idx1) then
            total_edges1 = total_edges1 + 1
            
            mapped_neighbor_idx1 = atomperm(neighbor_idx1)
            
            ! Check if edge (mapped_idx1, mapped_neighbor_idx1) exists in structure 2
            if (any(adjcs2(mapped_idx1)%adjlist(:) == mapped_neighbor_idx1)) then
               common_edges = common_edges + 1
            end if
         end if
      end do
   end do

   ! Calculate total edges in structure 2
   ! Use atomperm(atomset1) to get the corresponding atoms in molecule 2
   total_edges2 = 0
   do i = 1, size(atomset1)
      idx2 = atomperm(atomset1(i))
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

function adjacencydelta(adjcs1, adjmat2, atomperm, k, l) result(delta)
!------------------------------------------------------------------------------
! Efficiently compute the change in adjacency difference when swapping
! atoms k and l in the permutation. Uses adjacency lists for structure 1 and
! adjacency matrix for structure 2.
!------------------------------------------------------------------------------
   type(adjc_t), dimension(:), intent(in) :: adjcs1
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(in) :: atomperm
   integer, intent(in) :: k, l
   integer :: i, nkk, nkl, nll, nlk, delta, nadjs_k, nadjs_l

   nadjs_k = size(adjcs1(k)%adjlist)
   nadjs_l = size(adjcs1(l)%adjlist)

   nkk = 0
   nkl = 0

   do i = 1, nadjs_k
      if (adjcs1(k)%adjlist(i) /= l) then
         if (adjmat2(atomperm(k), atomperm(adjcs1(k)%adjlist(i)))) nkk = nkk + 1
         if (adjmat2(atomperm(l), atomperm(adjcs1(k)%adjlist(i)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, nadjs_l
      if (adjcs1(l)%adjlist(i) /= k) then
         if (adjmat2(atomperm(l), atomperm(adjcs1(l)%adjlist(i)))) nll = nll + 1
         if (adjmat2(atomperm(k), atomperm(adjcs1(l)%adjlist(i)))) nlk = nlk + 1
      end if
   end do

   ! The change in adjacency difference when swapping k and l:
   ! delta = (new_diff_kl + new_diff_lk) - (old_diff_kk + old_diff_ll)
   ! After simplification: delta = 2*(nkk + nll - nkl - nlk)
   delta = 2*(nkk + nll - nkl - nlk)
end function

subroutine minimize_adjdiff(atomset1, atomtypes, scnatypes, adjcs1, adjcs2, adjmat2, &
                  coords1, coords2, atomperm)
!------------------------------------------------------------------------------
! Find best correspondence between points of graphs
! Randomly selects starting atoms - fragment identification is implicit
!------------------------------------------------------------------------------

   integer, dimension(:), intent(in) :: atomset1
   type(partition_t), intent(in) :: atomtypes, scnatypes
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(inout) :: atomperm

   ! Local variables
   integer :: ntrack, moldiff
   integer, dimension(:), allocatable :: track
   logical, dimension(:), allocatable :: tracked
   integer, dimension(:), allocatable :: invperm
   real(rk) :: moldist
   integer :: num_atoms, i

   ! Variables for random selection
   integer, dimension(:), allocatable :: untracked_atoms
   integer :: num_untracked, random_idx, start_atom

   num_atoms = size(adjcs1)

   allocate(track(num_atoms))
   allocate(tracked(num_atoms))
   allocate(invperm(num_atoms))
   allocate(untracked_atoms(num_atoms))

   ! Set atoms block indices
   blkidx1 = atomtypes%itemdir1
   blkidx2 = atomtypes%itemdir2

   ! Set atoms equivalence indices
   eqvidx1 = scnatypes%itemdir1
   eqvidx2 = scnatypes%itemdir2

   ! Initialization
   ntrack = 0
   tracked(:) = .false.
   invperm = inverse_permutation(atomperm)
   moldiff = adjacencydiff(atomset1, atomperm, adjcs1, adjcs2)
!   moldist = sqdistsum(atomset1, atomperm, coords1, coords2)

   ! Process all atoms by randomly selecting untracked ones
   ! Each random selection implicitly starts a new fragment
   do while (ntrack < num_atoms)
      ! Build list of untracked atoms
      num_untracked = 0
      do i = 1, num_atoms
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
      call recurse_minimize_adjdiff(start_atom, adjcs1, adjcs2, adjmat2, atomperm, &
              invperm, tracked, moldiff, moldist, ntrack, track, coords1, coords2)
   end do

   if (DO_DEBUG_TESTS) then
      if (adjacencydiff(atomset1, atomperm, adjcs1, adjcs2) /= moldiff) then
         error stop 'incorrect edge difference'
      end if
   end if

   ! Deallocate arrays
   deallocate(untracked_atoms)
   deallocate(blkidx1, blkidx2, eqvidx1, eqvidx2)
end subroutine

subroutine match_neighbors(node, adjcs1, adjcs2, atomperm, tracked, nmatch, matches, &
                nmismatch1, mismatches1, nmismatch2, mismatches2)
! Classify the atoms connected to node as matches or unmatched
   integer, intent(in) :: node
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   integer, dimension(:), intent(in) :: atomperm
   logical, dimension(:), intent(in) :: tracked
   integer, intent(out) :: nmatch, nmismatch1, nmismatch2
   integer, dimension(:), intent(out) :: matches, mismatches1, mismatches2
   ! Local variables
   integer :: i, mapped_node

   mapped_node = atomperm(node)
   nmatch = 0
   nmismatch1 = 0
   nmismatch2 = 0

   ! Classify neighbors of node in structure 1
   do i = 1, size(adjcs1(node)%adjlist)
      if (any(adjcs2(mapped_node)%adjlist(:) == atomperm(adjcs1(node)%adjlist(i)))) then
         nmatch = nmatch + 1
         matches(nmatch) = adjcs1(node)%adjlist(i)
      else
         nmismatch1 = nmismatch1 + 1
         mismatches1(nmismatch1) = adjcs1(node)%adjlist(i)
      end if
   end do

   ! Find neighbors in structure 2 that don't match
   do i = 1, size(adjcs2(mapped_node)%adjlist)
      if (.not. any(atomperm(matches(:nmatch)) == adjcs2(mapped_node)%adjlist(i))) then
         nmismatch2 = nmismatch2 + 1
         mismatches2(nmismatch2) = adjcs2(mapped_node)%adjlist(i)
      end if
   end do
end subroutine

recursive subroutine recurse_minimize_adjdiff(node, adjcs1, adjcs2, adjmat2, atomperm, &
                invperm, tracked, moldiff, moldist, ntrack, track, coords1, coords2)
! Backtracks structure to find assignments that minimize moldiff
   integer, intent(in) :: node
   type(adjc_t), dimension(:), intent(in) :: adjcs1, adjcs2
   logical, dimension(:,:), intent(in) :: adjmat2
   integer, dimension(:), intent(inout) :: atomperm, invperm
   logical, dimension(:), intent(inout) :: tracked
   integer, intent(inout) :: moldiff, ntrack
   integer, dimension(:), intent(inout) :: track
   real(rk), intent(inout) :: moldist
   real(rk), dimension(:,:), intent(in) :: coords1, coords2

   ! Local variables
   integer :: nmatch, nmismatch1, nmismatch2
   integer, dimension(:), allocatable :: matches, mismatches1, mismatches2
   integer, dimension(:), allocatable :: mapping_branch, unmapping_branch, track_branch
   integer :: moldiff_branch, ntrack_branch
   logical, dimension(:), allocatable :: tracked_branch, matched1, matched2
   real(rk) :: moldist_branch
   integer :: i, j, num_atoms

   num_atoms = size(adjcs1)

   allocate(matches(num_atoms))
   allocate(mismatches1(num_atoms))
   allocate(mismatches2(num_atoms))
   allocate(mapping_branch(num_atoms))
   allocate(unmapping_branch(num_atoms))
   allocate(track_branch(num_atoms))
   allocate(tracked_branch(num_atoms))
   allocate(matched1(num_atoms))
   allocate(matched2(num_atoms))

   ! Reserve node as tracked
   ntrack = ntrack + 1
   track(ntrack) = node
   tracked(node) = .true.

   ! Classify neighbor atoms as matches or mismatched for coords1/coords2
   call match_neighbors(node, adjcs1, adjcs2, atomperm, tracked, nmatch, matches, &
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
         call recurse_minimize_adjdiff(matches(i), adjcs1, adjcs2, adjmat2, atomperm, &
                 invperm, tracked, moldiff, moldist, ntrack, track, coords1, coords2)
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
                  mapping_branch(:) = atomperm(:)
                  unmapping_branch(:) = invperm(:)

                  ! Apply swap to atomperm branch
                  mapping_branch(mismatches1(i)) = mismatches2(j)
                  mapping_branch(invperm(mismatches2(j))) = atomperm(mismatches1(i))

                  ! Apply swap to invperm branch
                  unmapping_branch(mismatches2(j)) = mismatches1(i)
                  unmapping_branch(atomperm(mismatches1(i))) = invperm(mismatches2(j))

                  ! Update adjd with swap
                  moldiff_branch = moldiff + adjacencydelta(adjcs1, adjmat2, &
                                atomperm, mismatches1(i), invperm(mismatches2(j)))

                  ! Update ssd with swap
!                  moldist_branch = moldist + ( &
!                     - sum((coords2(:, atomperm(mismatches1(i))) - coords1(:, mismatches1(i)))**2) &
!                     - sum((coords2(:, mismatches2(j)) - coords1(:, invperm(mismatches2(j))))**2) &
!                     + sum((coords2(:, mismatches2(j)) - coords1(:, mismatches1(i)))**2) &
!                     + sum((coords2(:, atomperm(mismatches1(i))) - coords1(:, invperm(mismatches2(j))))**2))

                  ! Backtrack swapped index
                  call recurse_minimize_adjdiff(mismatches1(i), adjcs1, adjcs2, adjmat2, mapping_branch, &
                     unmapping_branch, tracked_branch, moldiff_branch, moldist_branch, ntrack_branch, &
                     track_branch, coords1, coords2)

                  if ( &
                     moldiff_branch < moldiff &
                     .and. ( &
                        eqvidx1(mismatches1(i)) == eqvidx1(invperm(mismatches2(j))) &
                        .and. eqvidx2(atomperm(mismatches1(i))) == eqvidx2(mismatches2(j)) &
                     ) &
                  ) then
                     ntrack = ntrack_branch
                     track(:) = track_branch(:)
                     tracked(:) = tracked_branch(:)
                     atomperm(:) = mapping_branch(:)
                     invperm(:) = unmapping_branch(:)
                     moldiff = moldiff_branch
!                     moldist = moldist_branch
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
            call recurse_minimize_adjdiff(mismatches1(i), adjcs1, adjcs2, adjmat2, atomperm, &
                    invperm, tracked, moldiff, moldist, ntrack, track, coords1, coords2)
         end if
      end if
   end do
end subroutine

end module
