module tracking
! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragidcs)
use parameters
use random
use sorting
use permutation
use adjacency
use spatial
use lcrs_tree
use molecule
use basetypes
implicit none
private

type(atom_type), dimension(:), allocatable :: atoms

public minadjdiff
public find_molfrags

contains

subroutine find_molfrags( mol, eltypes, molfrags)
   type(mol_type), intent(in) :: mol
   type(semipartitionarray_t), intent(in) :: eltypes
   type(int_list), dimension(:), allocatable, intent(out) :: molfrags
   ! Local variables
   integer :: i, nfrag
   logical, dimension(:), allocatable :: tracked
   integer, dimension(:), allocatable :: fragszs
   integer, dimension(:,:), allocatable :: fragidcs
   integer, dimension(:), allocatable :: order

   allocate (fragszs(size(mol%atoms)))
   allocate (fragidcs(size(mol%atoms), size(mol%atoms)))
   allocate (tracked(size(mol%atoms)))
   atoms = mol%atoms

   ! initialization

   nfrag = 0
   fragszs(:) = 0
   tracked(:) = .false.

   ! detect fragments and populate frag arrays
   i = 1
   do while (i <= size(mol%atoms))
      if (tracked(i)) then
         i = i + 1
      else
         nfrag = nfrag + 1
         call recrun( tracked, i, nfrag, fragszs, fragidcs)
         i = 1
      end if
   end do

   ! Order molecular fragments
   do i = 1, nfrag
      order = sorted_order(eltypes%parts(eltypes%itemdir(fragidcs(:fragszs(i), i)))%num_items)
      fragidcs(:fragszs(i), i) = fragidcs(order, i)
!      write (stderr, *) fragidcs1(:fragszs1(i), i)
!      write (stderr, *)
   end do

   allocate (molfrags(nfrag))

   do i = 1, nfrag
      molfrags(i)%e = fragidcs(:fragszs(i), i)
   end do

end subroutine

recursive subroutine recrun( tracked, iatom, nfrag, fragszs, fragidcs)
! runs recursivelly over the structure and populates arrays
   logical, dimension(:), intent(inout) :: tracked
   integer, intent(in) :: iatom, nfrag
   integer, dimension(:), intent(inout) :: fragszs
   integer, dimension(:,:), intent(inout) :: fragidcs
   ! Local variables
   integer :: i

   if (tracked(iatom)) return

   tracked(iatom) = .true.
   fragszs(nfrag) = fragszs(nfrag) + 1
   fragidcs(fragszs(nfrag), nfrag) = iatom

   do i = 1, size(atoms(iatom)%adjlist)
      call recrun( tracked, atoms(iatom)%adjlist(i), nfrag, fragszs, fragidcs)
   end do

end subroutine

subroutine minadjdiff( eltypes, mnatypes, molfrags, mol1, mol2, coords1, &
                       coords2, atomperm)
!
! Find best correspondence between points of graphs
!

   ! Arguments
   type(partitionarray_t), intent(in) :: eltypes, mnatypes
   type(int_list), dimension(:), intent(in) :: molfrags
   type(mol_type), intent(in) :: mol1, mol2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(inout) :: atomperm

   ! Local variables
   integer :: num_atoms1, num_atoms2
   integer, dimension(:), allocatable :: nadjs1, nadjs2
   integer, dimension(:,:), allocatable :: adjlists1, adjlists2
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: blkidx1, blkidx2
   integer, dimension(:), allocatable :: eqvidx1, eqvidx2
   integer :: ntrack, moldiff
   integer, dimension(:), allocatable :: track
   logical, dimension(:), allocatable :: tracked
   integer, dimension(:), allocatable :: invatomperm
   logical, parameter :: print_info = .false.
   real(rk) :: moldist
   integer :: i

   num_atoms1 = size(mol1%atoms)
   num_atoms2 = size(mol2%atoms)

   allocate (nadjs1(num_atoms1))
   allocate (nadjs2(num_atoms2))
   allocate (adjlists1(num_atoms1, num_atoms1))
   allocate (adjlists2(num_atoms2, num_atoms1))
   allocate (track(num_atoms1))
   allocate (tracked(num_atoms1))
   allocate (invatomperm(num_atoms1))

   do i = 1, num_atoms1
      nadjs1(i) = size(mol1%atoms(i)%adjlist)
      adjlists1(:nadjs1(i), i) = mol1%atoms(i)%adjlist
   end do

   do i = 1, num_atoms2
      nadjs2(i) = size(mol2%atoms(i)%adjlist)
      adjlists2(:nadjs2(i), i) = mol2%atoms(i)%adjlist
   end do

   adjmat1 = get_adjmat(mol1)
   adjmat2 = get_adjmat(mol2)

   ! set atoms block indices

   blkidx1 = eltypes%itemdir1
   blkidx2 = eltypes%itemdir2

   ! set atoms equivalence indices

   eqvidx1 = mnatypes%itemdir1
   eqvidx2 = mnatypes%itemdir2

   !  initialization

   ntrack = 0
   tracked(:) = .false.
   invatomperm = inverse_perm( atomperm)
   moldiff = adjacencydiff( atomperm, adjmat1, adjmat2)
   moldist = totsqdist( atomperm, coords1, coords2)

   if ( print_info ) then
      print '(a,1x,i0)', "moldiff:", moldiff
      print '(a,1x,f0.4)', "moldist:", moldist
   end if

   do i = 1, size(molfrags)
      call recursive_backtrack( molfrags(i)%e(1), atomperm, invatomperm, tracked, &
         moldiff, moldist, ntrack, track)
!        print *, ntrack
   end do

   if ( print_info ) then
      print '(a,1x,i0)', "Fragments:", size(molfrags)
      print '(a,1x,i0,1x,i0)', "moldiff:", adjacencydiff( atomperm, adjmat1, adjmat2), moldiff
      print '(a,1x,f0.4,1x,f0.4)', "moldist:", totsqdist( atomperm, coords1, coords2), moldist
   end if

!    if (adjacencydiff( atomperm, adjmat1, adjmat2) /= moldiff) then
!        print '(a,x,i0,x,i0)', "moldiff:", adjacencydiff( atomperm, adjmat1, adjmat2), moldiff
!    end if

   contains    

! classify the atoms connected to node as matches or unmatched
   subroutine nodematch( node, atomperm, tracked, nmatch, matches, &
                   nmismatch1, mismatches1, nmismatch2, mismatches2)
      integer, intent(in) :: node
      integer, dimension(:), intent(in) :: atomperm
      logical, dimension(:), intent(in) :: tracked
      integer, intent(out) :: nmatch, nmismatch1, nmismatch2
      integer, dimension(:), intent(out) :: matches, mismatches1, mismatches2
      ! Local variables
      integer, allocatable, dimension(:) :: adj1, adj2
      integer :: i

      allocate (adj1(num_atoms1))
      allocate (adj2(num_atoms1))

      adj1(:nadjs1(node)) = adjlists1(:nadjs1(node), node)
      adj2(:nadjs2(atomperm(node))) = adjlists2(:nadjs2(atomperm(node)), atomperm(node))
      nmatch = 0
      nmismatch1 = 0
      nmismatch2 = 0

      do i = 1, nadjs1(node)
         if (find_index(atomperm(adj1(i)), nadjs2(atomperm(node)), adj2) /= 0) then
            call add_index(adj1(i), nmatch, matches)
         else
            call add_index(adj1(i), nmismatch1, mismatches1)
         end if
      end do

      do i = 1, nadjs2(atomperm(node))
         if (find_index(adj2(i), nmatch, atomperm(matches(:nmatch))) == 0) then
            call add_index(adj2(i), nmismatch2, mismatches2)
         end if
      end do
   end subroutine nodematch

! backtracks structure to find assignments that minimize moldiff
   recursive subroutine recursive_backtrack( node, atomperm, invatomperm, tracked, moldiff, moldist, &
                                  ntrack, track)
      integer, intent(in) :: node
      integer, dimension(:), intent(inout) :: atomperm, invatomperm
      logical, dimension(:), intent(inout) :: tracked
      integer, intent(inout) :: moldiff, ntrack
      integer, dimension(:), intent(inout) :: track
      real(rk), intent(inout) :: moldist

      logical, parameter :: print_info = .false.

      integer nmatch, matches(num_atoms1)
      integer nmismatch1, mismatches1(num_atoms1)
      integer nmismatch2, mismatches2(num_atoms1)
      integer mapping_branch(num_atoms1), moldiff_branch, unmapping_branch(num_atoms1)
      integer ntrack_branch, track_branch(num_atoms1)
      logical tracked_branch(num_atoms1)
      logical :: matched1(num_atoms1), matched2(num_atoms1)
      real(rk) :: moldist_branch
      integer :: i, j

      ! reserve node node as tracked
      ntrack = ntrack + 1
      track(ntrack) = node
      tracked(node) = .true.

      ! classify neighbor atoms as matches or mismatched for coords1/coords2
      call nodematch( node, atomperm, tracked, nmatch, matches, &
                  nmismatch1, mismatches1, nmismatch2, mismatches2)

      if (print_info) then
         print *, "node:", node
         print *, "matches:", matches(:nmatch)
         print *, "mismatches1:", mismatches1(:nmismatch1)
         print *, "mismatches2:", mismatches2(:nmismatch2)
      end if

      ! shuffle indices
      call shuffle (matches(:nmatch))
      call shuffle (mismatches1(:nmismatch1))
      call shuffle (mismatches2(:nmismatch2))

      ! run over matches neighbors
      do i = 1, nmatch
         if (.not. tracked(matches(i))) then
            call recursive_backtrack(matches(i), atomperm, invatomperm, tracked, moldiff, &
                            moldist, ntrack, track)
         end if
      end do

      matched1(:nmismatch1) = .false.
      matched2(:nmismatch2) = .false.

      ! run over mismatched neighbors
      do i = 1, nmismatch1
         if (.not. tracked(mismatches1(i))) then
            do j = 1, nmismatch2
               if (.not. matched2(j)) then
                  if (blkidx1(mismatches1(i)) == blkidx2(mismatches2(j))) then

                     ntrack_branch = ntrack
                     track_branch(:) = track(:)
                     tracked_branch(:) = tracked(:)
                     mapping_branch(:) = atomperm(:)
                     unmapping_branch(:) = invatomperm(:)

                     ! Apply swap to atomperm branch 
                     mapping_branch(mismatches1(i)) = mismatches2(j)
                     mapping_branch(invatomperm(mismatches2(j))) = atomperm(mismatches1(i))

                     ! Apply swap to invatomperm branch 
                     unmapping_branch(mismatches2(j)) = mismatches1(i)
                     unmapping_branch(atomperm(mismatches1(i))) = invatomperm(mismatches2(j))

                     ! Update ssd with swap
                     moldist_branch = moldist + ( &
                        - sum((coords2(:, atomperm(mismatches1(i))) - coords1(:, mismatches1(i)))**2) &
                        - sum((coords2(:, mismatches2(j)) - coords1(:, invatomperm(mismatches2(j))))**2) &
                        + sum((coords2(:, mismatches2(j)) - coords1(:, mismatches1(i)))**2) &
                        + sum((coords2(:, atomperm(mismatches1(i))) - coords1(:, invatomperm(mismatches2(j))))**2))

                     ! Update adjd with swap
                     moldiff_branch = moldiff + adjacencydelta(nadjs1, adjlists1, adjmat2, &
                                   atomperm, mismatches1(i), invatomperm(mismatches2(j)))

                     ! backtrack swapped index
                     call recursive_backtrack(mismatches1(i), mapping_branch, unmapping_branch, &
                        tracked_branch, moldiff_branch, moldist_branch, ntrack_branch, track_branch)

                     if ( &
                        moldiff_branch < moldiff &
                        .and. ( &
                           eqvidx1(mismatches1(i)) == eqvidx1(invatomperm(mismatches2(j))) &
                           .and. eqvidx2(atomperm(mismatches1(i))) == eqvidx2(mismatches2(j)) &
                        ) &
                     ) then
                        ntrack = ntrack_branch
                        track(:) = track_branch(:)
                        tracked(:) = tracked_branch(:)
                        atomperm(:) = mapping_branch(:)
                        invatomperm(:) = unmapping_branch(:)
                        moldiff = moldiff_branch
                        moldist = moldist_branch
                        matched1(i) = .true.
                        matched2(j) = .true.
                        exit   ! exits inner do loop
                     end if
                  end if
               end if
            end do
         end if
      end do

      ! run over non matches neighbors
      do i = 1, nmismatch1
         if (.not. matched1(i)) then
            if (.not. tracked(mismatches1(i))) then
               call recursive_backtrack(mismatches1(i), atomperm, invatomperm, &
                 tracked, moldiff, moldist, ntrack, track)
            end if
         end if
      end do

   end subroutine

end subroutine

! Adds a new element to a list of indices in growing order, no repeats
subroutine add_index( index, nind, indices)
   integer, intent(in) :: index
   integer, intent(inout) :: nind
   integer, dimension(:), intent(inout) :: indices
   integer :: p, e
   p = 1
   do while (p <= nind)   ! find position for growing order
      if (index == indices(p)) return
      if (index < indices(p)) then
         exit
      else
         p = p + 1
      end if
   end do
   nind = nind + 1
   e = nind
   do while (e > p)
      indices(e) = indices(e - 1)   ! shifts one position from the end
      e = e - 1
   end do
   indices(p) = index
end subroutine add_index

! Returns the position in indices where the element index is found or 0 otherwise
function find_index( index, nind, indices) result(pos)
   integer, intent(in) :: index, nind, indices(nind)
   integer :: pos
   integer :: i
   pos = 0
   do i = 1, nind
      if (index == indices(i)) then
         pos = i
         return
      end if
   end do
end function find_index

end module
