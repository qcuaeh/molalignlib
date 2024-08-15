! template for module backtracking in molalign program
module backtracking
use stdio
use sorting
use discrete
use adjacency
use alignment
use random
use subset

implicit none
private
public minadjdiff
public eqvatomperm

contains

subroutine minadjdiff (natom, weights, natomtype, atomtypelenlist, coords0, nadjs0, adjlists0, adjmat0, natomequiv0, &
   atomequivlenlist0, coords1, nadjs1, adjlists1, adjmat1, natomequiv1, atomequivlenlist1, mapping, nfrag0, fragroot0)
! Purpose: Find best correspondence between points of graphs

   integer, intent(in) :: natom, natomtype, natomequiv0, natomequiv1
   integer, dimension(:), intent(in) :: atomtypelenlist, nadjs0, nadjs1
   integer, dimension(:, :), intent(in) :: adjlists0, adjlists1
   integer, dimension(:), intent(in) :: atomequivlenlist0, atomequivlenlist1
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords0, coords1
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   integer, dimension(:), intent(inout) :: mapping
   integer, intent(in) :: nfrag0
   integer, dimension(:), intent(in) :: fragroot0

   logical, parameter :: printInfo = .false.

   integer blkidx(natom)
   integer h, i, offset, moldiff
   integer ntrack, track(natom)
   integer unmapping(natom)
   integer, dimension(natom) :: eqvidx0, eqvidx1
   logical tracked(natom)
   real(wp) moldist

   ! set atoms block indices

   offset = 0
   do h = 1, natomtype
      blkidx(offset+1:offset+atomtypelenlist(h)) = h
      offset = offset + atomtypelenlist(h)
   end do

   ! set atoms equivalence indices

   offset = 0
   do h = 1, natomequiv0
      eqvidx0(offset+1:offset+atomequivlenlist0(h)) = h
      offset = offset + atomequivlenlist0(h)
   end do

   offset = 0
   do h = 1, natomequiv1
      eqvidx1(offset+1:offset+atomequivlenlist1(h)) = h
      offset = offset + atomequivlenlist1(h)
   end do

   !  initialization

   ntrack = 0
   tracked(:) = .false.
   unmapping = inverse_mapping(mapping)
   moldiff = adjacencydiff (natom, adjmat0, adjmat1, mapping)
   moldist = squaredist (natom, weights, coords0, coords1, mapping)

   if ( printInfo ) then
      print '(a,1x,i0)', "moldiff:", moldiff
      print '(a,1x,f0.4)', "moldist:", moldist
   end if

   do i = 1, nfrag0
      call recursive_backtrack (fragroot0(i), mapping, unmapping, tracked, moldiff, moldist, ntrack, track)
!        print *, ntrack
   end do

   if ( printInfo ) then
      print '(a,1x,i0)', "countFrag:", nfrag0
      print '(a,1x,i0,1x,i0)', "moldiff:", adjacencydiff (natom, adjmat0, adjmat1, mapping), moldiff
      print '(a,1x,f0.4,1x,f0.4)', "moldist:", squaredist (natom, weights, coords0, coords1, mapping), moldist
   end if

!    if (adjacencydiff (natom, adjmat0, adjmat1, mapping) /= moldiff) then
!        print '(a,x,i0,x,i0)', "moldiff:", adjacencydiff (natom, adjmat0, adjmat1, mapping), moldiff
!    end if

   contains    

! classify the atoms connected to node as matches or unmatched
   subroutine nodematch (node, mapping, tracked, nmatch, matches, &
                   nmismatch0, mismatches0, nmismatch1, mismatches1)
      integer, intent(in) :: node, mapping(:)
      logical, intent(in) :: tracked(:)
      integer, intent(out) :: nmatch, matches(natom)
      integer, intent(out) :: nmismatch0, mismatches0(natom)
      integer, intent(out) :: nmismatch1, mismatches1(natom)
      integer :: i, adj0(natom), adj1(natom)

      adj0(:nadjs0(node)) = adjlists0(:nadjs0(node), node)
      adj1(:nadjs1(mapping(node))) = adjlists1(:nadjs1(mapping(node)), mapping(node))
      nmatch = 0
      nmismatch0 = 0
      nmismatch1 = 0

      do i = 1, nadjs0(node)
         if (findInd(mapping(adj0(i)), nadjs1(mapping(node)), adj1) /= 0) then
            call addInd(adj0(i), nmatch, matches)
         else
            call addInd(adj0(i), nmismatch0, mismatches0)
         end if
      end do

      do i = 1, nadjs1(mapping(node))
         if (findInd(adj1(i), nmatch, mapping(matches(:nmatch))) == 0) then
            call addInd(adj1(i), nmismatch1, mismatches1)
         end if
      end do
   end subroutine nodematch

! backtracks structure to find assignments that minimize moldiff
   recursive subroutine recursive_backtrack (node, mapping, unmapping, tracked, moldiff, moldist, &
                                  ntrack, track)
      integer, intent(in) :: node
      integer, intent(inout) :: mapping(natom), unmapping(natom)
      logical, intent(inout) :: tracked(natom)
      integer, intent(inout) :: moldiff, ntrack, track(natom)
      real(wp), intent(inout) :: moldist

      logical, parameter :: printInfo = .false.

      integer nmatch, matches(natom)
      integer nmismatch0, mismatches0(natom)
      integer nmismatch1, mismatches1(natom)
      integer mapping_branch(natom), moldiff_branch, unmapping_branch(natom)
      integer ntrack_branch, track_branch(natom)
      logical tracked_branch(natom)
      integer :: i, j
      logical :: matched0(natom), matched1(natom)
      real(wp) :: moldist_branch, mol0local(3,natom), mol1local(3,natom)

      ! reserve node node as tracked
      ntrack = ntrack + 1
      track(ntrack) = node
      tracked(node) = .true.

      ! classify neighbor atoms as matches or mismatched for coords0/coords1
      call nodematch (node, mapping, tracked, nmatch, matches, &
                  nmismatch0, mismatches0, nmismatch1, mismatches1)

      if (printInfo) then
         print *, "node:", node
         print *, "matches:", matches(:nmatch)
         print *, "mismatches0:", mismatches0(:nmismatch0)
         print *, "mismatches1:", mismatches1(:nmismatch1)
      end if

      ! shuffle indices
      call shuffle (matches, nmatch)
      call shuffle (mismatches0, nmismatch0)
      call shuffle (mismatches1, nmismatch1)

      ! run over matches neighbors
      do i = 1, nmatch
         if (.not. tracked(matches(i))) then
            call recursive_backtrack(matches(i), mapping, unmapping, tracked, moldiff, &
                            moldist, ntrack, track)
         end if
      end do

      matched0(:nmismatch0) = .false.
      matched1(:nmismatch1) = .false.

      ! run over mismatched neighbors
      do i = 1, nmismatch0
         if (.not. tracked(mismatches0(i))) then
            do j = 1, nmismatch1
               if (.not. matched1(j)) then
                  if (blkidx(mismatches0(i)) == blkidx(mismatches1(j))) then

                     ntrack_branch = ntrack
                     track_branch(:) = track(:)
                     tracked_branch(:) = tracked(:)
                     mapping_branch(:) = mapping(:)
                     unmapping_branch(:) = unmapping(:)

                     ! Apply swap to mapping branch 
                     mapping_branch(mismatches0(i)) = mismatches1(j)
                     mapping_branch(unmapping(mismatches1(j))) = mapping(mismatches0(i))

                     ! Apply swap to unmapping branch 
                     unmapping_branch(mismatches1(j)) = mismatches0(i)
                     unmapping_branch(mapping(mismatches0(i))) = unmapping(mismatches1(j))

                     ! Update ssd with swap
                     moldist_branch = moldist + weights(mismatches0(i))*( &
                        - sum((coords1(:, mapping(mismatches0(i))) - coords0(:, mismatches0(i)))**2) &
                        - sum((coords1(:, mismatches1(j)) - coords0(:, unmapping(mismatches1(j))))**2) &
                        + sum((coords1(:, mismatches1(j)) - coords0(:, mismatches0(i)))**2) &
                        + sum((coords1(:, mapping(mismatches0(i))) - coords0(:, unmapping(mismatches1(j))))**2))

                     ! Update adjd with swap
                     moldiff_branch = moldiff + adjacencydelta(nadjs0, adjlists0, adjmat1, &
                                   mapping, mismatches0(i), unmapping(mismatches1(j)))

                     ! backtrack swapped index
                     call recursive_backtrack(mismatches0(i), mapping_branch, unmapping_branch, &
                        tracked_branch, moldiff_branch, moldist_branch, ntrack_branch, track_branch)

                     if ( &
                        moldiff_branch < moldiff &
                        .and. ( &
                           eqvidx0(mismatches0(i)) == eqvidx0(unmapping(mismatches1(j))) &
                           .and. eqvidx1(mapping(mismatches0(i))) == eqvidx1(mismatches1(j)) &
                        ) &
                     ) then
                        ntrack = ntrack_branch
                        track(:) = track_branch(:)
                        tracked(:) = tracked_branch(:)
                        mapping(:) = mapping_branch(:)
                        unmapping(:) = unmapping_branch(:)
                        moldiff = moldiff_branch
                        moldist = moldist_branch
                        matched0(i) = .true.
                        matched1(j) = .true.
                        exit   ! exits inner do loop
                     end if
                  end if
               end if
            end do
         end if
      end do

      ! run over non matches neighbors
      do i = 1, nmismatch0
         if (.not. matched0(i)) then
            if (.not. tracked(mismatches0(i))) then
               call recursive_backtrack(mismatches0(i), mapping, unmapping, &
                 tracked, moldiff, moldist, ntrack, track)
            end if
         end if
      end do

   end subroutine

end subroutine

subroutine eqvatomperm (natom, weights, coords, adjmat, adjlists, natomequiv, atomequivlenlist, &
   nadjequivs, adjequivlenlists, refcoords, refadjmat, mapping, nfrag, fragroots)
! Purpose: Find best correspondence between points of graphs

   integer, intent(in) :: natom, natomequiv
   integer, dimension(:), intent(in) :: atomequivlenlist, nadjequivs
   integer, dimension(:, :), intent(in) :: adjlists, adjequivlenlists
   integer, dimension(:), intent(inout) :: mapping
   real(wp), dimension(:), intent(in) :: weights
   real(wp), dimension(:, :), intent(in) :: coords, refcoords
   logical, dimension(:, :), intent(in) :: adjmat, refadjmat
   integer, intent(in) :: nfrag
   integer, dimension(:), intent(in) :: fragroots

   integer :: unmap(natom), track(natom)
   integer, dimension(natom) :: eqvos, eqvidx
   integer :: permcount, h, i, n, diff, fragcount, ntrack
   logical :: tracked(natom), held(natom)
   real(wp) :: dist

   integer offset

   ! set equivalence group offsets

   eqvos(1) = 0
   do h = 1, natomequiv - 1
      eqvos(h+1) = eqvos(h) + atomequivlenlist(h)
   end do

   ! set atoms equivalence indices

   do h = 1, natomequiv
      eqvidx(eqvos(h)+1:eqvos(h)+atomequivlenlist(h)) = h
   end do

   ! print equiv types

!   do i = 1, natom
!      print *, i, eqvidx(i)
!   end do

   ! print equiv adjlists

!   do i = 1, natom
!      offset = 0
!      do h = 1, nadjequivs(i)
!         print *, i, ':', adjlists(offset+1:offset+adjequivlenlists(h, i), i)
!         offset = offset + adjequivlenlists(h, i)
!      end do
!   end do

   ! initialization

   permcount = 1
   fragcount = 0
   tracked(:) = .false.
   held(:) = .false.
   ntrack = 0

!    i = 1
!    do while ( i <= natom )
!        if ( tracked(i) ) then
!            i = i + 1
!        else
!            call recursive_permut (i, mapping, tracked, held, 0, 0)
!            i = 1   ! restart loop to find unconnected, untracked atoms
!            fragcount = fragcount + 1
!        end if
!    end do

   do i = 1, nfrag
      call recursive_permut (fragroots(i), mapping, tracked, held, ntrack, track)
!        print *, ntrack
   end do

!   print '(a,i0)', "Natoms: ", natom
!   print '(a,f8.4)', "dist: ", sqrt(leastsquaredist (natom, weights, coords, refcoords, mapping)/sum(weights))
!   print '(a,i0)', "permcount: ", permcount
!   print '(a,i0)', "fragcount: ", fragcount

   contains    

   recursive subroutine recursive_remap (nodea, nodeb, mapping_ref, mapping, held)
      integer, intent(in) :: nodea, nodeb, mapping_ref(natom)
      integer, intent(inout) :: mapping(natom)
      logical, dimension(natom), intent(inout) :: held
      logical, dimension(natom) :: locked_c
      integer meqvnei, equiva(natom), equivb(natom)
      integer h, i, offset, first, last

      first = eqvos(eqvidx(nodea)) + 1
      last = eqvos(eqvidx(nodea)) + atomequivlenlist(eqvidx(nodea))
      held(first:last) = .true.
      offset = 0
      do h = 1, nadjequivs(nodea)
         ! find not tracked adjlists in group
         meqvnei = 0
         do i = 1, adjequivlenlists(h, nodea)
            if (.not. held(adjlists(offset+i, nodea))) then 
               meqvnei = meqvnei + 1
               equiva(meqvnei) = adjlists(offset+meqvnei, nodea)
               equivb(meqvnei) = adjlists(offset+meqvnei, nodeb)
            end if
         end do
         locked_c(:) = held(:)
         do i = 1, meqvnei
            if ( equiva(i) /= equivb(i) ) then
               mapping(equiva(i)) = mapping_ref(equivb(i))
               call recursive_remap (equiva(i), equivb(i), mapping_ref, &
                  mapping, locked_c)
            end if
         end do
         offset = offset + adjequivlenlists(h, nodea)
      end do
   end subroutine recursive_remap

   recursive subroutine recursive_permut (node, mapping, tracked, held, ntrack, &
                                track)
      integer, intent(in) :: node
      integer, intent(inout) :: ntrack
      logical, dimension(natom), intent(inout) :: tracked, held
      integer, dimension(natom), intent(inout) :: mapping, track

      logical :: locked_c(natom)
      integer meqvnei, moldiff, track4ind(4), track4ind_c(4)
      integer, dimension(natom) :: mapping_p, mapping_min, equiv, perm, perm_min
      real(wp) :: moldist, moldist_p, moldist_min, dihed0(maxcoord), dihed1(maxcoord)
      logical :: nextperm, even, calcd, printInfo = .false.
      integer :: h, i, j, offset, first, last
      character(len=80) :: strfmt

      if ( printInfo ) then   ! print debugging info
         moldist = sqrt(leastsquaredist(natom, weights, coords, refcoords, mapping)/sum(weights))
         moldiff = adjacencydiff(natom, adjmat, refadjmat, mapping)
         write (strfmt, '(a,i0,a)') '(',1,'(2x),i0,a,i0,f8.4)'
         print strfmt, node,": ",moldiff,moldist
      end if

      ! reserves node as tracked
      tracked(node) = .true.
      ntrack = ntrack + 1
      track(ntrack) = node

      ! reserves atoms with the same type as n
      first = eqvos(eqvidx(node)) + 1
      last = eqvos(eqvidx(node)) + atomequivlenlist(eqvidx(node))
      held(first:last) = .true.

      offset = 0
      ! run over groups of atoms with equivalent type
      do h = 1, nadjequivs(node)
         ! find not tracked adjlists in group
         meqvnei = 0
         do i = 1, adjequivlenlists(h, node)
            if (.not. tracked(adjlists(offset+i, node)) &
               .and. .not. held(adjlists(offset+i, node))) then
               meqvnei = meqvnei + 1
               equiv(meqvnei) = adjlists(offset+meqvnei, node)
            end if
         end do
         ! check permutations
         if ( meqvnei > 0 ) then
            ! save initial state before permutations
            mapping_min(:) = mapping(:)
            ! run over all permutarions for the meqvnei equivalent atoms
            nextperm = .false.
            call perm1_next (meqvnei, perm, nextperm, even)
            perm_min(:) = perm(:)
            ! initialize state to test new permutation
            mapping_min(:) = mapping(:)
            moldist_min = 0.0
            ! apply permutation
            do i = 1, meqvnei
               mapping_min(equiv(i)) = mapping(equiv(perm_min(i)))
               moldist_min = moldist_min &
                        + sum((refcoords(:, mapping_min(equiv(i))) - coords(:, equiv(i)))**2)
            end do
            do while ( nextperm )
               call perm1_next (meqvnei, perm, nextperm, even)
               permcount = permcount + 1
               ! initialize state to test new permutation
               mapping_p(:) = mapping(:)
               moldist_p = 0.0
               ! apply permutation
               do i = 1, meqvnei
                  mapping_p(equiv(i)) = mapping(equiv(perm(i)))
                  moldist_p = moldist_p &
                         + sum((refcoords(:, mapping_p(equiv(i))) - coords(:, equiv(i)))**2) 
               end do
               ! save min dist permut
               if ( moldist_p < moldist_min ) then
                  mapping_min(:) = mapping_p(:)
                  moldist_min = moldist_p
                  perm_min(:) = perm(:)
               end if
            end do
            ! check chirality
!                if ( printInfo ) then 
!
!                    if ( meqvnei >= 3 ) then
!                        track4ind(1) = equiv(1)
!                        track4ind(2) = node
!                        track4ind(3) = track(ntrack)
!
!                        do j = 1, meqvnei
!                            track4ind(4) = equiv(j)
!                            call calc_dihedral (0, track4ind, coords, calcd, dihed0(j))
!                            track4ind_c(:) = mapping_min(track4ind(:))
!                            call calc_dihedral (0, track4ind_c, refcoords, calcd, dihed1(j))
!                        end do
!                        print '(a,3f10.3)', "coords: ", dihed0(:meqvnei)
!                        print '(a,3f10.3)', "refcoords: ", dihed1(:meqvnei)
!                    end if
!                end if
            ! remap connectivity of permuted branch path
            locked_c(:) = held(:)
            do i = 1, meqvnei
               if ( equiv(i) /= equiv(perm_min(i)) ) then
                  call recursive_remap (equiv(i), equiv(perm_min(i)), &
                               mapping, mapping_min, locked_c)
               end if
            end do
            ! update mapping
            mapping(:) = mapping_min(:)
            do i = 1, meqvnei
               call recursive_permut (equiv(i), mapping, tracked, &
                                locked_c, ntrack, track)
            end do
         end if
         offset = offset + adjequivlenlists(h, node)
      end do
   end subroutine

end subroutine

! Adds a new element to a list of indices in growing order, no repeats
subroutine addInd(ind, nInd, vecInd)
   integer, intent(in) :: ind
   integer, intent(inout) :: nInd, vecInd(:)
   integer :: p, e
   p = 1
   do while (p <= nInd)   ! find position for growing order
      if (ind == vecInd(p)) return
      if (ind < vecInd(p)) then
         exit
      else
         p = p + 1
      end if
   end do
   nInd = nInd + 1
   e = nInd
   do while (e > p)
      vecInd(e) = vecInd(e - 1)   ! shifts one position from the end
      e = e - 1
   end do
   vecInd(p) = ind
end subroutine addInd

! Returns the position in vecInd where the element ind is found or 0 otherwise
function findInd(ind, nInd, vecInd) result(pos)
   integer, intent(in) :: ind, nInd, vecInd(nInd)
   integer :: pos
   integer :: i
   pos = 0
   do i = 1, nInd
      if (ind == vecInd(i)) then
         pos = i
         return
      end if
   end do
end function findInd

end module
