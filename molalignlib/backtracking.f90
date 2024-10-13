! template for module backtracking in molalign program
module backtracking
use stdio
use sorting
use molecule
use partition
use permutation
use adjacency
use alignment
use random

implicit none
private
public minadjdiff
public eqvatomperm

logical, parameter :: printInfo = .false.

contains

! Find best correspondence between points of graphs
subroutine minadjdiff (mol0, mol1, mapping)
   type(molecule_type), intent(in) :: mol0, mol1
   integer, dimension(:), intent(inout) :: mapping

   ! Local variables

   integer :: eltypemap0(mol0%natom)
   integer :: h, i, moldiff
   integer :: ntrack, track(mol0%natom)
   integer :: unmapping(mol0%natom)
   integer, dimension(mol0%natom) :: mnatypemap0, mnatypemap1
   logical :: tracked(mol0%natom)
   real(rk) moldist

   integer :: natom
   integer :: neltype0
   integer :: nmnatype0, nmnatype1

   type(partition_type) :: eltypes0
   type(partition_type) :: mnatypes0, mnatypes1

   integer, allocatable, dimension(:) :: eltypepartsizes0
   integer, allocatable, dimension(:) :: mnatypepartsizes0, mnatypepartsizes1
   integer, allocatable, dimension(:) :: fragroots0, fragroots1

   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   real(rk), dimension(:, :), allocatable :: coords0, coords1
   real(rk), dimension(:), allocatable :: weights

   integer :: nadjs0(mol0%natom)
   integer :: nadjs1(mol1%natom)
   integer :: adjlists0(maxcoord, mol0%natom)
   integer :: adjlists1(maxcoord, mol1%natom)
   integer :: nadjmnatypes0(mol0%natom)
   integer :: nadjmnatypes1(mol1%natom)
   integer :: adjmnatypepartlens0(maxcoord, mol0%natom)
   integer :: adjmnatypepartlens1(maxcoord, mol1%natom)

   natom = mol0%get_natom()

   eltypes0 = mol0%gather_eltypes()
   eltypemap0 = eltypes0%idxmap

   neltype0 = eltypes0%size
   allocate (eltypepartsizes0(neltype0))
   do h = 1, neltype0
      eltypepartsizes0(h) = eltypes0%parts(h)%size
   end do

   mnatypes0 = mol0%gather_mnatypes()
   mnatypes1 = mol1%gather_mnatypes()
   mnatypemap0 = mnatypes0%idxmap
   mnatypemap1 = mnatypes1%idxmap

   nmnatype0 = size(mnatypes0%parts)
   nmnatype1 = size(mnatypes1%parts)
   allocate (mnatypepartsizes0(nmnatype0))
   allocate (mnatypepartsizes1(nmnatype1))
   do h = 1, nmnatype0
      mnatypepartsizes0(h) = mnatypes0%parts(h)%size
   end do
   do h = 1, nmnatype1
      mnatypepartsizes1(h) = mnatypes1%parts(h)%size
   end do

!   adjlists0 = mol0%gather_adjlists()
!   adjlists1 = mol1%gather_adjlists()
   nadjs0 = mol0%gather_nadjs()
   nadjs1 = mol1%gather_nadjs()
   adjlists0 = mol0%gather_adjlists()
   adjlists1 = mol1%gather_adjlists()

!   adjpartitions0 = mol0%gather_adjpartitions()
!   adjpartitions1 = mol1%gather_adjpartitions()
   nadjmnatypes0 = mol0%gather_nadjmnatypes()
   nadjmnatypes1 = mol1%gather_nadjmnatypes()
   adjmnatypepartlens0 = mol0%gather_adjmnatypepartlens()
   adjmnatypepartlens1 = mol1%gather_adjmnatypepartlens()

   coords0 = mol0%gather_coords()
   coords1 = mol1%gather_coords()
   weights = mol0%gather_weights()
   adjmat0 = mol0%gather_adjmatrix()
   adjmat1 = mol1%gather_adjmatrix()
   fragroots0 = mol0%gather_molfragroots()
   fragroots1 = mol1%gather_molfragroots()

   !  initialization

   ntrack = 0
   tracked(:) = .false.
   unmapping = inverse_permutation(mapping)
   moldiff = adjacencydiff (natom, adjmat0, adjmat1, mapping)
   moldist = squaredist (natom, weights, coords0, coords1, mapping)

   if ( printInfo ) then
      print '(a,1x,i0)', "moldiff:", moldiff
      print '(a,1x,f0.4)', "moldist:", moldist
   end if

   do i = 1, size(fragroots0)
      call recursive_backtrack (fragroots0(i), mapping, unmapping, tracked, moldiff, moldist, ntrack, track)
!        print *, ntrack
   end do

   if ( printInfo ) then
      print '(a,1x,i0)', "countFrag:", size(fragroots0)
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
      real(rk), intent(inout) :: moldist

      logical, parameter :: printInfo = .false.

      integer :: nmatch, matches(natom)
      integer :: nmismatch0, mismatches0(natom)
      integer :: nmismatch1, mismatches1(natom)
      integer :: mapping_branch(natom), moldiff_branch, unmapping_branch(natom)
      integer :: ntrack_branch, track_branch(natom)
      logical :: tracked_branch(natom)
      integer :: i, j
      logical :: matched0(natom), matched1(natom)
      real(rk) :: moldist_branch, mol0local(3,natom), mol1local(3,natom)

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
      call shuffle (matches(:nmatch))
      call shuffle (mismatches0(:nmismatch0))
      call shuffle (mismatches1(:nmismatch1))

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
                  if (eltypemap0(mismatches0(i)) == eltypemap0(mismatches1(j))) then

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
                           mnatypemap0(mismatches0(i)) == mnatypemap0(unmapping(mismatches1(j))) &
                           .and. mnatypemap1(mapping(mismatches0(i))) == mnatypemap1(mismatches1(j)) &
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

! Find best correspondence between points of graphs
subroutine eqvatomperm (mol0, mol1, workcoords1, mapping)
   type(molecule_type), intent(in) :: mol0, mol1
   real(rk), intent(in) :: workcoords1(:, :)
   integer, intent(inout) :: mapping(:)

   ! Local variables

   integer :: unmap(mol0%natom), track(mol0%natom)
   integer, dimension(mol0%natom) :: eqvos, eqvidx
   integer :: permcount, h, i, n, diff, fragcount, ntrack
   logical :: tracked(mol0%natom), held(mol0%natom)
   real(rk) :: dist

   integer :: natom
   integer :: neltype0
   integer :: nmnatype0, nmnatype1

   type(partition_type) :: eltypes0
   type(partition_type) :: mnatypes0, mnatypes1

   integer, allocatable, dimension(:) :: eltypepartsizes0
   integer, allocatable, dimension(:) :: mnatypepartsizes0, mnatypepartsizes1
   integer, allocatable, dimension(:) :: fragroots0, fragroots1

   logical, dimension(:, :), allocatable :: adjmat0, adjmat1
   real(rk), dimension(:, :), allocatable :: coords0, coords1
   real(rk), dimension(:), allocatable :: weights

   integer :: nadjs0(mol0%natom)
   integer :: nadjs1(mol1%natom)
   integer :: adjlists0(maxcoord, mol0%natom)
   integer :: adjlists1(maxcoord, mol1%natom)
   integer :: nadjmnatypes0(mol0%natom)
   integer :: nadjmnatypes1(mol1%natom)
   integer :: adjmnatypepartlens0(maxcoord, mol0%natom)
   integer :: adjmnatypepartlens1(maxcoord, mol1%natom)

   natom = mol0%get_natom()

   eltypes0 = mol0%gather_eltypes()

   neltype0 = eltypes0%size
   allocate (eltypepartsizes0(neltype0))
   do h = 1, neltype0
      eltypepartsizes0(h) = eltypes0%parts(h)%size
   end do

   mnatypes0 = mol0%gather_mnatypes()
   mnatypes1 = mol1%gather_mnatypes()

   nmnatype0 = size(mnatypes0%parts)
   nmnatype1 = size(mnatypes1%parts)
   allocate (mnatypepartsizes0(nmnatype0))
   allocate (mnatypepartsizes1(nmnatype1))
   do h = 1, nmnatype0
      mnatypepartsizes0(h) = mnatypes0%parts(h)%size
   end do
   do h = 1, nmnatype1
      mnatypepartsizes1(h) = mnatypes1%parts(h)%size
   end do

!   adjlists0 = mol0%gather_adjlists()
!   adjlists1 = mol1%gather_adjlists()
   nadjs0 = mol0%gather_nadjs()
   nadjs1 = mol1%gather_nadjs()
   adjlists0 = mol0%gather_adjlists()
   adjlists1 = mol1%gather_adjlists()

!   adjpartitions0 = mol0%gather_adjpartitions()
!   adjpartitions1 = mol1%gather_adjpartitions()
   nadjmnatypes0 = mol0%gather_nadjmnatypes()
   nadjmnatypes1 = mol1%gather_nadjmnatypes()
   adjmnatypepartlens0 = mol0%gather_adjmnatypepartlens()
   adjmnatypepartlens1 = mol1%gather_adjmnatypepartlens()

   coords0 = mol0%gather_coords()
   coords1 = mol1%gather_coords()
   weights = mol0%gather_weights()
   adjmat0 = mol0%gather_adjmatrix()
   adjmat1 = mol1%gather_adjmatrix()
   fragroots0 = mol0%gather_molfragroots()
   fragroots1 = mol1%gather_molfragroots()

   ! set equivalence group offsets

   eqvos(1) = 0
   do h = 1, nmnatype0 - 1
      eqvos(h+1) = eqvos(h) + mnatypepartsizes0(h)
   end do

   ! set atoms equivalence indices

   do h = 1, nmnatype0
      eqvidx(eqvos(h)+1:eqvos(h)+mnatypepartsizes0(h)) = h
   end do

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

   do i = 1, size(fragroots0)
      call recursive_permut (fragroots0(i), mapping, tracked, held, ntrack, track)
!        print *, ntrack
   end do

!   print '(a,i0)', "Natoms: ", natom
!   print '(a,f8.4)', "dist: ", sqrt(leastsquaredist (natom, weights, coords0, workcoords1, mapping)/sum(weights))
!   print '(a,i0)', "permcount: ", permcount
!   print '(a,i0)', "fragcount: ", fragcount

   contains    

   recursive subroutine recursive_remap (nodea, nodeb, mapping_ref, mapping, held)
      integer, intent(in) :: nodea, nodeb, mapping_ref(natom)
      integer, intent(inout) :: mapping(natom)
      logical, dimension(natom), intent(inout) :: held
      logical, dimension(natom) :: locked_c
      integer :: meqvnei, equiva(natom), equivb(natom)
      integer :: h, i, offset, first, last

      first = eqvos(eqvidx(nodea)) + 1
      last = eqvos(eqvidx(nodea)) + mnatypepartsizes0(eqvidx(nodea))
      held(first:last) = .true.
      offset = 0
      do h = 1, nadjmnatypes0(nodea)
         ! find not tracked adjlists0 in group
         meqvnei = 0
         do i = 1, adjmnatypepartlens0(h, nodea)
            if (.not. held(adjlists0(offset+i, nodea))) then 
               meqvnei = meqvnei + 1
               equiva(meqvnei) = adjlists0(offset+meqvnei, nodea)
               equivb(meqvnei) = adjlists0(offset+meqvnei, nodeb)
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
         offset = offset + adjmnatypepartlens0(h, nodea)
      end do
   end subroutine recursive_remap

   recursive subroutine recursive_permut (node, mapping, tracked, held, ntrack, track)
      integer, intent(in) :: node
      integer, intent(inout) :: ntrack
      logical, dimension(natom), intent(inout) :: tracked, held
      integer, dimension(natom), intent(inout) :: mapping, track

      logical :: locked_c(natom)
      integer :: meqvnei, moldiff, track4ind(4), track4ind_c(4)
      integer, dimension(natom) :: mapping_p, mapping_min, equiv, perm, perm_min
      real(rk) :: moldist, moldist_p, moldist_min, dihed0(maxcoord), dihed1(maxcoord)
      logical :: more, calcd, printInfo = .false.
      integer :: h, i, j, offset, first, last, rank
      character(len=80) :: strfmt

      if ( printInfo ) then   ! print debugging info
         moldist = sqrt(leastsquaredist(natom, weights, coords0, workcoords1, mapping)/sum(weights))
         moldiff = adjacencydiff(natom, adjmat0, adjmat1, mapping)
         write (strfmt, '(a,i0,a)') '(',1,'(2x),i0,a,i0,f8.4)'
         print strfmt, node,": ",moldiff,moldist
      end if

      ! reserves node as tracked
      tracked(node) = .true.
      ntrack = ntrack + 1
      track(ntrack) = node

      ! reserves atoms with the same type as n
      first = eqvos(eqvidx(node)) + 1
      last = eqvos(eqvidx(node)) + mnatypepartsizes0(eqvidx(node))
      held(first:last) = .true.

      offset = 0
      ! run over groups of atoms with equivalent type
      do h = 1, nadjmnatypes0(node)
         ! find not tracked adjlists0 in group
         meqvnei = 0
         do i = 1, adjmnatypepartlens0(h, node)
            if (.not. tracked(adjlists0(offset+i, node)) &
               .and. .not. held(adjlists0(offset+i, node))) then
               meqvnei = meqvnei + 1
               equiv(meqvnei) = adjlists0(offset+meqvnei, node)
            end if
         end do
         ! check permutations
         if ( meqvnei > 0 ) then
            ! save initial state before permutations
            mapping_min(:) = mapping(:)
            ! run over all permutarions for the meqvnei equivalent atoms
            more = .false.
            call perm1_next3 (meqvnei, perm, more, rank)
            perm_min(:) = perm(:)
            ! initialize state to test new permutation
            mapping_min(:) = mapping(:)
            moldist_min = 0.0
            ! apply permutation
            do i = 1, meqvnei
               mapping_min(equiv(i)) = mapping(equiv(perm_min(i)))
               moldist_min = moldist_min &
                        + sum((workcoords1(:, mapping_min(equiv(i))) - coords0(:, equiv(i)))**2)
            end do
            do while ( more )
               call perm1_next3 (meqvnei, perm, more, rank)
               permcount = permcount + 1
               ! initialize state to test new permutation
               mapping_p(:) = mapping(:)
               moldist_p = 0.0
               ! apply permutation
               do i = 1, meqvnei
                  mapping_p(equiv(i)) = mapping(equiv(perm(i)))
                  moldist_p = moldist_p &
                         + sum((workcoords1(:, mapping_p(equiv(i))) - coords0(:, equiv(i)))**2) 
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
!                            call calc_dihedral (0, track4ind, coords0, calcd, dihed0(j))
!                            track4ind_c(:) = mapping_min(track4ind(:))
!                            call calc_dihedral (0, track4ind_c, workcoords1, calcd, dihed1(j))
!                        end do
!                        print '(a,3f10.3)', "coords0: ", dihed0(:meqvnei)
!                        print '(a,3f10.3)', "workcoords1: ", dihed1(:meqvnei)
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
         offset = offset + adjmnatypepartlens0(h, node)
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
