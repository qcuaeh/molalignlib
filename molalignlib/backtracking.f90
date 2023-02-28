! template for module backtracking in molalign program
module backtracking

use iso_fortran_env, only: output_unit, error_unit
use sorting
use discrete
use adjacency
use alignment
use random
implicit none
private
public backtrack_bonds

contains

! subroutine matchbonds
! Purpose: Find best correspondence between points of graphs

! natom:  Number of atoms
! weights: atom masses
! coords0, coords1: matix of coordinates?
! nblock: number of blocks
! blksz: Number of atoms in each block
! blkid: atom mapping by blkid (atom type)
! nadj0, nadj1: coordination numbers for coords0/coords1
! adjlist0, adjlist1: list of coordinated neighbors for coords0/coords1
! adjmat0, adjmat1: Adjacency matrix for coords0/coords1
! atomperm: atom index mapping
subroutine backtrack_bonds (natom, weights, blkid, coords0, nadj0, adjlist0, adjmat0, &
               coords1, nadj1, adjlist1, adjmat1, atomperm, nfrag0, fragrt0)
! Purpose: Find best correspondence between points of graphs

    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: blkid, nadj0, nadj1
    integer, dimension(:, :), intent(in) :: adjlist0, adjlist1
    integer, dimension(:), intent(inout) :: atomperm
    real(wp), dimension(:), intent(in) :: weights
    real(wp), dimension(:, :), intent(in) :: coords0, coords1
    logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
    integer, intent(in) :: nfrag0
    integer, dimension(:), intent(in) :: fragrt0

    logical, parameter :: printInfo = .false.

    integer unmapping(natom), ntrack, track(natom)
    logical tracked(natom)
    integer moldiff
    real(wp) moldist
    integer i

    ntrack = 0
    tracked(:) = .false.
    unmapping = inverseperm(atomperm)
    moldiff = 0

!    moldist = squaredist (natom, weights, coords0, coords1, atomperm)
!    moldiff = adjacencydiff (natom, adjmat0, adjmat1, atomperm)

    if ( printInfo ) then
        print '(a,1x,i0)', "moldiff:", moldiff
        print '(a,1x,f0.4)', "moldist:", moldist
    end if

    do i = 1, nfrag0
        call recursive_backtrack (fragrt0(i), atomperm, unmapping, tracked, moldiff, moldist, ntrack, track)
!        print *, ntrack
    end do

    if ( printInfo ) then
        print '(a,1x,i0)', "countFrag:", nfrag0
        print '(a,1x,i0,1x,i0)', "moldiff:", adjacencydiff (natom, adjmat0, adjmat1, atomperm), moldiff
        print '(a,1x,f0.4,1x,f0.4)', "moldist:", squaredist (natom, weights, coords0, coords1, atomperm), moldist
    end if

!    if (adjacencydiff (natom, adjmat0, adjmat1, atomperm) /= moldiff) then
!        print '(a,x,i0,x,i0)', "moldiff:", adjacencydiff (natom, adjmat0, adjmat1, atomperm), moldiff
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
        integer :: i, nnadj0, nnadj1, adj0(natom), adj1(natom)

        nnadj0 = nadj0(node)
        adj0(:nnadj0) = adjlist0(:nnadj0, node)
        nnadj1 = nadj1(mapping(node))
        adj1(:nnadj1) = adjlist1(:nnadj1, mapping(node))
        nmatch = 0
        nmismatch0 = 0
        nmismatch1 = 0

        do i = 1, nnadj0
            if (findInd(mapping(adj0(i)), nnadj1, adj1) /= 0) then
                call addInd(adj0(i), nmatch, matches)
            else
                call addInd(adj0(i), nmismatch0, mismatches0)
            end if
        end do

        do i = 1, nnadj1
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
                        if (blkid(mismatches0(i)) == blkid(mismatches1(j))) then

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
!                            moldist_branch = moldist + weights(mismatches0(i))*( &
!                                - sum((coords1(:, mapping(mismatches0(i))) - coords0(:, mismatches0(i)))**2) &
!                                - sum((coords1(:, mismatches1(j)) - coords0(:, unmapping(mismatches1(j))))**2) &
!                                + sum((coords1(:, mismatches1(j)) - coords0(:, mismatches0(i)))**2) &
!                                + sum((coords1(:, mapping(mismatches0(i))) - coords0(:, unmapping(mismatches1(j))))**2))

                            ! Update adjd with swap
                            moldiff_branch = moldiff + adjacencydelta(nadj0, adjlist0, adjmat1, &
                                               mapping, mismatches0(i), unmapping(mismatches1(j)))

                            ! backtrack swapped index
                            call recursive_backtrack(mismatches0(i), mapping_branch, unmapping_branch, &
                                tracked_branch, moldiff_branch, moldist_branch, ntrack_branch, track_branch)

                            if (moldiff_branch < moldiff) then
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

