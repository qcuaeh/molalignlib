module remapping

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use utilities
use options
use random
use assignment
use translation
use alignment
use messages
use rotation
use biasing
use printing

implicit none

contains

logical function trial_stop_test(trials, matches) result(stop_test)
    integer, intent(in) :: trials, matches
    stop_test = trials < maxcount
end function

logical function match_stop_test(trials, matches) result(stop_test)
    integer, intent(in) :: trials, matches
    stop_test = matches < maxcount
end function

subroutine remapatoms(natom, nblock, blocksize, atoms0, atoms1, weights, mapcount, atomaplist, bias)

! natom: number of points to align
! atoms0: coordinates of first molecule
! atoms1: coordinates of second molecule
! nblock: number of block atoms
! blocksize: size of atom block (number of atoms in the block)
! newdist: total square distance between molecules
! mindist: minimum distance
! atomap: mapping between the corresponding atoms in both molecules
! atomaplist: mapping list of minimum distances or differences
! matches: number of equivalent orientations
! maxcount: maximum number of searching cycles
! ntrial: random rotation trial counter

    integer, intent(out) :: mapcount
    integer, intent(in) :: natom, nblock
    integer, dimension(:), intent(in) :: blocksize
    integer, dimension(:, :), intent(out) :: atomaplist
    real, dimension(:, :), intent(inout) :: atoms1
    real, dimension(:, :), intent(in) :: atoms0
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: bias

    logical map_found, overflow
    integer i, imap, jmap, ntrial, iteration, newdiff, matchcount
    integer, dimension(natom) :: atomap, auxmap, equivset0, equivset1
    integer, dimension(maxrecord) :: earliest, matches
    real u, olddist, newdist, meanrot
    real, dimension(4) :: rotquat, prodquat
    real, dimension(maxrecord) :: mindist, avgiter, avgmeanrot, avgtotrot
    integer :: nfrag0, nfrag1
    integer, dimension(natom) :: fragroot0, fragroot1, fragid0, fragid1
    real, dimension(3, natom) :: origmol1

! Print header and initial stats

    if (live) then
        write (output_unit, '(a)', advance='no') achar(27)//'[1H'//achar(27)//'[J'
        call print_header()
    end if


! Save original coordinates

    origmol1 = atoms1

! Initialize loop variables

    ntrial = 0
    mapcount = 0
    matchcount = 0
    overflow = .false.

! Loop for map searching

    do while (stop_test(ntrial, matchcount))

        ntrial = ntrial + 1

! Apply random rotation

        call rotate(natom, randrotquat(randvec3()), atoms1)

! Find the best mappping to minimize euclidean distance

        call assignatoms(natom, weights, atoms0, atoms1, nblock, blocksize, bias, atomap)
        rotquat = leastrotquat(natom, weights, atoms0, atoms1, atomap)
        prodquat = rotquat
        meanrot = rotangle(rotquat)
        call rotate(natom, rotquat, atoms1)
        iteration = 1

        do while (iterative)
            olddist = squaredist(natom, weights, atoms0, atoms1, atomap) &
                    + biasingdist(natom, weights, bias, atomap)
            call assignatoms(natom, weights, atoms0, atoms1, nblock, blocksize, bias, auxmap)
            if (all(auxmap == atomap)) exit
            newdist = squaredist(natom, weights, atoms0, atoms1, auxmap) &
                    + biasingdist(natom, weights, bias, auxmap)
            if (newdist > olddist) then
                call warning('new newdist is larger than old newdist!', 'mapatoms')
!                print *, olddist, newdist
            end if
            rotquat = leastrotquat(natom, weights, atoms0, atoms1, auxmap)
            prodquat = quatmul(rotquat, prodquat)
            call rotate(natom, rotquat, atoms1)
            iteration = iteration + 1
            meanrot = meanrot + (rotangle(rotquat) - meanrot)/iteration
            atomap = auxmap
        end do

        newdist = leastsquaredist(natom, weights, atoms0, atoms1, atomap)

! Check for new best mapping

        map_found = .false.

        do imap = 1, mapcount
            if (all(atomap == atomaplist(:, imap))) then
                if (imap == 1) matchcount = matchcount + 1
                matches(imap) = matches(imap) + 1
                avgiter(imap) = avgiter(imap) + (iteration - avgiter(imap))/matches(imap)
                avgtotrot(imap) = avgtotrot(imap) + (rotangle(prodquat) - avgtotrot(imap))/matches(imap)
                avgmeanrot(imap) = avgmeanrot(imap) + (meanrot - avgmeanrot(imap))/matches(imap)
                if (live) then
                    write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                    call print_stats(imap, earliest(imap), matches(imap), avgiter(imap), avgmeanrot(imap), &
                        avgtotrot(imap), mindist(imap))
                end if
                map_found = .true.
                exit
            end if
        end do

        if (.not. map_found) then
            if (mapcount >= maxrecord) then
                overflow = .true.
            end if
            do imap = 1, maxrecord
                if (imap > mapcount .or. newdist < mindist(imap)) then
                    if (imap == 1) matchcount = 1
                    if (mapcount < maxrecord) mapcount = mapcount + 1
                    do jmap = mapcount, imap + 1, -1
                        matches(jmap) = matches(jmap - 1)
                        avgiter(jmap) = avgiter(jmap - 1)
                        avgtotrot(jmap) = avgtotrot(jmap - 1)
                        avgmeanrot(jmap) = avgmeanrot(jmap - 1)
                        atomaplist(:, jmap) = atomaplist(:, jmap - 1)
                        mindist(jmap) = mindist(jmap - 1)
                        earliest(jmap) = earliest(jmap - 1)
                        if (live) then
                            write (output_unit, '(a)', advance='no') achar(27)//'['//str(jmap + 2)//'H'
                            call print_stats(jmap, earliest(jmap), matches(jmap), avgiter(jmap), avgmeanrot(jmap), &
                                avgtotrot(jmap), mindist(jmap))
                        end if
                    end do
                    matches(imap) = 1
                    avgiter(imap) = iteration
                    avgtotrot(imap) = rotangle(prodquat)
                    avgmeanrot(imap) = meanrot
                    atomaplist(:, imap) = atomap
                    mindist(imap) = newdist
                    earliest(imap) = ntrial
                    if (live) then
                        write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                        call print_stats(imap, earliest(imap), matches(imap), avgiter(imap), avgmeanrot(imap), &
                            avgtotrot(imap), mindist(imap))
                    end if
                    exit
                end if
            end do
        end if

        if (live) then
            write (output_unit, '(a)', advance='no') achar(27)//'['//str(mapcount + 3)//'H'
            call print_footer(.true., overflow, mapcount, ntrial)
        end if

! Restore original coordinates

        atoms1 = origmol1

    end do

    if (.not. live) then
        call print_header()
        do imap = 1, mapcount
            call print_stats(imap, earliest(imap), matches(imap), avgiter(imap), avgmeanrot(imap), &
                avgtotrot(imap), mindist(imap))
        end do
        call print_footer(.true., overflow, mapcount, ntrial)
    end if

end subroutine

end module
