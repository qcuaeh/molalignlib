module remapping

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use math
use random
use strutils
use assignment
use translation
use alignment
use messages
use rotation
use biasing
use printing

implicit none

abstract interface
    logical function generic_test(x, y)
        integer, intent(in) :: x, y
    end function
end interface

contains

logical function trueconst(counter, maxcount)
    integer, intent(in) :: counter, maxcount
    trueconst = .true.
end function

logical function lowerthan(counter, maxcount)
    integer, intent(in) :: counter, maxcount
    lowerthan = counter < maxcount
end function

subroutine remapatoms(natom, nblock, blocksize, weights, atoms0, atoms1, maxrecord, nrecord, &
                      atomaplist, trial_test, match_test)

    integer, intent(in) :: natom, nblock, maxrecord
    integer, dimension(:), intent(in) :: blocksize
    real(wp), dimension(:, :), intent(in) :: atoms0
    real(wp), dimension(:, :), intent(in) :: atoms1
    real(wp), dimension(:), intent(in) :: weights
    procedure (generic_test), pointer, intent(in) :: trial_test, match_test
    integer, intent(out) :: nrecord
    integer, dimension(:, :), intent(out) :: atomaplist

    logical map_found, overflow
    integer i, imap, jmap, ntrial, nmatch, iteration, newdiff
    integer, dimension(natom) :: atomap, auxmap, equivset0, equivset1
    integer, dimension(maxrecord) :: earliest, matches
    real(wp) u, dist, biased_dist, new_biased_dist, meanrot
    real(wp), dimension(natom, natom) :: bias
    real(wp), dimension(4) :: rotquat, prodquat
    real(wp), dimension(maxrecord) :: mindist, avgiter, avgmeanrot, avgtotrot
    real(wp), dimension(3, natom) :: atoms01

! Set bias for non equivalent atoms 

call setadjbias(natom, nblock, blocksize, atoms0, atoms1, bias)

! Print header and initial stats

    if (live) then
        write (output_unit, '(a)', advance='no') achar(27)//'[1H'//achar(27)//'[J'
        call print_header()
    end if

! Initialize loop variables

    nrecord = 0
    ntrial = 0
    nmatch = 0
    overflow = .false.

! Loop for map searching

    do while (trial_test(ntrial, maxtrial) .and. match_test(nmatch, maxmatch))

        ntrial = ntrial + 1

! Work with a copy of atoms1

        atoms01 = atoms1

! Apply random rotation

        call rotate(natom, getrotquat(randvec3()), atoms01)

! Minimize euclidean distance

        call assignatoms(natom, weights, atoms0, atoms01, nblock, blocksize, bias, atomap)
        rotquat = leastrotquat(natom, weights, atoms0, atoms01, atomap)
        prodquat = rotquat
        meanrot = rotangle(rotquat)
        call rotate(natom, rotquat, atoms01)
        iteration = 1

        do while (iterative)
            biased_dist = squaredist(natom, weights, atoms0, atoms01, atomap) &
                    + totalbias(natom, weights, bias, atomap)
            call assignatoms(natom, weights, atoms0, atoms01, nblock, blocksize, bias, auxmap)
            if (all(auxmap == atomap)) exit
            new_biased_dist = squaredist(natom, weights, atoms0, atoms01, auxmap) &
                    + totalbias(natom, weights, bias, auxmap)
            if (new_biased_dist > biased_dist) then
                call warning('new biased dist is larger than previous biased dist!', 'mapatoms')
!                print *, biased_dist, new_biased_dist
            end if
            rotquat = leastrotquat(natom, weights, atoms0, atoms01, auxmap)
            prodquat = quatmul(rotquat, prodquat)
            call rotate(natom, rotquat, atoms01)
            iteration = iteration + 1
            meanrot = meanrot + (rotangle(rotquat) - meanrot)/iteration
            atomap = auxmap
        end do

        dist = leastsquaredist(natom, weights, atoms0, atoms01, atomap)

! Check for new best mapping

        map_found = .false.

        do imap = 1, nrecord
            if (all(atomap == atomaplist(:, imap))) then
                if (imap == 1) nmatch = nmatch + 1
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
            if (nrecord >= maxrecord) then
                overflow = .true.
            end if
            do imap = 1, maxrecord
                if (imap > nrecord .or. dist < mindist(imap)) then
                    if (imap == 1) nmatch = 1
                    if (nrecord < maxrecord) nrecord = nrecord + 1
                    do jmap = nrecord, imap + 1, -1
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
                    mindist(imap) = dist
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
            write (output_unit, '(a)', advance='no') achar(27)//'['//str(nrecord + 3)//'H'
            call print_footer(.true., overflow, nrecord, ntrial)
        end if

    end do

    if (.not. live) then
        call print_header()
        do imap = 1, nrecord
            call print_stats(imap, earliest(imap), matches(imap), avgiter(imap), avgmeanrot(imap), &
                avgtotrot(imap), mindist(imap))
        end do
        call print_footer(.true., overflow, nrecord, ntrial)
    end if

end subroutine

end module