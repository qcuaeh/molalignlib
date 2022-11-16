module remapping

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use linear
use random
use strutils
use assignment
use translation
use alignment
use rotation
use biasing
use printing

implicit none

private
public f_logintint
public lessthan
public everless
public optimize

abstract interface
    logical function f_logintint(x, y)
        integer, intent(in) :: x, y
    end function
end interface

contains

logical function lessthan(a, b) result(res)
    integer, intent(in) :: a, b
    res = a < b
end function

logical function everless(a, b) result(res)
    integer, intent(in) :: a, b
    res = .true.
end function

subroutine optimize( &
    natom, nblock, blocksize, weights, coords0, coords1, records, nrec, &
    maplist, mapcount, mindist, trialtest &
)

    integer, intent(in) :: natom, nblock, records
    integer, dimension(:), intent(in) :: blocksize
    real, dimension(:, :), intent(in) :: coords0, coords1
    real, dimension(:), intent(in) :: weights
    integer, intent(out) :: nrec
    integer, intent(out) :: maplist(:, :), mapcount(:)
    real, intent(out) :: mindist(:)
    procedure (f_logintint), pointer, intent(in) :: trialtest

    logical found, overflow
    integer imap, jmap, ntrial, nmatch, cycles
    integer, dimension(natom) :: atomap, auxmap
    real :: dist, biased_dist, new_biased_dist, totalrot
    real, dimension(4) :: rotquat, prodquat
    real, dimension(records) :: avgiter, avgtotalrot, avgrealrot
    real bias(natom, natom)
    real auxcoords(3, natom)
    real randpos(3)

! Set bias for non equivalent atoms 

    call setadjbias(natom, nblock, blocksize, coords0, coords1, bias)

! Print header and initial stats

    if (live_flag) then
        write (output_unit, '(a)', advance='no') achar(27)//'[1H'//achar(27)//'[J'
        call print_header()
    end if

! Initialize loop variables

    nrec = 0
    ntrial = 0
    nmatch = 0
    overflow = .false.

! Loop for map searching

    do while (nmatch < maxcount .and. trialtest(ntrial, maxtrials))

        ntrial = ntrial + 1

! Work with a copy of coords1

        auxcoords = coords1

! Apply random rotation

        call getrandnum(randpos)
        call rotate(natom, auxcoords, getrotquat(randpos))

! Minimize euclidean distance

        call assignatoms(natom, coords0, auxcoords, nblock, blocksize, bias, atomap)
        rotquat = leastrotquat(natom, weights, coords0, auxcoords, atomap)
        prodquat = rotquat
        totalrot = rotangle(rotquat)
        call rotate(natom, auxcoords, rotquat)
        cycles = 1

        do while (conv_flag)
            biased_dist = squaredist(natom, weights, coords0, auxcoords, atomap) &
                    + totalbias(natom, weights, bias, atomap)
            call assignatoms(natom, coords0, auxcoords, nblock, blocksize, bias, auxmap)
            if (all(auxmap == atomap)) exit
            new_biased_dist = squaredist(natom, weights, coords0, auxcoords, auxmap) &
                    + totalbias(natom, weights, bias, auxmap)
            if (new_biased_dist > biased_dist) then
                write (error_unit, '(a)') 'new_biased_dist is larger than biased_dist!'
!                print *, biased_dist, new_biased_dist
            end if
            rotquat = leastrotquat(natom, weights, coords0, auxcoords, auxmap)
            prodquat = quatmul(rotquat, prodquat)
            call rotate(natom, auxcoords, rotquat)
            cycles = cycles + 1
            totalrot = totalrot + rotangle(rotquat)
            atomap = auxmap
        end do

        dist = squaredist(natom, weights, coords0, auxcoords, atomap)

! Check for new best mapping

        found = .false.

        do imap = 1, nrec
            if (all(atomap == maplist(:, imap))) then
                if (imap == 1) nmatch = nmatch + 1
                mapcount(imap) = mapcount(imap) + 1
                avgiter(imap) = avgiter(imap) + (cycles - avgiter(imap))/mapcount(imap)
                avgrealrot(imap) = avgrealrot(imap) + (rotangle(prodquat) - avgrealrot(imap))/mapcount(imap)
                avgtotalrot(imap) = avgtotalrot(imap) + (totalrot - avgtotalrot(imap))/mapcount(imap)
                if (live_flag) then
                    write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                    call print_stats(imap, mapcount(imap), avgiter(imap), avgtotalrot(imap), &
                        avgrealrot(imap), mindist(imap))
                end if
                found = .true.
                exit
            end if
        end do

        if (.not. found) then
            if (nrec >= records) then
                overflow = .true.
            end if
            do imap = 1, records
                if (imap > nrec .or. dist < mindist(imap)) then
                    if (imap == 1) nmatch = 1
                    if (nrec < records) nrec = nrec + 1
                    do jmap = nrec, imap + 1, -1
                        mapcount(jmap) = mapcount(jmap - 1)
                        avgiter(jmap) = avgiter(jmap - 1)
                        avgrealrot(jmap) = avgrealrot(jmap - 1)
                        avgtotalrot(jmap) = avgtotalrot(jmap - 1)
                        maplist(:, jmap) = maplist(:, jmap - 1)
                        mindist(jmap) = mindist(jmap - 1)
                        if (live_flag) then
                            write (output_unit, '(a)', advance='no') achar(27)//'['//str(jmap + 2)//'H'
                            call print_stats(jmap, mapcount(jmap), avgiter(jmap), avgtotalrot(jmap), &
                                avgrealrot(jmap), mindist(jmap))
                        end if
                    end do
                    mapcount(imap) = 1
                    avgiter(imap) = cycles
                    avgrealrot(imap) = rotangle(prodquat)
                    avgtotalrot(imap) = totalrot
                    maplist(:, imap) = atomap
                    mindist(imap) = dist
                    if (live_flag) then
                        write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                        call print_stats(imap, mapcount(imap), avgiter(imap), avgtotalrot(imap), &
                            avgrealrot(imap), mindist(imap))
                    end if
                    exit
                end if
            end do
        end if

        if (live_flag) then
            write (output_unit, '(a)', advance='no') achar(27)//'['//str(nrec + 3)//'H'
            call print_footer(overflow, nrec, ntrial)
        end if

    end do

    if (.not. live_flag) then
        call print_header()
        do imap = 1, nrec
            call print_stats(imap, mapcount(imap), avgiter(imap), avgtotalrot(imap), &
                avgrealrot(imap), mindist(imap))
        end do
        call print_footer(overflow, nrec, ntrial)
    end if

end subroutine

end module
