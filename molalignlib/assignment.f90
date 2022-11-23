module assignment

use iso_fortran_env, only: output_unit
use iso_fortran_env, only: error_unit

use options
use linear
use random
use strutils
use translation
use alignment
use rotation
use biasing
use printing
use hungarian

implicit none

private
public f_logintint
public lessthan
public everless
public optimize_assignment

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

subroutine optimize_assignment( &
    natom, &
    nblock, &
    blocksize, &
    weights, &
    coords0, &
    coords1, &
    nrec, &
    nmap, &
    maplist, &
    countlist, &
    rmsdlist, &
    stoptest &
)

    integer, intent(in) :: natom, nblock, nrec
    integer, dimension(:), intent(in) :: blocksize
    real, dimension(:, :), intent(in) :: coords0, coords1
    real, dimension(:), intent(in) :: weights
    integer, intent(out) :: nmap
    integer, intent(out) :: maplist(:, :)
    integer, intent(out) :: countlist(:)
    real, intent(out) :: rmsdlist(:)
    procedure (f_logintint), pointer, intent(in) :: stoptest

    logical found, overflow
    integer imap, jmap, ntrial, nmatch, cycles
    integer, dimension(natom) :: atomap, auxmap
    real :: sqdist, biased_dist, new_biased_dist, totalrot
    real, dimension(4) :: rotquat, prodquat
    real, dimension(nrec) :: avgiter, avgtotalrot, avgrealrot
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

    nmap = 0
    ntrial = 0
    nmatch = 0
    overflow = .false.

! Loop for map searching

    do while (nmatch < max_count .and. stoptest(ntrial, max_trials))

        ntrial = ntrial + 1

! Work with a copy of coords1

        auxcoords = coords1

! Apply random rotation

        call getrandnum(randpos)
        call rotate(natom, auxcoords, getrotquat(randpos))

! Minimize euclidean distance

        call localassign(natom, coords0, auxcoords, nblock, blocksize, bias, atomap)
        rotquat = leastrotquat(natom, weights, coords0, auxcoords, atomap)
        prodquat = rotquat
        totalrot = rotangle(rotquat)
        call rotate(natom, auxcoords, rotquat)
        cycles = 1

        do while (iter_flag)
            biased_dist = squaredist(natom, weights, coords0, auxcoords, atomap) &
                    + totalbias(natom, weights, bias, atomap)
            call localassign(natom, coords0, auxcoords, nblock, blocksize, bias, auxmap)
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

        sqdist = squaredist(natom, weights, coords0, auxcoords, atomap)

! Check for new best mapping

        found = .false.

        do imap = 1, nmap
            if (all(atomap == maplist(:, imap))) then
                if (imap == 1) nmatch = nmatch + 1
                countlist(imap) = countlist(imap) + 1
                avgiter(imap) = avgiter(imap) + (cycles - avgiter(imap))/countlist(imap)
                avgrealrot(imap) = avgrealrot(imap) + (rotangle(prodquat) - avgrealrot(imap))/countlist(imap)
                avgtotalrot(imap) = avgtotalrot(imap) + (totalrot - avgtotalrot(imap))/countlist(imap)
                if (live_flag) then
                    write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                    call print_stats(imap, countlist(imap), avgiter(imap), avgtotalrot(imap), &
                        avgrealrot(imap), rmsdlist(imap))
                end if
                found = .true.
                exit
            end if
        end do

        if (.not. found) then
            if (nmap >= nrec) then
                overflow = .true.
            end if
            do imap = 1, nrec
                if (imap > nmap .or. sqrt(sqdist) < rmsdlist(imap)) then
                    if (imap == 1) nmatch = 1
                    if (nmap < nrec) nmap = nmap + 1
                    do jmap = nmap, imap + 1, -1
                        countlist(jmap) = countlist(jmap - 1)
                        avgiter(jmap) = avgiter(jmap - 1)
                        avgrealrot(jmap) = avgrealrot(jmap - 1)
                        avgtotalrot(jmap) = avgtotalrot(jmap - 1)
                        maplist(:, jmap) = maplist(:, jmap - 1)
                        rmsdlist(jmap) = rmsdlist(jmap - 1)
                        if (live_flag) then
                            write (output_unit, '(a)', advance='no') achar(27)//'['//str(jmap + 2)//'H'
                            call print_stats(jmap, countlist(jmap), avgiter(jmap), avgtotalrot(jmap), &
                                avgrealrot(jmap), rmsdlist(jmap))
                        end if
                    end do
                    countlist(imap) = 1
                    avgiter(imap) = cycles
                    avgrealrot(imap) = rotangle(prodquat)
                    avgtotalrot(imap) = totalrot
                    maplist(:, imap) = atomap
                    rmsdlist(imap) = sqrt(sqdist)
                    if (live_flag) then
                        write (output_unit, '(a)', advance='no') achar(27)//'['//str(imap + 2)//'H'
                        call print_stats(imap, countlist(imap), avgiter(imap), avgtotalrot(imap), &
                            avgrealrot(imap), rmsdlist(imap))
                    end if
                    exit
                end if
            end do
        end if

        if (live_flag) then
            write (output_unit, '(a)', advance='no') achar(27)//'['//str(nmap + 3)//'H'
            call print_footer(overflow, nmap, ntrial)
        end if

    end do

    if (.not. live_flag) then
        call print_header()
        do imap = 1, nmap
            call print_stats(imap, countlist(imap), avgiter(imap), avgtotalrot(imap), &
                avgrealrot(imap), rmsdlist(imap))
        end do
        call print_footer(overflow, nmap, ntrial)
    end if

end subroutine

! Find best correspondence between points sets with fixed orientation
subroutine localassign(natom, coords0, coords1, nblock, blocksize, bias, atomap)

! nblock: Number of block atoms
! blocksize: Number of atoms in each block
! atomap: Map between correspondent points in the adjmat
! offset: First element of current block

    integer, intent(in) :: natom, nblock
    integer, dimension(:), intent(in) :: blocksize
    real, dimension(:, :), intent(in) :: coords0
    real, dimension(:, :), intent(in) :: bias
    real, dimension(:, :), intent(inout) :: coords1
    integer, dimension(:), intent(out) :: atomap

    integer h, i, j, offset
    integer, dimension(natom) :: blockmap
    real, dimension(natom, natom) :: costs
    real blocksum

! Fill distance matrix for each block

    offset = 0

    do h = 1, nblock

        do i = offset + 1, offset + blocksize(h)
            do j = offset + 1, offset + blocksize(h)
                costs(i - offset, j - offset) = sum((coords1(:, j) - coords0(:, i))**2) + bias(i, j)
            end do
        end do

! Find correspondence between points in the current block

        call assndx(1, costs, blocksize(h), blocksize(h), blockmap, blocksum)
        atomap(offset+1:offset+blocksize(h)) = blockmap(:blocksize(h)) + offset

        offset = offset + blocksize(h)

    end do

end subroutine

! Calculate total bias
real function totalbias(natom, weights, bias, mapping) result(total)

    integer, intent(in) :: natom
    real, dimension(:), intent(in) :: weights
    integer, dimension(:), intent(in) :: mapping
    real, dimension(:, :), intent(in) :: bias
    integer i

    total = 0.

    do i = 1, natom
        total = total + weights(i)*bias(i, mapping(i))
    end do

end function

end module
