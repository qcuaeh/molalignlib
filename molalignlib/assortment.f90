module assortment

use iso_fortran_env, only: error_unit

use options
use sorting
use maputils
use chemdata

implicit none

private
public getblocks

contains

subroutine getblocks(natom, znums, types, nblock, blocksize, blockindex, atomorder)
! Purpose: Group atoms by atomic numbers and types

    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: znums, types
    integer, intent(out) :: nblock
    integer, intent(out) :: atomorder(:)
    integer, dimension(:), intent(out) :: blocksize, blockindex

    integer i, j
    logical remaining(natom)
    integer blockznum(natom)
    integer blocktype(natom)
    integer blockorder(natom)

! Initialization

    nblock = 0
    remaining = .true.

! Create block list

    do i = 1, natom
        if (remaining(i)) then
            nblock = nblock + 1
            blocksize(nblock) = 1
            blockznum(nblock) = znums(i)
            blocktype(nblock) = types(i)
            blockindex(i) = nblock
            do j = i + 1, natom
                if (remaining(i)) then
                    if (znums(j) == znums(i) .and. types(j) == types(i)) then
                        blockindex(j) = nblock
                        blocksize(nblock) = blocksize(nblock) + 1
                        remaining(j) = .false.
                    end if
                end if
            end do
        end if
    end do

! Order blocks by atomic type

    blockorder(:nblock) = order(blocktype, nblock)
    blocksize(:nblock) = blocksize(blockorder(:nblock))
    blockorder(:nblock) = inversemap(blockorder(:nblock))
    blockindex = blockorder(blockindex)

! Order blocks by atomic number

    blockorder(:nblock) = order(blockznum, nblock)
    blocksize(:nblock) = blocksize(blockorder(:nblock))
    blockorder(:nblock) = inversemap(blockorder(:nblock))
    blockindex = blockorder(blockindex)

! Order blocks by block size

!    blockorder(:nblock) = order(blocksize, nblock)
!    blocksize(:nblock) = blocksize(blockorder(:nblock))
!    blockorder(:nblock) = inversemap(blockorder(:nblock))
!    blockindex = blockorder(blockindex)

! Save atom order

    atomorder = order(blockindex, natom)

end subroutine

end module
