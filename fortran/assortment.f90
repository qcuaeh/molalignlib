module assortment

use common
use options
use sorting
use maputils
use chemistry

implicit none

private
public grouplabels

contains

subroutine grouplabels(natom, labels, nblock, blocksize, blocktype)
! Group atoms by labels
    integer, intent(in) :: natom
    character(*), dimension(:), intent(in) :: labels
    integer, intent(out) :: nblock
    integer, dimension(:), intent(out) :: blocksize, blocktype

    integer i, j
    logical remaining(natom)
    integer blockorder(natom)
    character(32) blocklabel(natom)

! Initialization

    nblock = 0
    remaining = .true.

! Create block list

    do i = 1, natom
        if (remaining(i)) then
            nblock = nblock + 1
            blocksize(nblock) = 1
            blocklabel(nblock) = labels(i)
            blocktype(i) = nblock
            do j = i + 1, natom
                if (remaining(i)) then
                    if (labels(j) == labels(i)) then
                        blocktype(j) = nblock
                        blocksize(nblock) = blocksize(nblock) + 1
                        remaining(j) = .false.
                    end if
                end if
            end do
        end if
    end do

! Order blocks alphabetically

    blockorder(:nblock) = sortorder(blocklabel, nblock)
    blocksize(:nblock) = blocksize(blockorder(:nblock))
    blockorder(:nblock) = inversemap(blockorder(:nblock))
    blocktype = blockorder(blocktype)

! Order blocks by category size

    blockorder(:nblock) = sortorder(blocksize, nblock)
    blocksize(:nblock) = blocksize(blockorder(:nblock))
    blockorder(:nblock) = inversemap(blockorder(:nblock))
    blocktype = blockorder(blocktype)

!    print *, blocksize(:nblock)
!    print *, labels(blocksizeorder)

end subroutine

end module
