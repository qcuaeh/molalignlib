module assortment
use globals
use utilities
use sorting
use chemistry
implicit none
private
public grouplabels

contains

subroutine grouplabels(natom, label, nblock, blocksize, blocktype)
! Group atoms by label
    integer, intent(in) :: natom
    character(*), dimension(:), intent(in) :: label
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
            blocklabel(nblock) = label(i)
            blocktype(i) = nblock
            do j = i + 1, natom
                if (remaining(i)) then
                    if (label(j) == label(i)) then
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
!    print *, label(blocksizeorder)

end subroutine

end module
