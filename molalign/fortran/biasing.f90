module biasing

use options
use sorting

implicit none

contains

subroutine setadjbias(natom, nblock, blocksize, coords0, coords1, bias)
    integer, intent(in) :: natom, nblock
    integer, dimension(:), intent(in) :: blocksize
    real, dimension(:, :), intent(in) :: coords0, coords1
    real, dimension(:, :), intent(out) :: bias

    integer h, i, j, offset
    real d0(natom, natom), d1(natom, natom)

    do i = 1, natom
        offset = 0
        do h = 1, nblock
            do j = offset + 1, offset + blocksize(h)
                d0(j, i) = sqrt(sum((coords0(:, j) - coords0(:, i))**2))
            end do
            call sort(d0(:, i), offset + 1, offset + blocksize(h))
            offset = offset + blocksize(h)
        end do
    end do

    do i = 1, natom
        offset = 0
        do h = 1, nblock
            do j = offset + 1, offset + blocksize(h)
                d1(j, i) = sqrt(sum((coords1(:, j) - coords1(:, i))**2))
            end do
            call sort(d1(:, i), offset + 1, offset + blocksize(h))
            offset = offset + blocksize(h)
        end do
    end do

    offset = 0
    bias(:, :) = 0.

    if (bias_flag) then
        do h = 1, nblock
            do i = offset + 1, offset + blocksize(h)
                do j = offset + 1, offset + blocksize(h)
                    if (any(abs(d1(:, j) - d0(:, i)) > bias_tol)) then
                        bias(i, j) = bias_scale**2
                    end if
                end do
            end do
            offset = offset + blocksize(h)
        end do
    end if

end subroutine

end module
