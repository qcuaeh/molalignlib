module maputils
use sorting
implicit none

contains

function inversemap(mapping)
! Get the inverse mapping of mapping
    integer, intent(in), dimension(:) :: mapping
    integer, dimension(size(mapping)) :: inversemap
    integer i
    do i = 1, size(mapping)
        inversemap(mapping(i)) = i
    end do
end function

logical function propermap(mapping, n)
! Check if mapping is a proper mapping
    integer, intent(in) :: n, mapping(:)
    integer i
    propermap = .false.
    if (any(sorted(mapping, n) /= [(i, i=1, n)])) then
        propermap = .true.
    end if
end function

end module
