module maputils

use options

implicit none

private

public resetmap
public identitymap
public inversemap

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

function identitymap(n)
    integer, intent(in) :: n
    integer identitymap(n)
    integer i
    do i = 1, n
        identitymap(i) = i
    end do
end function

subroutine resetmap(mapping)
    integer, intent(out) :: mapping(:)
    integer i
    do i = 1, size(mapping)
        mapping(i) = i
    end do
end subroutine

end module
