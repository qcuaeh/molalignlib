module random

use options

implicit none

interface shuffle
   module procedure intshuffle
end interface

contains

function randvec3() result(x)
    real x(3)
    call random_number(x)
end function

subroutine intshuffle(array, n)
   integer, intent(in) :: n
   integer, dimension(:), intent(inout) :: array(:)
   integer i, j, k, temp
   real u

   do k = 1, 2
      do i = 1, n
         call random_number(u)
         j = floor(n*u) + 1
         ! switch values
         temp = array(j)
         array(j) = array(i)
         array(i) = temp
      enddo
   enddo

end subroutine

subroutine init_random_seed(seed)
    integer n
    integer, allocatable :: seed(:)
    if (testing) then
        return
    else
        ! Read /dev/urandom to seed the random number generator
        call random_seed(size=n)
        allocate(seed(n))
        open(unit=99, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old")
        read (99) seed
        close(99)
        call random_seed(put=seed)
    end if

end subroutine

end module

