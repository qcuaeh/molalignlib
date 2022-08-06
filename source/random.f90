module random

use options

implicit none

contains

function randvec3() result(x)
    real x(3)
    call random_number(x)
end function

subroutine shuffle(array, n)

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
      end do
   end do

end subroutine

subroutine init_random_seed(seed)

    integer n, u
    integer, allocatable :: seed(:)
    integer, parameter, dimension(12) :: test_seed = &
        [ 287027030, -719361131, 574274270, 292048305, &
          185733336, -1598963619, 572469522, 1446716853, &
          437591706, 1398099429, 570932571, -1177695979 ]

    if (testing) then
!        call random_seed(size=n)
!        allocate(seed(n))
!        call random_seed(get=seed)
        seed = test_seed
    else
        ! Read /dev/urandom to seed the random number generator
        call random_seed(size=n)
        allocate(seed(n))
        open(newunit=u, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old")
        read (u) seed
        close(u)
    end if
    call random_seed(put=seed)

end subroutine

end module

