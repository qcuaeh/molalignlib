module random

use iso_fortran_env, only: output_unit

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

    call random_seed(size=n)
    allocate(seed(n))

    if (.not. testing) then
        ! Read /dev/urandom to seed the random number generator
        open(newunit=u, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old")
        read (u) seed
        close(u)
        call random_seed(put=seed)
    end if

    if (debug) then
        call random_seed(get=seed)
        write (output_unit, '(a)') 'Random seed:'
        write (output_unit, *) n
        write (output_unit, *) seed
        write (output_unit, *)
    end if

end subroutine

end module

