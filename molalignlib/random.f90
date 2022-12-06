module random
use parameters
use settings
use rnglib

implicit none

private
public rand3
public shuffle
public initialize_random

contains

function rand3() result(x)
   real(wp) x(3)
!   call random_number(x)
   call rnglib_number(x)
end function

subroutine initialize_random()
!   integer i, n
!   integer, allocatable :: seed(:)
   integer :: u, stat, dt(8)
   integer :: seed1, seed2
   integer(int64) :: t

   call rnglib_init()
   if (.not. test_flag) then
!      call random_seed(size=n)
!      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=u, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=stat)
      if (stat == 0) then
!         read(u) seed
         read(u) seed1
         read(u) seed2
         close(u)
      else
         ! Fallback to XORing the current time
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
         end if
!         do i = 1, n
!            seed(i) = lcg(t)
!         end do
         seed1 = lcg(t)
         seed2 = lcg(t)
      end if
!      call random_seed(put=seed)
      call rnglib_seed(seed1, seed2)
   end if
end subroutine initialize_random

! This simple PRNG might not be good enough for real(wp) work, but is
! sufficient for seeding a better PRNG.
integer function lcg(s)
   integer(int64) :: s
   if (s == 0) then
      s = 104729
   else
      s = mod(s, 4294967296_int64)
   end if
   s = mod(s * 279470273_int64, 4294967291_int64)
   lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg

subroutine shuffle(array, n)
   integer, intent(in) :: n
   integer, dimension(:), intent(inout) :: array(:)

   integer :: i, j, k, temp
   real(wp) :: u

   do k = 1, 2
      do i = 1, n
!         call random_number(u)
         call rnglib_number(u)
         j = floor(n*u) + 1
         ! switch values
         temp = array(j)
         array(j) = array(i)
         array(i) = temp
      end do
   end do

end subroutine

end module
