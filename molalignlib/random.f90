! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module random
use stdio
use kinds
use flags
use rnglib

implicit none

abstract interface
   subroutine f_sp(x)
      use kinds
      real(sp) :: x
   end subroutine
end interface

abstract interface
   subroutine f_dp(x)
      use kinds
      real(dp) :: x
   end subroutine
end interface

interface randnum01
   procedure randnum01_sp
   procedure randnum01_dp
end interface

procedure(f_sp), pointer :: randnum01_sp
procedure(f_dp), pointer :: randnum01_dp

private
public randarray
public shuffle
public initialize_random

contains

subroutine random_number_sp(x)
   real(sp) :: x
   call random_number(x)
end subroutine

subroutine random_number_dp(x)
   real(dp) :: x
   call random_number(x)
end subroutine

subroutine initialize_random()
   integer n
   integer, allocatable :: seed(:)
   integer :: unit, stat

   if (test_flag) then
!      randnum01_sp => random_number_sp
!      randnum01_dp => random_number_dp
      randnum01_sp => sp_uni_01
      randnum01_dp => dp_uni_01
      call rnglib_init()
   else
      randnum01_sp => random_number_sp
      randnum01_dp => random_number_dp
      call random_seed(size=n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=unit, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=stat)
      if (stat == 0) then
         read(unit) seed
         close(unit)
      else
         write (stderr, '(a)') 'Error: can not read from /dev/urandom'
         stop
      end if
      call random_seed(put=seed)
   end if

end subroutine

function randarray(n) result(array)
   integer, intent(in) :: n
   real(wp) :: array(n)
   integer :: i
   do i = 1, n
      call randnum01(array(i))
   end do
end function

subroutine shuffle(array, n)
   integer, intent(in) :: n
   integer, intent(inout) :: array(:)
   integer :: i, j, k, temp
   real :: x
   do k = 1, 2
      do i = 1, n
         call randnum01(x)
         j = floor(n*x) + 1
         temp = array(j)
         array(j) = array(i)
         array(i) = temp
      end do
   end do
end subroutine

end module
