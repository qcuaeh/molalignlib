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
use parameters
use settings
use rnglib

implicit none

abstract interface
   subroutine f_sp(x)
      use parameters
      real(sp) :: x
   end subroutine
end interface

abstract interface
   subroutine f_dp(x)
      use parameters
      real(dp) :: x
   end subroutine
end interface

interface randnum
   procedure randnum_sp
   procedure randnum_dp
end interface

procedure(f_sp), pointer :: randnum_sp
procedure(f_dp), pointer :: randnum_dp

private
public rand3
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

   if (repro_flag) then
      randnum_sp => real_uni01_sp
      randnum_dp => real_uni01_dp
      call rnglib_init()
   else
      randnum_sp => random_number_sp
      randnum_dp => random_number_dp
      call random_seed(size=n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=unit, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=stat)
      if (stat == 0) then
         read(unit) seed
         close(unit)
      else
         write (error_unit, '(a)') 'Error: can not read from /dev/urandom'
         stop
      end if
      call random_seed(put=seed)
   end if

end subroutine initialize_random

function rand3() result(x)
   integer :: i
   real(wp) :: x(3)
   do i = 1, 3
      call randnum(x(i))
   end do
end function

subroutine shuffle(array, n)
   integer, intent(in) :: n
   integer, intent(inout) :: array(:)
   integer :: i, j, k, temp
   real :: x
   do k = 1, 2
      do i = 1, n
         call randnum(x)
         j = floor(n*x) + 1
         temp = array(j)
         array(j) = array(i)
         array(i) = temp
      end do
   end do
end subroutine

end module
