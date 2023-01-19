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

module chemutils
use parameters
use chemdata
use strutils

implicit none

private
public f_realint
public unity
public stdmass
public valency
public readlabel

abstract interface
   real(wp) function f_realint(z)
      use parameters
      integer, intent(in) :: z
   end function
end interface

contains

real(wp) function unity(z) result(res)
   integer, intent(in) :: z
   res = 1.
end function

real(wp) function stdmass(z) result(res)
   integer, intent(in) :: z
   res = stdmasses(z)
end function

real(wp) function valency(z) result(res)
   integer, intent(in) :: z
   res = real(valencies(z), wp)
end function

subroutine readlabel(label, znum, type)
   character(*), intent(in) :: label
   integer, intent(out) :: znum, type
   integer :: m, n

   n = len_trim(label)
   m = verify(uppercase(trim(label)), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

   if (m == 0) then
      type = 0
      znum = atomic_number(label)
   else if (verify(label(m:n), '1234567890') == 0) then
      read (label(m:n), *) type
      znum = atomic_number(label(1:m-1))
   else
      type = -1
      znum = atomic_number(label(1:m-1))
   end if

end subroutine

function atomic_number(symbol) result(z)
   character(*), intent(in) :: symbol
   integer :: z

   do z = 1, nelem
      if (uppercase(symbol) == uppercase(elsym(z))) then
         return
      end if
   end do

   select case (uppercase(symbol))
   case ('LJ')
      z = 1001
   case default
      z = -1
   end select

end function

end module
