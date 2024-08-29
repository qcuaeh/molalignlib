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
use kinds
use pointers
use chemdata
use strutils

implicit none

private
public unity
public stdmass
public valency
public readlabel

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

subroutine readlabel(element, atomelnum, atomlabel)
   character(*), intent(in) :: element
   integer, intent(out) :: atomelnum, atomlabel
   integer :: m, n

   n = len_trim(element)
   m = verify(uppercase(trim(element)), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

   if (m == 0) then
      atomlabel = 0
      atomelnum = atomic_number(element)
   else if (verify(element(m:n), '1234567890') == 0) then
      read (element(m:n), *) atomlabel
      atomelnum = atomic_number(element(1:m-1))
   else
      atomlabel = -1
      atomelnum = atomic_number(element(1:m-1))
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
