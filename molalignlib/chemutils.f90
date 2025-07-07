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
public split_symbol
public atomic_number

contains

function atomic_number(elsym) result(z)
   character(*), intent(in) :: elsym
   integer :: z

   do z = 1, num_elems
      if (uppercase(elsym) == uppercase(element_symbols(z))) then
         return
      end if
   end do

   select case (uppercase(elsym))
   case ('LJ')
      z = 1001
   case default
      z = -1
   end select

end function

subroutine split_symbol( elsym, elnum, label)
   character(*), intent(in) :: elsym
   integer, intent(out) :: elnum, label
   ! Local variables
   integer :: m, n

   n = len_trim(elsym)
   m = verify(uppercase(trim(elsym)), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

   if (m == 0) then
      label = 0
      elnum = atomic_number(elsym)
   else if (verify(elsym(m:n), '1234567890') == 0) then
      read (elsym(m:n), *) label
      elnum = atomic_number(elsym(1:m-1))
   else
      write (stderr, '(a,1x,a)') 'Invalid atomic elsym/label:', elsym
      stop
   end if

end subroutine

end module
