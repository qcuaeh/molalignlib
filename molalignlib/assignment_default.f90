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

module assignment_default
use molecule
use permutation
implicit none

contains

function default_atomperm(atoms1, atoms2) result(subperm)
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(subperm_t), target :: subperm
   integer, pointer :: n
   integer :: i

   call subperm_init( subperm, size(atoms1))
   n => subperm%count
   do i = 1, size(atoms1)
      if (atoms1(i)%mask) then
         if (atoms2(i)%mask) then
            n = n + 1
            subperm%subset(n) = i
            subperm%permut(n) = i
         else
            write (stderr, '(A)') "Error: Aligning atoms don't match"
            stop 1
         end if
      end if
   end do
end function

end module
