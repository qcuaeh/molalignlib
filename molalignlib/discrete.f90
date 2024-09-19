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

module discrete

implicit none

contains

! Get an identity permutation
function identity_mapping(n) result(mapping)
   integer, intent(in) :: n
   ! Local variables
   integer :: i
   integer :: mapping(n)

   do i = 1, n
      mapping(i) = i
   end do

end function

! Get the inverse mapping of mapping
function inverse_mapping(mapping)
   integer, intent(in) :: mapping(:)
   ! Local variables
   integer :: i
   integer, allocatable :: inverse_mapping(:)

   allocate(inverse_mapping(size(mapping)))

   do i = 1, size(mapping)
      inverse_mapping(mapping(i)) = i
   end do

end function

function intersection(list1, list2, hash_size)
   integer, intent(in) :: hash_size
   integer, intent(in) :: list1(:), list2(:)
   ! Local variables
   integer :: i, ni
   integer, allocatable :: intersection(:)
   integer, allocatable :: intersection_buffer(:)
   logical, allocatable :: hash_table(:)

   allocate(hash_table(hash_size))
   allocate(intersection_buffer(min(size(list1), size(list2))))

   hash_table = .false.  ! Inicializar todos los elementos a falso

   ! Añadir elementos de list1 al conjunto
   do i = 1, size(list1)
      hash_table(list1(i)) = .true.
   end do

   ! Verificar elementos de list2 y construir la intersección
   ni = 0
   do i = 1, size(list2)
      if (hash_table(list2(i))) then
         ni = ni + 1
         intersection_buffer(ni) = list2(i)
      end if
   end do

   intersection = intersection_buffer(:ni)

end function

end module
