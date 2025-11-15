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

module types_basic
use parameters
implicit none

type, public :: int_list
   integer, dimension(:), allocatable :: u
end type

type, public :: real_list
   real(rk), dimension(:), allocatable :: u
end type

type :: int_matrix
   integer, dimension(:,:), allocatable :: a
end type

type :: bool_matrix
   logical, dimension(:,:), allocatable :: a
end type

type :: real_matrix
   real(rk), dimension(:,:), allocatable :: a
end type

type, public :: int_listlist
   type(int_list), dimension(:), allocatable :: u
end type

type, public :: real_listlist
   type(real_list), dimension(:), allocatable :: u
end type

type, public :: int_listmatrix
   type(int_list), dimension(:,:), allocatable :: a
end type

! Part array
type, public :: partition_part_t
   integer :: elnum
   integer :: num_items1
   integer :: num_items2
   integer :: num_children
   integer, dimension(:), allocatable :: items1
   integer, dimension(:), allocatable :: items2
   integer, dimension(:), allocatable :: signature
   integer, dimension(:), allocatable :: children
end type

! Partition array
type, public :: partition_t
   integer :: num_parts
   type(partition_part_t), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir1
   integer, dimension(:), allocatable :: itemdir2
end type

end module
