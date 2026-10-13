! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module common_types
use parameters
implicit none

type, public :: int_list
   integer(IK), dimension(:), allocatable :: u
end type

type, public :: real_list
   real(RK), dimension(:), allocatable :: u
end type

type :: int_matrix
   integer(IK), dimension(:,:), allocatable :: a
end type

type :: bool_matrix
   logical(LK), dimension(:,:), allocatable :: a
end type

type :: real_matrix
   real(RK), dimension(:,:), allocatable :: a
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

! Part of a flat partition: atoms of each molecule in the part
type, public :: partition_part_t
   integer(IK) :: elnum
   integer(IK) :: n_items1
   integer(IK) :: n_items2
   integer(IK) :: n_children
   integer(IK), dimension(:), allocatable :: items1
   integer(IK), dimension(:), allocatable :: items2
   integer(IK), dimension(:), allocatable :: signature
   integer(IK), dimension(:), allocatable :: children
end type

! Flat partition; itemdir1/2 give the part of each atom
type, public :: partition_t
   integer(IK) :: n_parts
   type(partition_part_t), dimension(:), allocatable :: parts
   integer(IK), dimension(:), allocatable :: itemdir1
   integer(IK), dimension(:), allocatable :: itemdir2
end type

end module
