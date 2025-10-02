module derived_types
use parameters
implicit none

type, public :: int_list
   integer, dimension(:), allocatable :: e
end type

type, public :: real_list
   real(rk), dimension(:), allocatable :: e
end type

type :: int_matrix
   integer, dimension(:,:), allocatable :: ee
end type

type :: bool_matrix
   logical, dimension(:,:), allocatable :: ee
end type

type :: real_matrix
   real(rk), dimension(:,:), allocatable :: ee
end type

type, public :: int_listlist
   type(int_list), dimension(:), allocatable :: e
end type

type, public :: real_listlist
   type(real_list), dimension(:), allocatable :: e
end type

type, public :: int_listmatrix
   type(int_list), dimension(:,:), allocatable :: ee
end type

! Semipart array
type, public :: semipartition_part_t
   integer :: num_items
   integer, dimension(:), allocatable :: items
end type

! Semipartition array
type, public :: semipartition_t
   integer :: num_parts
   type(semipartition_part_t), dimension(:), allocatable :: parts
   integer, dimension(:), allocatable :: itemdir
end type

! Part array
type, public :: partition_part_t
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
