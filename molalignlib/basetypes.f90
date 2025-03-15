module basetypes
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

end module
