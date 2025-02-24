module common_types
use parameters

implicit none

type, public :: intlist_type
   integer, allocatable :: n(:)
end type

type, public :: reallist_type
   real(rk), allocatable :: x(:)
end type

type, public :: nested_intlist_type
   type(intlist_type), allocatable :: s(:)
end type

type, public :: nested_reallist_type
   type(reallist_type), allocatable :: s(:)
end type

type :: boolmatrix_type
   logical, allocatable :: b(:, :)
end type

type :: intmatrix_type
   integer, allocatable :: n(:, :)
end type

type :: realmatrix_type
   real(rk), allocatable :: x(:, :)
end type

end module
