module partition
use stdio

implicit none

type, public :: atomlist_type
   integer, allocatable :: atomidcs(:)
end type

type, public :: atompartition_type
   integer :: natom
   integer :: largest
   type(atomlist_type), allocatable :: subsets(:)
contains
   procedure :: get_atomtypes
   procedure :: mapped
end type

interface operator (==)
   module procedure is_equal
end interface

contains

function atompartition(ntype, atomtypes)
   integer, intent(in) :: ntype
   integer, intent(in) :: atomtypes(:)
   ! Result variable
   type(atompartition_type) :: atompartition
   ! Local variables
   integer :: h
   integer, allocatable :: n(:)

   allocate (n(ntype))
   allocate (atompartition%subsets(ntype))

   n(:) = 0
   do h = 1, size(atomtypes)
      n(atomtypes(h)) = n(atomtypes(h)) + 1
   end do

   atompartition%natom = sum(n)
   atompartition%largest = maxval(n)

   do h = 1, ntype
      allocate (atompartition%subsets(h)%atomidcs(n(h)))
   end do

   n(:) = 0
   do h = 1, size(atomtypes)
      n(atomtypes(h)) = n(atomtypes(h)) + 1
      atompartition%subsets(atomtypes(h))%atomidcs(n(atomtypes(h))) = h
   end do

end function

function get_atomtypes(self) result(atomtypes)
   class(atompartition_type), intent(in) :: self
   ! Result variable
   integer, allocatable :: atomtypes(:)
   ! Local variables
   integer :: h, i
   integer :: atomidx_i

   allocate (atomtypes(self%natom))

   do h = 1, size(self%subsets)
      do i = 1, size(self%subsets(h)%atomidcs)
         atomidx_i = self%subsets(h)%atomidcs(i)
         atomtypes(atomidx_i) = h
      end do
   end do

end function

function is_equal(self, other)
   type(atompartition_type), intent(in) :: self, other
   ! Result variable
   logical :: is_equal
   ! Local variables
   integer :: h

   if (size(self%subsets) /= size(other%subsets)) then
      is_equal = .false.
      return
   end if

   do h = 1, size(self%subsets)
      if (size(self%subsets(h)%atomidcs) /= size(other%subsets(h)%atomidcs)) then
         is_equal = .false.
         return
      end if
   end do

   is_equal = .true.

end function

function mapped(self, mapping) result(atompartition)
   class(atompartition_type), intent(in) :: self
   integer, intent(in) :: mapping(:)
   ! Result variable
   type(atompartition_type) :: atompartition
   ! Local variables
   integer :: h

   atompartition = self
   do h = 1, size(self%subsets)
      atompartition%subsets(h)%atomidcs = mapping(self%subsets(h)%atomidcs)
   end do

end function

end module
