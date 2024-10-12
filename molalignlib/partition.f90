module partition
use stdio

implicit none

type, public :: new_atomlist_type
   integer, pointer :: atomidcs(:) => null()
   integer, pointer :: atomidcs_allocation(:)
contains
   procedure :: append
end type

type, public :: new_atompartition_type
   integer :: subsetsum
   integer :: subsetmax
   integer :: maxsize  ! Store the maximum size
   type(new_atomlist_type), pointer :: subsets(:) => null()
   type(new_atomlist_type), pointer :: subsets_allocation(:)
contains
   procedure :: new_subset
   procedure :: allocate_partition
!   procedure :: get_atomtypes
end type

type, public :: atomlist_type
   integer, allocatable :: atomidcs(:)
end type

type, public :: atompartition_type
   integer :: subsetsum
   integer :: subsetmax
   type(atomlist_type), allocatable :: subsets(:)
contains
   procedure :: get_atomtypes
end type

interface operator (==)
   module procedure is_equal
end interface

contains

subroutine allocate_partition(self, maxsize)
   class(new_atompartition_type), intent(inout) :: self
   integer, intent(in) :: maxsize

   self%maxsize = maxsize
   allocate(self%subsets_allocation(maxsize))
   self%subsets => self%subsets_allocation(1:0)

end subroutine

subroutine append(self, element)
   class(new_atomlist_type), intent(inout) :: self
   integer, intent(in) :: element
   integer :: new_size

   new_size = size(self%atomidcs) + 1
   self%atomidcs_allocation(new_size) = element
   self%atomidcs => self%atomidcs_allocation(1:new_size)

end subroutine

function new_subset(self) result(index)
   class(new_atompartition_type), intent(inout) :: self
   integer :: index

   index = size(self%subsets) + 1
   allocate(self%subsets_allocation(index)%atomidcs_allocation(self%maxsize))
   self%subsets_allocation(index)%atomidcs => self%subsets_allocation(index)%atomidcs_allocation(1:0)
   self%subsets => self%subsets_allocation(1:index)

end function

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

   atompartition%subsetsum = sum(n)
   atompartition%subsetmax = maxval(n)

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

   allocate (atomtypes(self%subsetsum))

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

end module
