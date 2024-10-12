module partition
use stdio

implicit none

type, public :: atomlist_type
   integer, allocatable :: atomidcs(:)
end type

type, public :: atompartitionlist_type
   integer, pointer :: atomidcs(:) => null()
   integer, pointer :: atomidcs_allocation(:)
end type

type, public :: atompartition_type
   integer :: max_size
   integer :: tot_size
   integer :: maxsize  ! Store the maximum size
   type(atompartitionlist_type), pointer :: subsets(:) => null()
   type(atompartitionlist_type), pointer :: subsets_allocation(:)
contains
   procedure :: append
   procedure :: new_subset
   procedure :: add_subset
   procedure :: get_atomtypes
   procedure :: allocate => allocate_partition
end type

interface operator (==)
   module procedure is_equal
end interface

contains

subroutine allocate_partition(self, maxsize)
   class(atompartition_type), intent(inout) :: self
   integer, intent(in) :: maxsize

   self%maxsize = maxsize
   allocate (self%subsets_allocation(maxsize))
   self%subsets => self%subsets_allocation(:0)

end subroutine

subroutine append(self, index, element)
   class(atompartition_type), intent(inout) :: self
   integer, intent(in) :: index, element
   integer :: new_size

   new_size = size(self%subsets(index)%atomidcs) + 1
   self%subsets(index)%atomidcs_allocation(new_size) = element
   self%subsets(index)%atomidcs => self%subsets(index)%atomidcs_allocation(:new_size)

   if (new_size > self%max_size) then
      self%max_size = new_size
   end if

   self%tot_size = self%tot_size + 1

end subroutine

function new_subset(self) result(index)
   class(atompartition_type), intent(inout) :: self
   integer :: index

   index = size(self%subsets) + 1
   allocate (self%subsets_allocation(index)%atomidcs_allocation(self%maxsize))
   self%subsets_allocation(index)%atomidcs => self%subsets_allocation(index)%atomidcs_allocation(:0)
   self%subsets => self%subsets_allocation(:index)

end function

subroutine add_subset(self, index)
   class(atompartition_type), intent(inout) :: self
   integer, intent(in) :: index

   if (index /= size(self%subsets) + 1) then
      write (stderr, '(a)') 'Error: index does not match subset size'
   end if

   allocate (self%subsets_allocation(index)%atomidcs_allocation(self%maxsize))
   self%subsets_allocation(index)%atomidcs => self%subsets_allocation(index)%atomidcs_allocation(:0)
   self%subsets => self%subsets_allocation(:index)

end subroutine

function get_atomtypes(self) result(atomtypes)
   class(atompartition_type), intent(in) :: self
   ! Result variable
   integer, allocatable :: atomtypes(:)
   ! Local variables
   integer :: h, i
   integer :: atomidx_i

   allocate (atomtypes(self%tot_size))

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
