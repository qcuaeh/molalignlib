module partition
use stdio

implicit none

type, public :: countlist_type
   integer, allocatable :: counts(:)
end type

type, public :: indexlist_type
   integer, allocatable :: indices(:)
end type

type, public :: partitionlist_type
   integer, pointer :: indices(:) => null()
   integer, pointer :: indices_allocation(:)
contains
   procedure :: append
   procedure :: allocate => allocate_part
end type

type, public :: partition_type
   integer :: max_size
   integer :: tot_size
   integer :: allocation_size
   type(partitionlist_type), pointer :: parts(:) => null()
   type(partitionlist_type), pointer :: parts_allocation(:)
contains
   procedure :: new_part
   procedure :: add_part
   procedure :: get_item_types
   procedure :: append_to_part
   procedure :: allocate => allocate_partition
end type

interface operator (==)
   module procedure is_equal
end interface

contains

subroutine allocate_part(self, allocation_size)
   class(partitionlist_type), intent(inout) :: self
   integer, intent(in) :: allocation_size

   allocate (self%indices_allocation(allocation_size))
   self%indices => self%indices_allocation(:0)

end subroutine

subroutine allocate_partition(self, allocation_size)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: allocation_size

   self%allocation_size = allocation_size
   allocate (self%parts_allocation(allocation_size))
   self%parts => self%parts_allocation(:0)

end subroutine

subroutine append(self, element)
   class(partitionlist_type), intent(inout) :: self
   integer, intent(in) :: element
   integer :: new_size

   new_size = size(self%indices) + 1
   self%indices_allocation(new_size) = element
   self%indices => self%indices_allocation(:new_size)

end subroutine

subroutine append_to_part(self, index, element)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: index, element
   integer :: new_size

   call self%parts(index)%append(element)

   self%tot_size = self%tot_size + 1

   new_size = size(self%parts(index)%indices)
   if (new_size > self%max_size) then
      self%max_size = new_size
   end if

end subroutine

function new_part(self) result(index)
   class(partition_type), intent(inout) :: self
   integer :: index

   index = size(self%parts) + 1

   call self%parts_allocation(index)%allocate(self%allocation_size)
   self%parts => self%parts_allocation(:index)

end function

subroutine add_part(self, index)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: index

   if (index /= size(self%parts) + 1) then
      write (stderr, '(a)') 'Error: index does not match part size'
   end if

   call self%parts_allocation(index)%allocate(self%allocation_size)
   self%parts => self%parts_allocation(:index)

end subroutine

function get_item_types(self) result(types)
   class(partition_type), intent(in) :: self
   ! Result variable
   integer, allocatable :: types(:)
   ! Local variables
   integer :: h, i
   integer :: index

   allocate (types(self%tot_size))

   do h = 1, size(self%parts)
      do i = 1, size(self%parts(h)%indices)
         index = self%parts(h)%indices(i)
         types(index) = h
      end do
   end do

end function

function is_equal(self, other)
   type(partition_type), intent(in) :: self, other
   ! Result variable
   logical :: is_equal
   ! Local variables
   integer :: h

   if (size(self%parts) /= size(other%parts)) then
      is_equal = .false.
      return
   end if

   do h = 1, size(self%parts)
      if (size(self%parts(h)%indices) /= size(other%parts(h)%indices)) then
         is_equal = .false.
         return
      end if
   end do

   is_equal = .true.

end function

end module
