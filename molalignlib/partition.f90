module partition
use stdio

implicit none

type, public :: indexlist_type
   integer, allocatable :: indices(:)
end type

type, public :: pointertopart_type
   type(part_type), pointer :: ptr => null()
end type

type, public :: part_type
   integer :: size
   integer :: index
   integer, pointer :: total_size
   integer, pointer :: max_part_size
   integer, pointer :: index_part_map(:)
   integer, pointer :: indices(:)
   integer, pointer :: allocation(:)
contains
   procedure :: add => part_add
end type

type, public :: partition_type
   integer :: size
   integer :: allocation_size
   integer, pointer :: total_size
   integer, pointer :: max_part_size
   integer, pointer :: index_part_map(:)
   type(part_type), pointer :: parts(:)
contains
   procedure :: init => partition_init
   procedure :: get_new_part => partition_get_new_part
   procedure :: print_parts => partition_print_parts
end type

interface operator (==)
   module procedure equality
end interface

contains

function equality(self, other)
   type(partition_type), intent(in) :: self, other
   ! Result variable
   logical :: equality
   ! Local variables
   integer :: h

   if (self%size /= other%size) then
      equality = .false.
      return
   end if

   do h = 1, self%size
      if (self%parts(h)%size /= other%parts(h)%size) then
         equality = .false.
         return
      end if
   end do

   equality = .true.

end function

subroutine partition_init(self, allocation_size)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: allocation_size

   allocate (self%total_size)
   allocate (self%max_part_size)
   allocate (self%parts(allocation_size))
   allocate (self%index_part_map(allocation_size))

   self%size = 0
   self%total_size = 0
   self%max_part_size = 0
   self%allocation_size = allocation_size

end subroutine

subroutine part_add(self, element)
   class(part_type), intent(inout) :: self
   integer, intent(in) :: element

   self%size = self%size + 1
   self%indices => self%allocation(:self%size)
   self%indices(self%size) = element
   self%index_part_map(element) = self%index

   self%total_size = self%total_size + 1
   if (self%size > self%max_part_size) then
      self%max_part_size = self%size
   end if

end subroutine

function partition_get_new_part(self) result(part)
   class(partition_type), intent(inout) :: self
   ! Result variable
   type(part_type), pointer :: part

   self%size = self%size + 1

   self%parts(self%size)%size = 0
   allocate (self%parts(self%size)%allocation(self%allocation_size))
   self%parts(self%size)%indices => self%parts(self%size)%allocation(:0)

   self%parts(self%size)%index = self%size
   self%parts(self%size)%index_part_map => self%index_part_map
   self%parts(self%size)%total_size => self%total_size
   self%parts(self%size)%max_part_size => self%max_part_size

   part => self%parts(self%size)

end function

subroutine partition_print_parts(self)
   class(partition_type), intent(in) :: self
   integer :: h

   write (stderr, *)
   do h = 1, self%size
      write (stderr, *) h, self%parts(h)%indices
   end do

end subroutine

end module
