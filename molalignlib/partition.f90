module partition
use stdio

implicit none

type, public :: indexlist_type
   integer, allocatable :: indices(:)
end type

type, public :: pointer_to_part_type
   integer, allocatable :: signature(:)
   type(part_type), pointer :: ptr
end type

type, public :: part_type
   integer :: size
   integer :: index
   integer, pointer :: total_size
   integer, pointer :: maxpart_size
   integer, pointer :: idxmap(:)
   integer, pointer :: indices(:)
   integer, pointer :: allocation(:)
contains
   procedure :: append => part_append
end type

type, public :: partition_type
   integer :: size
   integer :: allocation_size
   integer, pointer :: total_size
   integer, pointer :: maxpart_size
   integer, pointer :: idxmap(:)
   type(part_type), pointer :: parts(:)
contains
   procedure :: print_partition
   procedure :: init => partition_init
   procedure :: new_part => partition_new_part
   procedure :: add_part => partition_add_part
end type

type, public :: subpartition_type
   integer :: size
   type(pointer_to_part_type), allocatable :: parts(:)
contains
   procedure :: init => subpartition_init
   procedure :: reset => subpartition_reset
   procedure :: add_part => subpartition_add_part
end type

interface operator (==)
   module procedure test_equality
end interface

contains

subroutine partition_init(self, allocation_size)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: allocation_size

   allocate (self%total_size)
   allocate (self%maxpart_size)
   allocate (self%parts(allocation_size))
   allocate (self%idxmap(allocation_size))

   self%size = 0
   self%total_size = 0
   self%maxpart_size = 0
   self%allocation_size = allocation_size

end subroutine

subroutine subpartition_init(self, allocation_size)
   class(subpartition_type), intent(inout) :: self
   integer, intent(in) :: allocation_size

   allocate (self%parts(allocation_size))

   self%size = 0

end subroutine

subroutine subpartition_reset(self)
   class(subpartition_type), intent(inout) :: self
   integer :: h

   do h = 1, self%size
      nullify (self%parts(h)%ptr)
   end do

   self%size = 0

end subroutine

subroutine part_append(self, element)
   class(part_type), intent(inout) :: self
   integer, intent(in) :: element

   self%size = self%size + 1
   self%indices(self%size) = element
   self%idxmap(element) = self%index
   self%indices => self%allocation(:self%size)

   self%total_size = self%total_size + 1
   if (self%size > self%maxpart_size) then
      self%maxpart_size = self%size
   end if

end subroutine

function partition_new_part(self) result(part)
   class(partition_type), intent(inout) :: self
   ! Result variable
   type(part_type), pointer :: part

   self%size = self%size + 1

   self%parts(self%size)%size = 0
   allocate (self%parts(self%size)%allocation(self%allocation_size))
   self%parts(self%size)%indices => self%parts(self%size)%allocation(:0)

   self%parts(self%size)%index = self%size
   self%parts(self%size)%idxmap => self%idxmap
   self%parts(self%size)%total_size => self%total_size
   self%parts(self%size)%maxpart_size => self%maxpart_size

   part => self%parts(self%size)

end function

subroutine partition_add_part(self, part)
   class(partition_type), intent(inout) :: self
   type(part_type), intent(in) :: part

   self%size = self%size + 1

   if (part%index /= self%size) then
      write (stderr, '(a,1x,i0,1x,a,1x,i0)') &
            'Error: part index', part%index, 'does not match partition size', self%size
      error stop
   end if

   self%parts(self%size)%size = 0
   allocate (self%parts(self%size)%allocation(self%allocation_size))
   self%parts(self%size)%indices => self%parts(self%size)%allocation(:0)

   self%parts(self%size)%index = part%index
   self%parts(self%size)%idxmap => self%idxmap
   self%parts(self%size)%total_size => self%total_size
   self%parts(self%size)%maxpart_size => self%maxpart_size

end subroutine

subroutine subpartition_add_part(self, part, signature)
   class(subpartition_type), intent(inout) :: self
   type(part_type), pointer, intent(in) :: part
   integer, intent(in) :: signature(:)

   self%size = self%size + 1
   self%parts(self%size)%ptr => part
   self%parts(self%size)%signature = signature

end subroutine

subroutine print_partition(self)
   class(partition_type), intent(in) :: self
   integer :: h

   write (stderr, *)
   do h = 1, self%size
      write (stderr, *) self%parts(h)%indices
   end do

end subroutine

function test_equality(self, other) result(equality)
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

end module
