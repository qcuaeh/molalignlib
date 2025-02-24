module metapartition
use parameters
use partition
implicit none

public metapart_type
public metapartition_type
public metapartpointer_type

type :: metapart_type
   integer :: index
   integer :: part_size
   type(partition_type) :: subpartition
   integer, pointer :: largest_part_size
   integer, pointer :: allocation(:)
   integer, pointer :: items(:)
contains
   procedure :: add => part_add
end type

type :: metapartition_type
   integer :: num_parts
   integer :: max_partition_size
   integer, pointer :: largest_part_size
   logical :: initialized = .false.
   type(metapart_type), pointer :: parts(:) => null()
contains
   final :: partition_finalize
   procedure :: initialize => partition_initialize
   procedure :: new_part => partition_new_part
   procedure :: add_part => partition_add_part
   procedure :: print_parts => partition_print_parts
end type

type :: metapartpointer_type
   type(metapart_type), pointer :: ptr => null()
end type

contains

subroutine partition_initialize(self, max_partition_size)
   class(metapartition_type), intent(inout) :: self
   integer, intent(in) :: max_partition_size

   if (associated(self%parts)) then
      error stop 'Memory leak'
   end if

   allocate (self%largest_part_size)
   allocate (self%parts(max_partition_size))

   self%num_parts = 0
   self%largest_part_size = 0
   self%max_partition_size = max_partition_size
   self%initialized = .true.

end subroutine

subroutine partition_finalize(self)
   type(metapartition_type), intent(inout) :: self
   integer :: h

   if (.not. self%initialized) then
      return
   end if

   do h = 1, self%num_parts
      deallocate (self%parts(h)%allocation)
   end do

   deallocate (self%parts)

end subroutine

function partition_new_part(self, max_size) result(part)
   class(metapartition_type), intent(inout) :: self
   integer, intent(in) :: max_size
   ! Result variable
   type(metapart_type), pointer :: part

   ! Increase partition size
   self%num_parts = self%num_parts + 1

   ! Initialize part size
   self%parts(self%num_parts)%part_size = 0

   ! Set index to current partition size
   self%parts(self%num_parts)%index = self%num_parts

   ! Allocate list allocation
   allocate (self%parts(self%num_parts)%allocation(max_size))

   ! Point list pointer to list allocation with null size
   self%parts(self%num_parts)%items => self%parts(self%num_parts)%allocation(:0)

   ! Point largest part size pointer to partition largest part size
   self%parts(self%num_parts)%largest_part_size => self%largest_part_size

   ! Return pointer to part
   part => self%parts(self%num_parts)

end function

subroutine partition_add_part(self, part)
   class(metapartition_type), intent(inout) :: self
   type(metapart_type), intent(in) :: part
   type(metapart_type), pointer :: newpart
   integer :: i

   newpart => self%new_part(part%part_size)

   do i = 1, part%part_size
      call newpart%add(part%items(i))
   end do

end subroutine

subroutine part_add(self, element)
   class(metapart_type), intent(inout) :: self
   integer, intent(in) :: element

   ! Increase part size
   self%part_size = self%part_size + 1

   ! Update list pointer
   self%items => self%allocation(:self%part_size)

   ! Add element to part
   self%items(self%part_size) = element

   ! Update largest part size
   if (self%part_size > self%largest_part_size) then
      self%largest_part_size = self%part_size
   end if

end subroutine

subroutine partition_print_parts(self)
   class(metapartition_type), intent(in) :: self
   integer :: h
   character(:), allocatable :: fmtstr

   write (stderr, *)
   do h = 1, self%num_parts
      fmtstr = "('{'" // repeat(',1x,i3', self%parts(h)%part_size) // ",' }')"
      write (stderr, fmtstr) self%parts(h)%items(:self%parts(h)%part_size)
   end do

end subroutine

end module
