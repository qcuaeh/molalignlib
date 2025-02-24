module bipartition
use parameters
use partition
implicit none
private

public bipart_type
public bipartition_type
public bipartpointer_type
public operator (==)
public assignment (=)

type :: bipart_type
   integer :: index
   integer :: part_size1
   integer :: part_size2
   integer, pointer :: idcs1(:)
   integer, pointer :: idcs2(:)
   integer, pointer :: largest_part_size
   integer, pointer :: allocation1(:)
   integer, pointer :: allocation2(:)
   integer, pointer :: items1(:)
   integer, pointer :: items2(:)
contains
   procedure :: add1 => part_add1
   procedure :: add2 => part_add2
end type

type :: bipartition_type
   integer :: num_parts
   integer :: num_items1
   integer :: num_items2
   integer, pointer :: idcs1(:)
   integer, pointer :: idcs2(:)
   integer, pointer :: largest_part_size
   logical :: initialized = .false.
   type(bipart_type), pointer :: parts(:) => null()
contains
   final :: partition_finalize
   procedure :: initialize => partition_initialize
   procedure :: new_part => partition_new_part
   procedure :: add_part => partition_add_part
   procedure :: partition1 => partition_first_partition
   procedure :: partition2 => partition_second_partition
   procedure :: print_parts => partition_print_parts
end type

type :: bipartpointer_type
   type(bipart_type), pointer :: ptr => null()
end type

interface operator (==)
   module procedure bipartition_equality
end interface

interface assignment (=)
   module procedure partition_assignment
end interface

contains

subroutine partition_assignment(left, right)
   class(bipartition_type), intent(out) :: left
   type(bipartition_type), intent(in) :: right
   ! Local variables
   integer :: h

   call left%initialize(right%num_items1, right%num_items2)

   do h = 1, right%num_parts
      call left%add_part(right%parts(h))
   end do

end subroutine

subroutine partition_initialize(self, num_items1, num_items2)
   class(bipartition_type), intent(inout) :: self
   integer, intent(in) :: num_items1, num_items2

   if (associated(self%parts)) then
      error stop 'Memory leak'
   end if

   allocate (self%largest_part_size)
   allocate (self%idcs1(num_items1))
   allocate (self%idcs2(num_items2))
   allocate (self%parts(num_items1 + num_items2))

   self%num_items1 = num_items1
   self%num_items2 = num_items2
   self%num_parts = 0
   self%largest_part_size = 0
   self%initialized = .true.

end subroutine

subroutine partition_finalize(self)
   type(bipartition_type), intent(inout) :: self
   integer :: h

   if (.not. self%initialized) then
      return
   end if

   do h = 1, self%num_parts
      deallocate (self%parts(h)%allocation1)
      deallocate (self%parts(h)%allocation2)
   end do

   deallocate (self%parts)
   deallocate (self%idcs1)
   deallocate (self%idcs2)
   deallocate (self%largest_part_size)

end subroutine

function partition_new_part(self, max_size1, max_size2) result(part)
   class(bipartition_type), intent(inout) :: self
   integer, intent(in) :: max_size1, max_size2
   ! Result variable
   type(bipart_type), pointer :: part

   ! Increase partition size
   self%num_parts = self%num_parts + 1

   ! Set part index to current part num
   self%parts(self%num_parts)%index = self%num_parts

   ! Initialize part size
   self%parts(self%num_parts)%part_size1 = 0
   self%parts(self%num_parts)%part_size2 = 0

   ! Allocate list allocation
   allocate (self%parts(self%num_parts)%allocation1(max_size1))
   allocate (self%parts(self%num_parts)%allocation2(max_size2))

   ! Point list pointer to list allocation with null size
   self%parts(self%num_parts)%items1 => self%parts(self%num_parts)%allocation1(:0)
   self%parts(self%num_parts)%items2 => self%parts(self%num_parts)%allocation2(:0)

   ! Point part max part size pointer to partition max part size
   self%parts(self%num_parts)%largest_part_size => self%largest_part_size

   ! Point part map pointer to partition maps
   self%parts(self%num_parts)%idcs1 => self%idcs1
   self%parts(self%num_parts)%idcs2 => self%idcs2

   ! Return pointer to part
   part => self%parts(self%num_parts)

end function

subroutine partition_add_part(self, part)
   class(bipartition_type), intent(inout) :: self
   type(bipart_type), intent(in) :: part
   type(bipart_type), pointer :: newpart
   integer :: i

   newpart => self%new_part(part%part_size1, part%part_size2)

   do i = 1, part%part_size1
      call newpart%add1(part%items1(i))
   end do

   do i = 1, part%part_size2
      call newpart%add2(part%items2(i))
   end do

end subroutine

subroutine part_add1(self, element)
   class(bipart_type), intent(inout) :: self
   integer, intent(in) :: element

   ! Increase part size
   self%part_size1 = self%part_size1 + 1

   ! Update list pointers
   self%items1 => self%allocation1(:self%part_size1)

   ! Add element to part
   self%items1(self%part_size1) = element

   ! Add index to part map
   self%idcs1(element) = self%index

   ! Update largest part size
   if (self%part_size1 > self%largest_part_size) then
      self%largest_part_size = self%part_size1
   end if

end subroutine

subroutine part_add2(self, element)
   class(bipart_type), intent(inout) :: self
   integer, intent(in) :: element

   ! Increase part size
   self%part_size2 = self%part_size2 + 1

   ! Update list pointers
   self%items2 => self%allocation2(:self%part_size2)

   ! Add element to part
   self%items2(self%part_size2) = element

   ! Add index to part map
   self%idcs2(element) = self%index

   ! Update largest part size
   if (self%part_size2 > self%largest_part_size) then
      self%largest_part_size = self%part_size2
   end if

end subroutine

subroutine partition_print_parts(self)
   class(bipartition_type), intent(in) :: self
   integer :: h
   character(:), allocatable :: fmtstr

   write (stderr, *)
   do h = 1, self%num_parts
      fmtstr = "(i3,':',2x,'{'" // repeat(',1x,i3', self%parts(h)%part_size1) &
            // ",1x,'}',2x,'{'" // repeat(',1x,i3', self%parts(h)%part_size2) // ",1x,'}')"
      write (stderr, fmtstr) h, self%parts(h)%items1(:self%parts(h)%part_size1), &
            self%parts(h)%items2(:self%parts(h)%part_size2)
   end do

end subroutine

function partition_first_partition(self) result(partition)
   class(bipartition_type), intent(in) :: self
   integer :: h, i
   type(partition_type) :: partition
   type(part_type), pointer :: newtype

   call partition%initialize(self%num_items1)

   do h = 1, self%num_parts
      if (self%parts(h)%part_size1 > 0) then
         newtype => partition%new_part(self%parts(h)%part_size1)
         do i = 1, self%parts(h)%part_size1
            call newtype%add(self%parts(h)%items1(i))
         end do
      end if
   end do

end function

function partition_second_partition(self) result(partition)
   class(bipartition_type), intent(in) :: self
   integer :: h, i
   type(partition_type) :: partition
   type(part_type), pointer :: newtype

   call partition%initialize(self%num_items2)

   do h = 1, self%num_parts
      if (self%parts(h)%part_size2 > 0) then
         newtype => partition%new_part(self%parts(h)%part_size2)
         do i = 1, self%parts(h)%part_size2
            call newtype%add(self%parts(h)%items2(i))
         end do
      end if
   end do

end function

function bipartition_equality(self, other) result(equality)
   type(bipartition_type), intent(in) :: self, other
   ! Result variable
   logical :: equality
   ! Local variables
   integer :: h

   if (self%num_parts /= other%num_parts) then
      equality = .false.
      return
   end if

   do h = 1, self%num_parts
      if (self%parts(h)%part_size1 /= other%parts(h)%part_size1) then
         equality = .false.
         return
      end if
      if (self%parts(h)%part_size2 /= other%parts(h)%part_size2) then
         equality = .false.
         return
      end if
   end do

   equality = .true.

end function

end module
