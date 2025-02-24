module partition
use parameters
implicit none
private

public part_type
public partition_type
public partpointer_type
public operator (==)
public operator (<)
public assignment (=)

type :: part_type
   integer :: index
   integer :: part_size
   integer, pointer :: idcs(:)
   integer, pointer :: largest_part_size
   integer, pointer :: allocation(:)
   integer, pointer :: items(:)
contains
   procedure :: add => part_add
end type

type :: partition_type
   integer :: num_parts
   integer :: num_items
   integer, pointer :: idcs(:)
   integer, pointer :: largest_part_size
   logical :: initialized = .false.
   type(part_type), pointer :: parts(:) => null()
contains
   final :: partition_finalize
   procedure :: initialize => partition_initialize
   procedure :: new_part => partition_new_part
   procedure :: add_part => partition_add_part
   procedure :: print_parts => partition_print_parts
end type

type :: partpointer_type
   type(part_type), pointer :: ptr => null()
end type

interface operator (==)
   module procedure partition_equality
end interface

interface operator (<)
   module procedure partition_precedence
end interface

interface assignment (=)
   module procedure partition_assignment
end interface

contains

subroutine partition_assignment(left, right)
   class(partition_type), intent(out) :: left
   type(partition_type), intent(in) :: right
   ! Local variables
   integer :: h

   call left%initialize(right%num_items)

   do h = 1, right%num_parts
      call left%add_part(right%parts(h))
   end do

end subroutine

subroutine partition_initialize(self, num_items)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: num_items

   if (associated(self%parts)) then
      error stop 'Memory leak'
   end if

   allocate (self%largest_part_size)
   allocate (self%idcs(num_items))
   allocate (self%parts(num_items))

   self%num_parts = 0
   self%largest_part_size = 0
   self%num_items = num_items
   self%initialized = .true.

end subroutine

subroutine partition_finalize(self)
   type(partition_type), intent(inout) :: self
   integer :: h

   if (.not. self%initialized) then
      return
   end if

   do h = 1, self%num_parts
      deallocate (self%parts(h)%allocation)
   end do

   deallocate (self%parts)
   deallocate (self%idcs)
   deallocate (self%largest_part_size)

end subroutine

function partition_new_part(self, max_size) result(part)
   class(partition_type), intent(inout) :: self
   integer, intent(in) :: max_size
   ! Result variable
   type(part_type), pointer :: part

   ! Increase partition size
   self%num_parts = self%num_parts + 1

   ! Set part index to current part num
   self%parts(self%num_parts)%index = self%num_parts

   ! Initialize part size
   self%parts(self%num_parts)%part_size = 0

   ! Allocate list allocation
   allocate (self%parts(self%num_parts)%allocation(max_size))

   ! Point list pointer to list allocation with null size
   self%parts(self%num_parts)%items => self%parts(self%num_parts)%allocation(:0)

   ! Point largest part size pointer to partition largest part size
   self%parts(self%num_parts)%largest_part_size => self%largest_part_size

   ! Point map pointer to partition map
   self%parts(self%num_parts)%idcs => self%idcs

   ! Return pointer to new part
   part => self%parts(self%num_parts)

end function

subroutine partition_add_part(self, part)
   class(partition_type), intent(inout) :: self
   type(part_type), intent(in) :: part
   type(part_type), pointer :: newpart
   integer :: i

   newpart => self%new_part(part%part_size)

   do i = 1, part%part_size
      call newpart%add(part%items(i))
   end do

end subroutine

subroutine part_add(self, element)
   class(part_type), intent(inout) :: self
   integer, intent(in) :: element

   ! Increase part size
   self%part_size = self%part_size + 1

   ! Update list pointer
   self%items => self%allocation(:self%part_size)

   ! Add element to part
   self%items(self%part_size) = element

   ! Add index to part map
   self%idcs(element) = self%index

   ! Update largest part size
   if (self%part_size > self%largest_part_size) then
      self%largest_part_size = self%part_size
   end if

end subroutine

subroutine partition_print_parts(self)
   class(partition_type), intent(in) :: self
   integer :: h
   character(:), allocatable :: fmtstr

   write (stderr, *)
   do h = 1, self%num_parts
      fmtstr = "(i3,':',2x,'{'" // repeat(',1x,i3', self%parts(h)%part_size) // ",1x,'}')"
      write (stderr, fmtstr) h, self%parts(h)%items(:self%parts(h)%part_size)
   end do

end subroutine

function partition_equality(left, right) result(equality)
   type(partition_type), intent(in) :: left, right
   ! Result variable
   logical :: equality
   ! Local variables
   integer :: h, i

   if (left%num_parts /= right%num_parts) then
      equality = .false.
      return
   end if

   do h = 1, left%num_parts

      if (left%parts(h)%part_size /= right%parts(h)%part_size) then
         equality = .false.
         return
      end if

      do i = 1, left%parts(h)%part_size
         if (left%parts(h)%items(i) /= right%parts(h)%items(i)) then
            equality = .false.
            return
         end if
      end do

   end do

   equality = .true.

end function

function partition_precedence(left, right) result(priority)
   type(partition_type), intent(in) :: left, right
   ! Result variable
   logical :: priority
   ! Local variables
   integer :: h, k, i, offset

   if (left%num_parts >= right%num_parts) then
!      write (stderr, *) 'left%num_parts >= right%num_parts'
      priority = .false.
      return
   end if

   k = 1
   do h = 1, left%num_parts

      if (left%parts(h)%part_size >= right%parts(k)%part_size) then

         offset = 0
         do while (offset < left%parts(h)%part_size)
            do i = 1, right%parts(k)%part_size
               if (left%parts(h)%items(offset+i) /= right%parts(k)%items(i)) then
!                  write (stderr, *) 'left%parts(h)%items(offset+i) /= right%parts(k)%items(i)'
!                  write (stderr, *) left%parts(h)%items(:left%parts(h)%part_size)
!                  write (stderr, *) right%parts(k)%items(:right%parts(k)%part_size)
                  priority = .false.
                  return
               end if
            end do
            offset = offset + right%parts(k)%part_size
            k = k + 1
         end do

         if (offset /= left%parts(h)%part_size) then
            error stop 'offset /= left%parts(h)%part_size'
         end if

      else

!         write (stderr, *) 'left%parts(h)%part_size < right%parts(k)%part_size'
         priority = .false.
         return

      end if

   end do

   priority = .true.

end function

end module
