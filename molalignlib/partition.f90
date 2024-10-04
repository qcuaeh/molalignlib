module partition

implicit none
private

type, public :: t_atomlist
   integer, allocatable :: atomidcs(:)
end type

type, public :: t_partition
   type(t_atomlist), allocatable :: parts(:)
   integer, allocatable :: atom_order(:)
   integer, allocatable :: atom_mapping(:)
contains
   procedure :: init => partition_init
   procedure :: get_lenlist => partition_lenlist
end type

contains

subroutine partition_init(self, npart, partidcs)
   class(t_partition), intent(out) :: self
   integer, intent(in) :: npart
   integer, intent(in) :: partidcs(:)
   ! Local variables
   integer :: h, i, offset
   integer, allocatable :: n(:)

   allocate (n(npart))
   allocate (self%parts(npart))
   allocate (self%atom_order(size(partidcs)))
   allocate (self%atom_mapping(size(partidcs)))

   n(:) = 0
   do i = 1, size(partidcs)
      n(partidcs(i)) = n(partidcs(i)) + 1
   end do

   do i = 1, npart
      allocate (self%parts(i)%atomidcs(n(i)))
   end do

   n(:) = 0
   do i = 1, size(partidcs)
      n(partidcs(i)) = n(partidcs(i)) + 1
      self%parts(partidcs(i))%atomidcs(n(partidcs(i))) = i
   end do

   offset = 0
   do h = 1, size(self%parts)
      do i = 1, size(self%parts(h)%atomidcs)
         self%atom_order(offset + i) = self%parts(h)%atomidcs(i)
         self%atom_mapping(self%parts(h)%atomidcs(i)) = offset + i
      end do
      offset = offset + size(self%parts(h)%atomidcs)
   end do

end subroutine

function partition_lenlist(self) result(lenlist)
   class(t_partition), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: lenlist(:)

   allocate (lenlist(size(self%parts)))

   do i = 1, size(self%parts)
      lenlist(i) = size(self%parts(i)%atomidcs)
   end do

end function

end module
