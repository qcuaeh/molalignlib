module partition

implicit none
private
public make_partition

type, public :: atomlist_type
   integer, allocatable :: atomidcs(:)
end type

type, public :: atompartition_type
   integer :: natom
   type(atomlist_type), allocatable :: subsets(:)
end type

type, public :: partition_type
   type(atomlist_type), allocatable :: parts(:)
   integer, allocatable :: atom_order(:)
   integer, allocatable :: atom_mapping(:)
contains
   procedure :: init => partition_init
   procedure :: get_lenlist => partition_lenlist
end type

contains

function make_partition(typepops, atomtypes) result(partition)
   integer, intent(in) :: typepops(:)
   integer, intent(in) :: atomtypes(:)
   ! Result variable
   type(atompartition_type) :: partition
   ! Local variables
   integer :: i
   integer, allocatable :: n(:)

   allocate (n(size(typepops)))
   allocate (partition%subsets(size(typepops)))

   do i = 1, size(typepops)
      allocate (partition%subsets(i)%atomidcs(typepops(i)))
   end do

   n(:) = 0
   do i = 1, size(atomtypes)
      n(atomtypes(i)) = n(atomtypes(i)) + 1
      partition%subsets(atomtypes(i))%atomidcs(n(atomtypes(i))) = i
   end do

   partition%natom = size(atomtypes)

end function

subroutine partition_init(self, npart, partidcs)
   class(partition_type), intent(out) :: self
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
   class(partition_type), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: lenlist(:)

   allocate (lenlist(size(self%parts)))

   do i = 1, size(self%parts)
      lenlist(i) = size(self%parts(i)%atomidcs)
   end do

end function

end module
