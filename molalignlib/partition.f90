module partition
use stdio

implicit none
private
public atompartition

type, public :: atomlist_type
   integer, allocatable :: atomidcs(:)
end type

type, public :: atompartition_type
   integer :: natom
   integer :: largest
   type(atomlist_type), allocatable :: subsets(:)
contains
   procedure :: mapped
end type

contains

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

   atompartition%natom = sum(n)
   atompartition%largest = maxval(n)

   do h = 1, ntype
      allocate (atompartition%subsets(h)%atomidcs(n(h)))
   end do

   n(:) = 0
   do h = 1, size(atomtypes)
      n(atomtypes(h)) = n(atomtypes(h)) + 1
      atompartition%subsets(atomtypes(h))%atomidcs(n(atomtypes(h))) = h
   end do

end function

function mapped(self, mapping) result(atompartition)
   class(atompartition_type), intent(in) :: self
   integer, intent(in) :: mapping(:)
   ! Result variable
   type(atompartition_type) :: atompartition
   ! Local variables
   integer :: h

   atompartition = self
   do h = 1, size(self%subsets)
      atompartition%subsets(h)%atomidcs = mapping(self%subsets(h)%atomidcs)
   end do

end function

end module
