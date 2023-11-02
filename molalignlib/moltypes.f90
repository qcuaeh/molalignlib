module moltypes
use kinds
use discrete
use rotation
use translation
use adjacency
use alignment
implicit none
private

type :: Atom
   character(:), allocatable :: label
   integer :: znum
   integer :: type
   real(wp) :: weight
   real(wp) :: coords(3)
end type

type :: Bond
   integer :: atom1
   integer :: atom2
   character(2) :: order
end type

type, public :: Block
   integer :: nblk
   integer, allocatable :: nblkidx
   integer, allocatable :: nblklen
   integer, allocatable :: nblkoff
end type

type, public :: Molecule
   integer :: natom
   integer :: nbond
   character(:), allocatable :: title
   type(Atom), allocatable :: atoms(:)
   type(Bond), allocatable :: bonds(:)
   logical, allocatable :: adjmat(:, :)
contains
   procedure :: get_znums => atoms_get_znums
   procedure :: get_types => atoms_get_types
   procedure :: get_weights => atoms_get_weights
   procedure :: get_coords => atoms_get_coords
   procedure :: get_labels => atoms_get_labels
   procedure :: set_coords => atoms_set_coords
   procedure :: rotate => atoms_rotate
   procedure :: translate => atoms_translate
   procedure :: permutate => atoms_permutate
end type

contains

subroutine atoms_permutate(self, order)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: order(:)

   self%atoms = self%atoms(order(1:self%natom))
   self%adjmat = self%adjmat(order(1:self%natom), order(1:self%natom))

end subroutine

subroutine atoms_rotate(self, rotmat)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_coords(rotated(self%natom, self%get_coords(), rotmat))

end subroutine

subroutine atoms_translate(self, travec)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: travec(3)

   call self%set_coords(translated(self%natom, self%get_coords(), travec))

end subroutine

function atoms_get_znums(self) result(znums)
   class(Molecule), intent(in) :: self
   integer, dimension(self%natom) :: znums
   integer i

   do i = 1, self%natom
      znums(i) = self%atoms(i)%znum
   end do

end function

function atoms_get_types(self) result(types)
   class(Molecule), intent(in) :: self
   integer, dimension(self%natom) :: types
   integer i

   do i = 1, self%natom
      types(i) = self%atoms(i)%type
   end do

end function

function atoms_get_weights(self) result(weights)
   class(Molecule), intent(in) :: self
   real(wp), dimension(self%natom) :: weights
   integer i

   do i = 1, self%natom
      weights(i) = self%atoms(i)%weight
   end do

end function

function atoms_get_coords(self) result(coords)
   class(Molecule), intent(in) :: self
   real(wp), dimension(3, self%natom) :: coords
   integer i

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords
   end do

end function

subroutine atoms_set_coords(self, coords)
   class(Molecule), intent(inout) :: self
   real(wp), dimension(3, self%natom), intent(in) :: coords
   integer i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine

function atoms_get_labels(self) result(labels)
   class(Molecule), intent(in) :: self
   character(wl), dimension(self%natom) :: labels
   integer i

   do i = 1, self%natom
      labels(i) = self%atoms(i)%label
   end do

end function

end module
