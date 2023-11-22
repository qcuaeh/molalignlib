module types
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
   procedure :: get_znums => get_znums
   procedure :: get_ztypes => get_ztypes
   procedure :: get_weights => get_weights
   procedure :: get_coords => get_coords
   procedure :: get_labels => get_labels
   procedure :: set_coords => set_coords
   procedure :: rotate_coords => rotate_coords
   procedure :: translate_coords => translate_coords
   procedure :: permutate_atoms => permutate_atoms
end type

contains

subroutine rotate_coords(self, rotmat)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_coords(rotated(self%natom, self%get_coords(), rotmat))

end subroutine

subroutine translate_coords(self, travec)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: travec(3)

   call self%set_coords(translated(self%natom, self%get_coords(), travec))

end subroutine

subroutine permutate_atoms(self, order)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: order(:)

   self%atoms = self%atoms(order(1:self%natom))
   self%adjmat = self%adjmat(order(1:self%natom), order(1:self%natom))

end subroutine

function get_znums(self) result(znums)
   class(Molecule), intent(in) :: self
   integer, dimension(self%natom) :: znums
   integer i

   do i = 1, self%natom
      znums(i) = self%atoms(i)%znum
   end do

end function

function get_ztypes(self) result(ztypes)
   class(Molecule), intent(in) :: self
   integer, dimension(self%natom) :: ztypes
   integer i

   do i = 1, self%natom
      ztypes(i) = self%atoms(i)%type
   end do

end function

function get_weights(self) result(weights)
   class(Molecule), intent(in) :: self
   real(wp), dimension(self%natom) :: weights
   integer i

   do i = 1, self%natom
      weights(i) = self%atoms(i)%weight
   end do

end function

function get_coords(self) result(coords)
   class(Molecule), intent(in) :: self
   real(wp), dimension(3, self%natom) :: coords
   integer i

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords
   end do

end function

subroutine set_coords(self, coords)
   class(Molecule), intent(inout) :: self
   real(wp), dimension(3, self%natom), intent(in) :: coords
   integer i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine

function get_labels(self) result(labels)
   class(Molecule), intent(in) :: self
   character(wl), dimension(self%natom) :: labels
   integer i

   do i = 1, self%natom
      labels(i) = self%atoms(i)%label
   end do

end function

end module
