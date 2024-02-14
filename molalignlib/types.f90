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
   integer :: ztype
   real(wp) :: weight
   real(wp) :: coords(3)
   integer :: nadj
   integer, allocatable :: adjidx(:)
end type

type :: Bond
   integer :: atom1
   integer :: atom2
   character(2) :: order
end type

type, public :: Equiv
   integer :: eqvlen
   integer, allocatable :: eqvidx(:)
end type

type, public :: Block
   integer :: blklen
   integer, allocatable :: blkidx(:)
end type

type, public :: Molecule
   integer :: natom
   integer :: nbond
   integer :: nblock
   integer :: nequiv
   character(:), allocatable :: title
   type(Atom), allocatable :: atoms(:)
   type(Bond), allocatable :: bonds(:)
   type(Block), allocatable :: blocks(:)
   type(Equiv), allocatable :: equivs(:)
   logical, allocatable :: adjmat(:, :)
contains
   procedure :: get_znums
   procedure :: get_ztypes
   procedure :: get_weights
   procedure :: get_coords
   procedure :: get_labels
   procedure :: set_coords
   procedure :: mirror_coords
   procedure, private :: matrix_rotate_coords
   procedure, private :: quater_rotate_coords
   generic, public :: rotate_coords => matrix_rotate_coords, quater_rotate_coords
   procedure :: translate_coords
   procedure :: permutate_atoms
end type

contains

subroutine mirror_coords(self)

   class(Molecule), intent(inout) :: self
   real(wp) :: coords(3, self%natom)

   coords = self%get_coords()
   coords(1, :) = -coords(1, :)
   call self%set_coords(coords)

end subroutine mirror_coords

subroutine matrix_rotate_coords(self, rotmat)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_coords(matrix_rotated(self%natom, self%get_coords(), rotmat))

end subroutine matrix_rotate_coords

subroutine quater_rotate_coords(self, rotquat)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: rotquat(4)

   call self%set_coords(quater_rotated(self%natom, self%get_coords(), rotquat))

end subroutine quater_rotate_coords

subroutine translate_coords(self, travec)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: travec(3)

   call self%set_coords(translated(self%natom, self%get_coords(), travec))

end subroutine translate_coords

subroutine permutate_atoms(self, order)

   class(Molecule), intent(inout) :: self
   integer, intent(in) :: order(:)

   self%atoms = self%atoms(order(1:self%natom))
   self%adjmat = self%adjmat(order(1:self%natom), order(1:self%natom))

end subroutine permutate_atoms

function get_znums(self) result(znums)

   class(Molecule), intent(in) :: self
   integer :: znums(self%natom)
   integer i

   do i = 1, self%natom
      znums(i) = self%atoms(i)%znum
   end do

end function get_znums

function get_ztypes(self) result(ztypes)

   class(Molecule), intent(in) :: self
   integer :: ztypes(self%natom)
   integer i

   do i = 1, self%natom
      ztypes(i) = self%atoms(i)%ztype
   end do

end function get_ztypes

function get_weights(self) result(weights)

   class(Molecule), intent(in) :: self
   real(wp) :: weights(self%natom)
   integer i

   do i = 1, self%natom
      weights(i) = self%atoms(i)%weight
   end do

end function get_weights

function get_coords(self) result(coords)

   class(Molecule), intent(in) :: self
   real(wp) :: coords(3, self%natom)
   integer i

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords
   end do

end function get_coords

subroutine set_coords(self, coords)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: coords(3, self%natom)
   integer i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine set_coords

function get_labels(self) result(labels)

   class(Molecule), intent(in) :: self
   character(wl) :: labels(self%natom)
   integer i

   do i = 1, self%natom
      labels(i) = self%atoms(i)%label
   end do

end function get_labels

function get_adjmat(self) result(adjmat)

   class(Molecule), intent(in) :: self
   type(Atom) :: iatom
   logical :: adjmat(self%natom, self%natom)
   integer i, k

   adjmat(:, :) = .false.

   do i = 1, self%natom
      iatom = self%atoms(i)
      do k = 1, iatom%nadj
         adjmat(i, iatom%adjidx(k)) = .true.
      end do
   end do

end function get_adjmat

end module
