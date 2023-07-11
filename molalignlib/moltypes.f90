module moltypes
use kinds
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
   procedure :: reorder => atoms_reorder
end type

contains

subroutine atoms_reorder(this, order)
   class(Molecule), intent(inout) :: this
   integer, intent(in) :: order(:)

   this%atoms = this%atoms(order(1:this%natom))
   this%adjmat = this%adjmat(order(1:this%natom), order(1:this%natom))

end subroutine

end module
