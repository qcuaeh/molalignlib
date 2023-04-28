module molecule
implicit none
private

type, public :: Block
   integer :: nblk
   integer, allocatable :: nblkidx
   integer, allocatable :: nblklen
   integer, allocatable :: nblkoff
end type

type :: Atom
   character(:), allocatable :: label
   integer, allocatable :: znum
   integer, allocatable :: type
   real(wp), allocatable :: weight
   real(wp), allocatable :: coords
end type

type, public :: Molecule
   integer :: natom
   character(:), allocatable :: title
   type(Atom), allocatable :: atoms(:)
   logical, allocatable :: adjmat(:)
contains
   procedure :: reorder => atoms_reorder
end type

contains

subroutine atoms_reorder(this, order)
   class(Atoms), intent(in) :: this
   integer, intent(in) :: order(:)

   this%atoms = this%atoms(order(1:this%natom))
   this%adjmat = this%adjmat(order(1:this%natom), order(1:this%natom))

end subroutine

end module
