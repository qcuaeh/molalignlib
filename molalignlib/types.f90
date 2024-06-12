module types
use kinds
use discrete
use rotation
use translation
use adjacency
use alignment
use bounds
implicit none
private

type :: MNA
   integer, allocatable :: lengths(:)
end type

type :: Atom
   character(:), allocatable, private :: label
   integer, private :: znum
   integer, private :: typeidx
   integer, private :: equividx
   integer, allocatable, private :: mnaid(:)
   real(wp), private :: weight
   real(wp), private :: coords(3)
   integer, private :: coonum
   integer, allocatable, private :: neighbors(:)
   integer :: nneieqv
   integer, allocatable :: neieqvlens(:)
   integer :: nneimna
   integer, allocatable :: neimnalens(:)
contains
   procedure :: print => print_atom
end type

type, public :: Molecule
   character(:), allocatable :: title
   integer :: natom
   type(Atom), allocatable :: atoms(:)
   integer, allocatable :: typeaggs(:)
   integer, allocatable :: equivaggs(:)
   integer :: nfrag
   integer, allocatable :: fragroots(:)
   integer :: nmna
   type(MNA), allocatable :: mnas(:)
contains
   procedure :: get_natom
   procedure :: get_ntype
   procedure :: get_nequiv
   procedure :: get_znums
   procedure :: set_znums
   procedure :: get_center
   procedure :: get_weights
   procedure :: set_weights
   procedure :: get_typeidcs
   procedure :: set_typeidcs
   procedure :: get_equividcs
   procedure :: set_equividcs
   procedure :: get_typeaggs
   procedure :: set_typeaggs
   procedure :: get_equivaggs
   procedure :: set_equivaggs
   procedure :: get_coords
   procedure :: set_coords
   procedure :: get_labels
   procedure :: set_labels
   procedure :: get_adjmat
   procedure :: get_coonums
   procedure :: get_neighbors
   procedure :: set_neighbors
   procedure :: get_fragroot
   procedure :: get_blocks
   procedure :: mirror_coords
   procedure, private :: matrix_rotate_coords
   procedure, private :: quater_rotate_coords
   generic, public :: rotate_coords => matrix_rotate_coords, quater_rotate_coords
   procedure :: translate_coords
   procedure :: permutate_atoms
   procedure :: print => print_molecule
   procedure :: bonded
   procedure :: remove_bond
   procedure :: add_bond
end type

type, public :: Block
   integer, allocatable :: atomidx(:)
end type

contains

integer function get_natom(self) result(natom)
   class(Molecule), intent(in) :: self

   natom = self%natom

end function get_natom

integer function get_ntype(self) result(ntype)
   class(Molecule), intent(in) :: self

   ntype = size(self%typeaggs)

end function get_ntype

integer function get_nequiv(self) result(nequiv)
   class(Molecule), intent(in) :: self

   nequiv = size(self%equivaggs)

end function get_nequiv

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
   integer :: i, k
   integer :: invorder(self%natom)

   invorder = inverse_permut(order)
   self%atoms = self%atoms(order(:))
   do i = 1, self%natom
      do k = 1, self%atoms(i)%coonum
         self%atoms(i)%neighbors(k) = invorder(self%atoms(i)%neighbors(k))
      end do
   end do

end subroutine permutate_atoms

function get_znums(self) result(znums)
   class(Molecule), intent(in) :: self
   integer :: znums(self%natom)

   znums(:) = self%atoms(:)%znum

end function get_znums

subroutine set_znums(self, znums)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: znums(self%natom)

   self%atoms(:)%znum = znums(:)

end subroutine set_znums

subroutine set_equividcs(self, equividcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: equividcs(self%natom)

   self%atoms(:)%equividx = equividcs(:)

end subroutine set_equividcs

subroutine set_typeidcs(self, typeidcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: typeidcs(self%natom)

   self%atoms(:)%typeidx = typeidcs(:)

end subroutine set_typeidcs

function get_center(self) result(cntrcoords)
! Purpose: Get the centroid coordinates
   class(Molecule), intent(in) :: self
   real(wp) :: cntrcoords(3)
   integer :: i

! Calculate the coordinates of the center of mass

   cntrcoords(:) = 0

   do i = 1, self%natom
      cntrcoords(:) = cntrcoords(:) + self%atoms(i)%weight*self%atoms(i)%coords(:)
   end do

   cntrcoords(:) = cntrcoords(:)/sum(self%get_weights())

end function get_center

function get_blocks(self) result(blocks)
   class(Molecule), intent(in) :: self
   type(Block) :: blocks(size(self%typeaggs))
   integer :: h, i, k(size(self%typeaggs))

   k(:) = 0

   do h = 1, size(self%typeaggs)
      allocate(blocks(h)%atomidx(self%typeaggs(h)))
   end do

   do i = 1, self%natom
      k(self%atoms(i)%typeidx) = k(self%atoms(i)%typeidx) + 1
      blocks(self%atoms(i)%typeidx)%atomidx(k(self%atoms(i)%typeidx)) = i
   end do

end function get_blocks

function get_typeidcs(self) result(typeidcs)
   class(Molecule), intent(in) :: self
   integer :: typeidcs(self%natom)

   typeidcs(:) = self%atoms(:)%typeidx

end function get_typeidcs

function get_equividcs(self) result(equividcs)
   class(Molecule), intent(in) :: self
   integer :: equividcs(self%natom)

   equividcs(:) = self%atoms(:)%equividx

end function get_equividcs

function get_weights(self) result(weights)
   class(Molecule), intent(in) :: self
   real(wp) :: weights(self%natom)

   weights(:) = self%atoms(:)%weight

end function get_weights

subroutine set_weights(self, weights)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: weights(self%natom)

   self%atoms(:)%weight = weights(:)

end subroutine set_weights

function get_coords(self) result(coords)
   class(Molecule), intent(in) :: self
   real(wp) :: coords(3, self%natom)
   integer :: i

!   coords(:, :) = self%atoms(:)%coords(:)
   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords(:)
   end do

end function get_coords

subroutine set_coords(self, coords)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: coords(3, self%natom)
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine set_coords

function get_labels(self) result(labels)
   class(Molecule), intent(in) :: self
   character(wl) :: labels(self%natom)
   integer :: i

   do i = 1, self%natom
      labels(i) = self%atoms(i)%label
   end do

end function get_labels

subroutine set_labels(self, labels)
   class(Molecule), intent(inout) :: self
   character(wl), intent(in) :: labels(self%natom)
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%label = labels(i)
   end do

end subroutine set_labels

function get_adjmat(self) result(adjmat)
   class(Molecule), intent(in) :: self
   type(Atom) :: iatom
   logical :: adjmat(self%natom, self%natom)
   integer :: i, k

   adjmat(:, :) = .false.

   do i = 1, self%natom
      iatom = self%atoms(i)
      do k = 1, iatom%coonum
         adjmat(i, iatom%neighbors(k)) = .true.
      end do
   end do

end function get_adjmat

function get_coonums(self) result(coonums)
   class(Molecule), intent(in) :: self
   integer, allocatable :: coonums(:)

   coonums = self%atoms(:)%coonum

end function get_coonums

subroutine set_neighbors(self, coonums, neighbors)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: coonums(:)
   integer, intent(in) :: neighbors(:, :)
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%coonum = coonums(i)
      self%atoms(i)%neighbors = neighbors(:coonums(i), i)
   end do

end subroutine set_neighbors

function get_neighbors(self) result(neighbors)
   class(Molecule), intent(in) :: self
   integer, allocatable :: neighbors(:, :)
   integer :: i

   allocate(neighbors(maxcoord, self%natom))

   do i = 1, self%natom
      neighbors(:self%atoms(i)%coonum, i) = self%atoms(i)%neighbors
   end do

end function get_neighbors

function get_fragroot(self) result(fragroots)
   class(Molecule), intent(in) :: self
   integer :: fragroots(self%natom)

   fragroots = self%fragroots

end function get_fragroot

subroutine print_atom(self, ind, outLvl)
   class(Atom), intent(in) :: self
   integer, intent(in) :: ind
   integer, intent(in), optional :: outLvl
   character(255) :: frmt
   character(2) :: num
   integer :: outLevel

! *** code to manage unitlbl pending ***

   outLevel = 1
   if (present(outLvl)) outLevel = outLvl
   
   select case (outLevel)

   case (1)
      if (self%coonum == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%coonum
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
               self%weight, " {", self%coords(:), " } [", &
               self%neighbors(:self%coonum), " ]"
      end if

!  case (2)
!
   case default
      if (self%coonum == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%coonum
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
               self%weight, " {", self%coords(:), " } [", &
               self%neighbors(:self%coonum), " ]"
      end if
   end select

end subroutine print_atom

subroutine print_molecule(self)
   class(Molecule), intent(in) :: self
   integer :: i

! *** code to manage unitlbl pending ***
   
   write (stderr, '(a,i0,a)') "Contents of molecule structure:   (", &
                                         self%natom, " atoms)"
   write (stderr, '(2a)') 'Title: ', self%title
   write (stderr, '(a,a4,a5,a6,a7,2a17)') "ind:", "lbl", "znum", "typeidx", &
                                          "weight", "{ coords }", "[ neighbors ]"

   do i = 1, self%natom
      call self%atoms(i)%print(i)
   end do

end subroutine print_molecule

function bonded(self, ind1, ind2) result(isbond)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: ind1, ind2
   integer :: i
   logical :: isbond, found1, found2

! initialization
   found1 = .false.
   found2 = .false.

! check cross reference of ind1 and ind2 in both neighbors
   do i = 1, self%atoms(ind1)%coonum
      if (ind2 == self%atoms(ind1)%neighbors(i)) then
         found1 = .true.
         exit
      end if
   end do
   do i = 1, self%atoms(ind2)%coonum
      if (ind1 == self%atoms(ind2)%neighbors(i)) then
         found2 = .true.
         exit
      end if
   end do

! report or stop program
   if (found1 .and. found2) then
      isbond = .true.
   else if (found1 .or. found2) then
      write (stderr, '(a,i0,2x,i0)') 'Inconsistent bond for atoms: ', ind1, ind2
      stop
   else
      isbond = .false.
   end if

end function bonded

subroutine remove_bond(self, ind1, ind2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: ind1, ind2
   integer :: i, pos1, pos2

! initialization
   pos1 = 0   ! position of ind2 in neighbors of atom 1
   pos2 = 0   ! position of ind1 in neighbors of atom 2

! find position of ind2 and ind1 in neighbors of atoms ind1 and ind2, resp.
   do i = 1, self%atoms(ind1)%coonum
      if (ind2 == self%atoms(ind1)%neighbors(i)) then
         pos1 = i
         exit
      end if
   end do
   do i = 1, self%atoms(ind2)%coonum
      if (ind1 == self%atoms(ind2)%neighbors(i)) then
         pos2 = i
         exit
      end if
   end do

! delete ind2 and ind1 from the ajdlists where they appear
   if ((pos1 /= 0) .and. (pos2 /= 0)) then
      self%atoms(ind1)%coonum = self%atoms(ind1)%coonum - 1
      do i = pos1, self%atoms(ind1)%coonum
         self%atoms(ind1)%neighbors(i) = self%atoms(ind1)%neighbors(i+1)
      end do
      self%atoms(ind2)%coonum = self%atoms(ind2)%coonum - 1
      do i = pos2, self%atoms(ind2)%coonum
         self%atoms(ind2)%neighbors(i) = self%atoms(ind2)%neighbors(i+1)
      end do
   else
      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', ind1, ind2
   end if

end subroutine remove_bond

subroutine add_bond(self, ind1, ind2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: ind1, ind2
   integer :: pos1, pos2

! initialization
   pos1 = self%atoms(ind1)%coonum
   pos2 = self%atoms(ind2)%coonum

   if (.not. self%bonded(ind1, ind2)) then
! indices in neighbors are supposed to be sorted; inserting new indices
      self%atoms(ind1)%coonum = self%atoms(ind1)%coonum + 1
! find position to insert ind2 and shift indices greater than ind2
      do while ((pos1 >= 1) .and. (ind2 < self%atoms(ind1)%neighbors(pos1)))
         self%atoms(ind1)%neighbors(pos1+1) = self%atoms(ind1)%neighbors(pos1)
         pos1 = pos1 - 1
      end do
      self%atoms(ind1)%neighbors(pos1+1) = ind2
      
      self%atoms(ind2)%coonum = self%atoms(ind2)%coonum + 1
! find position to insert ind1 and shift indices greater than ind1
      do while ((pos2 >= 1) .and. (ind1 < self%atoms(ind2)%neighbors(pos2)))
         self%atoms(ind2)%neighbors(pos2+1) = self%atoms(ind2)%neighbors(pos2)
         pos2 = pos2 - 1
      end do
      self%atoms(ind2)%neighbors(pos2+1) = ind1
   else
      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", ind1, ind2
   end if

end subroutine add_bond

function get_typeaggs(self) result(typeaggs)
   class(Molecule), intent(in) :: self
   integer, allocatable :: typeaggs(:)

   typeaggs = self%typeaggs

end function get_typeaggs

subroutine set_typeaggs(self, ntype, typeaggs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: ntype
   integer, intent(in) :: typeaggs(:)

!   allocate(self%typeaggs(ntype))
   self%typeaggs = typeaggs(:ntype)

end subroutine set_typeaggs

function get_equivaggs(self) result(equivaggs)
   class(Molecule), intent(in) :: self
   integer, allocatable :: equivaggs(:)

   equivaggs = self%equivaggs

end function get_equivaggs

subroutine set_equivaggs(self, nequiv, equivaggs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nequiv
   integer, intent(in) :: equivaggs(:)

!   allocate(self%equivaggs(nequiv))
   self%equivaggs = equivaggs(:nequiv)

end subroutine set_equivaggs

end module
