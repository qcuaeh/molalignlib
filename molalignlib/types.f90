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

type :: Atom
   character(:), allocatable :: label
   integer :: znum
   integer :: blkid
   integer :: eqvid
   real(wp) :: weight
   real(wp) :: coords(3)
   integer :: nadj
   integer, allocatable :: adjlist(:)
   integer :: nadjeqv
   integer, allocatable :: adjeqvlens(:)
contains
   procedure :: print => print_atom
end type

type, public :: Molecule
   character(:), allocatable :: title
   integer :: natom
   type(Atom), allocatable :: atoms(:)
   integer :: nblock
   integer, allocatable :: blklens(:)
   integer :: nequiv
   integer, allocatable :: eqvlens(:)
   integer :: nfrag
   integer, allocatable :: fragroot(:)
contains
   procedure :: get_natom
   procedure :: get_nblock
   procedure :: get_znums
   procedure :: get_center
   procedure :: get_weights
   procedure :: get_blkids
   procedure :: get_eqvids
   procedure :: get_blklen
   procedure :: get_blklens
   procedure :: get_eqvlen
   procedure :: get_coords
   procedure :: set_coords
   procedure :: get_labels
   procedure :: get_adjmat
   procedure :: get_nadj
   procedure :: get_adjlist
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
   integer, allocatable :: atmid(:)
end type

contains

integer function get_natom(self) result(natom)
   class(Molecule), intent(in) :: self

   natom = self%natom

end function get_natom

integer function get_nblock(self) result(nblock)
   class(Molecule), intent(in) :: self

   nblock = self%nblock

end function get_nblock

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
   integer i, k, invorder(self%natom)

   invorder = inverse_permut(order)
   self%atoms = self%atoms(order(:))
   do i = 1, self%natom
      do k = 1, self%atoms(i)%nadj
         self%atoms(i)%adjlist(k) = invorder(self%atoms(i)%adjlist(k))
      end do
   end do

end subroutine permutate_atoms

function get_znums(self) result(znums)
   class(Molecule), intent(in) :: self
   integer :: znums(self%natom)
   integer i

   znums(:) = self%atoms(:)%znum

end function get_znums

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
   type(Block) :: blocks(self%nblock)
   integer :: h, i, k(self%nblock)

   k(:) = 0

   do h = 1, self%nblock
      allocate(blocks(h)%atmid(self%blklens(h)))
   end do

   do i = 1, self%natom
      k(self%atoms(i)%blkid) = k(self%atoms(i)%blkid) + 1
      blocks(self%atoms(i)%blkid)%atmid(k(self%atoms(i)%blkid)) = i
   end do

end function get_blocks

function get_blkids(self) result(blkids)
   class(Molecule), intent(in) :: self
   integer :: blkids(self%natom)
   integer i

   blkids(:) = self%atoms(:)%blkid

end function get_blkids

function get_eqvids(self) result(eqvids)
   class(Molecule), intent(in) :: self
   integer :: eqvids(self%natom)
   integer i

   eqvids(:) = self%atoms(:)%eqvid

end function get_eqvids

function get_weights(self) result(weights)
   class(Molecule), intent(in) :: self
   real(wp) :: weights(self%natom)
   integer i

   weights(:) = self%atoms(:)%weight

end function get_weights

function get_coords(self) result(coords)
   class(Molecule), intent(in) :: self
   real(wp) :: coords(3, self%natom)
   integer i

!   coords(:, :) = self%atoms(:)%coords(:)
   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords(:)
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
         adjmat(i, iatom%adjlist(k)) = .true.
      end do
   end do

end function get_adjmat

function get_nadj(self) result(nadj)
   class(Molecule), intent(in) :: self
   integer :: nadj(self%natom)
   integer i

   do i = 1, self%natom
      nadj(i) = self%atoms(i)%nadj
   end do

end function get_nadj

function get_adjlist(self) result(adjlist)
   class(Molecule), intent(in) :: self
   integer :: adjlist(maxcoord, self%natom)
   integer i

   do i = 1, self%natom
      adjlist(:self%atoms(i)%nadj, i) = self%atoms(i)%adjlist(:self%atoms(i)%nadj)
   end do

end function get_adjlist

function get_fragroot(self) result(fragroot)
   class(Molecule), intent(in) :: self
   integer :: fragroot(self%natom)

   fragroot = self%fragroot

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
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%blkid, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%blkid, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:self%nadj), " ]"
      end if

!  case (2)
!
   case default
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%blkid, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%blkid, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:self%nadj), " ]"
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
   write (stderr, '(a,a4,a5,a6,a7,2a17)') "ind:", "lbl", "znum", "blkid", &
                                          "weight", "{ coords }", "[ adjlist ]"

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

! check cross reference of ind1 and ind2 in both adjlists
   do i = 1, self%atoms(ind1)%nadj
      if (ind2 == self%atoms(ind1)%adjlist(i)) then
         found1 = .true.
         exit
      end if
   end do
   do i = 1, self%atoms(ind2)%nadj
      if (ind1 == self%atoms(ind2)%adjlist(i)) then
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
   integer i, pos1, pos2

! initialization
   pos1 = 0   ! position of ind2 in adjlist of atom 1
   pos2 = 0   ! position of ind1 in adjlist of atom 2

! find position of ind2 and ind1 in adjlists of atoms ind1 and ind2, resp.
   do i = 1, self%atoms(ind1)%nadj
      if (ind2 == self%atoms(ind1)%adjlist(i)) then
         pos1 = i
         exit
      end if
   end do
   do i = 1, self%atoms(ind2)%nadj
      if (ind1 == self%atoms(ind2)%adjlist(i)) then
         pos2 = i
         exit
      end if
   end do

! delete ind2 and ind1 from the ajdlists where they appear
   if ((pos1 /= 0) .and. (pos2 /= 0)) then
      self%atoms(ind1)%nadj = self%atoms(ind1)%nadj - 1
      do i = pos1, self%atoms(ind1)%nadj
         self%atoms(ind1)%adjlist(i) = self%atoms(ind1)%adjlist(i+1)
      end do
      self%atoms(ind2)%nadj = self%atoms(ind2)%nadj - 1
      do i = pos2, self%atoms(ind2)%nadj
         self%atoms(ind2)%adjlist(i) = self%atoms(ind2)%adjlist(i+1)
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
   pos1 = self%atoms(ind1)%nadj
   pos2 = self%atoms(ind2)%nadj

   if (.not. self%bonded(ind1, ind2)) then
! indices in adjlists are supposed to be sorted; inserting new indices
      self%atoms(ind1)%nadj = self%atoms(ind1)%nadj + 1
! find position to insert ind2 and shift indices greater than ind2
      do while ((pos1 >= 1) .and. (ind2 < self%atoms(ind1)%adjlist(pos1)))
         self%atoms(ind1)%adjlist(pos1+1) = self%atoms(ind1)%adjlist(pos1)
         pos1 = pos1 - 1
      end do
      self%atoms(ind1)%adjlist(pos1+1) = ind2
      
      self%atoms(ind2)%nadj = self%atoms(ind2)%nadj + 1
! find position to insert ind1 and shift indices greater than ind1
      do while ((pos2 >= 1) .and. (ind1 < self%atoms(ind2)%adjlist(pos2)))
         self%atoms(ind2)%adjlist(pos2+1) = self%atoms(ind2)%adjlist(pos2)
         pos2 = pos2 - 1
      end do
      self%atoms(ind2)%adjlist(pos2+1) = ind1
   else
      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", ind1, ind2
   end if

end subroutine add_bond

integer function get_blklen(self, blockid) result(blklen)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: blockid
   integer :: i

   blklen = 0

   do i = 1, self%natom
      if (self%atoms(i)%blkid == blockid) then
         blklen = blklen + 1
      end if
   end do

end function get_blklen

function get_blklens(self) result(blklens)
   class(Molecule), intent(in) :: self
   integer, allocatable :: blklens(:)

   blklens = self%blklens

end function get_blklens

integer function get_eqvlen(self, equivid) result(eqvlen)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: equivid
   integer :: i

   eqvlen = 0

   do i = 1, self%natom
      if (self%atoms(i)%eqvid == equivid) then
         eqvlen = eqvlen + 1
      end if
   end do

end function get_eqvlen



end module
