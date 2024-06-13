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
   integer, private :: nadj
   integer, allocatable, private :: adjlist(:)
   integer, allocatable :: adeqlenlist(:)
   integer, allocatable :: mnalenlist(:)
contains
   procedure :: print => print_atom
   procedure :: get_adjlist
   procedure :: set_adjlist
end type

type, public :: Molecule
   character(:), allocatable :: title
   integer :: natom
   type(Atom), allocatable :: atoms(:)
   integer, allocatable :: typelenlist(:)
   integer, allocatable :: equivlenlist(:)
   integer :: nfrag
   integer, allocatable :: fragroots(:)
   type(MNA), allocatable :: mnas(:)
contains
   procedure :: get_natom
   procedure :: get_ntype
   procedure :: get_nequiv
   procedure :: get_nadjs
   procedure :: get_nadeqs
   procedure :: get_typeidcs
   procedure :: set_typeidcs
   procedure :: get_equividcs
   procedure :: set_equividcs
   procedure :: get_typelenlist
   procedure :: set_typelenlist
   procedure :: get_equivlenlist
   procedure :: set_equivlenlist
   procedure :: get_adeqlenlists
   procedure :: set_adeqlenlists
   procedure :: get_coords
   procedure :: set_coords
   procedure :: get_labels
   procedure :: set_labels
   procedure :: get_adjmat
   procedure :: get_adjlists
   procedure :: set_adjlists
   procedure :: get_fragroot
   procedure :: get_blocks
   procedure :: get_znums
   procedure :: set_znums
   procedure :: get_weights
   procedure :: set_weights
   procedure :: get_center
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

   ntype = size(self%typelenlist)

end function get_ntype

integer function get_nequiv(self) result(nequiv)
   class(Molecule), intent(in) :: self

   nequiv = size(self%equivlenlist)

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
      do k = 1, self%atoms(i)%nadj
         self%atoms(i)%adjlist(k) = invorder(self%atoms(i)%adjlist(k))
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
   type(Block) :: blocks(size(self%typelenlist))
   integer :: h, i, k(size(self%typelenlist))

   k(:) = 0

   do h = 1, size(self%typelenlist)
      allocate(blocks(h)%atomidx(self%typelenlist(h)))
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
      do k = 1, size(iatom%adjlist)
         adjmat(i, iatom%adjlist(k)) = .true.
      end do
   end do

end function get_adjmat

function get_nadjs(self) result(nadjs)
   class(Molecule), intent(in) :: self
   integer, allocatable :: nadjs(:)
   integer :: i

   allocate(nadjs(size(self%atoms)))
   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(i)%adjlist)
   end do

end function get_nadjs

subroutine set_adjlist(self, nadj, adjlist)
   class(Atom), intent(inout) :: self
   integer, intent(in) :: nadj
   integer, dimension(nadj), intent(in) :: adjlist

   self%adjlist = adjlist
! provisional:
   self%nadj = nadj

end subroutine set_adjlist

function get_adjlist(self) result(adjlist)
   class(Atom), intent(in) :: self
   integer, dimension(:), allocatable :: adjlist

   adjlist = self%adjlist

end function get_adjlist

subroutine set_adjlists(self, nadjs, adjlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nadjs(:)
   integer, intent(in) :: adjlists(:, :)
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%nadj = nadjs(i)
      self%atoms(i)%adjlist = adjlists(:nadjs(i), i)
   end do

end subroutine set_adjlists

function get_adjlists(self) result(adjlists)
   class(Molecule), intent(in) :: self
   integer, allocatable :: adjlists(:, :)
   integer :: i

   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      adjlists(:self%atoms(i)%nadj, i) = self%atoms(i)%adjlist
   end do

end function get_adjlists

subroutine set_adeqlenlists(self, nadeqs, adeqlenlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nadeqs(:)
   integer, intent(in) :: adeqlenlists(:, :)
   integer :: i

   do i = 1, size(self%atoms)
      if (.not. allocated(self%atoms(i)%adeqlenlist)) then
         allocate(self%atoms(i)%adeqlenlist(nadeqs(i)))
      end if
      self%atoms(i)%adeqlenlist = adeqlenlists(:nadeqs(i), i)
   end do

end subroutine set_adeqlenlists

function get_adeqlenlists(self) result(adeqlenlists)
   class(Molecule), intent(inout) :: self
   integer, allocatable :: adeqlenlists(:, :)
   integer :: i

   allocate(adeqlenlists(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      adeqlenlists(:size(self%atoms(i)%adeqlenlist), i) = self%atoms(i)%adeqlenlist
   end do

end function get_adeqlenlists

function get_nadeqs(self) result(nadeqs)
   class(Molecule), intent(inout) :: self
   integer, allocatable :: nadeqs(:)
   integer :: i

   allocate(nadeqs(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadeqs(i) = size(self%atoms(i)%adeqlenlist)
   end do

end function get_nadeqs

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
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:self%nadj), " ]"
      end if

!  case (2)
!
   case default
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%typeidx, &
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
   write (stderr, '(a,a4,a5,a6,a7,2a17)') "ind:", "lbl", "znum", "typeidx", &
                                          "weight", "{ coords }", "[ adjlist ]"

   do i = 1, self%natom
      call self%atoms(i)%print(i)
   end do

end subroutine print_molecule

function bonded(self, idx1, idx2) result(isbond)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: idx1, idx2
   integer, dimension(maxcoord) :: adjlist1, adjlist2
   integer :: i
   logical :: isbond, found1, found2

! copy arrays of adjlist
   adjlist1 = self%atoms(idx1)%get_adjlist()
   adjlist2 = self%atoms(idx2)%get_adjlist()

! initialization
   found1 = .false.
   found2 = .false.

! check cross reference of idx1 and idx2 in both adjlists
   do i = 1, size(self%atoms(idx1)%adjlist)
      if (idx2 == adjlist1(i)) then
         found1 = .true.
         exit
      end if
   end do
   do i = 1, size(self%atoms(idx2)%adjlist)
      if (idx1 == adjlist2(i)) then
         found2 = .true.
         exit
      end if
   end do

! report or stop program
   if (found1 .and. found2) then
      isbond = .true.
   else if (found1 .or. found2) then
      write (stderr, '(a,i0,2x,i0)') 'Inconsistent bond for atoms: ', idx1, idx2
      stop
   else
      isbond = .false.
   end if

end function bonded

subroutine remove_bond(self, idx1, idx2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   integer, dimension(maxcoord) :: adjlist1, adjlist2
   integer :: i, pos1, pos2, nadj1, nadj2

! copy adjlist arrays
   nadj1 = size(self%atoms(idx1)%adjlist)
   nadj2 = size(self%atoms(idx2)%adjlist)
   adjlist1 = self%atoms(idx1)%get_adjlist()
   adjlist2 = self%atoms(idx2)%get_adjlist()

! initialization
   pos1 = 0   ! position of idx2 in adjlist of atom 1
   pos2 = 0   ! position of idx1 in adjlist of atom 2

! find position of idx2 and idx1 in adjlist of atoms idx1 and idx2, resp.
   do i = 1, nadj1
      if (idx2 == adjlist1(i)) then
         pos1 = i
         exit
      end if
   end do
   do i = 1, nadj2
      if (idx1 == adjlist2(i)) then
         pos2 = i
         exit
      end if
   end do

! delete idx2 and idx1 from the ajdlists where they appear
   if ((pos1 /= 0) .and. (pos2 /= 0)) then
      nadj1 = nadj1 - 1
      do i = pos1, nadj1
         adjlist1(i) = adjlist1(i+1)
      end do
      nadj2 = nadj2 - 1
      do i = pos2, nadj2
         adjlist2(i) = adjlist2(i+1)
      end do
! update neighbor arrays for atoms idx1 and idx2
      call self%atoms(idx1)%set_adjlist(nadj1, adjlist1(:nadj1))
      call self%atoms(idx2)%set_adjlist(nadj2, adjlist2(:nadj2))
   else
      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', idx1, idx2
   end if

end subroutine remove_bond

subroutine add_bond(self, idx1, idx2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   integer, dimension(maxcoord) :: adjlist1, adjlist2
   integer :: pos1, pos2, nadj1, nadj2

! copy array of adjlist
   nadj1 = size(self%atoms(idx1)%adjlist)
   nadj2 = size(self%atoms(idx2)%adjlist)
   adjlist1 = self%atoms(idx1)%get_adjlist()
   adjlist2 = self%atoms(idx2)%get_adjlist()

! initialization
   pos1 = nadj1
   pos2 = nadj2

   if (.not. self%bonded(idx1, idx2)) then
! indices in adjlist are supposed to be sorted; inserting new indices
!      self%atoms(idx1)%nadj = self%atoms(idx1)%nadj + 1
      nadj1 = nadj1 + 1
! find position to insert idx2 and shift indices greater than idx2
      do while ((pos1 >= 1) .and. (idx2 < adjlist1(pos1)))
         adjlist1(pos1+1) = adjlist1(pos1)
         pos1 = pos1 - 1
      end do
      adjlist1(pos1+1) = idx2
      
      nadj2 = nadj2 + 1
! find position to insert idx1 and shift indices greater than idx1
      do while ((pos2 >= 1) .and. (idx1 < adjlist2(pos2)))
         adjlist2(pos2+1) = adjlist2(pos2)
         pos2 = pos2 - 1
      end do
      adjlist2(pos2+1) = idx1
! update neighbor arrays for atoms in idx1 and idx2
      call self%atoms(idx1)%set_adjlist(nadj1, adjlist1(:nadj1))
      call self%atoms(idx2)%set_adjlist(nadj2, adjlist2(:nadj2))
   else
      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", idx1, idx2
   end if

end subroutine add_bond

function get_typelenlist(self) result(typelenlist)
   class(Molecule), intent(in) :: self
   integer, allocatable :: typelenlist(:)

   typelenlist = self%typelenlist

end function get_typelenlist

subroutine set_typelenlist(self, ntype, typelenlist)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: ntype
   integer, intent(in) :: typelenlist(:)

   allocate(self%typelenlist(ntype))
   self%typelenlist = typelenlist(:ntype)

end subroutine set_typelenlist

function get_equivlenlist(self) result(equivlenlist)
   class(Molecule), intent(in) :: self
   integer, allocatable :: equivlenlist(:)

   equivlenlist = self%equivlenlist

end function get_equivlenlist

subroutine set_equivlenlist(self, nequiv, equivlenlist)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nequiv
   integer, intent(in) :: equivlenlist(:)

   allocate(self%equivlenlist(nequiv))
   self%equivlenlist = equivlenlist(:nequiv)

end subroutine set_equivlenlist

end module
