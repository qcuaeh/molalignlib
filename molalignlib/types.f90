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
   integer, private :: atomtypeidx
   integer, private :: atomequividx
   integer, allocatable, private :: mnaid(:)
   real(wp), private :: weight
   real(wp), private :: coords(3)
   integer, private :: nadj
   integer, allocatable, private :: adjlist(:)
   integer, allocatable :: adjequivlenlist(:)
contains
   procedure :: print => print_atom
end type

type, public :: Molecule
   character(:), allocatable :: title
   integer :: natom
   type(Atom), allocatable :: atoms(:)
   integer, allocatable :: atomtypelenlist(:)
   integer, allocatable :: atomequivlenlist(:)
   integer :: nfrag
   integer, allocatable :: fragroots(:)
   type(MNA), allocatable :: mnas(:)
contains
   procedure :: get_natom
   procedure :: get_natomtype
   procedure :: get_natomequiv
   procedure :: get_nadjs
   procedure :: get_nadjequivs
   procedure :: get_atomtypeidcs
   procedure :: set_atomtypeidcs
   procedure :: get_atomequividcs
   procedure :: set_atomequividcs
   procedure :: get_atomtypelenlist
   procedure :: set_atomtypelenlist
   procedure :: get_atomequivlenlist
   procedure :: set_atomequivlenlist
   procedure :: get_adjequivlenlists
   procedure :: set_adjequivlenlists
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

integer function get_natomtype(self) result(natomtype)
   class(Molecule), intent(in) :: self

   natomtype = size(self%atomtypelenlist)

end function get_natomtype

integer function get_natomequiv(self) result(natomequiv)
   class(Molecule), intent(in) :: self

   natomequiv = size(self%atomequivlenlist)

end function get_natomequiv

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

subroutine set_atomequividcs(self, atomequividcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomequividcs(self%natom)

   self%atoms(:)%atomequividx = atomequividcs(:)

end subroutine set_atomequividcs

subroutine set_atomtypeidcs(self, atomtypeidcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomtypeidcs(self%natom)

   self%atoms(:)%atomtypeidx = atomtypeidcs(:)

end subroutine set_atomtypeidcs

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
   type(Block) :: blocks(size(self%atomtypelenlist))
   integer :: h, i, k(size(self%atomtypelenlist))

   k(:) = 0

   do h = 1, size(self%atomtypelenlist)
      allocate(blocks(h)%atomidx(self%atomtypelenlist(h)))
   end do

   do i = 1, self%natom
      k(self%atoms(i)%atomtypeidx) = k(self%atoms(i)%atomtypeidx) + 1
      blocks(self%atoms(i)%atomtypeidx)%atomidx(k(self%atoms(i)%atomtypeidx)) = i
   end do

end function get_blocks

function get_atomtypeidcs(self) result(atomtypeidcs)
   class(Molecule), intent(in) :: self
   integer :: atomtypeidcs(self%natom)

   atomtypeidcs(:) = self%atoms(:)%atomtypeidx

end function get_atomtypeidcs

function get_atomequividcs(self) result(atomequividcs)
   class(Molecule), intent(in) :: self
   integer :: atomequividcs(self%natom)

   atomequividcs(:) = self%atoms(:)%atomequividx

end function get_atomequividcs

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

subroutine set_adjequivlenlists(self, nadjequivs, adjequivlenlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nadjequivs(:)
   integer, intent(in) :: adjequivlenlists(:, :)
   integer :: i

   do i = 1, size(self%atoms)
      if (.not. allocated(self%atoms(i)%adjequivlenlist)) then
         allocate(self%atoms(i)%adjequivlenlist(nadjequivs(i)))
      end if
      self%atoms(i)%adjequivlenlist = adjequivlenlists(:nadjequivs(i), i)
   end do

end subroutine set_adjequivlenlists

function get_adjequivlenlists(self) result(adjequivlenlists)
   class(Molecule), intent(inout) :: self
   integer, allocatable :: adjequivlenlists(:, :)
   integer :: i

   allocate(adjequivlenlists(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      adjequivlenlists(:size(self%atoms(i)%adjequivlenlist), i) = self%atoms(i)%adjequivlenlist
   end do

end function get_adjequivlenlists

function get_nadjequivs(self) result(nadjequivs)
   class(Molecule), intent(inout) :: self
   integer, allocatable :: nadjequivs(:)
   integer :: i

   allocate(nadjequivs(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadjequivs(i) = size(self%atoms(i)%adjequivlenlist)
   end do

end function get_nadjequivs

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
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%atomtypeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%atomtypeidx, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:self%nadj), " ]"
      end if

!  case (2)
!
   case default
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%atomtypeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%label(:2), self%znum, self%atomtypeidx, &
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
   write (stderr, '(a,a4,a5,a6,a7,2a17)') "ind:", "lbl", "znum", "atomtypeidx", &
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
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

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
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

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
      self%atoms(idx1)%adjlist = adjlist1(:nadj1)
      self%atoms(idx2)%adjlist = adjlist2(:nadj2)
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
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

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
      self%atoms(idx1)%adjlist = adjlist1(:nadj1)
      self%atoms(idx2)%adjlist = adjlist2(:nadj2)
   else
      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", idx1, idx2
   end if

end subroutine add_bond

function get_atomtypelenlist(self) result(atomtypelenlist)
   class(Molecule), intent(in) :: self
   integer, allocatable :: atomtypelenlist(:)

   atomtypelenlist = self%atomtypelenlist

end function get_atomtypelenlist

subroutine set_atomtypelenlist(self, natomtype, atomtypelenlist)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: natomtype
   integer, intent(in) :: atomtypelenlist(:)

   allocate(self%atomtypelenlist(natomtype))
   self%atomtypelenlist = atomtypelenlist(:natomtype)

end subroutine set_atomtypelenlist

function get_atomequivlenlist(self) result(atomequivlenlist)
   class(Molecule), intent(in) :: self
   integer, allocatable :: atomequivlenlist(:)

   atomequivlenlist = self%atomequivlenlist

end function get_atomequivlenlist

subroutine set_atomequivlenlist(self, natomequiv, atomequivlenlist)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: natomequiv
   integer, intent(in) :: atomequivlenlist(:)

   allocate(self%atomequivlenlist(natomequiv))
   self%atomequivlenlist = atomequivlenlist(:natomequiv)

end subroutine set_atomequivlenlist

end module
