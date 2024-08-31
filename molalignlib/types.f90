module types
use kinds
use bounds
use discrete
use rotation
use translation
use adjacency
use alignment
use strutils
implicit none
private

type, public :: Part
   integer, allocatable :: atomidcs(:)
end type

type, public :: Partition
   type(Part), allocatable :: parts(:)
contains
   procedure :: set_partition
   procedure :: get_lenlist
end type

type :: Atom
   integer, private :: elnum
   integer, private :: label
   integer, private :: typeidx
   integer, private :: equividx
   integer, private :: fragidx
   integer, allocatable, private :: mnaid(:)
   real(wp), private :: weight
   real(wp), private :: coords(3)
   integer, allocatable, private :: adjlist(:)
   integer, allocatable :: adjequivlenlist(:)
contains
   procedure :: print => print_atom
end type

type, public :: Molecule
   integer :: natom
   integer :: nfrag
   character(:), allocatable :: title
   type(Atom), allocatable :: atoms(:)
   type(Partition) :: molfragpartition
   type(Partition) :: atomtypepartition
   type(Partition) :: atomequivpartition
contains
   procedure :: set_atomelnums
   procedure :: set_atomlabels
   procedure :: set_atomcoords
   procedure :: set_weights
   procedure :: set_adjlists
   procedure :: set_molfrags
   procedure :: set_atomtypes
   procedure :: set_atomequivs
   procedure :: set_adjequivlenlists
   procedure :: get_natom
   procedure :: get_natomtype
   procedure :: get_natomequiv
   procedure :: get_atomtypelenlist
   procedure :: get_atomequivlenlist
   procedure :: mirror_coords
   procedure :: translate_coords
   procedure :: permutate_atoms
   procedure :: bonded
   procedure :: remove_bond
   procedure :: add_bond
   procedure :: get_fragroots
   procedure :: print => print_molecule
   procedure, private :: get_inner_nadjs
   procedure, private :: get_sorted_nadjs
   procedure, private :: get_inner_nadjequivs
   procedure, private :: get_sorted_nadjequivs
   procedure, private :: get_inner_atomtypeidcs
   procedure, private :: get_sorted_atomtypeidcs
   procedure, private :: get_inner_atomequividcs
   procedure, private :: get_sorted_atomequividcs
   procedure, private :: get_inner_adjequivlenlists
   procedure, private :: get_sorted_adjequivlenlists
   procedure, private :: get_inner_atomcoords
   procedure, private :: get_sorted_atomcoords
   procedure, private :: get_inner_adjlists
   procedure, private :: get_sorted_adjlists
   procedure, private :: get_inner_adjmat
   procedure, private :: get_sorted_adjmat
   procedure, private :: get_inner_fragparts
   procedure, private :: get_sorted_fragparts
   procedure, private :: get_inner_atomelnums
   procedure, private :: get_sorted_atomelnums
   procedure, private :: get_inner_atomlabels
   procedure, private :: get_sorted_atomlabels
   procedure, private :: get_inner_atomweights
   procedure, private :: get_sorted_atomweights
   procedure, private :: get_inner_atomtypeparts
   procedure, private :: get_sorted_atomtypeparts
   procedure, private :: matrix_rotate_coords
   procedure, private :: quater_rotate_coords
   generic :: get_nadjs => get_inner_nadjs, get_sorted_nadjs
   generic :: get_nadjequivs => get_inner_nadjequivs, get_sorted_nadjequivs
   generic :: get_atomtypeidcs => get_inner_atomtypeidcs, get_sorted_atomtypeidcs
   generic :: get_atomequividcs => get_inner_atomequividcs, get_sorted_atomequividcs
   generic :: get_adjequivlenlists => get_inner_adjequivlenlists, get_sorted_adjequivlenlists
   generic :: get_atomcoords => get_inner_atomcoords, get_sorted_atomcoords
   generic :: get_adjlists => get_inner_adjlists, get_sorted_adjlists
   generic :: get_adjmat => get_inner_adjmat, get_sorted_adjmat
   generic :: get_fragparts => get_inner_fragparts, get_sorted_fragparts
   generic :: get_atomtypeparts => get_inner_atomtypeparts, get_sorted_atomtypeparts
   generic :: get_atomelnums => get_inner_atomelnums, get_sorted_atomelnums
   generic :: get_atomlabels => get_inner_atomlabels, get_sorted_atomlabels
   generic :: get_atomweights => get_inner_atomweights, get_sorted_atomweights
   generic :: rotate_coords => matrix_rotate_coords, quater_rotate_coords
end type

contains

integer function get_natom(self) result(natom)
   class(Molecule), intent(in) :: self

   natom = size(self%atoms)

end function

integer function get_natomtype(self) result(natomtype)
   class(Molecule), intent(in) :: self

   natomtype = size(self%atomtypepartition%parts)

end function

integer function get_natomequiv(self) result(natomequiv)
   class(Molecule), intent(in) :: self

   natomequiv = size(self%atomequivpartition%parts)

end function

subroutine mirror_coords(self)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   coords = self%get_atomcoords()
   coords(1, :) = -coords(1, :)
   call self%set_atomcoords(coords)

end subroutine

subroutine matrix_rotate_coords(self, rotmat)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_atomcoords(matrix_rotated(self%natom, self%get_atomcoords(), rotmat))

end subroutine

subroutine quater_rotate_coords(self, rotquat)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: rotquat(4)

   call self%set_atomcoords(quater_rotated(self%natom, self%get_atomcoords(), rotquat))

end subroutine

subroutine translate_coords(self, travec)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: travec(3)

   call self%set_atomcoords(translated(self%natom, self%get_atomcoords(), travec))

end subroutine

subroutine permutate_atoms(self, atomorder)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i, k
   integer, allocatable :: atomunorder(:)

   allocate(atomunorder(self%natom))

   atomunorder = inverse_mapping(atomorder)
   self%atoms = self%atoms(atomorder(:))
   do i = 1, self%natom
      do k = 1, size(self%atoms(i)%adjlist)
         self%atoms(i)%adjlist(k) = atomunorder(self%atoms(i)%adjlist(k))
      end do
   end do

end subroutine

subroutine set_atomelnums(self, atomelnums)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomelnums(:)

   self%atoms(:)%elnum = atomelnums(:)

end subroutine

function get_inner_atomelnums(self) result(atomelnums)
   class(Molecule), intent(in) :: self
   integer, allocatable :: atomelnums(:)

   atomelnums = self%atoms(:)%elnum

end function

function get_sorted_atomelnums(self, atomorder) result(atomelnums)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomelnums(:)

   atomelnums = self%atoms(atomorder)%elnum

end function

subroutine set_atomlabels(self, atomlabels)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomlabels(:)

   self%atoms(:)%label = atomlabels(:)

end subroutine

function get_inner_atomlabels(self) result(atomlabels)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomlabels(:)

   atomlabels = self%atoms(:)%label

end function

function get_sorted_atomlabels(self, atomorder) result(atomlabels)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomlabels(:)

   atomlabels = self%atoms(atomorder)%label

end function

subroutine set_atomequivs(self, natomequiv, atomequividcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: natomequiv
   integer, intent(in) :: atomequividcs(:)

   self%atoms(:)%equividx = atomequividcs(:)
   call self%atomequivpartition%set_partition(natomequiv, atomequividcs)

end subroutine

subroutine set_molfrags(self, nfrag, fragidcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nfrag
   integer, intent(in) :: fragidcs(:)

   self%atoms(:)%fragidx = fragidcs(:)
   call self%molfragpartition%set_partition(nfrag, fragidcs)

end subroutine

subroutine set_partition(self, npart, partidcs)
   class(Partition), intent(inout) :: self
   integer, intent(in) :: npart
   integer, intent(in) :: partidcs(:)
   ! Local variables
   integer :: i
   integer, allocatable :: n(:)

   allocate(n(npart))
   allocate(self%parts(npart))

   n(:) = 0
   do i = 1, size(partidcs)
      n(partidcs(i)) = n(partidcs(i)) + 1
   end do

   do i = 1, npart
      allocate(self%parts(i)%atomidcs(n(i)))
   end do

   n(:) = 0
   do i = 1, size(partidcs)
      n(partidcs(i)) = n(partidcs(i)) + 1
      self%parts(partidcs(i))%atomidcs(n(partidcs(i))) = i
   end do

end subroutine

function get_lenlist(self) result(lenlist)
   class(Partition), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: lenlist(:)

   allocate(lenlist(size(self%parts)))

   do i = 1, size(self%parts)
      lenlist(i) = size(self%parts(i)%atomidcs)
   end do

end function

function get_atomequivlenlist(self) result(atomequivlenlist)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomequivlenlist(:)

   atomequivlenlist = self%atomequivpartition%get_lenlist()

end function

subroutine set_atomtypes(self, natomtype, atomtypeidcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: natomtype
   integer, intent(in) :: atomtypeidcs(:)

   self%atoms(:)%typeidx = atomtypeidcs(:)
   call self%atomtypepartition%set_partition(natomtype, atomtypeidcs)

end subroutine

function get_atomtypelenlist(self) result(atomtypelenlist)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomtypelenlist(:)

   atomtypelenlist = self%atomtypepartition%get_lenlist()

end function

function get_inner_atomtypeparts(self) result(parts)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   type(Part), allocatable :: parts(:)

   allocate(parts(size(self%atomtypepartition%parts)))

   do i = 1, size(self%atomtypepartition%parts)
      parts(i)%atomidcs = self%atomtypepartition%parts(i)%atomidcs(:)
   end do

end function

function get_sorted_atomtypeparts(self, atomorder) result(parts)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: atomunorder(:)
   type(Part), allocatable :: parts(:)

   atomunorder = inverse_mapping(atomorder)
   allocate(parts(size(self%atomtypepartition%parts)))

   do i = 1, size(self%atomtypepartition%parts)
      parts(i)%atomidcs = atomunorder(self%atomtypepartition%parts(i)%atomidcs(:))
   end do

end function

function get_inner_fragparts(self) result(parts)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   type(Part), allocatable :: parts(:)

   allocate(parts(size(self%molfragpartition%parts)))

   do i = 1, size(self%molfragpartition%parts)
      parts(i)%atomidcs = self%molfragpartition%parts(i)%atomidcs(:)
   end do

end function

function get_sorted_fragparts(self, atomorder) result(parts)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: atomunorder(:)
   type(Part), allocatable :: parts(:)

   atomunorder = inverse_mapping(atomorder)
   allocate(parts(size(self%molfragpartition%parts)))

   do i = 1, size(self%molfragpartition%parts)
      parts(i)%atomidcs = atomunorder(self%molfragpartition%parts(i)%atomidcs(:))
   end do

end function

function get_fragroots(self, atomorder) result(fragroots)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: order(:)
   integer, allocatable :: fragroots(:)
   integer, allocatable :: atomequividcs(:)
   integer, allocatable :: atomequivlenlist(:)
   type(Part), allocatable :: fragparts(:)

   fragparts = self%get_fragparts(atomorder)
   atomequividcs = self%get_atomequividcs(atomorder)
   atomequivlenlist = self%get_atomequivlenlist()
   allocate(fragroots(size(fragparts)))

   do i = 1, size(fragparts)
      order = sorted_order(atomequivlenlist(atomequividcs(fragparts(i)%atomidcs)))
      fragroots(i) = fragparts(i)%atomidcs(order(1))
   end do

end function

function get_inner_atomtypeidcs(self) result(atomtypeidcs)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomtypeidcs(:)

   atomtypeidcs = self%atoms(:)%typeidx

end function

function get_sorted_atomtypeidcs(self, atomorder) result(atomtypeidcs)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomtypeidcs(:)

   atomtypeidcs = self%atoms(atomorder)%typeidx

end function

function get_inner_atomequividcs(self) result(atomequividcs)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomequividcs(:)

   atomequividcs = self%atoms(:)%equividx

end function

function get_sorted_atomequividcs(self, atomorder) result(atomequividcs)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomequividcs(:)

   atomequividcs = self%atoms(atomorder)%equividx

end function

subroutine set_weights(self, weights)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: weights(self%natom)

   self%atoms(:)%weight = weights(:)

end subroutine

function get_inner_atomweights(self) result(weights)
   class(Molecule), intent(in) :: self
   ! Local variables
   real(wp), allocatable :: weights(:)

   weights = self%atoms(:)%weight

end function

function get_sorted_atomweights(self, atomorder) result(weights)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   real(wp), allocatable :: weights(:)

   weights = self%atoms(atomorder)%weight

end function

subroutine set_atomcoords(self, coords)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: coords(3, self%natom)
   ! Local variables
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine

function get_inner_atomcoords(self) result(coords)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords(:)
   end do

end function

function get_sorted_atomcoords(self, atomorder) result(coords)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   do i = 1, self%natom
      coords(:, i) = self%atoms(atomorder(i))%coords(:)
   end do

end function

function get_inner_adjmat(self) result(adjmat)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i, k
   type(Atom) :: iatom
   logical, allocatable :: adjmat(:, :)

   allocate(adjmat(self%natom, self%natom))

   adjmat(:, :) = .false.

   do i = 1, self%natom
      iatom = self%atoms(i)
      do k = 1, size(iatom%adjlist)
         adjmat(i, iatom%adjlist(k)) = .true.
      end do
   end do

end function

function get_sorted_adjmat(self, atomorder) result(adjmat)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i, k
   type(Atom) :: atom_i
   integer, allocatable :: atomunorder(:)
   logical, allocatable :: adjmat(:, :)

   atomunorder = inverse_mapping(atomorder)
   allocate(adjmat(self%natom, self%natom))
   adjmat(:, :) = .false.

   do i = 1, self%natom
      atom_i = self%atoms(atomorder(i))
      do k = 1, size(atom_i%adjlist)
         adjmat(i, atomunorder(atom_i%adjlist(k))) = .true.
      end do
   end do

end function

function get_inner_nadjs(self) result(nadjs)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: nadjs(:)

   allocate(nadjs(size(self%atoms)))
   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(i)%adjlist)
   end do

end function

function get_sorted_nadjs(self, atomorder) result(nadjs)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: nadjs(:)

   allocate(nadjs(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(atomorder(i))%adjlist)
   end do

end function

subroutine set_adjlists(self, nadjs, adjlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nadjs(:)
   integer, intent(in) :: adjlists(:, :)
   ! Local variables
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%adjlist = adjlists(:nadjs(i), i)
   end do

end subroutine

function get_inner_adjlists(self) result(adjlists)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: adjlists(:, :)

   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      adjlists(:size(self%atoms(i)%adjlist), i) = self%atoms(i)%adjlist
   end do

end function

function get_sorted_adjlists(self, atomorder) result(adjlists)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   type(Atom) :: atom_i
   integer, allocatable :: atomunorder(:)
   integer, allocatable :: adjlists(:, :)

   atomunorder = inverse_mapping(atomorder)
   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      atom_i = self%atoms(atomorder(i))
      adjlists(:size(atom_i%adjlist), i) = atomunorder(atom_i%adjlist)
   end do

end function

subroutine set_adjequivlenlists(self, nadjequivs, adjequivlenlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nadjequivs(:)
   integer, intent(in) :: adjequivlenlists(:, :)
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      if (.not. allocated(self%atoms(i)%adjequivlenlist)) then
         allocate(self%atoms(i)%adjequivlenlist(nadjequivs(i)))
      end if
      self%atoms(i)%adjequivlenlist = adjequivlenlists(:nadjequivs(i), i)
   end do

end subroutine

function get_inner_adjequivlenlists(self) result(adjequivlenlists)
   class(Molecule), intent(inout) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: adjequivlenlists(:, :)

   allocate(adjequivlenlists(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      adjequivlenlists(:size(self%atoms(i)%adjequivlenlist), i) = self%atoms(i)%adjequivlenlist
   end do

end function

function get_sorted_adjequivlenlists(self, atomorder) result(adjequivlenlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   type(Atom) :: atom_i
   integer, allocatable :: adjequivlenlists(:, :)

   allocate(adjequivlenlists(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      atom_i = self%atoms(atomorder(i))
      adjequivlenlists(:size(atom_i%adjequivlenlist), i) = atom_i%adjequivlenlist
   end do

end function

function get_inner_nadjequivs(self) result(nadjequivs)
   class(Molecule), intent(inout) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: nadjequivs(:)

   allocate(nadjequivs(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadjequivs(i) = size(self%atoms(i)%adjequivlenlist)
   end do

end function

function get_sorted_nadjequivs(self, atomorder) result(nadjequivs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   type(Atom) :: atom_i
   integer, allocatable :: nadjequivs(:)

   allocate(nadjequivs(size(self%atoms)))

   do i = 1, size(self%atoms)
      atom_i = self%atoms(atomorder(i))
      nadjequivs(i) = size(atom_i%adjequivlenlist)
   end do

end function

subroutine print_atom(self, ind, outLvl)
   class(Atom), intent(in) :: self
   integer, intent(in) :: ind
   integer, intent(in), optional :: outLvl
   ! Local variables
   character(255) :: frmt
   character(2) :: num
   integer :: outLevel

! *** code to manage unitlbl pending ***

   outLevel = 1
   if (present(outLvl)) outLevel = outLvl
   
   select case (outLevel)

   case (1)
      if (size(self%adjlist) == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%typeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') size(self%adjlist)
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%typeidx, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:size(self%adjlist)), " ]"
      end if

!  case (2)
!
   case default
      if (size(self%adjlist) == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%typeidx, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') size(self%adjlist)
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%typeidx, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:size(self%adjlist)), " ]"
      end if
   end select

end subroutine

subroutine print_molecule(self)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i

! *** code to manage unitlbl pending ***
   
   write (stderr, '(a,i0,a)') "Contents of molecule structure:   (", &
                                         self%natom, " atoms)"
   write (stderr, '(2a)') 'Title: ', self%title
   write (stderr, '(a,a4,a5,a12,a7,2a14)') "ind:", "lbl", "elnum", "typeidx", &
                                          "weight", "{ coords }", "[ adjlist ]"

   do i = 1, self%natom
      call self%atoms(i)%print(i)
   end do

end subroutine

function bonded(self, idx1, idx2) result(isbond)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i
   logical :: isbond, found1, found2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate(adjlist1(maxcoord), adjlist2(maxcoord))

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

end function

subroutine remove_bond(self, idx1, idx2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i, pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate(adjlist1(maxcoord), adjlist2(maxcoord))

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

end subroutine

subroutine add_bond(self, idx1, idx2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate(adjlist1(maxcoord), adjlist2(maxcoord))

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
!      size(self%atoms(idx1)%adjlist) = size(self%atoms(idx1)%adjlist) + 1
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

end subroutine

end module
