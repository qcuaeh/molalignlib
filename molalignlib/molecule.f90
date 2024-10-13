module molecule
use kinds
use bounds
use setutils
use partition
use permutation
use translation
use rotation
use adjacency
use alignment
use strutils
use chemdata

implicit none
private

type, public :: atom_type
   integer :: elnum
   integer :: label
   real(rk) :: weight
   real(rk) :: coords(3)
   integer, allocatable :: adjlist(:)
   integer, allocatable :: adjmnatypepartlens(:)
end type

type, public :: bond_type
   integer :: atomidx1
   integer :: atomidx2
end type

type, public :: molecule_type
   integer :: natom
   character(:), allocatable :: title
   type(atom_type), allocatable :: atoms(:)
   type(partition_type) :: eltypes
   type(partition_type) :: mnatypes
   type(partition_type) :: molfrags
contains
   procedure :: set_elnums
   procedure :: set_labels
   procedure :: set_coords
   procedure :: set_weights
   procedure :: set_adjlists
   procedure :: set_eltypes
   procedure :: set_mnatypes
!   procedure :: set_molfrags
   procedure :: set_adjmnatypepartlens
   procedure :: get_natom
   procedure :: get_atoms
   procedure :: get_elnums
   procedure :: get_labels
   procedure :: get_title
   procedure :: get_bonds
   procedure :: get_adjlists
   procedure :: get_adjpartitions
   procedure :: get_eltypes
   procedure :: get_mnatypes
   procedure :: get_molfrags
   procedure :: get_weights
   procedure :: get_coords
   procedure :: get_molfragroots
   procedure :: get_adjmatrix
   procedure :: get_nadjmnatypes
   procedure :: get_adjmnatypepartlens
   procedure :: permutate_atoms
   procedure :: mirror_coords
   procedure :: translate_coords
   procedure :: rotate_coords
   procedure :: bonded
   procedure :: add_bond
   procedure :: remove_bond
   procedure :: print_atoms
   procedure :: print_bonds
end type

contains

function get_title(self) result(title)
   class(molecule_type), intent(in) :: self
   character(:), allocatable :: title

   title = self%title

end function

function get_natom(self) result(natom)
   class(molecule_type), intent(in) :: self
   integer :: natom

   natom = size(self%atoms)

end function

function get_atoms(self) result(atoms)
   class(molecule_type), intent(in) :: self
   type(atom_type), allocatable :: atoms(:)

   atoms = self%atoms

end function

subroutine mirror_coords(self)
   class(molecule_type), intent(inout) :: self
   ! Local variables
   real(rk), allocatable :: coords(:, :)

   allocate (coords(3, size(self%atoms)))

   coords = self%get_coords()
   coords(1, :) = -coords(1, :)
   call self%set_coords(coords)

end subroutine

subroutine rotate_coords(self, rotquat)
   class(molecule_type), intent(inout) :: self
   ! Local variables
   real(rk), intent(in) :: rotquat(4)

   call self%set_coords(rotated(size(self%atoms), self%get_coords(), rotquat))

end subroutine

subroutine translate_coords(self, travec)
   class(molecule_type), intent(inout) :: self
   ! Local variables
   real(rk), intent(in) :: travec(3)

   call self%set_coords(translated(size(self%atoms), self%get_coords(), travec))

end subroutine

subroutine permutate_atoms(self, atom_order)
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i, k
   integer, allocatable :: atom_mapping(:)

   allocate (atom_mapping(size(self%atoms)))

   atom_mapping = inverse_permutation(atom_order)
   self%atoms = self%atoms(atom_order)
   do i = 1, size(self%atoms)
      do k = 1, size(self%atoms(i)%adjlist)
         self%atoms(i)%adjlist(k) = atom_mapping(self%atoms(i)%adjlist(k))
      end do
   end do

end subroutine

subroutine set_elnums(self, elnums)
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: elnums(:)

   self%atoms%elnum = elnums

end subroutine

function get_elnums(self) result(elnums)
   class(molecule_type), intent(in) :: self
   integer, allocatable :: elnums(:)

   elnums = self%atoms%elnum

end function

subroutine set_labels(self, labels)
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: labels(:)

   self%atoms%label = labels

end subroutine

function get_labels(self) result(labels)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer, allocatable :: labels(:)

   labels = self%atoms%label

end function

subroutine set_eltypes(self, eltypes)
   class(molecule_type), intent(inout) :: self
   type(partition_type), intent(in) :: eltypes

   self%eltypes = eltypes

end subroutine

function get_eltypes(self) result(eltypes)
   class(molecule_type), intent(in) :: self
   type(partition_type) :: eltypes

   eltypes = self%eltypes

end function

subroutine set_mnatypes(self, mnatypes)
   class(molecule_type), intent(inout) :: self
   type(partition_type), intent(in) :: mnatypes

   self%mnatypes = mnatypes

end subroutine

function get_mnatypes(self) result(mnatypes)
   class(molecule_type), intent(in) :: self
   type(partition_type) :: mnatypes

   mnatypes = self%mnatypes

end function

!subroutine set_molfrags(self, nfrag, fragidcs)
!   class(molecule_type), intent(inout) :: self
!   integer, intent(in) :: nfrag
!   integer, intent(in) :: fragidcs(:)
!
!   self%molfrags = atompartition(nfrag, fragidcs)
!
!end subroutine

function get_molfrags(self) result(parts)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   type(indexlist_type), allocatable :: parts(:)

   allocate (parts(size(self%molfrags%parts)))

   do i = 1, size(self%molfrags%parts)
      parts(i)%indices = self%molfrags%parts(i)%indices
   end do

end function

function get_molfragroots(self) result(fragroots)
   class(molecule_type), intent(in) :: self
   ! Result variable
   integer, allocatable :: fragroots(:)
   ! Local variables
   integer :: h, i
   integer :: atomidx_i, atomidx_min
   integer :: eltypepop_i, eltypepop_min
   integer, allocatable :: atomeltypes(:)

   allocate (fragroots(size(self%molfrags%parts)))

   atomeltypes = self%eltypes%get_item_types()
   do h = 1, size(self%molfrags%parts)
      eltypepop_min = huge(eltypepop_min)
      do i = 1, size(self%molfrags%parts(h)%indices)
         atomidx_i = self%molfrags%parts(h)%indices(i)
         eltypepop_i = size(self%eltypes%parts(atomeltypes(atomidx_i))%indices)
         if (eltypepop_i < eltypepop_min) then
            atomidx_min = atomidx_i
            eltypepop_min = eltypepop_i
         end if
         fragroots(h) = atomidx_min
      end do
   end do

end function

subroutine set_weights(self, weights)
   class(molecule_type), intent(inout) :: self
   real(rk), intent(in) :: weights(size(self%atoms))

   self%atoms%weight = weights

end subroutine

function get_weights(self) result(weights)
   class(molecule_type), intent(in) :: self
   ! Local variables
   real(rk), allocatable :: weights(:)

   weights = self%atoms%weight

end function

subroutine set_coords(self, coords)
   class(molecule_type), intent(inout) :: self
   real(rk), intent(in) :: coords(3, size(self%atoms))
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine

function get_coords(self) result(coords)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   real(rk), allocatable :: coords(:, :)

   allocate (coords(3, size(self%atoms)))

   do i = 1, size(self%atoms)
      coords(:, i) = self%atoms(i)%coords
   end do

end function

subroutine set_adjlists(self, nadjs, adjlists)
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: nadjs(:)
   integer, intent(in) :: adjlists(:, :)
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      self%atoms(i)%adjlist = adjlists(:nadjs(i), i)
   end do

end subroutine

function get_adjlists(self) result(adjlists)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   type(indexlist_type), allocatable :: adjlists(:)

   allocate (adjlists(size(self%atoms)))

   do i = 1, size(self%atoms)
      adjlists(i)%indices = self%atoms(i)%adjlist
   end do

end function

function get_adjpartitions(self, partition) result(adjpartitions)
   class(molecule_type), intent(in) :: self
   type(partition_type), intent(in) :: partition
   ! Result variable
   type(partition_type), allocatable :: adjpartitions(:)
   ! Local variables
   integer :: h, i
   type(atom_type) :: atom
   integer, allocatable :: adjlistsubset(:)

   do i = 1, size(self%atoms)
      atom = self%atoms(i)
      allocate (adjpartitions(i)%parts(size(partition%parts)))
      do h = 1, size(partition%parts)
         adjlistsubset = intersection(atom%adjlist, partition%parts(h)%indices, size(self%atoms))
         allocate (adjpartitions(i)%parts(h)%indices(size(adjlistsubset)))
         adjpartitions(i)%parts(h)%indices = adjlistsubset
      end do
   end do

end function

function get_adjmatrix(self) result(adjmat)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i, k
   type(atom_type) :: atom
   logical, allocatable :: adjmat(:, :)

   allocate (adjmat(size(self%atoms), size(self%atoms)))

   adjmat(:, :) = .false.

   do i = 1, size(self%atoms)
      atom = self%atoms(i)
      do k = 1, size(atom%adjlist)
         adjmat(i, atom%adjlist(k)) = .true.
      end do
   end do

end function

subroutine set_adjmnatypepartlens(self, nadjmnatypes, adjmnatypepartlens)
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: nadjmnatypes(:)
   integer, intent(in) :: adjmnatypepartlens(:, :)
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      if (.not. allocated(self%atoms(i)%adjmnatypepartlens)) then
         allocate (self%atoms(i)%adjmnatypepartlens(nadjmnatypes(i)))
      end if
      self%atoms(i)%adjmnatypepartlens = adjmnatypepartlens(:nadjmnatypes(i), i)
   end do

end subroutine

function get_adjmnatypepartlens(self) result(adjmnatypepartlens)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: adjmnatypepartlens(:, :)

   allocate (adjmnatypepartlens(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      adjmnatypepartlens(:size(self%atoms(i)%adjmnatypepartlens), i) = self%atoms(i)%adjmnatypepartlens
   end do

end function

function get_nadjmnatypes(self) result(nadjmnatypes)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: nadjmnatypes(:)

   allocate (nadjmnatypes(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadjmnatypes(i) = size(self%atoms(i)%adjmnatypepartlens)
   end do

end function

function get_bonds(self) result(bonds)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i, j, nbond
   logical, allocatable :: adjmat(:, :)
   type(bond_type), allocatable :: bonds(:)

   adjmat = self%get_adjmatrix()
   allocate (bonds(count(adjmat)/2))

   nbond = 0
   do i = 1, size(self%atoms)
      do j = i + 1, size(self%atoms)
         if (adjmat(i, j)) then
            nbond = nbond + 1
            bonds(nbond)%atomidx1 = i
            bonds(nbond)%atomidx2 = j
         end if
      end do
   end do

end function

subroutine print_atoms(self)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   character(ll) :: fmtstr
   type(atom_type) :: atom

   write (stderr, '(a,1x,i0)') 'Atoms:', size(self%atoms)
   write (stderr, '(a,a4,a12,a7,2a14)') "ind:", "elsym", "weight", &
         "{ coords }", "[ adjlist ]"

   do i = 1, size(self%atoms)
      atom = self%atoms(i)
      fmtstr = '(i3,": ",a2,1x,f6.2," {",3(1x,f8.4)," } ["' // &
            repeat(',1x,i3', size(atom%adjlist)) // '," ]")'
      write (stderr, fmtstr) i, elsyms(atom%elnum), atom%weight, &
            atom%coords, atom%adjlist
   end do

end subroutine

subroutine print_bonds(self)
   class(molecule_type), intent(in) :: self
   ! Local variables
   integer :: i
   type(bond_type), allocatable :: bonds(:)

   bonds = self%get_bonds()

   write (stderr, '(a,1x,i0)') 'Bonds:', size(bonds)
   write (stderr, '(a)') "idx1 idx2"

   do i = 1, size(bonds)
      write (stderr, '(i3,2x,i3)') bonds(i)%atomidx1, bonds(i)%atomidx2 
   end do

end subroutine

function bonded(self, idx1, idx2) result(isbond)
   class(molecule_type), intent(in) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i
   logical :: isbond, found1, found2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate (adjlist1(maxcoord), adjlist2(maxcoord))

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
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i, pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate (adjlist1(maxcoord), adjlist2(maxcoord))

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
!   else
!      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', idx1, idx2
   end if

end subroutine

subroutine add_bond(self, idx1, idx2)
   class(molecule_type), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate (adjlist1(maxcoord), adjlist2(maxcoord))

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
