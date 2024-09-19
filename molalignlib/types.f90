module types
use kinds
use bounds
use discrete
use rotation
use translation
use adjacency
use alignment
use strutils
use chemdata
implicit none
private

type, public :: cPart
   integer, allocatable :: atomidcs(:)
end type

type, public :: cPartition
   type(cPart), allocatable :: parts(:)
   integer, allocatable :: atom_mapping(:)
contains
   procedure :: init => partition_init
   procedure :: get_lenlist => partition_lenlist
end type

type, public :: cBond
   integer :: atomidx1
   integer :: atomidx2
end type

type, public :: cAtom
   integer :: elnum
   integer :: label
   integer :: eltype
   integer :: mnatype
   integer :: fragidx
   integer, allocatable :: mnaid(:)
   real(wp) :: weight
   real(wp) :: coords(3)
   integer, allocatable :: adjlist(:)
   integer, allocatable :: atomneimnatypepartlens(:)
end type

type, public :: cMol
   integer :: natom
   integer :: nfrag
   character(:), allocatable :: title
   type(cAtom), allocatable :: atoms(:)
   type(cPartition) :: eltypepartition
   type(cPartition) :: mnatypepartition
   type(cPartition) :: molfragpartition
contains
   procedure :: print_atoms
   procedure :: print_bonds
   procedure :: set_atomelnums
   procedure :: set_atomlabels
   procedure :: set_atomcoords
   procedure :: set_weights
   procedure :: set_adjlists
   procedure :: set_eltypes
   procedure :: set_mnatypes
   procedure :: set_molfrags
   procedure :: set_atomneimnatypepartlens
   procedure :: get_title
   procedure :: get_neltype
   procedure :: get_nmnatype
   procedure :: get_eltypepartlens
   procedure :: get_mnatypepartlens
   procedure :: mirror_atomcoords
   procedure :: translate_atomcoords
   procedure :: permutate_atoms
   procedure :: get_bonds
   procedure :: bonded
   procedure :: add_bond
   procedure :: remove_bond
   procedure :: get_atoms
   procedure :: get_sorted_atoms
   procedure :: get_mnatypeparts
   procedure :: get_molfragparts
   procedure :: get_sorted_molfragparts
   procedure :: get_sorted_fragroots
   procedure, private :: get_default_nadjs
   procedure, private :: get_default_adjlists
   procedure, private :: get_default_adjmat
   procedure, private :: get_default_atomcoords
   procedure, private :: get_default_atomeltypes
   procedure, private :: get_default_atomweights
   procedure, private :: get_default_natomneimnatypes
   procedure, private :: get_default_atomneimnatypepartlens
   procedure, private :: get_default_atommnatypes
   procedure, private :: get_sorted_nadjs
   procedure, private :: get_sorted_adjlists
   procedure, private :: get_sorted_adjmat
   procedure, private :: get_sorted_atomcoords
   procedure, private :: get_sorted_atomeltypes
   procedure, private :: get_sorted_atomweights
   procedure, private :: get_sorted_natomneimnatypes
   procedure, private :: get_sorted_atomneimnatypepartlens
   procedure, private :: get_sorted_atommnatypes
   procedure, private :: matrix_rotate_atomcoords
   procedure, private :: quater_rotate_atomcoords
   generic :: get_nadjs => get_default_nadjs, get_sorted_nadjs
   generic :: get_adjlists => get_default_adjlists, get_sorted_adjlists
   generic :: get_adjmat => get_default_adjmat, get_sorted_adjmat
   generic :: get_atomcoords => get_default_atomcoords, get_sorted_atomcoords
   generic :: get_atomeltypes => get_default_atomeltypes, get_sorted_atomeltypes
   generic :: get_atommnatypes => get_default_atommnatypes, get_sorted_atommnatypes
   generic :: get_atomneimnatypepartlens => get_default_atomneimnatypepartlens, get_sorted_atomneimnatypepartlens
   generic :: get_atomweights => get_default_atomweights, get_sorted_atomweights
   generic :: get_natomneimnatypes => get_default_natomneimnatypes, get_sorted_natomneimnatypes
   generic :: rotate_atomcoords => matrix_rotate_atomcoords, quater_rotate_atomcoords
end type

contains

function get_title(self) result(title)
   class(cMol), intent(in) :: self
   character(:), allocatable :: title

   title = self%title

end function

integer function get_neltype(self) result(neltype)
   class(cMol), intent(in) :: self

   neltype = size(self%eltypepartition%parts)

end function

integer function get_nmnatype(self) result(nmnatype)
   class(cMol), intent(in) :: self

   nmnatype = size(self%mnatypepartition%parts)

end function

subroutine mirror_atomcoords(self)
   class(cMol), intent(inout) :: self
   ! Local variables
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   coords = self%get_atomcoords()
   coords(1, :) = -coords(1, :)
   call self%set_atomcoords(coords)

end subroutine

subroutine matrix_rotate_atomcoords(self, rotmat)
   class(cMol), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_atomcoords(matrix_rotated(self%natom, self%get_atomcoords(), rotmat))

end subroutine

subroutine quater_rotate_atomcoords(self, rotquat)
   class(cMol), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: rotquat(4)

   call self%set_atomcoords(quater_rotated(self%natom, self%get_atomcoords(), rotquat))

end subroutine

subroutine translate_atomcoords(self, travec)
   class(cMol), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: travec(3)

   call self%set_atomcoords(translated(self%natom, self%get_atomcoords(), travec))

end subroutine

subroutine permutate_atoms(self, atom_order)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i, k
   integer, allocatable :: atom_mapping(:)

   allocate(atom_mapping(self%natom))

   atom_mapping = inverse_mapping(atom_order)
   self%atoms = self%atoms(atom_order(:))
   do i = 1, self%natom
      do k = 1, size(self%atoms(i)%adjlist)
         self%atoms(i)%adjlist(k) = atom_mapping(self%atoms(i)%adjlist(k))
      end do
   end do

end subroutine

subroutine set_atomelnums(self, atomelnums)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: atomelnums(:)

   self%atoms(:)%elnum = atomelnums(:)

end subroutine

function get_atoms(self) result(atoms)
   class(cMol), intent(in) :: self
   ! Local variables
   type(cAtom), allocatable :: atoms(:)

   atoms = self%atoms

end function

function get_sorted_atoms(self) result(atoms)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: p, offset, partsize
   integer, allocatable :: atom_order(:)
   type(cAtom), allocatable :: atoms(:)

   allocate(atoms(size(self%atoms)))

   offset = 0
   do p = 1, size(self%mnatypepartition%parts)
      partsize = size(self%mnatypepartition%parts(p)%atomidcs)
      atoms(offset+1:offset+partsize) = self%atoms(self%mnatypepartition%parts(p)%atomidcs)
      offset = offset + partsize
   end do

end function

subroutine set_atomlabels(self, atomlabels)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: atomlabels(:)

   self%atoms(:)%label = atomlabels(:)

end subroutine

subroutine set_eltypes(self, neltype, atomeltypes)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: neltype
   integer, intent(in) :: atomeltypes(:)

   self%atoms(:)%eltype = atomeltypes(:)
   call self%eltypepartition%init(neltype, atomeltypes)

end subroutine

subroutine set_mnatypes(self, nmnatype, atommnatypes)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: nmnatype
   integer, intent(in) :: atommnatypes(:)

   self%atoms(:)%mnatype = atommnatypes(:)
   call self%mnatypepartition%init(nmnatype, atommnatypes)

end subroutine

subroutine set_molfrags(self, nfrag, fragidcs)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: nfrag
   integer, intent(in) :: fragidcs(:)

   self%atoms(:)%fragidx = fragidcs(:)
   call self%molfragpartition%init(nfrag, fragidcs)

end subroutine

subroutine partition_init(self, npart, partidcs)
   class(cPartition), intent(out) :: self
   integer, intent(in) :: npart
   integer, intent(in) :: partidcs(:)
   ! Local variables
   integer :: i, p, offset
   integer, allocatable :: n(:)

   allocate(n(npart))
   allocate(self%parts(npart))
   allocate(self%atom_mapping(size(partidcs)))

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

   offset = 0
   do p = 1, size(self%parts)
      do i = 1, size(self%parts(p)%atomidcs)
         self%atom_mapping(self%parts(p)%atomidcs(i)) = offset + i
      end do
      offset = offset + size(self%parts(p)%atomidcs)
   end do

end subroutine

function partition_lenlist(self) result(lenlist)
   class(cPartition), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: lenlist(:)

   allocate(lenlist(size(self%parts)))

   do i = 1, size(self%parts)
      lenlist(i) = size(self%parts(i)%atomidcs)
   end do

end function

function get_mnatypeparts(self) result(parts)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   type(cPart), allocatable :: parts(:)

   allocate(parts(size(self%mnatypepartition%parts)))

   do i = 1, size(self%mnatypepartition%parts)
      parts(i)%atomidcs = self%mnatypepartition%parts(i)%atomidcs(:)
   end do

end function

function get_mnatypepartlens(self) result(mnatypepartlens)
   class(cMol), intent(in) :: self
   ! Local variables
   integer, allocatable :: mnatypepartlens(:)

   mnatypepartlens = self%mnatypepartition%get_lenlist()

end function

function get_eltypepartlens(self) result(eltypepartlens)
   class(cMol), intent(in) :: self
   ! Local variables
   integer, allocatable :: eltypepartlens(:)

   eltypepartlens = self%eltypepartition%get_lenlist()

end function

function get_molfragparts(self) result(parts)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   type(cPart), allocatable :: parts(:)

   allocate(parts(size(self%molfragpartition%parts)))

   do i = 1, size(self%molfragpartition%parts)
      parts(i)%atomidcs = self%molfragpartition%parts(i)%atomidcs(:)
   end do

end function

function get_sorted_molfragparts(self) result(parts)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   type(cPart), allocatable :: parts(:)

   allocate(parts(size(self%molfragpartition%parts)))

   do i = 1, size(self%molfragpartition%parts)
      parts(i)%atomidcs = self%mnatypepartition%atom_mapping(self%molfragpartition%parts(i)%atomidcs(:))
   end do

end function

function get_sorted_fragroots(self) result(fragroots)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable, dimension(:) :: order, order1, order2
   integer, allocatable :: fragroots(:)
   integer, allocatable :: atommnatypes(:)
   integer, allocatable :: eltypepartlens(:)
   integer, allocatable :: mnatypepartlens(:)
   type(cPart), allocatable :: fragparts(:)
   type(cAtom), allocatable :: atoms(:)

   atoms = self%get_sorted_atoms()
   fragparts = self%get_sorted_molfragparts()
   eltypepartlens = self%get_eltypepartlens()
   mnatypepartlens = self%get_mnatypepartlens()
   allocate(fragroots(size(fragparts)))

   do i = 1, size(fragparts)
      order1 = sorted_order(mnatypepartlens(atoms(fragparts(i)%atomidcs)%mnatype))
      order2 = sorted_order(eltypepartlens(atoms(fragparts(i)%atomidcs(order1))%eltype))
      order = order1(order2)
      fragroots(i) = fragparts(i)%atomidcs(order(1))
   end do

end function

function get_default_atomeltypes(self) result(atomeltypes)
   class(cMol), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomeltypes(:)

   atomeltypes = self%atoms(:)%eltype

end function

function get_sorted_atomeltypes(self, atom_order) result(atomeltypes)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer, allocatable :: atomeltypes(:)

   atomeltypes = self%atoms(atom_order)%eltype

end function

function get_default_atommnatypes(self) result(atommnatypes)
   class(cMol), intent(in) :: self
   ! Local variables
   integer, allocatable :: atommnatypes(:)

   atommnatypes = self%atoms(:)%mnatype

end function

function get_sorted_atommnatypes(self, atom_order) result(atommnatypes)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer, allocatable :: atommnatypes(:)

   atommnatypes = self%atoms(atom_order)%mnatype

end function

subroutine set_weights(self, weights)
   class(cMol), intent(inout) :: self
   real(wp), intent(in) :: weights(self%natom)

   self%atoms(:)%weight = weights(:)

end subroutine

function get_default_atomweights(self) result(weights)
   class(cMol), intent(in) :: self
   ! Local variables
   real(wp), allocatable :: weights(:)

   weights = self%atoms(:)%weight

end function

function get_sorted_atomweights(self, atom_order) result(weights)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   real(wp), allocatable :: weights(:)

   weights = self%atoms(atom_order)%weight

end function

subroutine set_atomcoords(self, coords)
   class(cMol), intent(inout) :: self
   real(wp), intent(in) :: coords(3, self%natom)
   ! Local variables
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine

function get_default_atomcoords(self) result(coords)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords(:)
   end do

end function

function get_sorted_atomcoords(self, atom_order) result(coords)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   do i = 1, self%natom
      coords(:, i) = self%atoms(atom_order(i))%coords(:)
   end do

end function

function get_default_adjmat(self) result(adjmat)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i, k
   type(cAtom) :: atom
   logical, allocatable :: adjmat(:, :)

   allocate(adjmat(self%natom, self%natom))

   adjmat(:, :) = .false.

   do i = 1, self%natom
      atom = self%atoms(i)
      do k = 1, size(atom%adjlist)
         adjmat(i, atom%adjlist(k)) = .true.
      end do
   end do

end function

function get_sorted_adjmat(self, atom_order) result(adjmat)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i, k
   type(cAtom) :: atom
   integer, allocatable :: atom_mapping(:)
   logical, allocatable :: adjmat(:, :)

   atom_mapping = inverse_mapping(atom_order)
   allocate(adjmat(self%natom, self%natom))
   adjmat(:, :) = .false.

   do i = 1, self%natom
      atom = self%atoms(atom_order(i))
      do k = 1, size(atom%adjlist)
         adjmat(i, atom_mapping(atom%adjlist(k))) = .true.
      end do
   end do

end function

function get_default_nadjs(self) result(nadjs)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: nadjs(:)

   allocate(nadjs(size(self%atoms)))
   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(i)%adjlist)
   end do

end function

function get_sorted_nadjs(self, atom_order) result(nadjs)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i
   integer, allocatable :: nadjs(:)

   allocate(nadjs(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(atom_order(i))%adjlist)
   end do

end function

subroutine set_adjlists(self, nadjs, adjlists)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: nadjs(:)
   integer, intent(in) :: adjlists(:, :)
   ! Local variables
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%adjlist = adjlists(:nadjs(i), i)
   end do

end subroutine

function get_default_adjlists(self) result(adjlists)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: adjlists(:, :)

   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      adjlists(:size(self%atoms(i)%adjlist), i) = self%atoms(i)%adjlist
   end do

end function

function get_sorted_adjlists(self, atom_order) result(adjlists)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i
   type(cAtom) :: atom
   integer, allocatable :: atom_mapping(:)
   integer, allocatable :: adjlists(:, :)

   atom_mapping = inverse_mapping(atom_order)
   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      atom = self%atoms(atom_order(i))
      adjlists(:size(atom%adjlist), i) = atom_mapping(atom%adjlist)
   end do

end function

subroutine set_atomneimnatypepartlens(self, natomneimnatypes, atomneimnatypepartlens)
   class(cMol), intent(inout) :: self
   integer, intent(in) :: natomneimnatypes(:)
   integer, intent(in) :: atomneimnatypepartlens(:, :)
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      if (.not. allocated(self%atoms(i)%atomneimnatypepartlens)) then
         allocate(self%atoms(i)%atomneimnatypepartlens(natomneimnatypes(i)))
      end if
      self%atoms(i)%atomneimnatypepartlens = atomneimnatypepartlens(:natomneimnatypes(i), i)
   end do

end subroutine

function get_default_atomneimnatypepartlens(self) result(atomneimnatypepartlens)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: atomneimnatypepartlens(:, :)

   allocate(atomneimnatypepartlens(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      atomneimnatypepartlens(:size(self%atoms(i)%atomneimnatypepartlens), i) = self%atoms(i)%atomneimnatypepartlens
   end do

end function

function get_sorted_atomneimnatypepartlens(self, atom_order) result(atomneimnatypepartlens)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i
   type(cAtom) :: atom
   integer, allocatable :: atomneimnatypepartlens(:, :)

   allocate(atomneimnatypepartlens(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      atom = self%atoms(atom_order(i))
      atomneimnatypepartlens(:size(atom%atomneimnatypepartlens), i) = atom%atomneimnatypepartlens
   end do

end function

function get_default_natomneimnatypes(self) result(natomneimnatypes)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: natomneimnatypes(:)

   allocate(natomneimnatypes(size(self%atoms)))

   do i = 1, size(self%atoms)
      natomneimnatypes(i) = size(self%atoms(i)%atomneimnatypepartlens)
   end do

end function

function get_sorted_natomneimnatypes(self, atom_order) result(natomneimnatypes)
   class(cMol), intent(in) :: self
   integer, intent(in) :: atom_order(:)
   ! Local variables
   integer :: i
   type(cAtom) :: atom
   integer, allocatable :: natomneimnatypes(:)

   allocate(natomneimnatypes(size(self%atoms)))

   do i = 1, size(self%atoms)
      atom = self%atoms(atom_order(i))
      natomneimnatypes(i) = size(atom%atomneimnatypepartlens)
   end do

end function

function get_bonds(self) result(bonds)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i, j, nbond
   logical, allocatable :: adjmat(:, :)
   type(cBond), allocatable :: bonds(:)
   type(cBond), allocatable :: bonds_buffer(:)

   adjmat = self%get_adjmat()
   allocate(bonds_buffer(size(self%atoms)*maxcoord))

   nbond = 0
   do i = 1, size(self%atoms)
      do j = i + 1, size(self%atoms)
         if (adjmat(i, j)) then
            nbond = nbond + 1
            bonds_buffer(nbond)%atomidx1 = i
            bonds_buffer(nbond)%atomidx2 = j
         end if
      end do
   end do

   bonds = bonds_buffer(:nbond)

end function

subroutine print_atoms(self)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   character(ll) :: fmtstr
   type(cAtom) :: atom

   write (stderr, '(a,1x,i0)') 'Atoms:', size(self%atoms)
   write (stderr, '(a,a4,a5,a12,a7,2a14)') "ind:", "elsym", "eltype", "weight", &
         "{ coords }", "[ adjlist ]"

   do i = 1, size(self%atoms)
      atom = self%atoms(i)
      fmtstr = '(i3,": ",a2,1x,i3,1x,f6.2," {",3(1x,f8.4)," } ["' // &
            repeat(',1x,i3', size(atom%adjlist)) // '," ]")'
      write (stderr, fmtstr) i, elsym(atom%elnum), atom%eltype, atom%weight, &
            atom%coords, atom%adjlist
   end do

end subroutine

subroutine print_bonds(self)
   class(cMol), intent(in) :: self
   ! Local variables
   integer :: i
   type(cBond), allocatable :: bonds(:)

   bonds = self%get_bonds()

   write (stderr, '(a,1x,i0)') 'Bonds:', size(bonds)
   write (stderr, '(a)') "idx1 idx2"

   do i = 1, size(bonds)
      write (stderr, '(i3,2x,i3)') bonds(i)%atomidx1, bonds(i)%atomidx2 
   end do

end subroutine

function bonded(self, idx1, idx2) result(isbond)
   class(cMol), intent(in) :: self
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
   class(cMol), intent(inout) :: self
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
!   else
!      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', idx1, idx2
   end if

end subroutine

subroutine add_bond(self, idx1, idx2)
   class(cMol), intent(inout) :: self
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
