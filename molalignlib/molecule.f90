module molecule
use parameters
use chemdata

implicit none
private

type, public :: atom_type
   integer :: elnum
   integer :: label
   integer :: weight
   real(rk) :: coords(3)
   integer, allocatable :: adjlist(:)
end type

type, public :: bond_type
   integer :: atomidx1
   integer :: atomidx2
end type

type, public :: mol_type
   character(:), allocatable :: title
   type(atom_type), allocatable :: atoms(:)
end type

public set_coords
public set_adjlists
public get_bonds
public get_coords
public get_adjmat
public add_bond
public remove_bond
public print_atoms
public print_bonds

contains

subroutine set_coords(self, coords)
   class(mol_type), intent(inout) :: self
   real(rk), intent(in) :: coords(:, :)
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      self%atoms(i)%coords = coords(:, i)
   end do
end subroutine

subroutine set_adjlists(self, nadjs, adjlists)
   class(mol_type), intent(inout) :: self
   integer, intent(in) :: nadjs(:)
   integer, intent(in) :: adjlists(:, :)
   ! Local variables
   integer :: i, num_atoms

   num_atoms = size(self%atoms)

   do i = 1, num_atoms
      self%atoms(i)%adjlist = adjlists(:nadjs(i), i)
   end do
end subroutine

function get_coords(self) result(coords)
   class(mol_type), intent(in) :: self
   ! Local variables
   real(rk), allocatable :: coords(:, :)
   integer :: i

   allocate (coords(3, size(self%atoms)))

   do i = 1, size(self%atoms)
      coords(:, i) = self%atoms(i)%coords
   end do
end function

function get_adjmat(self) result(adjmat)
   class(mol_type), intent(in) :: self
   ! Local variables
   logical, allocatable :: adjmat(:,:)
   integer :: i, k, num_atoms

   num_atoms = size(self%atoms)
   allocate (adjmat(num_atoms, num_atoms))
   adjmat(:, :) = .false.

   do i = 1, num_atoms
      do k = 1, size(self%atoms(i)%adjlist)
         adjmat(i, self%atoms(i)%adjlist(k)) = .true.
      end do
   end do
end function

function get_bonds(self) result(bonds)
   class(mol_type), intent(in) :: self
   ! Local variables
   logical, allocatable :: adjmat(:,:)
   type(bond_type), allocatable :: bonds(:)
   integer :: i, j, nbond

   adjmat = get_adjmat(self)
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

subroutine remove_bond(self, idx1, idx2)
   class(mol_type), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i, pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate (adjlist1(max_coord), adjlist2(max_coord))

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
   class(mol_type), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate (adjlist1(max_coord), adjlist2(max_coord))

! copy array of adjlist
   nadj1 = size(self%atoms(idx1)%adjlist)
   nadj2 = size(self%atoms(idx2)%adjlist)
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

! initialization
   pos1 = nadj1
   pos2 = nadj2

!   if (.not. adjmat(idx1, idx2)) then
!      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", idx1, idx2
!   end if

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
end subroutine

subroutine print_atoms(self)
   class(mol_type), intent(in) :: self
   ! Local variables
   integer :: i
   character(:), allocatable :: fmtstr
   type(atom_type) :: atom

   write (stderr, '(a,1x,i0)') 'Atoms:', size(self%atoms)
   write (stderr, '(a,a4,a12,a7,2a14)') "ind:", "elsym", "label", &
         "{ coords }", "[ adjlist ]"

   do i = 1, size(self%atoms)
      atom = self%atoms(i)
      fmtstr = '(i3,": ",a2,1x,i3," {",3(1x,f8.4)," } ["' // &
            repeat(',1x,i3', size(atom%adjlist)) // '," ]")'
      write (stderr, fmtstr) i, element_symbols(atom%elnum), atom%label, &
            atom%coords, atom%adjlist
   end do
end subroutine

subroutine print_bonds(self)
   class(mol_type), intent(in) :: self
   ! Local variables
   integer :: i
   type(bond_type), allocatable :: bonds(:)

   bonds = get_bonds(self)

   write (stderr, '(a,1x,i0)') 'Bonds:', size(bonds)
   write (stderr, '(a)') "idx1 idx2"

   do i = 1, size(bonds)
      write (stderr, '(i3,2x,i3)') bonds(i)%atomidx1, bonds(i)%atomidx2 
   end do
end subroutine

end module
