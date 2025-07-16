module molecule
use parameters
use chemistry
implicit none
private

public set_coords
public add_bond
public remove_bond
public bond_atoms
public get_coords
public get_adjmat
public get_bonds
public print_atoms
public print_bonds

type, public :: atom_t
   integer :: elnum
   integer :: label
   integer :: weight
   real(rk) :: coords(3)
   integer :: adjlist_allocation(MAX_COORD)
   integer, pointer :: adjlist(:)
end type

type, public :: bond_t
   integer :: idx1
   integer :: idx2
end type

contains

subroutine set_coords(atoms, coords)
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(rk), dimension(:,:), intent(in) :: coords
   ! Local variables
   integer :: i

   do i = 1, size(atoms)
      atoms(i)%coords = coords(:, i)
   end do
end subroutine

function get_coords(atoms) result(coords)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   real(rk), dimension(:,:), allocatable :: coords
   integer :: i

   allocate (coords(3, size(atoms)))

   do i = 1, size(atoms)
      coords(:, i) = atoms(i)%coords
   end do
end function

subroutine bond_atoms(atoms)
   type(atom_t), dimension(:), target, intent(inout) :: atoms
   ! Local variables
   integer :: i, j, num_atoms
   integer, dimension(:), allocatable :: nadjs
   real(rk), dimension(:), allocatable :: atom_radii
   real(rk) :: atom_dist

   num_atoms = size(atoms)
   allocate (nadjs(num_atoms))

   ! Set atom radii
   atom_radii = 0.75*covalent_radii(atoms%elnum) + 0.25*vdw_radii(atoms%elnum)

   ! Register adjacency matrix i,j if atoms i and j are closer
   ! than the sum of their adjacency radius
   nadjs = 0
   do i = 1, num_atoms
      do j = i + 1, num_atoms
         atom_dist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
         if (atom_dist < atom_radii(i) + atom_radii(j)) then
            nadjs(i) = nadjs(i) + 1
            nadjs(j) = nadjs(j) + 1
            if (nadjs(i) > MAX_COORD .or. nadjs(j) > MAX_COORD) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            atoms(i)%adjlist_allocation(nadjs(i)) = j
            atoms(j)%adjlist_allocation(nadjs(j)) = i
         end if
      end do
   end do

   do i = 1, num_atoms
      atoms(i)%adjlist => atoms(i)%adjlist_allocation(1:nadjs(i))
   end do
end subroutine

function get_adjmat(atoms) result(adjmat)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   logical, dimension(:,:), allocatable :: adjmat
   integer :: i, k, num_atoms

   num_atoms = size(atoms)
   allocate (adjmat(num_atoms, num_atoms))
   adjmat(:, :) = .false.

   do i = 1, num_atoms
      do k = 1, size(atoms(i)%adjlist)
         adjmat(i, atoms(i)%adjlist(k)) = .true.
      end do
   end do
end function

function get_bonds(atoms) result(bonds)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   logical, dimension(:,:), allocatable :: adjmat
   type(bond_t), dimension(:), allocatable :: bonds
   integer :: i, j, num_bonds

   adjmat = get_adjmat(atoms)
   allocate (bonds(count(adjmat)/2))

   num_bonds = 0
   do i = 1, size(atoms)
      do j = i + 1, size(atoms)
         if (adjmat(i, j)) then
            num_bonds = num_bonds + 1
            bonds(num_bonds)%idx1 = i
            bonds(num_bonds)%idx2 = j
         end if
      end do
   end do
end function

subroutine remove_bond(atoms, idx1, idx2)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i, pos1, pos2, nadj1, nadj2
   integer, dimension(:), allocatable :: adjlist1, adjlist2

   allocate (adjlist1(MAX_COORD), adjlist2(MAX_COORD))

! copy adjlist arrays
   nadj1 = size(atoms(idx1)%adjlist)
   nadj2 = size(atoms(idx2)%adjlist)
   adjlist1 = atoms(idx1)%adjlist
   adjlist2 = atoms(idx2)%adjlist

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
      atoms(idx1)%adjlist = adjlist1(:nadj1)
      atoms(idx2)%adjlist = adjlist2(:nadj2)
!   else
!      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', idx1, idx2
   end if
end subroutine

subroutine add_bond(atoms, idx1, idx2)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: pos1, pos2, nadj1, nadj2
   integer, dimension(:), allocatable :: adjlist1, adjlist2

   allocate (adjlist1(MAX_COORD), adjlist2(MAX_COORD))

! copy array of adjlist
   nadj1 = size(atoms(idx1)%adjlist)
   nadj2 = size(atoms(idx2)%adjlist)
   adjlist1 = atoms(idx1)%adjlist
   adjlist2 = atoms(idx2)%adjlist

! initialization
   pos1 = nadj1
   pos2 = nadj2

!   if (.not. adjmat(idx1, idx2)) then
!      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", idx1, idx2
!   end if

! indices in adjlist are supposed to be sorted; inserting new indices
!      size(atoms(idx1)%adjlist) = size(atoms(idx1)%adjlist) + 1
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
   atoms(idx1)%adjlist = adjlist1(:nadj1)
   atoms(idx2)%adjlist = adjlist2(:nadj2)
end subroutine

subroutine print_atoms(atoms)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   integer :: i
   character(:), allocatable :: fmtstr
   type(atom_t) :: atom

   write (stderr, '(A,2X,A,1X,A,4X,A,8X,A,8X,A,4X,A)') "idx", "sym", "label", &
         "X","Y","Z", "adjlist"

   do i = 1, size(atoms)
      atom = atoms(i)
      fmtstr = '(I3,3X,A2,1X,I3,3(1X,f8.4),2X,"["' // &
            repeat(',1X,I0', size(atom%adjlist)) // ',1X,"]")'
      write (stderr, fmtstr) i, element_symbols(atom%elnum), atom%label, &
            atom%coords, atom%adjlist
   end do
end subroutine

subroutine print_bonds(atoms)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   integer :: i
   type(bond_t), dimension(:), allocatable :: bonds

   bonds = get_bonds(atoms)

   write (stderr, '(a)') "idx1 idx2"

   do i = 1, size(bonds)
      write (stderr, '(I3,2X,I3)') bonds(i)%idx1, bonds(i)%idx2 
   end do
end subroutine

end module
