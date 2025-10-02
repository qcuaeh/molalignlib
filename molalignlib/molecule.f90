module molecule
use parameters
use chemistry
implicit none
private
public set_coords
public include_all_atoms
public include_heavy_atoms
public adjacency_from_bonds
public adjacency_from_distance
public get_coords
public get_mirrored_coords
public get_weighted_coords
public get_centroid
public print_atoms
public print_bonds
!public get_adjmat
!public add_bond
!public remove_bond

type, public :: atom_t
   integer :: elnum
   integer :: typeid
   real(rk) :: coords(3)
end type

type, public :: bond_t
   integer :: typeid
   integer :: atomidx1
   integer :: atomidx2
end type

type, public :: adjc_t
   integer, allocatable :: adjlist(:)
end type

interface get_weighted_coords
   module procedure get_weighted_coords_base
   module procedure get_weighted_coords_center
end interface

contains

subroutine include_all_atoms(atoms, atomset, subset_alloc)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer, dimension(:), pointer, intent(out) :: atomset
   integer, dimension(:), allocatable, target, intent(out) :: subset_alloc
   ! Local variables
   integer :: i

   allocate (subset_alloc(size(atoms)))
   atomset => subset_alloc

   do i = 1, size(atoms)
      subset_alloc(i) = i
   end do
end subroutine

subroutine include_heavy_atoms(atoms, atomset, subset_alloc)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer, dimension(:), pointer, intent(out) :: atomset
   integer, dimension(:), allocatable, target, intent(out) :: subset_alloc
   ! Local variables
   integer :: atomidx, n

   allocate (subset_alloc(size(atoms)))

   n = 0
   do atomidx = 1, size(atoms)
      if (atoms(n)%elnum > 1) then
         n = n + 1
         subset_alloc(n) = atomidx
      end if
   end do

   atomset => subset_alloc(1:n)
end subroutine

subroutine set_coords(atoms, coords)
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(rk), dimension(:,:), intent(in) :: coords
   ! Local variables
   integer :: i

   do i = 1, size(atoms)
      atoms(i)%coords = coords(:, i)
   end do
end subroutine

subroutine adjacency_from_bonds(atomset, atoms, bonds, adjcs)
   integer, dimension(:), intent(in) :: atomset
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   ! Local variables
   integer, allocatable ::  nadjs(:)
   integer, allocatable :: adjlist(:,:)
   integer :: i, idx1, idx2

   allocate (adjcs(size(atoms)))
   allocate (nadjs(size(atoms)))
   allocate (adjlist(MAX_COORD, size(atoms)))

   nadjs = 0
   do i = 1, size(bonds)
      idx1 = bonds(i)%atomidx1
      idx2 = bonds(i)%atomidx2
      if (any(atomset == idx1) .and. any(atomset == idx2)) then
         nadjs(idx1) = nadjs(idx1) + 1
         nadjs(idx2) = nadjs(idx2) + 1
         adjlist(nadjs(idx1), idx1) = idx2
         adjlist(nadjs(idx2), idx2) = idx1
      end if
   end do

   do i = 1, size(adjcs)
      adjcs(i)%adjlist = adjlist(1:nadjs(i), i)
   end do
end subroutine

subroutine adjacency_from_distance(atomset, atoms, adjcs)
   integer, dimension(:), intent(in) :: atomset
   type(atom_t), dimension(:), intent(in) :: atoms
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   ! Local variables
   integer, allocatable :: nadjs(:)
   integer, allocatable :: adjlist(:,:)
   integer :: i, j
   real(rk), allocatable :: atom_radii(:)
   real(rk) :: atom_dist

   allocate (adjcs(size(atoms)))
   allocate (nadjs(size(atoms)))
   allocate (adjlist(MAX_COORD, size(atoms)))

   ! Set atom radii
   atom_radii = 0.75*covalent_radii(atoms%elnum) + 0.25*vdw_radii(atoms%elnum)

   ! Register adjacency matrix i,j if atoms i and j are closer
   ! than the sum of their adjacency radius
   nadjs = 0
   do i = 1, size(atoms)
      if (any(atomset == i)) then
         do j = i + 1, size(atoms)
            if (any(atomset == j)) then
               atom_dist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
               if (atom_dist < atom_radii(i) + atom_radii(j)) then
                  nadjs(i) = nadjs(i) + 1
                  nadjs(j) = nadjs(j) + 1
                  if (nadjs(i) > MAX_COORD .or. nadjs(j) > MAX_COORD) then
                     write (stderr, '(A)') 'Error: Maximum coordination number exceeded'
                     stop 1
                  end if
                  adjlist(nadjs(i), i) = j
                  adjlist(nadjs(j), j) = i
               end if
            end if
         end do
      end if
   end do

   do i = 1, size(adjcs)
      adjcs(i)%adjlist = adjlist(1:nadjs(i), i)
   end do
end subroutine

function get_coords(atoms) result(coords)
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   integer :: i

   allocate (coords(3, size(atoms)))
   do i = 1, size(atoms)
      coords(:, i) = atoms(i)%coords
   end do
end function

function get_mirrored_coords(atoms) result(coords)
! Reflect coordinates on the YZ plane
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   integer :: i

   allocate (coords(3, size(atoms)))
   do i = 1, size(atoms)
      coords(1,i) = -atoms(i)%coords(1)
      coords(2:3,i) = atoms(i)%coords(2:3)
   end do
end function

function get_weighted_coords_base(atoms, weights) result(coords)
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:), intent(in) :: weights
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   real(rk) :: total_weight
   integer :: i

   allocate (coords(3, size(atoms)))
   total_weight = sum(weights)

   do i = 1, size(atoms)
      coords(:, i) = sqrt(weights(i)/total_weight)*atoms(i)%coords
   end do
end function

function get_weighted_coords_center(atoms, weights, center) result(coords)
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:), intent(in) :: weights
   real(rk), intent(in) :: center(3)
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   real(rk) :: total_weight
   integer :: i

   allocate (coords(3, size(atoms)))
   total_weight = sum(weights)

   do i = 1, size(atoms)
      coords(:, i) = sqrt(weights(i)/total_weight)*(atoms(i)%coords - center)
   end do
end function

function get_centroid(atomset, atoms, weights) result(centroid)
! Calculate the coordinates of the center of mass
   integer, dimension(:), intent(in) :: atomset
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:), intent(in) :: weights
   ! Local variables
   real(rk) :: centroid(3)
   real(rk) :: total_weight, total_coords(3)
   integer :: i

   total_weight = 0
   total_coords = 0
   do i = 1, size(atoms)
      if (any(atomset == i)) then
         total_weight = total_weight + weights(i)
         total_coords = total_coords + weights(i)*atoms(i)%coords
      end if
   end do
   centroid = total_coords/total_weight
end function

!function get_adjmat(atoms) result(adjmat)
!   type(atom_t), dimension(:), intent(in) :: atoms
!   ! Local variables
!   logical, dimension(:,:), allocatable :: adjmat
!   integer :: i, k, num_atoms
!
!   num_atoms = size(atoms)
!   allocate (adjmat(num_atoms, num_atoms))
!   adjmat(:, :) = .false.
!
!   do i = 1, num_atoms
!      do k = 1, size(atoms(i)%adjlist)
!         adjmat(i, atoms(i)%adjlist(k)) = .true.
!      end do
!   end do
!end function
!
!subroutine remove_bond(atoms, idx1, idx2)
!   type(atom_t), dimension(:), intent(inout) :: atoms
!   integer, intent(in) :: idx1, idx2
!   ! Local variables
!   integer :: i, pos1, pos2, nadj1, nadj2
!   integer, dimension(:), allocatable :: adjlist1, adjlist2
!
!   allocate (adjlist1(MAX_COORD), adjlist2(MAX_COORD))
!
!! copy adjlist arrays
!   nadj1 = size(atoms(idx1)%adjlist)
!   nadj2 = size(atoms(idx2)%adjlist)
!   adjlist1 = atoms(idx1)%adjlist
!   adjlist2 = atoms(idx2)%adjlist
!
!! initialization
!   pos1 = 0   ! position of idx2 in adjlist of atom 1
!   pos2 = 0   ! position of idx1 in adjlist of atom 2
!
!! find position of idx2 and idx1 in adjlist of atoms idx1 and idx2, resp.
!   do i = 1, nadj1
!      if (idx2 == adjlist1(i)) then
!         pos1 = i
!         exit
!      end if
!   end do
!   do i = 1, nadj2
!      if (idx1 == adjlist2(i)) then
!         pos2 = i
!         exit
!      end if
!   end do
!
!! delete idx2 and idx1 from the ajdlists where they appear
!   if ((pos1 /= 0) .and. (pos2 /= 0)) then
!      nadj1 = nadj1 - 1
!      do i = pos1, nadj1
!         adjlist1(i) = adjlist1(i+1)
!      end do
!      nadj2 = nadj2 - 1
!      do i = pos2, nadj2
!         adjlist2(i) = adjlist2(i+1)
!      end do
!! update neighbor arrays for atoms idx1 and idx2
!      atoms(idx1)%adjlist = adjlist1(:nadj1)
!      atoms(idx2)%adjlist = adjlist2(:nadj2)
!!   else
!!      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', idx1, idx2
!   end if
!end subroutine
!
!subroutine add_bond(atoms, idx1, idx2)
!   type(atom_t), dimension(:), intent(inout) :: atoms
!   integer, intent(in) :: idx1, idx2
!   ! Local variables
!   integer :: pos1, pos2, nadj1, nadj2
!   integer, dimension(:), allocatable :: adjlist1, adjlist2
!
!   allocate (adjlist1(MAX_COORD), adjlist2(MAX_COORD))
!
!! copy array of adjlist
!   nadj1 = size(atoms(idx1)%adjlist)
!   nadj2 = size(atoms(idx2)%adjlist)
!   adjlist1 = atoms(idx1)%adjlist
!   adjlist2 = atoms(idx2)%adjlist
!
!! initialization
!   pos1 = nadj1
!   pos2 = nadj2
!
!!   if (.not. adjmat(idx1, idx2)) then
!!      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", idx1, idx2
!!   end if
!
!! indices in adjlist are supposed to be sorted; inserting new indices
!!      size(atoms(idx1)%adjlist) = size(atoms(idx1)%adjlist) + 1
!   nadj1 = nadj1 + 1
!! find position to insert idx2 and shift indices greater than idx2
!   do while ((pos1 >= 1) .and. (idx2 < adjlist1(pos1)))
!      adjlist1(pos1+1) = adjlist1(pos1)
!      pos1 = pos1 - 1
!   end do
!   adjlist1(pos1+1) = idx2
!
!   nadj2 = nadj2 + 1
!! find position to insert idx1 and shift indices greater than idx1
!   do while ((pos2 >= 1) .and. (idx1 < adjlist2(pos2)))
!      adjlist2(pos2+1) = adjlist2(pos2)
!      pos2 = pos2 - 1
!   end do
!   adjlist2(pos2+1) = idx1
!! update neighbor arrays for atoms in idx1 and idx2
!   atoms(idx1)%adjlist = adjlist1(:nadj1)
!   atoms(idx2)%adjlist = adjlist2(:nadj2)
!end subroutine

subroutine print_atoms(atoms)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   integer :: i
   character(:), allocatable :: fmtstr
   type(atom_t) :: atom

   write (stderr, '(A,2X,A,1X,A,4X,A,8X,A,8X,A,4X,A)') "idx", "sym", "label", &
         "X","Y","Z"

   do i = 1, size(atoms)
      atom = atoms(i)
      fmtstr = '(I3,3X,A2,1X,I3,3(1X,f8.4),2X)'
      write (stderr, fmtstr) i, element_symbols(atom%elnum), atom%typeid, atom%coords
   end do
end subroutine

subroutine print_bonds(bonds)
   type(bond_t), dimension(:), intent(in) :: bonds
   ! Local variables
   integer :: i

   write (stderr, '(a)') "idx1 idx2"

   do i = 1, size(bonds)
      write (stderr, '(I3,2X,I3)') bonds(i)%atomidx1, bonds(i)%atomidx2 
   end do
end subroutine

end module
