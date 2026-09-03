! MolAlignLib
! Copyright (C) 2025 José M. Vásquez, Carlos Z. Gómez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module molecule
use parameters
use chemdata
use adjacency
implicit none
private
public set_coords
public get_coords
public get_mirrored_coords
public get_weighted_coords
public get_centroid
public include_all_atoms
public include_heavy_atoms
public bonds_from_atoms
public adjacency_from_bonds
public print_atoms
public print_bonds
!public get_adjmat
!public add_bond
!public remove_bond

type, public :: atom_t
   integer(ik) :: elnum
   integer(ik) :: group
   real(rk) :: coords(3)
end type

type, public :: bond_t
   integer(ik) :: atomidx1
   integer(ik) :: atomidx2
   integer(ik) :: bondtype
end type

interface get_weighted_coords
   module procedure get_weighted_coords_base
   module procedure get_weighted_coords_center
end interface

contains

subroutine include_all_atoms(atoms, atomset)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: i

   allocate (atomset(size(atoms)))

   do i = 1, size(atoms)
      atomset(i) = i
   end do
end subroutine

subroutine include_heavy_atoms(atoms, atomset)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: nel, atomidx

   allocate (atomset(count(atoms%elnum > 1)))

   nel = 0
   do atomidx = 1, size(atoms)
      if (atoms(atomidx)%elnum > 1) then
         nel = nel + 1
         atomset(nel) = atomidx
      end if
   end do
end subroutine

subroutine set_coords(atoms, coords)
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(rk), dimension(:,:), intent(in) :: coords
   ! Local variables
   integer(ik) :: i

   do i = 1, size(atoms)
      atoms(i)%coords = coords(:, i)
   end do
end subroutine

subroutine bonds_from_atoms(atoms, bonds)
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   ! Local variables
   integer(ik) :: i, j, n_atoms, n_bonds
   real(rk), allocatable :: atom_radii(:)
   real(rk) :: atom_dist

   n_atoms = size(atoms)

   ! Set atom radii
   atom_radii = bond_tol*covalent_radii(atoms%elnum)

   ! First pass: count bonds
   n_bonds = 0
   do i = 1, n_atoms
      do j = i + 1, n_atoms
         atom_dist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
         if (atom_dist < atom_radii(i) + atom_radii(j)) then
            n_bonds = n_bonds + 1
         end if
      end do
   end do

   ! Allocate bonds array with exact size
   allocate(bonds(n_bonds))

   ! Second pass: populate bonds
   n_bonds = 0
   do i = 1, n_atoms
      do j = i + 1, n_atoms
         atom_dist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
         if (atom_dist < atom_radii(i) + atom_radii(j)) then
            n_bonds = n_bonds + 1
            bonds(n_bonds)%atomidx1 = i
            bonds(n_bonds)%atomidx2 = j
            bonds(n_bonds)%bondtype = 1
         end if
      end do
   end do
end subroutine

subroutine adjacency_from_bonds(atomset, atoms, bonds, adjcs)
   integer(ik), dimension(:), intent(in) :: atomset
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(adjc_t), dimension(:), allocatable, intent(out) :: adjcs
   ! Local variables
   logical(lk), dimension(:,:), allocatable :: adjmat
   integer(ik) :: n_atoms, atomidx1, atomidx2, i

   n_atoms = size(atoms)

   allocate (adjmat(n_atoms, n_atoms))
   adjmat = .FALSE.

   ! Build adjacency matrix from bonds
   do i = 1, size(bonds)
      atomidx1 = bonds(i)%atomidx1
      atomidx2 = bonds(i)%atomidx2
      adjmat(atomidx1, atomidx2) = .TRUE.
      adjmat(atomidx2, atomidx1) = .TRUE.
   end do

   ! Convert to adjacency lists
   call adjmat_to_adjcs(atomset, adjmat, adjcs)

   deallocate (adjmat)
end subroutine

function get_coords(atoms) result(coords)
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:,:), allocatable :: coords
   ! Local variables
   integer(ik) :: i

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
   integer(ik) :: i

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
   integer(ik) :: i

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
   integer(ik) :: i

   allocate (coords(3, size(atoms)))
   total_weight = sum(weights)

   do i = 1, size(atoms)
      coords(:, i) = sqrt(weights(i)/total_weight)*(atoms(i)%coords - center)
   end do
end function

function get_centroid(atomset, atoms, weights) result(centroid)
! Calculate the coordinates of the center of mass
   integer(ik), dimension(:), intent(in) :: atomset
   type(atom_t), dimension(:), intent(in) :: atoms
   real(rk), dimension(:), intent(in) :: weights
   ! Local variables
   real(rk) :: centroid(3)
   real(rk) :: total_weight, total_coords(3)
   integer(ik) :: atomidx, i

   total_weight = 0
   total_coords = 0
   do i = 1, size(atomset)
      atomidx = atomset(i)
      total_weight = total_weight + weights(atomidx)
      total_coords = total_coords + weights(atomidx)*atoms(atomidx)%coords
   end do
   centroid = total_coords/total_weight
end function

subroutine print_atoms(atoms)
   type(atom_t), dimension(:), intent(in) :: atoms
   ! Local variables
   integer(ik) :: i
   character(:), allocatable :: fmtstr
   type(atom_t) :: atom

   write (stderr, '(A,2X,A,1X,A,4X,A,8X,A,8X,A,4X,A)') "idx", "sym", "type", &
         "X","Y","Z"

   do i = 1, size(atoms)
      atom = atoms(i)
      fmtstr = '(I3,3X,A2,1X,I3,3(1X,f8.4),2X)'
      write (stderr, fmtstr) i, atomic_symbols(atom%elnum), atom%group, atom%coords
   end do
end subroutine

subroutine print_bonds(bonds)
   type(bond_t), dimension(:), intent(in) :: bonds
   ! Local variables
   integer(ik) :: i

   write (stderr, '(a)') "atomidx1 atomidx2"

   do i = 1, size(bonds)
      write (stderr, '(I3,2X,I3)') bonds(i)%atomidx1, bonds(i)%atomidx2
   end do
end subroutine

end module
