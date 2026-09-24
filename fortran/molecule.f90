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
public pad_atoms
public complete_atomperm
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
   integer(ik) :: nel, atomidx

   ! Dummy atoms (elnum = 0) are never included
   allocate (atomset(count(atoms%elnum > 0)))

   nel = 0
   do atomidx = 1, size(atoms)
      if (atoms(atomidx)%elnum > 0) then
         nel = nel + 1
         atomset(nel) = atomidx
      end if
   end do
end subroutine

subroutine include_heavy_atoms(atoms, atomset)
   type(atom_t), dimension(:), intent(inout) :: atoms
   integer(ik), dimension(:), allocatable, intent(out) :: atomset
   ! Local variables
   integer(ik) :: nel, atomidx

   ! elnum > 1 excludes both hydrogens and dummy atoms (elnum = 0)
   allocate (atomset(count(atoms%elnum > 1)))

   nel = 0
   do atomidx = 1, size(atoms)
      if (atoms(atomidx)%elnum > 1) then
         nel = nel + 1
         atomset(nel) = atomidx
      end if
   end do
end subroutine

subroutine pad_atoms(atoms, n_pad)
! Append dummy atoms (elnum = 0) until the molecule has n_pad atoms. The
! original atoms keep their indices, so bond tables and file line numbers
! stay valid. Padding atoms are recognised by index (i > original size)
! everywhere else.
   type(atom_t), dimension(:), allocatable, intent(inout) :: atoms
   integer(ik), intent(in) :: n_pad
   ! Local variables
   type(atom_t), dimension(:), allocatable :: padded
   integer(ik) :: n_real, i

   n_real = size(atoms)
   if (n_real >= n_pad) return

   allocate (padded(n_pad))
   padded(1:n_real) = atoms
   do i = n_real + 1, n_pad
      padded(i)%elnum = 0
      padded(i)%group = 0
      padded(i)%coords = 0.0_rk
   end do

   call move_alloc(padded, atoms)
end subroutine

subroutine neighbor_table(n_atoms, bonds, nbr_start, nbr_list)
! Compressed neighbour lists: the neighbours of atom i are
! nbr_list(nbr_start(i) : nbr_start(i+1)-1).
   integer(ik), intent(in) :: n_atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(ik), dimension(:), allocatable, intent(out) :: nbr_start, nbr_list
   ! Local variables
   integer(ik), dimension(:), allocatable :: fill
   integer(ik) :: i, a1, a2

   allocate (nbr_start(n_atoms + 1))
   allocate (fill(n_atoms))
   fill = 0

   do i = 1, size(bonds)
      a1 = bonds(i)%atomidx1
      a2 = bonds(i)%atomidx2
      fill(a1) = fill(a1) + 1
      fill(a2) = fill(a2) + 1
   end do

   nbr_start(1) = 1
   do i = 1, n_atoms
      nbr_start(i + 1) = nbr_start(i) + fill(i)
   end do

   allocate (nbr_list(nbr_start(n_atoms + 1) - 1))
   fill = nbr_start(1:n_atoms)

   do i = 1, size(bonds)
      a1 = bonds(i)%atomidx1
      a2 = bonds(i)%atomidx2
      nbr_list(fill(a1)) = a2
      fill(a1) = fill(a1) + 1
      nbr_list(fill(a2)) = a1
      fill(a2) = fill(a2) + 1
   end do
end subroutine

subroutine complete_atomperm(atomset1, atoms1, atoms2, n_real1, n_real2, &
      bonds1, bonds2, coords1, coords2, atomperm1)
! Turn a partial atom mapping (defined only on atomset1) into a full
! permutation of 1..size(atoms1). Both molecules must already be padded to
! the same size. Excluded atoms are paired in this order of preference:
!   1. Same element, bonded to the image of one of its included neighbours
!      (e.g. an H follows its heavy atom), closest first.
!   2. Same element, closest remaining real atom.
!   3. Whatever is left, in ascending index order. Because padding atoms are
!      appended at the end, real atoms are used up before padding atoms.
! coords1 and coords2 must be in the same (aligned) frame.
   integer(ik), dimension(:), intent(in) :: atomset1
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   integer(ik), intent(in) :: n_real1, n_real2
   type(bond_t), dimension(:), intent(in) :: bonds1, bonds2
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer(ik), dimension(:), intent(inout) :: atomperm1
   ! Local variables
   logical(lk), dimension(:), allocatable :: in_set1, used2
   integer(ik), dimension(:), allocatable :: nbr1_start, nbr1_list, nbr2_start, nbr2_list
   integer(ik) :: n_atoms, i, j, k, p, q, h, best
   real(rk) :: dist, best_dist

   n_atoms = size(atoms1)

   if (size(atoms2) /= n_atoms .or. size(atomperm1) /= n_atoms) then
      error stop 'complete_atomperm: molecules and permutation must have the padded size'
   end if

   allocate (in_set1(n_atoms), used2(n_atoms))
   in_set1 = .FALSE.
   used2 = .FALSE.

   ! Validate the partial mapping: every included atom must be assigned,
   ! and no target may be used twice
   do i = 1, size(atomset1)
      j = atomperm1(atomset1(i))
      if (j < 1 .or. j > n_atoms) then
         error stop 'complete_atomperm: included atom was left unassigned'
      end if
      if (used2(j)) then
         error stop 'complete_atomperm: atom of molecule 2 assigned twice'
      end if
      in_set1(atomset1(i)) = .TRUE.
      used2(j) = .TRUE.
   end do

   ! Clear excluded slots so the result does not depend on how the
   ! permutation was initialised
   do i = 1, n_atoms
      if (.not. in_set1(i)) atomperm1(i) = 0
   end do

   ! Pass 1: follow bonds from included neighbours
   call neighbor_table(n_atoms, bonds1, nbr1_start, nbr1_list)
   call neighbor_table(n_atoms, bonds2, nbr2_start, nbr2_list)

   do i = 1, n_real1
      if (atomperm1(i) /= 0) cycle
      best = 0
      best_dist = huge(best_dist)
      do p = nbr1_start(i), nbr1_start(i + 1) - 1
         h = nbr1_list(p)
         if (.not. in_set1(h)) cycle
         j = atomperm1(h)
         do q = nbr2_start(j), nbr2_start(j + 1) - 1
            k = nbr2_list(q)
            if (k > n_real2) cycle
            if (used2(k)) cycle
            if (atoms2(k)%elnum /= atoms1(i)%elnum) cycle
            dist = sum((coords1(:, i) - coords2(:, k))**2)
            if (dist < best_dist) then
               best_dist = dist
               best = k
            end if
         end do
      end do
      if (best > 0) then
         atomperm1(i) = best
         used2(best) = .TRUE.
      end if
   end do

   ! Pass 2: closest remaining real atom of the same element
   do i = 1, n_real1
      if (atomperm1(i) /= 0) cycle
      best = 0
      best_dist = huge(best_dist)
      do k = 1, n_real2
         if (used2(k)) cycle
         if (atoms2(k)%elnum /= atoms1(i)%elnum) cycle
         dist = sum((coords1(:, i) - coords2(:, k))**2)
         if (dist < best_dist) then
            best_dist = dist
            best = k
         end if
      end do
      if (best > 0) then
         atomperm1(i) = best
         used2(best) = .TRUE.
      end if
   end do

   ! Pass 3: fill the remaining slots in ascending order. The counts of
   ! free slots and free targets are equal, so k never runs past n_atoms.
   k = 0
   do i = 1, n_atoms
      if (atomperm1(i) /= 0) cycle
      do
         k = k + 1
         if (.not. used2(k)) exit
      end do
      atomperm1(i) = k
      used2(k) = .TRUE.
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
   real(rk) :: atom_dist
   real(rk), dimension(:), allocatable :: atom_radii
   logical(lk), dimension(:,:), allocatable :: is_bonded

   n_atoms = size(atoms)
   allocate (is_bonded(n_atoms, n_atoms))
   is_bonded = .FALSE.

   ! Set atom radii
   atom_radii = covalent_radii(atoms%elnum)

   ! Single pass: compute distances once, cache result in matrix
   do i = 1, n_atoms
      do j = i + 1, n_atoms
         atom_dist = sqrt(sum((atoms(i)%coords - atoms(j)%coords)**2))
         is_bonded(i, j) = atom_dist < atom_radii(i) + atom_radii(j) + bond_tol
      end do
   end do

   ! Count bonds from cached matrix (cheap, no distance calc)
   n_bonds = count(is_bonded)

   ! Allocate bonds array with exact size
   allocate (bonds(n_bonds))

   ! Populate bonds from cached matrix
   n_bonds = 0
   do i = 1, n_atoms
      do j = i + 1, n_atoms
         if (is_bonded(i, j)) then
            n_bonds = n_bonds + 1
            bonds(n_bonds)%atomidx1 = i
            bonds(n_bonds)%atomidx2 = j
            bonds(n_bonds)%bondtype = 1
         end if
      end do
   end do

   deallocate (is_bonded)
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
