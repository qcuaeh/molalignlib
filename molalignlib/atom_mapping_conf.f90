! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

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

module atom_mapping_conf
use parameters
use globals
use random
use molecule
use strutils
use chemdata
use permutation
use rigid_body
use rotation
use alignment
use lap_driver
use lap_solvers
use adjacency
use biasing
use pruning
use printing
use lcrs_tree
use partitioning
!use backtracking

implicit none

contains

subroutine remap_conformations(mol1, mol2, eltypes, results)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_container), intent(in) :: eltypes
   type(registry_type), target, intent(out) :: results

   ! Local variables
   type(tree_node), pointer :: mnatypes_tree
   integer, dimension(:), allocatable :: atomperm, auxperm
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: rmsd, dist
   real(rk) :: eigquat(4), totquat(4)
   integer :: num_atoms1, num_trials, num_steps
   integer, pointer :: lead_count

   num_atoms1 = size(mol1%atoms)
   coords1 = mol1%get_weightcoords()
   coords2 = mol2%get_weightcoords()
   call results%initialize(max_records)

   allocate (atomperm(num_atoms1))
   allocate (auxperm(num_atoms1))

   ! Compute MNA types
   call tree_from_partition(eltypes, mnatypes_tree)
   call compute_consistent_mnatypes(mol1, mol2, mnatypes_tree)

   ! Reflect atoms
   if (mirror_flag) then
      call reflect_coords(coords2)
   end if

   ! Translate atoms to their centroids
   call translate_coords(coords1, -centroid(coords1))
   call translate_coords(coords2, -centroid(coords2))

   ! Initialize random number generator
   call random_initialize()

   ! Optimize atom permutation

   num_trials = 0
   lead_count => results%records(1)%count

   do while (lead_count < max_count .and. num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat())

      ! Assign atoms with current orientation
      call assign_atoms_conf(mnatypes_tree, mol1, mol2, coords1, coords2, atomperm, dist)
      totquat = leasteigquat(atomperm, coords1, coords2)
      call rotate_coords(coords2, totquat)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_conf(mnatypes_tree, mol1, mol2, coords1, coords2, auxperm, dist)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         eigquat = leasteigquat(atomperm, coords1, coords2)
         call rotate_coords(coords2, eigquat)
         totquat = quatmul(eigquat, totquat)
         num_steps = num_steps + 1
      end do

      ! Update results
      rmsd = sqrt(sqdistsum(atomperm, coords1, coords2))
      call results%push_rmsd(atomperm, num_steps, angle(totquat), rmsd)

   end do

   results%num_trials = num_trials

end subroutine

subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node), intent(in) :: mnatypes
   real(rk), dimension(:, :), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(tree_node), pointer :: submnatypes
   logical :: remaining

   write (stderr, *)
   write (stderr, '(A)') repeat('submnatypes ', 8)
   submnatypes = mnatypes
   call print_tree(submnatypes)

   do
      call assign_remaining_items(submnatypes, remaining)
      call compute_consistent_mnatypes(mol1, mol2, submnatypes)
      call print_tree(submnatypes)
      if (.not. remaining) exit
   end do
   stop

end subroutine

!subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   type(mol_type), intent(in) :: mol1, mol2
!   type(bipartition_container), intent(in) :: mnatypes
!   real(rk), dimension(:, :), intent(in) :: coords1, coords2
!   integer, dimension(:), intent(out) :: atomperm
!   real(rk), intent(out) :: dist
!   ! Local variables
!   type(bipartition_container) :: submnatypes
!
!   submnatypes = mnatypes
!   call assign_atoms_conf_rec(submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)
!
!end subroutine
!
!recursive subroutine assign_atoms_conf_rec( submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   type(mol_type), intent(in) :: mol1, mol2
!   type(bipartition_container), intent(inout) :: submnatypes
!   real(rk), dimension(:, :), intent(in) :: coords1, coords2
!   integer, dimension(:), intent(out) :: atomperm
!   real(rk), intent(out) :: dist
!   ! Local variables
!   type(metapartition_type) :: metatypes
!   integer :: h, i
!
!   call collect_mnatypes(mol1, submnatypes%partition1(), metatypes)
!!   call metatypes%print_parts()
!!   call submnatypes%print_parts()
!
!   do i = 1, metatypes%num_parts
!!      write (stderr, *) 'loop:', i
!      h = random_element(metatypes%parts(i)%items)
!      call minperm(submnatypes%parts(h), coords1, coords2, atomperm, dist)
!      call split_crossmnatypes(h, atomperm, submnatypes)
!      call compute_crossmnatypes(mol1, mol2, submnatypes)
!      call assign_atoms_conf_rec(submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   end do
!
!end subroutine
!
!subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   type(mol_type), intent(in) :: mol1, mol2
!   type(bipartition_container), intent(in) :: mnatypes
!   real(rk), dimension(:, :), intent(in) :: coords1, coords2
!   integer, dimension(:), intent(out) :: atomperm
!   real(rk), intent(out) :: dist
!   ! Local variables
!   type(bipartition_container) :: submnatypes
!   type(metapartition_type) :: metatypes
!   integer :: h, i, j, k
!
!   write (stderr, *) repeat('*', 80)
!   call mnatypes%print_parts()
!
!   call collect_mnatypes(mol1, mnatypes%partition1(), metatypes)
!   call metatypes%print_parts()
!   submnatypes = mnatypes
!
!   do while (metatypes%num_parts > 0)
!      write (stderr, *) repeat('+', 80)
!      do i = 1, metatypes%num_parts
!         h = random_element(metatypes%parts(i)%items)
!         do j = 1, metatypes%num_parts
!            do k = 1, metatypes%parts(j)%num_items
!               if (metatypes%parts(j)%items(k) > h) then
!                  metatypes%parts(j)%items(k) = metatypes%parts(j)%items(k) + submnatypes%parts(h)%num_items1 - 1
!               end if
!            end do
!         end do
!         call minperm(submnatypes%parts(h), coords1, coords2, atomperm, dist)
!         call split_crossmnatypes(h, atomperm, submnatypes)
!         write (stderr, *)
!         write (stderr, *) repeat(str(h)//'   ', 8)
!         call submnatypes%print_parts()
!         call metatypes%print_parts()
!      end do
!      write (stderr, *) repeat('-', 80)
!      call compute_crossmnatypes(mol1, mol2, submnatypes)
!      call submnatypes%print_parts()
!      call collect_mnatypes(mol1, submnatypes%partition1(), metatypes)
!      call metatypes%print_parts()
!   end do
!
!   call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)
!
!end subroutine

end module
