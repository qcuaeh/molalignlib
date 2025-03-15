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

module remapping_bonded
use parameters
use globals
use random
use molecule
use strutils
use chemdata
use permutation
use spatial
use assignment
use adjacency
use biasing
use pruning
use lcrs_tree
use partitioning
use reactivity
use registry

implicit none

contains

subroutine remap_bonded_atoms(mol1, mol2, results)
   type(mol_type), intent(inout) :: mol1, mol2
   type(adjd_registry), target, intent(out) :: results

   ! Local variables
   type(tree_node), pointer :: eltree, mnatree
   type(bipartition_container) :: eltypes
   integer, dimension(:), allocatable :: atomperm, auxperm
   integer :: num_trials, num_steps
   integer, dimension(:), allocatable :: elnums1, elnums2
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   real(rk) :: step_rotation(4), total_rotation(4)
   real(rk) :: center1(3), center2(3)
   real(rk) :: dist

   ! Abort if molecules have different number of atoms
   if (size(mol1%atoms) /= size(mol2%atoms)) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Abort if molecules are not isomers
   if (any(sorted(mol1%atoms%elnum) /= sorted(mol2%atoms%elnum))) then
      write (stderr, '(a)') 'Error: These molecules are not isomers'
      stop
   end if

   ! Compute atomic types
   call compute_eltypes(mol1, mol2, eltree)
   call partition_from_tree(eltree, eltypes)

   ! Abort if there are conflicting atomic types
!   if (any(sorted(eltypes%itemdir1) /= sorted(eltypes%itemdir2))) then
   if (any(eltypes%parts%num_items1 /= eltypes%parts%num_items2)) then
      write (stderr, '(a)') 'Error: There are conflicting atomic types'
      stop
   end if

   allocate (atomperm(size(mol1%atoms)))
   allocate (auxperm(size(mol1%atoms)))

   elnums1 = mol1%atoms%elnum
   elnums2 = mol2%atoms%elnum
   coords1 = get_coords(mol1)
   coords2 = get_coords(mol2)
   adjmat1 = get_adjmat(mol1)
   adjmat2 = get_adjmat(mol2)

   if (reac_flag) then
      write (stderr, *) 'adjd before', adjacencydiff( adjmat1, adjmat2)
      call remove_reactive_bonds( mol1, mol2, eltypes, atomperm)
      adjmat1 = get_adjmat(mol1)
      adjmat2 = get_adjmat(mol2)
      write (stderr, *) 'adjd after', adjacencydiff( atomperm, adjmat1, adjmat2)
   end if

   ! Recompute MNA types
   call tree_from_partition( eltypes, mnatree)
   call compute_consistent_mnatypes( mol1, mol2, mnatree)
!   call print_tree( mnatree)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( coords2)
   end if

   ! Mass weight coordinates
   call weight_coords(coords1, atomic_weights(elnums1))
   call weight_coords(coords2, atomic_weights(elnums2))

   ! Translate atoms to their centroids
   center1 = centroid(coords1)
   center2 = centroid(coords2)
   call translate_coords(coords2, -center2)
   call translate_coords(coords2, center1)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call registry_init(results, max_records, coords1, adjmat1)

   ! Optimize atom permutation
   do while (results%records(1)%count < max_count .and. results%num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat(), center1)

      ! Assign atoms with current orientation
      call assign_atoms_conf(mnatree, mol1, mol2, coords1, coords2, atomperm, dist)
      total_rotation = optimal_rotation(atomperm, coords1, coords2, center1)
      call rotate_coords(coords2, total_rotation, center1)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_conf(mnatree, mol1, mol2, coords1, coords2, auxperm, dist)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         step_rotation = optimal_rotation(atomperm, coords1, coords2, center1)
         call rotate_coords(coords2, step_rotation, center1)
         total_rotation = quatmul(step_rotation, total_rotation)
         num_steps = num_steps + 1
      end do

      ! Update results
      call registry_push(results, coords2, adjmat2, atomperm, num_steps, total_rotation)

   end do

end subroutine

subroutine assign_atoms_conf( mnatree, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(tree_node), intent(in) :: mnatree
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(tree_node), pointer :: submnatree
   logical :: assigned

   write (stderr, *)
   write (stderr, *) repeat('assign atoms   ', 3)

   submnatree = mnatree
   call print_tree(submnatree)
   do
      assigned = .false.
      call reduce_partial_matches(submnatree, assigned)
      if (.not. assigned) exit
      call compute_consistent_mnatypes(mol1, mol2, submnatree)
      call flatten_tree(submnatree)
      call print_tree(submnatree)
   end do
   stop

end subroutine

!subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
!   type(mol_type), intent(in) :: mol1, mol2
!   type(bipartition_container), intent(in) :: mnatypes
!   real(rk), dimension(:,:), intent(in) :: coords1, coords2
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
!   real(rk), dimension(:,:), intent(in) :: coords1, coords2
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
!      call solve_lap(submnatypes%parts(h), coords1, coords2, atomperm, dist)
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
!   real(rk), dimension(:,:), intent(in) :: coords1, coords2
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
!         call solve_lap(submnatypes%parts(h), coords1, coords2, atomperm, dist)
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
