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

module remapping_free
use parameters
use globals
use random
use molecule
use chemdata
use spatial
use assignment
use lcrs_tree
use partitioning
use pruning
use registry

implicit none

contains

subroutine remap_free_atoms(mol1, mol2, results)
   type(mol_type), intent(in) :: mol1, mol2
   type(atomperm_registry), intent(out) :: results

   ! Local variables
   type(tree_node), pointer :: eltree
   type(bipartition_container) :: eltypes
   type(bool_matrix), dimension(:), allocatable :: prunes
   integer :: num_steps
   integer, dimension(:), allocatable :: atomperm, auxperm
   integer, dimension(:), allocatable :: elnums1, elnums2
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: step_rotation(4), total_rotation(4)
   real(rk) :: center1(3), center2(3)

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

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords(coords2)
   end if

   ! Mass weight coordinates
   call weight_coords(coords1, atomic_weights(elnums1))
   call weight_coords(coords2, atomic_weights(elnums2))

   ! Translate atoms to their centroids
   center1 = centroid(coords1)
   center2 = centroid(coords2)
   call translate_coords(coords2, -center2)
   call translate_coords(coords2, center1)

   ! Find unfeasible assignments
   call prune_procedure(eltypes, mol1, mol2, prunes)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call registry_init(results, max_records, coords1)

   ! Optimize atom permutation
   do while (results%records(1)%count < max_count .and. results%num_trials < max_trials)

      ! Aply a random rotation to coords2
      call rotate_coords(coords2, randrotquat(), center1)

      ! Assign atoms with current orientation
      call assign_atoms_pruned(eltypes, coords1, coords2, prunes, atomperm)
      total_rotation = optimal_rotation(atomperm, coords1, coords2, center1)
      call rotate_coords(coords2, total_rotation, center1)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_pruned(eltypes, coords1, coords2, prunes, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         step_rotation = optimal_rotation(atomperm, coords1, coords2, center1)
         call rotate_coords(coords2, step_rotation, center1)
         total_rotation = quatmul(step_rotation, total_rotation)
         num_steps = num_steps + 1
      end do

      ! Push local minimum to registry
      call registry_push(results, coords2, atomperm, num_steps, total_rotation)

   end do

end subroutine

end module
