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

module assignment_conformer
use parameters
use globals
use random
use molecule
use strutils
use chemdata
use permutation
use spatial_transforms
use assignment
use adjacency
use biasing
use pruning
use lcrs_tree
use array_trees
use atom_mnas
use assigntree_precompute
use assigntree_recompute
use assigntree_distribute
use registration
use fileio

implicit none

contains

subroutine optimize_atomperm_conformer( mol1, mol2, registry)
   type(mol_type), intent(inout) :: mol1, mol2
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   type(partition_t) :: atomtypes
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: atomperm, auxperm
   integer, dimension(:), allocatable :: elnums1, elnums2
   type(partree_node_t), pointer :: part_tree, temp_part
   type(assigntree_node_t), pointer :: mnachain, assignment_tree
   type(array_trees_t) :: array_trees
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: step_rotation(4), total_rotation(4)
   real(rk) :: center1(3), center2(3)
   integer :: num_trials, num_steps

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
   call set_eltypes( mol1%atoms, mol2%atoms, atomtypes)

   ! Abort if there are conflicting atomic types
   if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
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

   call init_chain_from_partition( atomtypes, mnachain, temp_part)
   call compute_consistent_mnas( mol1, mol2, mnachain)
   call precompute_assignment_tree( mol1, mol2, part_tree, mnachain%last_link, assignment_tree)
   call print_chain_tree( assignment_tree)
   call convert_trees_to_arrays( part_tree, assignment_tree, array_trees, mol1, mol2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( coords2)
   end if

   ! Mass weight coordinates
   call weight_coords( coords1, atomic_weights(elnums1))
   call weight_coords( coords2, atomic_weights(elnums2))

   ! Translate atoms to their centroids
   center1 = centroid( coords1)
   center2 = centroid( coords2)
   call translate_coords( coords2, -center2)
   call translate_coords( coords2, center1)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
!   call init_rmsd_registry( registry, max_records, coords1)
   call init_dual_registry( registry, max_records, coords1, adjmat1)

   ! Optimize atom permutation
   do while (registry%records(1)%count < max_count .and. registry%num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords( coords2, randrotquat(), center1)

      ! Assign atoms with current orientation
      call distribute_items_dfs( coords1, coords2, array_trees, atomperm)
      call optimize_rotation( atomperm, coords1, coords2, center1, total_rotation)
      call rotate_coords( coords2, total_rotation, center1)
      num_steps = 1

      do while (iter_flag)
         call distribute_items_dfs( coords1, coords2, array_trees, auxperm)
!         write (stderr, '(F8.4,F8.4,F8.4)') sqrt(total_sqdist(auxperm, coords1, coords2))
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         call optimize_rotation( atomperm, coords1, coords2, center1, step_rotation)
         call rotate_coords( coords2, step_rotation, center1)
         total_rotation = quatmul( step_rotation, total_rotation)
         num_steps = num_steps + 1
      end do

      ! Update results
!      call push_record( registry, atomperm, coords2=coords2, num_steps=num_steps, rotation=total_rotation)
      call push_record( registry, atomperm, coords2=coords2, adjmat2=adjmat2, num_steps=num_steps, rotation=total_rotation)

   end do
end subroutine

end module
