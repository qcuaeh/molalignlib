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

subroutine optimize_atomperm_conformer( mol1, mol2, atomtypes, registry)
   type(mol_type), intent(inout) :: mol1, mol2
   type(partition_t), intent(in) :: atomtypes
   type(registry_t), target, intent(out) :: registry

   ! Local variables
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: atomperm, auxperm
   type(partree_node_t), pointer :: part_tree, temp_part
   type(assigntree_node_t), pointer :: mnachain, assignment_tree
   type(array_trees_t) :: array_trees
   real(rk), dimension(:), allocatable :: weights1, weights2
   real(rk), dimension(:,:), allocatable :: wcoords1, wcoords2, rcoords2
   real(rk) :: rmsd, center1(3), center2(3), rotation_step(4), rotation(4)
   integer :: num_trials, num_steps

   allocate (atomperm(size(mol1%atoms)))
   allocate (auxperm(size(mol1%atoms)))
   weights1 = atomic_weights(mol1%atoms%elnum)
   weights2 = atomic_weights(mol2%atoms%elnum)
   wcoords1 = get_coords( mol1)
   wcoords2 = get_coords( mol2)
   adjmat1 = get_adjmat( mol1)
   adjmat2 = get_adjmat( mol2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords( wcoords2)
   end if

   ! Calculate centroids
   center1 = centroid( wcoords1, weights1)
   center2 = centroid( wcoords2, weights2)

   ! Translate atoms to their centroids
   call translate_coords( wcoords1, -center1)
   call translate_coords( wcoords2, -center2)

   ! Weight coordinates
   call weight_coords( wcoords1, weights1)
   call weight_coords( wcoords2, weights2)

   ! Pre-compute assignment tree
   call init_chain_from_partition( atomtypes, mnachain, temp_part)
   call compute_consistent_mnas( mol1, mol2, mnachain)
   call precompute_assignment_tree( mol1, mol2, part_tree, mnachain%last_link, assignment_tree)
   call print_chain_tree( assignment_tree)
   call convert_trees_to_arrays( part_tree, assignment_tree, array_trees, mol1, mol2)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, max_records)

   ! Optimize atom permutation
   do while (registry%records(1)%count < max_count .and. registry%num_trials < max_trials)

      num_trials = num_trials + 1

      ! Get randomly rotated wcoords2
      rotation = randrotquat()
      rcoords2 = rotated_coords( wcoords2, rotation)

      ! Assign atoms with current orientation
      call distribute_items_dfs( wcoords1, rcoords2, array_trees, atomperm)
      call align_coords( atomperm, wcoords1, rcoords2, rotation_step)
      rotation = quatmul( rotation, rotation_step)
      num_steps = 1

      do while (iter_flag)
         call distribute_items_dfs( wcoords1, rcoords2, array_trees, auxperm)
!         write (stderr,'(F8.4)') sqrt( total_sqdist( auxperm, wcoords1, rcoords2))
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         call align_coords( atomperm, wcoords1, rcoords2, rotation_step)
         rotation = quatmul( rotation, rotation_step)
         num_steps = num_steps + 1
      end do

      ! Update results
      rmsd = sqrt( total_sqdist( atomperm, wcoords1, rcoords2))
      call push_record( registry, atomperm, num_steps, rmsd=rmsd, rotation=rotation)
!      write (stderr,'(I0)') adjacencydiff( atomperm, adjmat1, adjmat2)

   end do
end subroutine

end module
