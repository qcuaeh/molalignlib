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

module assignment_cluster
use parameters
use options
use random
use molecule
use chemistry
use spatial_transforms
use assignment_atoms
use lcrs_tree
use atom_mnas
use atom_types
use pruning
use registration

implicit none

contains

subroutine optimize_atomperm_cluster(atoms1, atoms2, atomtypes, registry)
   type(atom_t), dimension(:), intent(in) :: atoms1, atoms2
   type(partition_t), intent(in) :: atomtypes
   type(registry_t), intent(out) :: registry

   ! Local variables
   type(bool_matrix), dimension(:), allocatable :: prunes
   integer :: num_steps
   integer, dimension(:), allocatable :: atomperm, auxperm
   real(rk), dimension(:), allocatable :: weights1, weights2
   real(rk), dimension(:,:), allocatable :: wcoords1, wcoords2, rcoords2
   real(rk) :: rmsd, center1(3), center2(3), rotation_step(4), rotation(4)

   allocate (atomperm(size(atoms1)))
   allocate (auxperm(size(atoms1)))
   weights1 = atomic_weights(atoms1%elnum)
   weights2 = atomic_weights(atoms2%elnum)
   wcoords1 = get_coords( atoms1)
   wcoords2 = get_coords( atoms2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords(wcoords2)
   end if

   ! Find unfeasible assignments
   call prune_procedure( atomtypes, atoms1, atoms2, prunes)

   ! Calculate centroids
   center1 = centroid( wcoords1, weights1)
   center2 = centroid( wcoords2, weights2)

   ! Translate atoms to their centroids
   call translate_coords( wcoords1, -center1)
   call translate_coords( wcoords2, -center2)

   ! Weight coordinates
   call weight_coords( wcoords1, weights1)
   call weight_coords( wcoords2, weights2)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, max_records)

   ! Optimize atom permutation
   do while (registry%records(1)%count < max_count .and. registry%num_trials < max_trials)

      ! Aply a random rotation to wcoords2
      rotation = randrotquat()
      rcoords2 = rotated_coords( wcoords2, rotation)

      ! Assign atoms with current orientation
      call assign_atoms_pruned( atomtypes, wcoords1, rcoords2, prunes, atomperm)
      call align_coords( atomperm, wcoords1, rcoords2, rotation_step)
      rotation = quatmul( rotation, rotation_step)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_pruned( atomtypes, wcoords1, rcoords2, prunes, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         call align_coords( atomperm, wcoords1, rcoords2, rotation_step)
         rotation = quatmul( rotation, rotation_step)
         num_steps = num_steps + 1
      end do

      ! Push local minimum to registry
      rmsd = sqrt( total_sqdist( atomperm, wcoords1, rcoords2))
      call push_record( registry, atomperm, num_steps, rmsd=rmsd, rotation=rotation)

   end do

end subroutine

end module
