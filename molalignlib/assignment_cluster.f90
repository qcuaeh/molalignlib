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
use globals
use random
use molecule
use chemdata
use spatial_transforms
use assignment
use lcrs_tree
use atom_mnas
use atom_types
use pruning
use registration

implicit none

contains

subroutine optimize_atomperm_atoms(mol1, mol2, registry)
   type(mol_type), intent(in) :: mol1, mol2
   type(registry_t), intent(out) :: registry

   ! Local variables
   type(partition_t) :: atomtypes
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
   coords1 = get_coords( mol1)
   coords2 = get_coords( mol2)

   ! Mirror coordinates
   if (mirror_flag) then
      call mirror_coords(coords2)
   end if

   ! Mass weight coordinates
   call weight_coords( coords1, atomic_weights(elnums1))
   call weight_coords( coords2, atomic_weights(elnums2))

   ! Translate atoms to their centroids
   center1 = centroid( coords1)
   center2 = centroid( coords2)
   call translate_coords( coords2, -center2)
   call translate_coords( coords2, center1)

   ! Find unfeasible assignments
   call prune_procedure( atomtypes, mol1, mol2, prunes)

   ! Initialize random number generator
   call random_initialize()

   ! Initialize local minima registry
   call init_rmsd_registry( registry, max_records, coords1)

   ! Optimize atom permutation
   do while (registry%records(1)%count < max_count .and. registry%num_trials < max_trials)

      ! Aply a random rotation to coords2
      call rotate_coords( coords2, randrotquat(), center1)

      ! Assign atoms with current orientation
      call assign_atoms_pruned( atomtypes, coords1, coords2, prunes, atomperm)
      call optimize_rotation( atomperm, coords1, coords2, center1, total_rotation)
      call rotate_coords( coords2, total_rotation, center1)
      num_steps = 1

      do while (iter_flag)
         call assign_atoms_pruned( atomtypes, coords1, coords2, prunes, auxperm)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         call optimize_rotation( atomperm, coords1, coords2, center1, step_rotation)
         call rotate_coords( coords2, step_rotation, center1)
         total_rotation = quatmul( step_rotation, total_rotation)
         num_steps = num_steps + 1
      end do

      ! Push local minimum to registry
      call push_record( registry, atomperm, coords2=coords2, num_steps=num_steps, rotation=total_rotation)

   end do

end subroutine

end module
