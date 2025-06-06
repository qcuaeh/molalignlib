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
use mna_compute
use mna_precompute
use mna_recompute
use reactivity
use registry
use fileio

implicit none

contains

subroutine remap_bonded_atoms(mol1, mol2, results)
   type(mol_type), intent(inout) :: mol1, mol2
   type(bondatomperm_registry), target, intent(out) :: results

   ! Local variables
   type(partitionarray_t) :: eltypes
   logical, dimension(:,:), allocatable :: adjmat1, adjmat2
   integer, dimension(:), allocatable :: atomperm, auxperm
   integer, dimension(:), allocatable :: elnums1, elnums2
   integer :: num_trials, num_steps
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   real(rk) :: step_rotation(4), total_rotation(4)
   real(rk) :: center1(3), center2(3)
!   real(rk) :: dist

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
   call set_eltypes( mol1%atoms, mol2%atoms, eltypes)

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
      call find_reactive_bonds( mol1, mol2, eltypes, atomperm)
      write (stderr, *) 'adjd before', adjacencydiff( atomperm, adjmat1, adjmat2)
      call remove_reactive_bonds( mol1, mol2, eltypes, atomperm)
      adjmat1 = get_adjmat( mol1)
      adjmat2 = get_adjmat( mol2)
      write (stderr, *) 'adjd after ', adjacencydiff( atomperm, adjmat1, adjmat2)
   end if

   call assign_conform_atoms( mol1, mol2, eltypes)
   stop

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
   call registry_init( results, max_records, coords1, adjmat1)

   ! Optimize atom permutation
   do while (results%records(1)%count < max_count .and. results%num_trials < max_trials)

      num_trials = num_trials + 1

      ! Aply a random rotation to coords2
      call rotate_coords( coords2, randrotquat(), center1)

      ! Assign atoms with current orientation
!      call minimize_conformation_distance( split_branch, coords1, coords2, atomperm, dist)
      total_rotation = optimal_rotation( atomperm, coords1, coords2, center1)
      call rotate_coords( coords2, total_rotation, center1)
      num_steps = 1

      do while (iter_flag)
!         call minimize_conformation_distance( split_branch, coords1, coords2, auxperm, dist)
         if (all(auxperm == atomperm)) exit
         atomperm = auxperm
         step_rotation = optimal_rotation( atomperm, coords1, coords2, center1)
         call rotate_coords( coords2, step_rotation, center1)
         total_rotation = quatmul( step_rotation, total_rotation)
         num_steps = num_steps + 1
      end do

      ! Update results
      call registry_push( results, coords2, adjmat2, atomperm, num_steps, total_rotation)

   end do
end subroutine

subroutine assign_conform_atoms( mol1, mol2, eltypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(partitionarray_t), intent(in) :: eltypes
   ! Local variables
   type(part_node_t), pointer :: root_part, temp_part
   type(split_node_t), pointer :: mnachain, temp_chain, root_branch
   type(link_node_t), pointer :: branch_parts
   type(part_node_t), pointer :: child_part
   integer unit1, unit2

   open(newunit=unit1, file='molec1.mol2', action='write', status='replace')
   open(newunit=unit2, file='molec2.mol2', action='write', status='replace')
   call writefile( unit1, 'mol2', mol1)
   call writefile( unit2, 'mol2', mol2)

!   call print_atoms( mol1)
!   call print_atoms( mol2)

!   call init_chain_from_partarray( eltypes, mnachain, root_part)
!   root_branch => new_root_branch( mnachain%tot_items1, mnachain%tot_items2)
!   leaf_link => new_generic_link( root_branch)
!   call update_itemdir_children( leaf_link, root_part)
!   branch_parts => new_root_link()
!   call update_branch_parts(branch_parts, root_part)
!   call precompute_consistent_mnas(mol1, mol2, mnachain, root_branch, branch_parts)

   call init_chain_from_partarray( eltypes, temp_chain, temp_part)
   call compute_consistent_mnas( mol1, mol2, temp_chain)
   call init_chain_from_link( temp_chain%last_link, mnachain, root_part)
   call print_tree_items( root_part)

   root_branch => new_root_branch( mnachain%tot_items1, mnachain%tot_items2)
   branch_parts => new_root_link()
   child_part => root_part%first_child_part
   do while (associated(child_part))
      call add_branch_part(branch_parts, child_part)
      child_part => child_part%next_sibling_part
   end do

   call split_independent_parts( mol1, mol2, mnachain, root_branch, branch_parts)
   call print_part_tree( root_part)
!   call print_tree_items( root_part)

!   call print_tree_signatures( root_part)
   call print_split_tree( root_branch)
!   call print_chain( mnachain)
!   call print_chain( root_branch%first_child_branch)

   call random_init(.false., .true.)
   call redistribute_items( mol1, mol2, root_branch)
   call print_leaf_items( root_part)
end subroutine

end module
