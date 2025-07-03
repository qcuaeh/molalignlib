module mna_testing_arrays
use parameters
use derived_types
use molecule
use lcrs_tree
use array_trees
use atom_mnas
use assigntree_precompute
use assigntree_distribute
implicit none

contains

subroutine redistribute_items_random(coords1, coords2, array_trees, random_permutation)
   ! Random exploration wrapper - generates one random assignment
   ! Similar to distribute_items_dfs but generates random assignment instead of optimal
   real(rk), intent(in) :: coords1(:,:), coords2(:,:)
   type(array_trees_t), intent(inout) :: array_trees
   integer, allocatable, intent(out) :: random_permutation(:)
   ! Local variables
   type(assignment_t) :: random_assignment
   integer :: num_atoms, assigned_count
   real(rk) :: total_distance

   num_atoms = array_trees%num_atoms1

   ! Initialize random assignment
   call init_assignment(random_assignment, num_atoms)

   ! Initialize assignment with preassigned pairs
   call collect_leaf_assignments(array_trees, 1, random_assignment)

   ! Perform random exploration to generate one assignment (starting from root chain at index 1)
   call redistribute_items_random_recursive(coords1, coords2, array_trees, 1, random_assignment)

   ! Copy final permutation array from assignment
   random_permutation = random_assignment%permutation

   ! Calculate distance from the random permutation array
   total_distance = total_sqdist(random_assignment%assigned_indices(1:random_assignment%num_assigned), &
      random_assignment%permutation, coords1, coords2)

   ! Count assigned atoms
   assigned_count = random_assignment%num_assigned

   ! Report final results
   write(stderr, '(A)') repeat("=", 60)
   write(stderr, '(A,I0,A,I0,A)') "Atoms assigned: ", assigned_count, " out of ", num_atoms, " total atoms"
   write(stderr, '(A,F10.4)') "Random assignment total squared distance: ", total_distance
   write(stderr, '(A)') repeat("=", 60)

   ! Validate permutation consistency
   call validate_perm(random_permutation)
end subroutine

subroutine assign_conform_atoms( mol1, mol2, atomtypes)
   type(mol_type), intent(in) :: mol1, mol2
   type(partition_t), intent(in) :: atomtypes
   ! Local variables
   type(partree_node_t), pointer :: part_tree, temp_part
   type(assigntree_node_t), pointer :: mnachain, assignment_tree
   type(array_trees_t) :: array_trees
!   type(chain_node_t), pointer :: branch_parts
!   type(partree_node_t), pointer :: child_part
   real(rk), dimension(:,:), allocatable :: coords1, coords2
   integer, allocatable :: atomperm(:)
!   integer :: unit1, unit2

!   open(newunit=unit1, file='molec1.mol2', action='write', status='replace')
!   open(newunit=unit2, file='molec2.mol2', action='write', status='replace')
!   call write_file( unit1, 'mol2', mol1)
!   call write_file( unit2, 'mol2', mol2)
!   call print_atoms( mol1)
!   call print_atoms( mol2)

   coords1 = get_coords(mol1)
   coords2 = get_coords(mol2)

!   call init_chain_from_partition( atomtypes, mnachain, part_tree)
!   assignment_tree => new_root_chain( mnachain%tot_items1, mnachain%tot_items2)
!   leaf_link => new_chain_link( assignment_tree)
!   call update_itemdir_children( leaf_link, part_tree)
!   branch_parts => new_bare_link()
!   call update_branch_parts(branch_parts, part_tree)
!   call precompute_consistent_mnas(mol1, mol2, mnachain, assignment_tree, branch_parts)

   call init_chain_from_partition( atomtypes, mnachain, temp_part)
   call compute_consistent_mnas( mol1, mol2, mnachain)
   call precompute_assignment_tree( mol1, mol2, part_tree, mnachain%last_link, assignment_tree)

!   call print_part_tree( part_tree)
!   call print_part_indices( part_tree)
!   call print_chain_indices( assignment_tree)

!   call print_tree_signatures( part_tree)
!   call print_chain_tree( assignment_tree)

!   call random_init(.true., .true.)
!   call recompute_assignment_tree( mol1, mol2, assignment_tree)
!   call print_tree_items( part_tree)

   call convert_trees_to_arrays( part_tree, assignment_tree, array_trees, mol1, mol2)
   call validate_conversion(part_tree, assignment_tree, array_trees)
!   call print_part_tree_array( array_trees)
!   call print_first_level_items_array( array_trees)
   call print_chain_tree_array( array_trees)
!   call random_init(.true., .true.)
!   call redistribute_items_random( coords1, coords2, array_trees, atomperm)
!   call print_leaf_items_array( array_trees)
   call distribute_items_dfs( coords1, coords2, array_trees, atomperm)
end subroutine

end module
