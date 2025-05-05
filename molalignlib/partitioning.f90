module partitioning
use parameters
use sorting
use chemdata
use molecule
use lcrs_tree

implicit none

type :: atomtype_item
   integer :: elnum
   integer :: label
   type(leaf_node), pointer :: node
end type

type :: atomtype_table
   integer :: num_items
   type(atomtype_item), dimension(:), allocatable :: items
end type

contains

subroutine add_atomtype(atomtypetable, elnum, label, leaf)
   type(atomtype_table), intent(inout) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(leaf_node), pointer, intent(in) :: leaf

   atomtypetable%num_items = atomtypetable%num_items + 1
   atomtypetable%items(atomtypetable%num_items)%elnum = elnum
   atomtypetable%items(atomtypetable%num_items)%label = label
   atomtypetable%items(atomtypetable%num_items)%node => leaf
end subroutine

function find_atomtype(atomtypetable, elnum, label) result(leaf)
   type(atomtype_table), intent(in) :: atomtypetable
   integer, intent(in) :: elnum
   integer, intent(in) :: label
   type(leaf_node), pointer :: leaf
   integer :: i

   do i = 1, atomtypetable%num_items
      if (atomtypetable%items(i)%elnum == elnum .and. &
          atomtypetable%items(i)%label == label) then
         leaf => atomtypetable%items(i)%node
         return
      end if
   end do

   leaf => null()
end function

! Partition atoms by atomic number and label
subroutine compute_eltypes(mol1, mol2, eltypetree)
   type(mol_type), intent(in) :: mol1, mol2
   type(root_node), pointer, intent(out) :: eltypetree
   ! Local variables
   type(leaf_node), pointer :: node
   type(atomtype_table) :: atomtypetable
   integer :: i, elnum, label, num_atoms1, num_atoms2
   type(leaf_node_ptr) :: typehood(0)

   num_atoms1 = size(mol1%atoms)
   num_atoms2 = size(mol2%atoms)

   eltypetree => make_new_root(num_atoms1, num_atoms2)
   allocate (atomtypetable%items(num_atoms1 + num_atoms2))
   atomtypetable%num_items = 0

   ! First molecule
   do i = 1, num_atoms1
      elnum = mol1%atoms(i)%elnum
      label = mol1%atoms(i)%label
      node => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(node)) then
         node => add_new_leaf(eltypetree, typehood)
         call add_atomtype(atomtypetable, elnum, label, node)
      end if
      call add_new_item1(node, i)
   end do

   ! Second molecule
   do i = 1, num_atoms2
      elnum = mol2%atoms(i)%elnum
      label = mol2%atoms(i)%label
      node => find_atomtype(atomtypetable, elnum, label)
      if (.not. associated(node)) then
         node => add_new_leaf(eltypetree, typehood)
         call add_atomtype(atomtypetable, elnum, label, node)
      end if
      call add_new_item2(node, i)
   end do
end subroutine

subroutine compute_nextlevel_children(mol1, mol2, next_root, leaf)
   type(mol_type), intent(in) :: mol1, mol2
   type(root_node), intent(inout) :: next_root
   type(leaf_node), intent(inout) :: leaf
   ! Local variables
   type(leaf_node), pointer :: heir_leaf
   type(item_node), pointer :: item
   type(leaf_node_ptr), dimension(:), allocatable :: typehood

   ! First molecule
   item => leaf%first_item1
   do while (associated(item))
      typehood = leaf%tree_root%itemdir1(mol1%atoms(item%index)%adjlist)
      heir_leaf => find_heir(leaf, typehood)
      if (.not. associated(heir_leaf)) then
         heir_leaf => add_new_leaf(next_root, typehood)
         call add_heir(leaf, heir_leaf)
      end if
      call add_new_item1(heir_leaf, item%index)
      item => item%next_item
   end do

   ! Second molecule
   item => leaf%first_item2
   do while (associated(item))
      typehood = leaf%tree_root%itemdir2(mol2%atoms(item%index)%adjlist)
      heir_leaf => find_heir(leaf, typehood)
      if (.not. associated(heir_leaf)) then
         heir_leaf => add_new_leaf(next_root, typehood)
         call add_heir(leaf, heir_leaf)
      end if
      call add_new_item2(heir_leaf, item%index)
      item => item%next_item
   end do
end subroutine

subroutine compute_nextlevel_mnas(mol1, mol2, mnapolytree)
! Compute next level MNA types
   type(mol_type), intent(in) :: mol1, mol2
   type(poly_node), intent(inout) :: mnapolytree

   ! Local variables
   type(leaf_node), pointer :: leaf
   type(root_node), pointer :: next_root

   leaf => mnapolytree%last_root%first_leaf
   next_root => add_new_root(mnapolytree)

   do while (associated(leaf))
      call compute_nextlevel_children(mol1, mol2, next_root, leaf)
      leaf => leaf%next_leaf
   end do
end subroutine

subroutine compute_consistent_mnas(mol1, mol2, mnapolytree)
! Iteratively compute MNA types
   type(mol_type), intent(in) :: mol1, mol2
   type(poly_node), intent(inout) :: mnapolytree
   integer :: prev_num_leaves

   do
      ! Compute MNA upper level types
      prev_num_leaves = mnapolytree%last_root%num_leaves
      call compute_nextlevel_mnas(mol1, mol2, mnapolytree)

      ! Exit loop if types did not change
      if (mnapolytree%last_root%num_leaves == prev_num_leaves) exit

!      call print_tree(mnapolytree%last_root)
   end do
end subroutine

subroutine polytree_to_stack(poly, stack)
   type(poly_node), intent(in) :: poly
   type(partition_stack) :: stack
   type(root_node), pointer :: root
   type(leaf_node), pointer :: leaf
   type(leaf_node), pointer :: heir
   type(item_node), pointer :: item
   integer :: i, j, k, m, n

   allocate (stack%partitions(poly%num_roots))
   stack%num_partitions = poly%num_roots

   stack%partitions(1) = partition_from_tree(poly%first_root)

   ! Copy trees
   root => poly%first_root
   do i = 2, poly%num_roots

      allocate (stack%partitions(i)%itemdir1(poly%size_itemdir1))
      allocate (stack%partitions(i)%itemdir2(poly%size_itemdir2))
      allocate (stack%partitions(i)%parts(root%next_root%num_leaves))
      stack%partitions(i)%num_parts = root%next_root%num_leaves

      j = 1
      leaf => root%first_leaf
      do m = 1, root%num_leaves

         allocate (stack%partitions(i-1)%parts(m)%heirs(leaf%num_heirs))
         stack%partitions(i-1)%parts(m)%num_heirs = leaf%num_heirs

         ! Copy heirs
         heir => leaf%first_heir
         do n = 1, leaf%num_heirs

            stack%partitions(i-1)%parts(m)%heirs(n) = j
            stack%partitions(i)%parts(j)%num_items1 = heir%num_items1
            stack%partitions(i)%parts(j)%num_items2 = heir%num_items2

            allocate (stack%partitions(i)%parts(j)%items1(heir%num_items1))
            allocate (stack%partitions(i)%parts(j)%items2(heir%num_items2))

            item => heir%first_item1
            do k = 1, heir%num_items1
               stack%partitions(i)%parts(j)%items1(k) = item%index
               stack%partitions(i)%itemdir1(item%index) = j
               item => item%next_item
            end do

            item => heir%first_item2
            do k = 1, heir%num_items2
               stack%partitions(i)%parts(j)%items2(k) = item%index
               stack%partitions(i)%itemdir2(item%index) = j
               item => item%next_item
            end do

            heir => heir%next_heir
            j = j + 1
         end do

         leaf => leaf%next_leaf
      end do

      root => root%next_root
   end do

   stack%partitions(poly%num_roots)%parts(:)%num_heirs = 0
end subroutine

end module
