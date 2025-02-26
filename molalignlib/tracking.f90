! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragidcs)
module tracking
use parameters
use sorting
use molecule
use common_types
use lcrs_tree

implicit none
private
public find_molfrags

type(atom_type), allocatable :: atoms(:)

contains

subroutine find_molfrags( mol, eltypes, molfrags)
   type(mol_type), intent(in) :: mol
   type(partition_container), intent(in) :: eltypes
   type(intlist_type), allocatable, intent(out) :: molfrags(:)
   ! Local variables
   integer :: i, nfrag
   logical, allocatable :: tracked(:)
   integer, allocatable :: fragszs(:), fragidcs(:,:)
   integer, allocatable :: order(:)

   allocate (fragszs(size(mol%atoms)))
   allocate (fragidcs(size(mol%atoms), size(mol%atoms)))
   allocate (tracked(size(mol%atoms)))
   atoms = mol%atoms

   ! initialization

   nfrag = 0
   fragszs(:) = 0
   tracked(:) = .false.

   ! detect fragments and populate frag arrays
   i = 1
   do while (i <= size(mol%atoms))
      if (tracked(i)) then
         i = i + 1
      else
         nfrag = nfrag + 1
         call recrun( tracked, i, nfrag, fragszs, fragidcs)
         i = 1
      end if
   end do

   ! Order molecular fragments
   do i = 1, nfrag
      order = sorted_order(eltypes%parts(eltypes%itemdir(fragidcs(:fragszs(i), i)))%num_items)
      fragidcs(:fragszs(i), i) = fragidcs(order, i)
!      write (stderr, *) fragidcs1(:fragszs1(i), i)
!      write (stderr, *)
   end do

   allocate (molfrags(nfrag))

   do i = 1, nfrag
      molfrags(i)%n = fragidcs(:fragszs(i), i)
   end do

end subroutine

recursive subroutine recrun( tracked, iatom, nfrag, fragszs, fragidcs)
! runs recursivelly over the structure and populates arrays
   logical, intent(inout) :: tracked(:)
   integer, intent(in) :: iatom, nfrag
   integer, intent(inout) :: fragidcs(:,:), fragszs(:)
   ! Local variables
   integer :: i

   if (tracked(iatom)) return

   tracked(iatom) = .true.
   fragszs(nfrag) = fragszs(nfrag) + 1
   fragidcs(fragszs(nfrag), nfrag) = iatom

   do i = 1, size(atoms(iatom)%adjlist)
      call recrun( tracked, atoms(iatom)%adjlist(i), nfrag, fragszs, fragidcs)
   end do

end subroutine

end module
