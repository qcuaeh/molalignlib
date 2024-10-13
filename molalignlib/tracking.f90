! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragidcs)
module tracking

use stdio
use molecule
use bounds
use sorting
use partition

implicit none
private
public find_molfrags

type(indexlist_type), allocatable :: adjlists(:)

contains

subroutine find_molfrags (mol)
   type(molecule_type), intent(inout) :: mol
   ! Local variables
   integer :: n, nfrag
   integer, allocatable :: fragidcs(:), fragsize(:)
   logical, allocatable :: tracked(:)

   allocate (fragidcs(size(mol%atoms)))
   allocate (fragsize(size(mol%atoms)))
   allocate (tracked(size(mol%atoms)))

   adjlists = mol%get_adjlists()

   ! initialization

   nfrag = 0
   tracked(:) = .false.
   fragidcs(:) = 0
   fragsize(:) = 0

   ! detect fragments and populate frag arrays
   n = 1
   do while (n <= size(mol%atoms))
      if (tracked(n)) then
         n = n + 1
      else
         nfrag = nfrag + 1
         call recrun(tracked, n, nfrag, fragidcs, fragsize)
         n = 1
      end if
   end do

   ! register nfrag and fragidcs in the mol strucutre
   call mol%set_molfrags(nfrag, fragidcs)

end subroutine

recursive subroutine recrun (tracked, n, nfrag, fragidcs, fragsize)
! runs recursivelly over the structure and populates arrays
   logical, intent(inout) :: tracked(:)
   integer, intent(in) :: n, nfrag
   integer, intent(inout) :: fragidcs(:), fragsize(:)
   ! Local variables
   integer :: i

   if (tracked(n)) return

   tracked(n) = .true.
   fragidcs(n) = nfrag
   fragsize(nfrag) = fragsize(nfrag) + 1
   
   do i = 1, size(adjlists(n)%indices)
      call recrun(tracked, adjlists(n)%indices(i), nfrag, fragidcs, fragsize)
   end do

end subroutine

end module
