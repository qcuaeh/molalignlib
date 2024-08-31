! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragidcs)
module tracking

use stdio
use types
use bounds
use sorting
implicit none
private
public find_molfrags

contains

subroutine find_molfrags (mol)
   type(Molecule), intent(inout) :: mol

   integer :: i, n
   integer :: natom, nfrag
   integer, dimension(mol%natom) :: nadjs
   integer, dimension(maxcoord, mol%natom) :: adjlists
   integer, dimension(mol%natom) :: fragidcs
   integer, dimension(mol%natom) :: fragsize
   logical :: tracked(mol%natom)

   natom = mol%get_natom()
   nadjs = mol%get_nadjs()
   adjlists = mol%get_adjlists()

   ! initialization

   nfrag = 0
   tracked(:) = .false.
   fragidcs(:) = 0
   fragsize(:) = 0

   ! detect fragments and populate frag arrays
   n = 1
   do while (n <= natom)
      if (tracked(n)) then
         n = n + 1
      else
         nfrag = nfrag + 1
         call recrun(n, tracked)
         n = 1
      end if
   end do

   ! register nfrag and fragidcs in the mol strucutre
   call mol%set_molfrags(nfrag, fragidcs)

contains

   ! runs recursivelly over the structure and populates arrays
   recursive subroutine recrun (n, tracked)
      integer, intent(in) :: n
      logical, dimension(:), intent(inout) :: tracked
      integer i

      if (tracked(n)) return

      tracked(n) = .true.
      fragidcs(n) = nfrag
      fragsize(nfrag) = fragsize(nfrag) + 1
      
      do i = 1, nadjs(n)
         call recrun(adjlists(i, n), tracked)
      end do

   end subroutine

end subroutine

end module
