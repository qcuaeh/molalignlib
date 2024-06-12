module adjacency
use stdio
use kinds
use flags
use bounds
use sorting
use chemdata

implicit none

contains

subroutine adjmat2list(natom, adjmat, coonums, neighbors)
   integer, intent(in) :: natom
   logical, dimension(:, :), intent(in) :: adjmat
   integer, dimension(:), intent(out) :: coonums
   integer, dimension(:, :), intent(out) :: neighbors
   integer :: i, j

   coonums(:) = 0

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            coonums(i) = coonums(i) + 1
            coonums(j) = coonums(j) + 1
            if (coonums(i) > maxcoord .or. coonums(j) > maxcoord) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            neighbors(coonums(i), i) = j
            neighbors(coonums(j), j) = i
         end if
      end do
   end do

end subroutine

subroutine adjmat2bonds(natom, adjmat, nbond, bonds)
   integer, intent(in) :: natom
   logical, dimension(:, :), intent(in) :: adjmat
   integer, intent(out) :: nbond
   integer, dimension(:, :), intent(out) :: bonds
   integer :: i, j

   nbond = 0

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            nbond = nbond + 1
            bonds(1, nbond) = i
            bonds(2, nbond) = j
         end if
      end do
   end do

end subroutine

subroutine adjlist2bonds(natom, coonums, neighbors, nbond, bonds)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: coonums
   integer, dimension(:, :), intent(in) :: neighbors
   integer, intent(out) :: nbond
   integer, dimension(:, :), intent(out) :: bonds
   integer :: i, j

   nbond = 0

   do i = 1, natom
      do j =1, coonums(i)
         nbond = nbond + 1
         bonds(1, nbond) = i
         bonds(2, nbond) = neighbors(j, i)
      end do
   end do

end subroutine

function adjacencydiff(natom, adjmat0, adjmat1, mapping) result(diff)
! Purpose: Check if two graphs are equal.
! Return the number of differences between graphs.
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: mapping
   logical, dimension(:, :), intent(in) :: adjmat0, adjmat1
   integer diff

   integer i, j

   diff = 0

! Check differences element by element

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat0(i, j) .neqv. adjmat1(mapping(i), mapping(j))) then
            diff = diff + 1
!                print *, i, j, adjmat0(i, j), adjmat1(mapping(i), mapping(j))
         end if
      end do
   end do

end function

function adjacencydelta(coonums0, neighbors0, adjmat1, mapping, k, l) result(delta)
   integer, intent(in) :: k, l
   integer, dimension(:), intent(in) :: mapping, coonums0
   integer, dimension(:, :), intent(in) :: neighbors0
   logical, dimension(:, :), intent(in) :: adjmat1
   integer i, nkk, nkl, nll, nlk, delta

   nkk = 0
   nkl = 0

   do i = 1, coonums0(k)
      if (neighbors0(i, k) /= l) then
         if (adjmat1(mapping(k), mapping(neighbors0(i, k)))) nkk = nkk + 1
         if (adjmat1(mapping(l), mapping(neighbors0(i, k)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, coonums0(l)
      if (neighbors0(i, l) /= k) then
         if (adjmat1(mapping(l), mapping(neighbors0(i, l)))) nll = nll + 1
         if (adjmat1(mapping(k), mapping(neighbors0(i, l)))) nlk = nlk + 1
      end if
   end do

!        dkk = coonums0(k) + coonums1(mapping(k)) - 2*nkk
!        dll = coonums0(l) + coonums1(mapping(l)) - 2*nll
!        dkl = coonums0(k) + coonums1(mapping(l)) - 2*nkl
!        dlk = coonums0(l) + coonums1(mapping(k)) - 2*nlk
!        delta = dkl + dlk - dkk - dll

   ! Notice that dkl + dlk - dkk - dll == 2*(nkk + nll - nkl - nlk)
   delta = 2*(nkk + nll - nkl - nlk)

end function

end module
