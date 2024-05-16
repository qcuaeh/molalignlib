module adjacency
use stdio
use kinds
use flags
use bounds
use sorting
use chemdata

implicit none

contains

subroutine adjmat2list(natom, adjmat, nadjs, adjlists)
   integer, intent(in) :: natom
   logical, dimension(:, :), intent(in) :: adjmat
   integer, dimension(:), intent(out) :: nadjs
   integer, dimension(:, :), intent(out) :: adjlists
   integer :: i, j

   nadjs(:) = 0

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            nadjs(i) = nadjs(i) + 1
            nadjs(j) = nadjs(j) + 1
            if (nadjs(i) > maxcoord .or. nadjs(j) > maxcoord) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlists(nadjs(i), i) = j
            adjlists(nadjs(j), j) = i
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

subroutine adjlist2bonds(natom, nadjs, adjlists, nbond, bonds)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: nadjs
   integer, dimension(:, :), intent(in) :: adjlists
   integer, intent(out) :: nbond
   integer, dimension(:, :), intent(out) :: bonds
   integer :: i, j

   nbond = 0

   do i = 1, natom
      do j =1, nadjs(i)
         nbond = nbond + 1
         bonds(1, nbond) = i
         bonds(2, nbond) = adjlists(j, i)
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

function adjacencydelta(nadjs0, adjlists0, adjmat1, mapping, k, l) result(delta)
   integer, intent(in) :: k, l
   integer, dimension(:), intent(in) :: mapping, nadjs0
   integer, dimension(:, :), intent(in) :: adjlists0
   logical, dimension(:, :), intent(in) :: adjmat1
   integer i, nkk, nkl, nll, nlk, delta

   nkk = 0
   nkl = 0

   do i = 1, nadjs0(k)
      if (adjlists0(i, k) /= l) then
         if (adjmat1(mapping(k), mapping(adjlists0(i, k)))) nkk = nkk + 1
         if (adjmat1(mapping(l), mapping(adjlists0(i, k)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, nadjs0(l)
      if (adjlists0(i, l) /= k) then
         if (adjmat1(mapping(l), mapping(adjlists0(i, l)))) nll = nll + 1
         if (adjmat1(mapping(k), mapping(adjlists0(i, l)))) nlk = nlk + 1
      end if
   end do

!        dkk = nadjs0(k) + nadjs1(mapping(k)) - 2*nkk
!        dll = nadjs0(l) + nadjs1(mapping(l)) - 2*nll
!        dkl = nadjs0(k) + nadjs1(mapping(l)) - 2*nkl
!        dlk = nadjs0(l) + nadjs1(mapping(k)) - 2*nlk
!        delta = dkl + dlk - dkk - dll

   ! Notice that dkl + dlk - dkk - dll == 2*(nkk + nll - nkl - nlk)
   delta = 2*(nkk + nll - nkl - nlk)

end function

end module
