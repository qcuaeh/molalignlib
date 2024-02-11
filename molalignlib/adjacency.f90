module adjacency
use stdio
use kinds
use flags
use bounds
use sorting
use chemdata

implicit none

contains

subroutine getadjmat(natom, coords, znums, adjmat)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: znums
   real(wp), dimension(:, :), intent(in) :: coords
   logical, dimension(:, :), intent(out) :: adjmat

   integer :: i, j
   real(wp) :: atomdist
   real(wp) :: adjrad(natom)

   ! Initialization

   adjmat(:, :) = .false.

   ! Set adjacency radii

   adjrad(:) = covrad(znums) + 0.25*(vdwrad(znums) - covrad(znums))

   ! Change adjacency matrix i,j to true if atoms i and j are closer
   ! than the sum of their adjacency radius

   do i = 1, natom
      do j = i + 1, natom
         atomdist = sqrt(sum((coords(:, i) - coords(:, j))**2))
         if (atomdist < adjrad(i) + adjrad(j)) then
            adjmat(i, j) = .true.
            adjmat(j, i) = .true.
         end if
      end do
   end do

end subroutine

subroutine adjmat2list(natom, adjmat, nadj, adjlist)
   integer, intent(in) :: natom
   logical, dimension(:, :), intent(in) :: adjmat
   integer, dimension(:), intent(out) :: nadj
   integer, dimension(:, :), intent(out) :: adjlist
   integer :: i, j

   nadj(:) = 0

   do i = 1, natom
      do j = i + 1, natom
         if (adjmat(i, j)) then
            nadj(i) = nadj(i) + 1
            nadj(j) = nadj(j) + 1
            if (nadj(i) > maxcoord .or. nadj(j) > maxcoord) then
               write (stderr, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlist(nadj(i), i) = j
            adjlist(nadj(j), j) = i
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

subroutine adjlist2bonds(natom, nadj, adjlist, nbond, bonds)
   integer, intent(in) :: natom
   integer, dimension(:), intent(in) :: nadj
   integer, dimension(:, :), intent(in) :: adjlist
   integer, intent(out) :: nbond
   integer, dimension(:, :), intent(out) :: bonds
   integer :: i, j

   nbond = 0

   do i = 1, natom
      do j =1, nadj(i)
         nbond = nbond + 1
         bonds(1, nbond) = i
         bonds(2, nbond) = adjlist(j, i)
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

function adjacencydelta(nadj0, adjlist0, adjmat1, mapping, k, l) result(delta)
   integer, intent(in) :: k, l
   integer, dimension(:), intent(in) :: mapping, nadj0
   integer, dimension(:, :), intent(in) :: adjlist0
   logical, dimension(:, :), intent(in) :: adjmat1
   integer i, nkk, nkl, nll, nlk, delta

   nkk = 0
   nkl = 0

   do i = 1, nadj0(k)
      if (adjlist0(i, k) /= l) then
         if (adjmat1(mapping(k), mapping(adjlist0(i, k)))) nkk = nkk + 1
         if (adjmat1(mapping(l), mapping(adjlist0(i, k)))) nkl = nkl + 1
      end if
   end do

   nll = 0
   nlk = 0

   do i = 1, nadj0(l)
      if (adjlist0(i, l) /= k) then
         if (adjmat1(mapping(l), mapping(adjlist0(i, l)))) nll = nll + 1
         if (adjmat1(mapping(k), mapping(adjlist0(i, l)))) nlk = nlk + 1
      end if
   end do

!        dkk = nadj0(k) + nadj1(mapping(k)) - 2*nkk
!        dll = nadj0(l) + nadj1(mapping(l)) - 2*nll
!        dkl = nadj0(k) + nadj1(mapping(l)) - 2*nkl
!        dlk = nadj0(l) + nadj1(mapping(k)) - 2*nlk
!        delta = dkl + dlk - dkk - dll

   ! Notice that dkl + dlk - dkk - dll == 2*(nkk + nll - nkl - nlk)
   delta = 2*(nkk + nll - nkl - nlk)

end function

end module
