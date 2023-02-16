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

   ! Quick return

   if (.not. bond_flag) return

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
               write (error_unit, '("Maximum coordination number exceeded!")')
               stop
            end if
            adjlist(nadj(i), i) = j
            adjlist(nadj(j), j) = i
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
      
end module
