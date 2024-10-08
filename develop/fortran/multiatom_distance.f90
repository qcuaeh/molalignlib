real(rk) function equivdist(natom, i, j, maxlevel, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights)
   integer, intent(in) :: natom, i, j, maxlevel
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights

   real(rk) :: totdist, totweight
   logical, dimension(natom) :: mapped0, mapped1

   totdist = 0
   totweight = 0
   mapped0(:) = .false.
   mapped1(:) = .false.

   call recursivemap(i, j, 1, maxlevel, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, totdist, totweight, mapped0, mapped1)
   equivdist = totdist/totweight

end function

recursive subroutine recursivemap(i, j, level, maxlevel, nadjmna0, adjmnalen0, adjmnalist0, coords0, &
      nadjmna1, adjmnalen1, adjmnalist1, coords1, weights, totdist, totweight, mapped0, mapped1)
   integer, intent(in) :: i, j, level, maxlevel
   integer, dimension(:, :), intent(in) :: nadjmna0, nadjmna1
   integer, dimension(:, :, :), intent(in) :: adjmnalen0, adjmnalen1
   integer, dimension(:, :, :), intent(in) :: adjmnalist0, adjmnalist1
   real(rk), dimension(:, :), intent(in) :: coords0, coords1
   real(rk), dimension(:), intent(in) :: weights
   real(rk), intent(out) :: totdist, totweight
   logical, dimension(:), intent(inout) :: mapped0, mapped1

   integer :: h, k, l, m, n, offset
   integer, dimension(maxcoord) :: mapping, idx0, idx1
   real(rk) :: distmat(maxcoord, maxcoord)
   real(rk) :: dummy

!   print *, '>>>', i, j, ':', level, maxlevel - level
!!   print *, i, ':', adjmnalist0(:4, i, maxlevel - level)
!!   print *, j, ':', adjmnalist1(:4, j, maxlevel - level)
   mapped0(i) = .true.
   mapped1(j) = .true.
   totdist = totdist + weights(i)*sum((coords0(:, i) - coords1(:, j))**2)
   totweight = totweight + weights(i)
   if (level < maxlevel) then
      offset = 0
!      print *, 'num:', nadjmna0(i, maxlevel - level), &
!         nadjmna0(i, maxlevel - level) == nadjmna1(j, maxlevel - level)
      do h = 1, nadjmna0(i, maxlevel - level)
!         print *, 'len:', adjmnalen0(h, i, maxlevel - level), &
!            adjmnalen0(h, i, maxlevel - level) == adjmnalen1(h, j, maxlevel - level)
         m = 0
         do k = offset + 1, offset + adjmnalen0(h, i, maxlevel - level)
            if (.not. mapped0(adjmnalist0(k, i, maxlevel - level))) then
!!               print *, 'not mapped0:', adjmnalist0(k, i, maxlevel - level)
               m = m + 1
               idx0(m) = adjmnalist0(k, i, maxlevel - level)
               n = 0
               do l = offset + 1, offset + adjmnalen1(h, j, maxlevel - level)
                  if (.not. mapped1(adjmnalist1(l, j, maxlevel - level))) then
!!                     print *, 'not mapped1:', adjmnalist1(l, j, maxlevel - level)
                     n = n + 1
                     idx1(n) = adjmnalist1(l, j, maxlevel - level)
                     distmat(m, n) = sum((coords0(:, idx0(m)) - coords1(:, idx1(n)))**2)
!!                  else
!!                     print *, 'mapped1:', adjmnalist1(l, j, maxlevel - level)
                  end if
               end do
!!            else
!!               print *, 'mapped0:', adjmnalist0(k, i, maxlevel - level)
            end if
         end do
         if (m > 0) then
            call assndx(1, distmat, m, m, mapping, dummy)
            do n = 1, m
               call recursivemap(idx0(n), idx1(mapping(n)), level + 1, maxlevel, &
                  nadjmna0, adjmnalen0, adjmnalist0, coords0, nadjmna1, adjmnalen1, adjmnalist1, &
                  coords1, weights, totdist, totweight, mapped0, mapped1)
            end do
         end if
         offset = offset + adjmnalen0(h, i, maxlevel - level)
      end do
   end if
!   print *, '<<<'

end subroutine