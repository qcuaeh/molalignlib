! MolAlignLib
! Copyright (C) 2022 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module lap_jvc_dense
use parameters
implicit none
contains

subroutine jvc_dense(costs, n, rowsol, lapcost)
!---------------------------------------------------------------------------------------
! Solves the Linear Assignment Problem using the JVC algorithm for dense square matrices
! Finds the assignment that minimizes the total costs
!
! Input:
!   costs(n,n) - Cost matrix
!   n          - Problem dimension
!
! Output:
!   rowsol(n)  - Solution vector where rowsol(i) = j means row i is assigned to column j
!   lapcost    - Total costs of the optimal assignment
!---------------------------------------------------------------------------------------

   integer, intent(in) :: n
   real(rk), intent(in) :: costs(:,:)
   integer, intent(out) :: rowsol(:)
   real(rk), intent(out) :: lapcost
   
   ! Local variables
   integer :: i, j, j1, j2, k, f
   integer :: imin, numfree, prvnumfree, i0, freerow
   integer :: last, low, up, endofpath, loopcnt
   real(rk) :: dmin, h, umin, usubmin, v2
   real(rk), parameter :: resolution = 1.0e-10_rk
   logical :: unassignedfound
   
   ! Working arrays
   real(rk), allocatable :: v(:), d(:), u(:)
   integer, allocatable :: free(:), collist(:), matches(:)
   integer, allocatable :: pred(:), colsol(:)
   
   ! Allocate working arrays
   allocate(v(n), d(n), u(n))
   allocate(free(n), collist(n), matches(n), pred(n), colsol(n))
   
   ! Initialize
   v = 0.0_rk
   u = 0.0_rk
   matches = 0
   rowsol = -1
   colsol = -1
   
   !---------------------------------------------------
   ! COLUMN REDUCTION
   !---------------------------------------------------
   do j = n, 1, -1  ! Reverse order gives better results
      ! Find minimum costs over rows
      dmin = costs(1, j)
      imin = 1
      do i = 2, n
         if (costs(i, j) < dmin) then
            dmin = costs(i, j)
            imin = i
         end if
      end do
      
      v(j) = dmin
      matches(imin) = matches(imin) + 1
      
      if (matches(imin) == 1) then
         ! Initialize assignment if minimum row assigned for first time
         rowsol(imin) = j
         colsol(j) = imin
      else if (v(j) < v(rowsol(imin))) then
         j1 = rowsol(imin)
         rowsol(imin) = j
         colsol(j) = imin
         colsol(j1) = -1
      else
         colsol(j) = -1  ! Row already assigned, column not assigned
      end if
   end do
   
   !---------------------------------------------------
   ! REDUCTION TRANSFER
   !---------------------------------------------------
   numfree = 0
   do i = 1, n
      if (matches(i) == 0) then
         ! Fill list of unassigned 'free' rows
         numfree = numfree + 1
         free(numfree) = i
      else if (matches(i) == 1) then
         ! Transfer reduction from rows that are assigned once
         j1 = rowsol(i)
         dmin = huge(1.0_rk)
         do j = 1, n
            if (j /= j1) then
               h = costs(i, j) - v(j)
               if (h < dmin) dmin = h
            end if
         end do
         v(j1) = v(j1) - dmin
      end if
   end do
   
   !---------------------------------------------------
   ! AUGMENTING ROW REDUCTION (done exactly twice)
   !---------------------------------------------------
   loopcnt = 0
   do while (loopcnt < 2)
      loopcnt = loopcnt + 1
      k = 1
      prvnumfree = numfree
      numfree = 0
      
      do while (k <= prvnumfree)
         i = free(k)
         k = k + 1
         
         ! Find minimum and second minimum reduced costs over columns
         umin = costs(i, 1) - v(1)
         j1 = 1
         usubmin = huge(1.0_rk)
         
         do j = 2, n
            h = costs(i, j) - v(j)
            if (h < usubmin) then
               if (h < umin) then
                  usubmin = umin
                  j2 = j1
                  umin = h
                  j1 = j
               else
                  usubmin = h
                  j2 = j
               end if
            end if
         end do
         
         i0 = colsol(j1)
         
         if ((usubmin - umin) > resolution) then
            ! Increase minimum reduced costs to subminimum
            v(j1) = v(j1) - (usubmin - umin)
         else if (i0 > 0) then
            ! Minimum and subminimum equal, swap columns if j1 is assigned
            j1 = j2
            i0 = colsol(j2)
         end if
         
         ! (Re-)assign i to j1, possibly de-assigning i0
         rowsol(i) = j1
         colsol(j1) = i
         
         if (i0 > 0) then
            if ((usubmin - umin) > resolution) then
               ! Continue augmenting path with i0
               k = k - 1
               free(k) = i0
            else
               ! No further augmenting reduction possible
               numfree = numfree + 1
               free(numfree) = i0
            end if
         end if
      end do
   end do
   
   !---------------------------------------------------
   ! AUGMENT SOLUTION FOR EACH FREE ROW
   !---------------------------------------------------
   do f = 1, numfree
      freerow = free(f)
      
      ! Dijkstra shortest path algorithm
      do j = 1, n
         d(j) = costs(freerow, j) - v(j)
         pred(j) = freerow
         collist(j) = j
      end do
      
      low = 1
      up = 1
      unassignedfound = .FALSE.
      
      do while (.not. unassignedfound)
         if (up == low) then
            ! No more columns to scan, find new minimum
            last = low - 1
            dmin = d(collist(up))
            up = up + 1
            
            do k = up, n
               j = collist(k)
               h = d(j)
               if (h < dmin .or. abs(h - dmin) < resolution) then
                  if (h < dmin) then
                     up = low
                     dmin = h
                  end if
                  ! Swap collist elements
                  collist(k) = collist(up)
                  collist(up) = j
                  up = up + 1
               end if
            end do
            
            ! Check for unassigned column
            do k = low, up - 1
               if (colsol(collist(k)) < 0) then
                  endofpath = collist(k)
                  unassignedfound = .TRUE.
                  exit
               end if
            end do
         end if
         
         if (.not. unassignedfound) then
            ! Update distances via next scanned column
            j1 = collist(low)
            low = low + 1
            i = colsol(j1)
            h = costs(i, j1) - v(j1) - dmin
            
            do k = up, n
               j = collist(k)
               v2 = costs(i, j) - v(j) - h
               if (v2 < d(j)) then
                  pred(j) = i
                  if (abs(v2 - dmin) < resolution) then
                     if (colsol(j) < 0) then
                        endofpath = j
                        unassignedfound = .TRUE.
                        exit
                     else
                        ! Add to scan list
                        collist(k) = collist(up)
                        collist(up) = j
                        up = up + 1
                     end if
                  end if
                  d(j) = v2
               end if
            end do
         end if
      end do
      
      ! Update column prices
      do k = 1, last
         j1 = collist(k)
         v(j1) = v(j1) + d(j1) - dmin
      end do
      
      ! Update assignments along alternating path
      do
         i = pred(endofpath)
         colsol(endofpath) = i
         j1 = endofpath
         endofpath = rowsol(i)
         rowsol(i) = j1
         if (i == freerow) exit
      end do
   end do
   
   !---------------------------------------------------
   ! Calculate total costs
   !---------------------------------------------------
   lapcost = 0.0_rk
   do i = 1, n
      j = rowsol(i)
      u(i) = costs(i, j) - v(j)
      lapcost = lapcost + costs(i, j)
   end do
   
   ! Cleanup
   deallocate(v, d, u, free, collist, matches, pred, colsol)
end subroutine jvc_dense

end module lap_jvc_dense
