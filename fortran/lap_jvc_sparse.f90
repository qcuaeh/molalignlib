! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

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

module lap_jvc_sparse
use parameters
implicit none
contains

subroutine jvc_sparse(n, sz, cc, kk, first, rowsol, lapcost, ierr)
!-----------------------------------------------------------------------------
! Solves sparse Linear Assignment Problem using JVC algorithm
! Works with real costs (no integer scaling needed)
!
! Input:
!   n         - Number of rows and columns (assumes square matrix)
!   sz        - Total number of nonzero elements
!   cc(sz)    - Cost values (real)
!   kk(sz)    - Column indices (1-based)
!   first(n+1)- Row pointers: row i has elements from first(i) to first(i+1)-1
!
! Output:
!   rowsol(n) - Assignment: rowsol(i) = j means row i assigned to column j
!   lapcost   - Total cost of assignment
!   ierr      - Error flag: 0=success, 1=failure
!-----------------------------------------------------------------------------

   integer(ik), intent(in) :: n, sz
   real(rk), intent(in) :: cc(sz)
   integer(ik), intent(in) :: kk(sz), first(n+1)
   integer(ik), intent(out) :: rowsol(n)
   real(rk), intent(out) :: lapcost
   integer(ik), intent(out) :: ierr
   
   ! Parameters
   real(rk), parameter :: infValue = 1.0e30_rk
   real(rk), parameter :: resolution = 1.0e-10_rk
   
   ! Working arrays
   integer(ik) :: x(n), y(n), free(n), todo(n), lab(n), ok(n)
   integer(ik) :: xinv(n)
   real(rk) :: u(n), v(n), d(n)
   
   ! Local variables
   integer(ik) :: i, jp, j0p, j1p, i0, t, tp, h, lp, l0p, tel
   integer(ik) :: l, l0, td1
   real(rk) :: min_val, v0, vj, dj
   logical(lk) :: fail
   
   ! Initialize (0 means unassigned in Fortran)
   ierr = 0
   x = 0
   y = 0
   free = 0
   todo = 0
   xinv = 0
   v = infValue
   u = 0.0_rk
   
   !---------------------------------------------------
   ! INITIALIZATION (only for square matrices)
   !---------------------------------------------------
   
   ! Find column minima
   do i = 1, n
      do t = first(i), first(i+1) - 1
         jp = kk(t)
         if (cc(t) < v(jp)) then
            v(jp) = cc(t)
            y(jp) = i
         end if
      end do
   end do
   
   ! Build initial assignment
   do jp = n, 1, -1
      i = y(jp)
      if (i == 0) then
         ! Column has no entries - check if problem is feasible
         ! This can happen with sparse pruned matrices
         cycle  ! Skip this column for now, will be handled in augmentation
      end if
      if (x(i) == 0) then
         x(i) = jp
      else
         y(jp) = 0
         xinv(i) = 1
      end if
   end do
   
   ! Build free list and compute dual variables
   lp = 0
   do i = 1, n
      if (xinv(i) /= 0) cycle
      
      if (x(i) /= 0) then
         min_val = infValue
         j0p = x(i)
         do t = first(i), first(i+1) - 1
            jp = kk(t)
            if (jp /= j0p) then
               if (cc(t) - v(jp) < min_val) then
                  min_val = cc(t) - v(jp)
               end if
            end if
         end do
         u(i) = min_val
         
         ! Find position of j0p in sparse structure
         tp = first(i)
         do while (kk(tp) /= j0p)
            tp = tp + 1
         end do
         v(j0p) = cc(tp) - min_val
      else
         lp = lp + 1
         free(lp) = i
      end if
   end do
   
   !---------------------------------------------------
   ! AUGMENTING ROW REDUCTION (done twice)
   !---------------------------------------------------
   do tel = 1, 2
      h = 1
      l0p = lp
      lp = 0
      
      do while (h <= l0p)
         i = free(h)
         h = h + 1
         
         ! Find minimum and second minimum
         j0p = -1
         j1p = -1
         v0 = infValue
         vj = infValue
         
         do t = first(i), first(i+1) - 1
            jp = kk(t)
            dj = cc(t) - v(jp)
            if (dj < vj) then
               if (dj >= v0) then
                  vj = dj
                  j1p = jp
               else
                  vj = v0
                  v0 = dj
                  j1p = j0p
                  j0p = jp
               end if
            end if
         end do
         
         if (j0p < 0) then
            ierr = 1
            return
         end if
         
         i0 = y(j0p)
         u(i) = vj
         
         if (v0 < vj) then
            v(j0p) = v(j0p) + (v0 - vj)
         else if (i0 /= 0) then
            j0p = j1p
            i0 = y(j0p)
         end if
         
         x(i) = j0p
         y(j0p) = i
         
         if (i0 /= 0) then
            if (v0 < vj) then
               h = h - 1
               free(h) = i0
            else
               lp = lp + 1
               free(lp) = i0
            end if
         end if
      end do
   end do
   
   l0 = lp
   
   !---------------------------------------------------
   ! AUGMENT SOLUTION FOR EACH FREE ROW
   !---------------------------------------------------
   td1 = -1
   do l = 1, l0
      call solve_for_one_row(n, sz, cc, kk, first, l, free, &
                             x, y, u, v, d, ok, lab, todo, &
                             td1, infValue, resolution, fail)
      if (fail) then
         ierr = 1
         return
      end if
   end do
   
   !---------------------------------------------------
   ! Prepare output
   !---------------------------------------------------
   lapcost = 0.0_rk
   do i = 1, n
      rowsol(i) = x(i)
      jp = rowsol(i)
      
      if (jp == 0) then
         ! Row could not be assigned
         ierr = 1
         return
      end if
      
      ! Find cost in sparse matrix
      do t = first(i), first(i+1) - 1
         if (kk(t) == jp) then
            lapcost = lapcost + cc(t)
            exit
         end if
      end do
   end do
end subroutine jvc_sparse

subroutine solve_for_one_row(n, sz, cc, kk, first, l, free, &
                             x, y, u, v, d, ok, lab, todo, &
                             td1, infValue, resolution, fail)
!-----------------------------------------------------------------------------
! Solves augmenting path for one free row (Dijkstra's algorithm)
!-----------------------------------------------------------------------------
   integer(ik), intent(in) :: n, sz, l
   real(rk), intent(in) :: cc(sz), infValue, resolution
   integer(ik), intent(in) :: kk(sz), first(n+1), free(n)
   integer(ik), intent(inout) :: x(n), y(n), lab(n), todo(n), ok(n), td1
   real(rk), intent(inout) :: u(n), v(n), d(n)
   logical(lk), intent(out) :: fail
   
   integer(ik) :: i, i0, j, jp, t, tp, last, td2, hp
   real(rk) :: min_val, h, v2, dj
   
   fail = .FALSE.
   
   ! Initialize
   d = infValue
   ok = 0
   min_val = infValue
   i0 = free(l)
   
   ! Build initial todo list
   td1 = 0
   do t = first(i0), first(i0+1) - 1
      j = kk(t)
      dj = cc(t) - v(j)
      d(j) = dj
      lab(j) = i0
      if (dj <= min_val) then
         if (dj < min_val) then
            td1 = 0
            min_val = dj
         end if
         td1 = td1 + 1
         todo(td1) = j
      end if
   end do
   
   ! Check if any minimum column is unassigned
   do hp = 1, td1
      j = todo(hp)
      if (y(j) == 0) then
         call update_assignments(lab, y, x, j, i0)
         return
      end if
      ok(j) = 1
   end do
   
   td2 = n
   last = n
   
   ! Main loop
   do
      if (td1 < 1) then
         fail = .TRUE.
         return
      end if
      
      j = todo(td1)
      td1 = td1 - 1
      i = y(j)
      todo(td2) = j
      td2 = td2 - 1
      
      ! Find j in row i
      tp = first(i)
      do while (kk(tp) /= j)
         tp = tp + 1
      end do
      
      h = cc(tp) - v(j) - min_val
      
      ! Update distances
      do t = first(i), first(i+1) - 1
         j = kk(t)
         if (ok(j) == 0) then
            v2 = cc(t) - v(j) - h
            if (v2 < d(j)) then
               d(j) = v2
               lab(j) = i
               if (abs(v2 - min_val) < resolution) then
                  if (y(j) == 0) then
                     call update_dual(n, d, v, todo, last, min_val)
                     call update_assignments(lab, y, x, j, i0)
                     return
                  end if
                  td1 = td1 + 1
                  todo(td1) = j
                  ok(j) = 1
               end if
            end if
         end if
      end do
      
      ! Find new minimum if needed
      if (td1 == 0) then
         min_val = infValue
         last = td2 + 1
         do jp = 1, n
            if (ok(jp) == 0 .and. d(jp) < infValue) then
               if (d(jp) <= min_val) then
                  if (d(jp) < min_val) then
                     td1 = 0
                     min_val = d(jp)
                  end if
                  td1 = td1 + 1
                  todo(td1) = jp
               end if
            end if
         end do
         
         do hp = 1, td1
            j = todo(hp)
            if (y(j) == 0) then
               call update_dual(n, d, v, todo, last, min_val)
               call update_assignments(lab, y, x, j, i0)
               return
            end if
            ok(j) = 1
         end do
      end if
   end do
end subroutine solve_for_one_row

subroutine update_dual(n, d, v, todo, last, min_val)
   integer(ik), intent(in) :: n, last
   integer(ik), intent(in) :: todo(n)
   real(rk), intent(in) :: min_val
   real(rk), intent(inout) :: d(n), v(n)
   integer(ik) :: k, j
   
   do k = last, n
      j = todo(k)
      v(j) = v(j) + d(j) - min_val
   end do
end subroutine update_dual

subroutine update_assignments(lab, y, x, j, i0)
   integer(ik), intent(inout) :: lab(:), y(:), x(:)
   integer(ik), intent(inout) :: j
   integer(ik), intent(in) :: i0
   integer(ik) :: i, tmp
   
   i = lab(j)
   do
      y(j) = i
      tmp = j
      j = x(i)
      x(i) = tmp
      if (i == i0) exit
      i = lab(j)
   end do
end subroutine update_assignments

end module lap_jvc_sparse
