module math
use parameters

implicit none

private
public trace
public det33
public matmul2
public union
public intersection

contains

real(wp) function trace(natom, matrix)
   integer, intent(in) :: natom
   real(wp), dimension(:, :), intent(in) :: matrix
   integer :: i

   trace = 0

   do i = 1, natom
      trace = trace + matrix(i, i)
   end do

end function

real(wp) function det33(a) result(det)
   real(wp), dimension(3,3), intent(in) :: a

   det = a(1,1)*a(2,2)*a(3,3)  &
       - a(1,1)*a(2,3)*a(3,2)  &
       - a(1,2)*a(2,1)*a(3,3)  &
       + a(1,2)*a(2,3)*a(3,1)  &
       + a(1,3)*a(2,1)*a(3,2)  &
       - a(1,3)*a(2,2)*a(3,1)

end function

function matmul2(a, b, m, o, n) result(ab)
   integer, intent(in) :: m, o ,n
   real(wp), intent (in) :: a(m, o)
   real(wp), intent (in) :: b(o, n)
   real(wp) :: ab(m, n)
   integer :: i, j

   do i = 1, n
      ab(:, i) = 0.0
      do j = 1, o
         ab(:, i) = ab(:, i) + a(:, j)*b(j, i)
      end do
   end do

end function

! Find the union of two sorted lists without repeated elements
subroutine union(n1, list1, n2, list2, n3, list3)
   integer, intent(in) :: n1, n2, list1(:), list2(:)
   integer, intent(out) :: n3, list3(:)
   integer :: i1, i2

   i1 = 1
   i2 = 1
   n3 = 0

   do while (i1 <= n1 .and. i2 <= n2)
      n3 = n3 + 1
      if (list1(i1) < list2(i2)) then
         list3(n3) = list1(i1)
         i1 = i1 + 1
      else if (list1(i1) > list2(i2)) then
         list3(n3) = list2(i2)
         i2 = i2 + 1
      else
         list3(n3) = list2(i2)
         i1 = i1 + 1
         i2 = i2 + 1
      end if
   end do

   do while (i1 <= n1)
      n3 = n3 + 1
      list3(n3) = list1(i1)
      i1 = i1 + 1
   end do

   do while (i2 <= n2)
      n3 = n3 + 1
      list3(n3) = list2(i2)
      i2 = i2 + 1
   end do

end subroutine

! Find the intersection of two sorted lists without repeated elements
subroutine intersection(n1, list1, n2, list2, n3, list3)
   integer, intent(in) :: n1, n2, list1(:), list2(:)
   integer, intent(out) :: n3, list3(:)
   integer :: i1, i2

   i1 = 1
   i2 = 1
   n3 = 0

   do while (i1 <= n1 .and. i2 <= n2)
      if (list1(i1) < list2(i2)) then
         i1 = i1 + 1
      else if (list1(i1) > list2(i2)) then
         i2 = i2 + 1
      else
         n3 = n3 + 1
         list3(n3) = list2(i2)
         i1 = i1 + 1
         i2 = i2 + 1
      end if
   end do

end subroutine

end module
