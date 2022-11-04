module discrete

implicit none

contains

subroutine union(n1, list1, n2, list2, n3, list3)
! Find the union of two sorted lists without repeated elements

    integer, intent(in) :: n1, n2, list1(:), list2(:)
    integer, intent(out) :: n3, list3(:)

    integer i1, i2

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

subroutine intersection(n1, list1, n2, list2, n3, list3)
! Find the intersection of two sorted lists without repeated elements

    integer, intent(in) :: n1, n2, list1(:), list2(:)
    integer, intent(out) :: n3, list3(:)

    integer i1, i2

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
