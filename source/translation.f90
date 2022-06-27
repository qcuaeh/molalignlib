module translation

use options

implicit none

private
public translate
public translated
public centroid

contains

function centroid(natom, weights, coords)
! Purpose: Get coordinates of centroid

    integer, intent(in) :: natom
    real(wp), intent(in) :: weights(natom), coords(3, natom)
    real(wp) centroid(3)

    integer i

! Calculate the coordinates of the center of mass

    centroid(:) = 0

    do i = 1, natom
        centroid(:) = centroid(:) + weights(i)*coords(:, i)
    end do

    centroid(:) = centroid(:)/sum(weights)

end function

function translated(natom, vector, coords)
! Purpose: Translate atomic coordinates to its vector of geometry

    integer, intent(in) :: natom
    real(wp), intent(in) :: coords(3, natom), vector(3)
    real(wp) translated(3, natom)

    integer i

    do i = 1, natom
        translated(:, i) = coords(:, i) + vector(:)
    end do

end function

subroutine translate(natom, vector, coords)
! Purpose: Move the coords to the vector of mass

    integer, intent(in) :: natom
    real(wp), intent(in) :: vector(3)
    real(wp), intent(inout) :: coords(3, natom)

    integer i

    do i = 1, natom
        coords(:, i) = coords(:, i) + vector(:)
    end do

end subroutine

end module
