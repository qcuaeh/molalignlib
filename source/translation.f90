module translation

use options

implicit none

private
public translate
public translated
public centroid

contains

function centroid(natom, weights, coords) result(vector)
! Purpose: Get coordinates of centroid
! natom: Number of atoms
! znum: Atomic numbers
! coords: Atomic coordinates

    integer, intent(in) :: natom
    real(wp), intent(in) :: weights(natom), coords(3, natom)
    real(wp) vector(3)

    integer i

! Calculate the coordinates of the center of mass

    vector(:) = 0

    do i = 1, natom
        vector = vector + weights(i)*coords(:, i)
    end do

    vector = vector/sum(weights)

end function

function translated(natom, vector, coords)
! Purpose: Translate atomic coordinates to its vector of geometry

! natom: Number of atoms
! vector: Transalte vector
! coords: Atomic coordinates

    integer, intent(in) :: natom
    real(wp), intent(in) :: coords(3, natom), vector(3)
    real(wp) translated(3, natom)

    integer i

    do i = 1, natom
        translated(:, i) = coords(:, i) - vector(:)
    end do

end function

subroutine translate(natom, vector, coords)
! Purpose: Move the coords to the vector of mass
! natom: Number of atoms
! vector: Transalte vector
! coords: Atomic coordinates

    integer, intent(in) :: natom
    real(wp), intent(in) :: vector(3)
    real(wp), intent(inout) :: coords(3, natom)

    integer i

    do i = 1, natom
        coords(:, i) = coords(:, i) - vector(:)
    end do

end subroutine

end module
