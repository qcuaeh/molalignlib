module alignment

use linear
use rotation
use options

implicit none

private
public aligned
public squaredist
public leastsquaredist
public leastrotquat

contains

function squaredist(natom, weights, coords0, coords1, atomap) result(dist)
    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomap
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: coords0, coords1

    integer i
    real dist

    dist = 0

    do i = 1, natom
        dist = dist + weights(i)*sum((coords1(:, atomap(i)) - coords0(:, i))**2)
    end do

end function

!function leastsquaredist(natom, weights, coords0, coords1, atomap) result(dist)
! Calculate least square distance from eigenvalues
!    integer, intent(in) :: natom
!    integer, dimension(:), intent(in) :: atomap
!    real, dimension(:), intent(in) :: weights
!    real, dimension(:, :), intent(in) :: coords0, coords1
!
!    real dist, kearsleymat(4, 4), eigval(4)
!
!    call buildkearsleymat(natom, weights, coords0, coords1, atomap, kearsleymat)
!    call syeval4(kearsleymat, eigval)
!    dist = eigval(1)
!
!end function

function leastsquaredist(natom, weights, coords0, coords1, atomap) result(dist)
! Calculate least square distance from aligned coordinates
    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomap
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: coords0, coords1

    real dist

    dist = squaredist(natom, weights, coords0, aligned(natom, weights, coords0, coords1, atomap), atomap)

end function

function leastrotquat(natom, weights, coords0, coords1, atomap) result(quat)
    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomap
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: coords0, coords1

    real quat(4), kearsleymat(4, 4), eigval(4)

    call buildkearsleymat(natom, weights, coords0, coords1, atomap, kearsleymat)
    call syevec4(kearsleymat, eigval)
    quat = kearsleymat(:, 1)

end function

function aligned(natom, weights, coords0, coords1, atomap)
    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomap
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: coords0, coords1

    real aligned(3, natom), kearsleymat(4, 4), eigval(4)

    call buildkearsleymat(natom, weights, coords0, coords1, atomap, kearsleymat)
    call syevec4(kearsleymat, eigval)
    aligned = rotated(natom, coords1, kearsleymat(:, 1))

end function

subroutine buildkearsleymat(natom, weights, coords0, coords1, atomap, kearsleymat)
! Purpose: Find the best orientation by least squares minimization
! Reference: Acta Cryst. (1989). A45, 208-210

! weights: Atomic weights
! coords0: Reference atomic coordinates
! coords1: Atomic coordinates
! rotquat: Rotation quaternion vector
! ssd: Root mean square distance from least squares

    integer, intent(in) :: natom
    integer, dimension(:), intent(in) :: atomap
    real, dimension(:), intent(in) :: weights
    real, dimension(:, :), intent(in) :: coords0, coords1

    integer i
    real kearsleymat(4, 4), p(3, natom), q(3, natom), auxmat(4, 4)

    kearsleymat = 0.0

    do i = 1, natom
        p(:, i) = coords0(:, i) + coords1(:, atomap(i))
        q(:, i) = coords0(:, i) - coords1(:, atomap(i))
    end do

! Calculate upper diagonal elements of the matrix

    do i = 1, natom
        auxmat(1, 1) = q(1, i)**2 + q(2, i)**2 + q(3, i)**2
        auxmat(1, 2) = p(2, i)*q(3, i) - q(2, i)*p(3, i)
        auxmat(1, 3) = q(1, i)*p(3, i) - p(1, i)*q(3, i)
        auxmat(1, 4) = p(1, i)*q(2, i) - q(1, i)*p(2, i)
        auxmat(2, 2) = p(2, i)**2 + p(3, i)**2 + q(1, i)**2
        auxmat(2, 3) = q(1, i)*q(2, i) - p(1, i)*p(2, i)
        auxmat(2, 4) = q(1, i)*q(3, i) - p(1, i)*p(3, i)
        auxmat(3, 3) = p(1, i)**2 + p(3, i)**2 + q(2, i)**2
        auxmat(3, 4) = q(2, i)*q(3, i) - p(2, i)*p(3, i)
        auxmat(4, 4) = p(1, i)**2 + p(2, i)**2 + q(3, i)**2
        kearsleymat = kearsleymat + weights(i)*auxmat
    end do

end subroutine

end module
