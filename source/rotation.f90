module rotation

use options

implicit none

private
public rotate
public rotated
public quat2mat
public quatmul
public rotangle
public getrotmat
public getrotquat

interface rotate
    module procedure matrotate
    module procedure quatrotate
end interface

interface rotated
    module procedure matrotated
    module procedure quatrotated
end interface

contains

function quat2mat(rotquat) result(rotmat)
! Purpose: Apply rotation quaternion to atomic coordinates

! rotquat: Rotation quaternion
! rotmat: Rotation matrix

    real(wp), intent(in) :: rotquat(4)
    real(wp) rotmat(3, 3)

! Calculate the rotation matrix

    rotmat(1, 1) = rotquat(1)**2 + rotquat(2)**2 - rotquat(3)**2 - rotquat(4)**2
    rotmat(1, 2) = 2*(rotquat(2)*rotquat(3) + rotquat(1)*rotquat(4))
    rotmat(1, 3) = 2*(rotquat(2)*rotquat(4) - rotquat(1)*rotquat(3))
    rotmat(2, 1) = 2*(rotquat(2)*rotquat(3) - rotquat(1)*rotquat(4))
    rotmat(2, 2) = rotquat(1)**2 + rotquat(3)**2 - rotquat(2)**2 - rotquat(4)**2
    rotmat(2, 3) = 2*(rotquat(3)*rotquat(4) + rotquat(1)*rotquat(2))
    rotmat(3, 1) = 2*(rotquat(2)*rotquat(4) + rotquat(1)*rotquat(3))
    rotmat(3, 2) = 2*(rotquat(3)*rotquat(4) - rotquat(1)*rotquat(2))
    rotmat(3, 3) = rotquat(1)**2 + rotquat(4)**2 - rotquat(2)**2 - rotquat(3)**2

end function

function quatrotated(natom, rotquat, coords)
! Purpose: Apply rotation quaternion to atomic coordinates
    integer, intent(in) :: natom
    real(wp), intent(in) :: rotquat(4)
    real(wp), intent(in) :: coords(3, natom)
    real(wp) quatrotated(3, natom)

    quatrotated = matrotated(natom, quat2mat(rotquat), coords)

end function

function matrotated(natom, rotmat, coords) result(rotated)
! Purpose: Apply rotation matrix to atomic coordinates
    integer, intent(in) :: natom
    real(wp), intent (in) :: rotmat(3, 3)
    real(wp), intent (in) :: coords(3, natom)
    real(wp) rotated(3, natom)

    rotated = matmul(rotmat, coords)

end function

subroutine matrotate(natom, rotmat, coords)
! Purpose: Apply rotation matrix to atomic coordinates
    integer, intent(in) :: natom
    real(wp), intent (in) :: rotmat(3, 3)
    real(wp), intent (inout) :: coords(3, natom)

    integer i, k
    real(wp) vecrot(3)

    coords = matmul(rotmat, coords)

end subroutine

subroutine quatrotate(natom, rotquat, coords)
! Purpose: Apply rotation quaternion to atomic coordinates
    integer, intent(in) :: natom
    real(wp), intent(in) :: rotquat(4)
    real(wp), intent(inout) :: coords(3, natom)

! Rotate atomic coordinates

!    call matrotate(natom, quat2mat(rotquat), coords)
     coords = matmul(quat2mat(rotquat), coords)

end subroutine

function quatmul(p, q)
    real(wp), dimension(4), intent(in) :: p, q
    real(wp) quatmul(4)

    real(wp) left(4, 4)

    left(:, 1) = [ p(1), -p(2), -p(3), -p(4) ]
    left(:, 2) = [ p(2),  p(1), -p(4),  p(3) ]
    left(:, 3) = [ p(3),  p(4),  p(1), -p(2) ]
    left(:, 4) = [ p(4), -p(3),  p(2),  p(1) ]

    quatmul = matmul(left, q)

end function

function rotangle(q)
    real(wp), dimension(4), intent(in) :: q
    real(wp) rotangle

!    rotangle = 2*acos(q(1))
!    rotangle = 2*asin(sqrt(sum(q(2:4)**2)))
    rotangle = 2*abs(atan(sqrt(sum(q(2:4)**2))/q(1)))
!    rotangle = 2*atan2(sqrt(sum(q(2:4)**2)), q(1))

!    print *, sum(q**2), &
!        90./asin(1.)*2*acos(q(1)), &
!        90./asin(1.)*2*asin(sqrt(sum(q(2:4)**2))), &
!        90./asin(1.)*2*atan(sqrt(sum(q(2:4)**2))/q(1)), &
!        90./asin(1.)*2*atan2(sqrt(sum(q(2:4)**2)), q(1))
end function

function getrotquat(x) result(rotquat)
! Purpose: Generate a random rotation quaternion.
! Reference: Academic Press Graphics Gems Series archive Graphics
!            Gems III archive. Pages: 129 - 132.

    real(wp), intent(in) :: x(3)
    real(wp) rotquat(4)

    real(wp) pi, a1, a2, r1, r2, s1, s2, c1, c2

! Calculate auxiliar vectors and constants

    pi = 2*asin(1.)
    a1 = 2*pi*x(1)
    a2 = 2*pi*x(2)
    r1 = sqrt(1.0 - x(3))
    r2 = sqrt(x(3))

    s1 = sin(a1)
    c1 = cos(a1)
    s2 = sin(a2)
    c2 = cos(a2)

!   Quaternion (w, x, y, z), w must be first as required by quat2mat
    rotquat = [ c2*r2, s1*r1, c1*r1, s2*r2 ]

!    print *, sum(rotquat**2)

end function

function getrotmat(x) result(rotmat)
! Purpose: Generate a random rotation matrix.
! Reference: Academic Press Graphics Gems Series archive Graphics
!            Gems III archive. Pages: 117 - 120.

    real(wp), intent(in) :: x(3)
    real(wp) rotmat(3, 3)

    integer i, j, k
    real(wp) pi, a1, a2, r1, r2, s1, c1, s2, c2
    real(wp) v(3), right(3, 3), left(3, 3)

! Calculate auxiliar vectors and constants

    pi = 2*asin(1.)
    a1 = 2*pi*x(1)
    a2 = 2*pi*x(2)
    r1 = sqrt(1.0 - x(3))
    r2 = sqrt(x(3))

    s1 = sin(a1)
    c1 = cos(a1)
    s2 = sin(a2)
    c2 = cos(a2)

    right(1, :) = [ c1, s1, 0.0_wp ]
    right(2, :) = [ -s1, c1, 0.0_wp ]
    right(3, :) = [ 0.0_wp, 0.0_wp, 1.0_wp ]

    v = [ c2*r2, s2*r2, r1 ]

    do j = 1, 3
        left(:, j) = v(:)*v(j)
    end do

    left = 2*left

    do i = 1, 3
        left(i, i) = left(i, i) - 1.0_wp
    end do

! Calculate the rotation matrix

    rotmat = matmul(left, right)

!    print *, det33(rotmat)

end function

end module