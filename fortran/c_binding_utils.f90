module c_binding_utils
   use parameters   ! ik == c_int, rk == c_double - no iso_c_binding needed
   use molecule
   use euclidean
   use str_utils
   use chemdata
   use permutation
   use assorting
   use pruning_atoms
   use recording
   use alignment_atoms
   use options
   implicit none

contains

! ---------------------------------------------------------------------------
! Build an atom_t array from a packed C atomdata array and a coords array.
!
! atomdata layout (length n*2):
!   [elnum0, label0, elnum1, label1, ...]
!
! coords_in is row-major from C: [x0,y0,z0, x1,y1,z1, ...] (length n*3).
! Because ik==c_int and rk==c_double, assignments are direct - no real/kind
! conversion is required.
! ---------------------------------------------------------------------------
subroutine build_atoms(n, atomdata, coords_in, atoms)
   integer(ik), intent(in), value :: n
   integer(ik), dimension(n*2), intent(in) :: atomdata
   real(rk),    dimension(n*3), intent(in) :: coords_in
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   integer(ik) :: i, abase, cbase

   allocate(atoms(n))
   do i = 1, n
      abase = (i - 1)*2
      cbase = (i - 1)*3
      atoms(i)%elnum     = atomdata(abase + 1)
      atoms(i)%group     = atomdata(abase + 2)
      atoms(i)%coords(1) = coords_in(cbase + 1)
      atoms(i)%coords(2) = coords_in(cbase + 2)
      atoms(i)%coords(3) = coords_in(cbase + 3)
   end do
end subroutine build_atoms

! ---------------------------------------------------------------------------
! Build a bond_t array from a flat C array.
!
! bonddata layout (length n_bonds*3):
!   [atom1_0, atom2_0, type_0,  atom1_1, atom2_1, type_1, ...]
! Atom indices are 1-based (as stored in mol/sdf files).
! ---------------------------------------------------------------------------
subroutine build_bonds(n, bonddata, bonds)
   integer(ik), intent(in), value :: n
   integer(ik), dimension(n*3), intent(in) :: bonddata
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   integer(ik) :: i, base

   allocate(bonds(n))
   do i = 1, n
      base = (i - 1)*3
      bonds(i)%atomidx1 = bonddata(base + 1)
      bonds(i)%atomidx2 = bonddata(base + 2)
      bonds(i)%bondtype = bonddata(base + 3)
   end do
end subroutine build_bonds

! ---------------------------------------------------------------------------
! Build a row-major 4x4 homogeneous transformation matrix.
! Maps a point in molecule-2 frame to molecule-1 frame:
!   p_out = R * p_in + t,   t = center1 - R * center2
!
! Layout (row-major, 1-based Fortran indexing):
!   [1..4]  -> row 0: R(1,*), tx
!   [5..8]  -> row 1: R(2,*), ty
!   [9..12] -> row 2: R(3,*), tz
!   [13..16]-> 0 0 0 1
! ---------------------------------------------------------------------------
subroutine build_homogeneous_transform(rotquat, center1, center2, htrans)
   real(rk), intent(in) :: rotquat(4)
   real(rk), intent(in) :: center1(3), center2(3)
   real(rk), intent(out) :: htrans(16)
   real(rk) :: R(3,3), t(3)

   R = quatrotmat(rotquat)
   t = center1 - matmul(R, center2)

   ! ik==c_int, rk==c_double - direct assignment, no real() conversion needed
   htrans(1)  = R(1,1);  htrans(2)  = R(1,2);  htrans(3)  = R(1,3);  htrans(4)  = t(1)
   htrans(5)  = R(2,1);  htrans(6)  = R(2,2);  htrans(7)  = R(2,3);  htrans(8)  = t(2)
   htrans(9)  = R(3,1);  htrans(10) = R(3,2);  htrans(11) = R(3,3);  htrans(12) = t(3)
   htrans(13) = 0.0_rk;  htrans(14) = 0.0_rk;  htrans(15) = 0.0_rk;  htrans(16) = 1.0_rk
end subroutine build_homogeneous_transform

! ---------------------------------------------------------------------------
! Return the 4x4 identity matrix.
! ---------------------------------------------------------------------------
subroutine set_identity_transform(htrans)
   real(rk), intent(out) :: htrans(16)
   htrans     = 0.0_rk
   htrans(1)  = 1.0_rk;  htrans(6)  = 1.0_rk
   htrans(11) = 1.0_rk;  htrans(16) = 1.0_rk
end subroutine set_identity_transform

end module c_binding_utils