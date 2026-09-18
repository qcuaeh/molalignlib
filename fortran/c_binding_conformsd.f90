module c_binding_conformsd
   use parameters
   use str_utils
   use molecule
   use euclidean
   use chemdata
   use permutation
   use assorting
   use recording
   use assignment_conformer
   use alignment_conformer
   use c_binding_utils
   use flags
   implicit none

contains

! C callable wrapper for conformsd functionality.
!
! Caller supplies pre-read molecular data as flat C arrays:
!   c_atom_data1/2 : packed atom data, length n_atoms*2:
!                    [elnum0, label0, elnum1, label1, ...]
!                    label = 0 means unlabelled.
!   coords1/2      : XYZ coordinates, row-major (n_atoms x 3), length n_atoms*3
!   n_bonds1/2     : number of bonds (may be 0 when bond_flag=1, i.e. derive from geometry)
!   c_bond_data1/2 : flat bond array, length n_bonds*3, layout: [atom1, atom2, type, ...]
!                    (1-based atom indices as in the original file)
!   c_bond_tol     : bond detection tolerance (only used, and required, when
!                    c_bond_flag is true; no default)
!
! Multiple ranked candidate solutions:
!   c_n_records requests up to that many ranked candidate solutions. Records
!   beyond the first are only ever produced when both c_align_flag and
!   c_remap_flag are true (the optimize_atomperm_conformer search); in every
!   other case exactly one record is written regardless of c_n_records.
!   c_occ_records reports how many were actually written; only the first
!   c_occ_records entries of c_rmsd_list, c_atomperm_list, and
!   c_transform_list are meaningful. All output arrays are flattened and
!   must be allocated by the caller with at least c_n_records elements per
!   record (natoms for c_atomperm_list, 16 for c_transform_list).
!
! c_error_code values (numbered to match atormsd_calculate where applicable):
!   0 = success
!   1 = not isomers
!   2 = atom type mismatch
!   3 = missing bonds
!   4 = bond connectivity mismatch (only possible when c_remap_flag is false)
subroutine conformsd_calculate(                                              &
      n_atoms1,  c_atom_data1,  c_coords1,                                    &
      n_bonds1,  c_bond_data1,                                                  &
      n_atoms2,  c_atom_data2,  c_coords2,                                    &
      n_bonds2,  c_bond_data2,                                                  &
      c_align_flag, c_remap_flag, c_heavy_flag, c_mass_flag,                &
      c_mirror_flag, c_label_flag, c_bond_flag, c_bond_tol,                 &
      c_print_stats, c_print_assigntree, c_random_flag,                      &
      c_conv_freq, c_max_trials,                                           &
      c_n_records,                                                          &
      c_rmsd_list, c_natoms, c_atomperm_list,                                &
      c_transform_list, c_occ_records, c_error_code)                         &
      bind(C, name="conformsd_calculate")

   ! Molecule 1
   integer(ik), intent(in), value :: n_atoms1
   integer(ik), dimension(n_atoms1*2), intent(in) :: c_atom_data1
   real(rk),    dimension(n_atoms1*3), intent(in) :: c_coords1
   integer(ik), intent(in), value :: n_bonds1
   integer(ik), dimension(n_bonds1*3), intent(in) :: c_bond_data1

   ! Molecule 2
   integer(ik), intent(in), value :: n_atoms2
   integer(ik), dimension(n_atoms2*2), intent(in) :: c_atom_data2
   real(rk),    dimension(n_atoms2*3), intent(in) :: c_coords2
   integer(ik), intent(in), value :: n_bonds2
   integer(ik), dimension(n_bonds2*3), intent(in) :: c_bond_data2

   ! Flags
   logical(lk), intent(in), value :: c_align_flag, c_remap_flag, c_heavy_flag, c_mass_flag
   logical(lk), intent(in), value :: c_mirror_flag, c_label_flag, c_bond_flag
   real(rk),    intent(in), value :: c_bond_tol
   logical(lk), intent(in), value :: c_print_stats, c_print_assigntree, c_random_flag
   integer(ik), intent(in), value :: c_conv_freq, c_max_trials

   ! Requested number of ranked records
   integer(ik), intent(in), value :: c_n_records

   ! Outputs
   real(rk),    dimension(*), intent(out) :: c_rmsd_list
   integer(ik),                intent(out) :: c_natoms
   integer(ik), dimension(*), intent(out) :: c_atomperm_list
   real(rk),    dimension(*), intent(out) :: c_transform_list
   integer(ik),                intent(out) :: c_occ_records
   integer(ik),                intent(out) :: c_error_code

   ! Local variables
   type(atom_t), dimension(:), allocatable :: atoms1, atoms2
   type(bond_t), dimension(:), allocatable :: bonds1, bonds2
   type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
   type(partition_t) :: atomtypes
   type(registry_t)  :: registry
   real(rk) :: rmsd, center1(3), center2(3), rotquat(4)
   real(rk), dimension(:),   allocatable :: weights1, weights2
   real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
   integer(ik), dimension(:), allocatable :: atomset1, atomset2, atomperm1
   integer(ik) :: num_records, max_trials, conv_freq
   integer(ik) :: i, j, natoms_local, base

   c_error_code = 0;  c_occ_records = 0;  c_natoms = 0

   ! Requested record count (must be at least 1)
   num_records = max(1_ik, c_n_records)

   ! Initialise all requested transform slots to identity so that, even on
   ! an early error return, every slot the caller allocated is well-defined.
   do i = 1, num_records
      call set_identity_transform(c_transform_list((i-1)*16+1:i*16))
   end do

   ! Set options
   align_flag       = c_align_flag
   remap_flag       = c_remap_flag
   heavy_flag       = c_heavy_flag
   mass_flag        = c_mass_flag
   mirror_flag      = c_mirror_flag
   label_flag       = c_label_flag
   bond_flag        = c_bond_flag
   bond_tol         = c_bond_tol
   random_flag      = c_random_flag
   print_stats      = c_print_stats
   print_assigntree = c_print_assigntree
   stochastic_flag  = .TRUE.
   adaptive_flag    = .TRUE.

   conv_freq   = c_conv_freq
   max_trials  = c_max_trials

   ! Build atom_t and bond_t arrays from flat C arrays
   call build_atoms(n_atoms1, c_atom_data1, c_coords1, atoms1)
   call build_atoms(n_atoms2, c_atom_data2, c_coords2, atoms2)
   call build_bonds(n_bonds1, c_bond_data1, bonds1)
   call build_bonds(n_bonds2, c_bond_data2, bonds2)

   if (heavy_flag) then
      ! Include only heavy atoms
      call include_heavy_atoms(atoms1, atomset1)
      call include_heavy_atoms(atoms2, atomset2)
   else
      ! Include all atoms
      call include_all_atoms(atoms1, atomset1)
      call include_all_atoms(atoms2, atomset2)
   end if

   ! Collect atom types in a partition
   call collect_atomtypes(atomset1, atomset2, atoms1, atoms2, atomtypes)

   ! Abort if molecules are not isomers
   if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
      c_error_code = 1
      return
   end if

   ! Reset bonds
   if (bond_flag) then
      call bonds_from_atoms(atoms1, bonds1)
      call bonds_from_atoms(atoms2, bonds2)
   end if

   ! Abort if either molecule has no bonds
   if (size(bonds1) < 1 .or. size(bonds2) < 1) then
      c_error_code = 3
      return
   end if

   ! Set adjacency lists
   call adjacency_from_bonds(atomset1, atoms1, bonds1, adjcs1)
   call adjacency_from_bonds(atomset2, atoms2, bonds2, adjcs2)

   ! Get user defined atom weights
   if (mass_flag) then
      weights1 = atomic_masses(atoms1%elnum)
      weights2 = atomic_masses(atoms2%elnum)
   else
      allocate (weights1(size(atoms1)), source=1.0_rk)
      allocate (weights2(size(atoms2)), source=1.0_rk)
   end if

   ! Get mol1 coordinates
   coords1 = get_coords(atoms1)

   ! Get mol2 coordinates
   if (mirror_flag) then
      coords2 = get_mirrored_coords(atoms2)
   else
      coords2 = get_coords(atoms2)
   end if

   if (align_flag) then

      center1 = get_centroid(atomset1, atoms1, weights1)
      center2 = get_centroid(atomset2, atoms2, weights2)
      call translate_coords(coords2, center1 - center2)

      coords1w = get_weighted_coords(atoms1, weights1, center1)
      coords2w = get_weighted_coords(atoms2, weights2, center2)

   else

      coords1w = get_weighted_coords(atoms1, weights1)
      coords2w = get_weighted_coords(atoms2, weights2)

   end if

   if (remap_flag) then

      if (align_flag) then

         call allocate_registry(registry, num_records)
         call optimize_atomperm_conformer(atomset1, atomset2, adjcs1, adjcs2, atomtypes, &
               coords1w, coords2w, conv_freq, max_trials, registry)

         if (print_stats) call print_records(registry)

         c_occ_records = registry%occ_records
         do i = 1, registry%occ_records
            atomperm1 = registry%records(i)%atomperm1
            natoms_local = size(atomperm1)
            rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
            coords2r = rotated_coords(coords2, rotquat, center1)
            rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))
            call build_homogeneous_transform(rotquat, center1, center2, &
                  c_transform_list((i-1)*16+1:i*16))

            c_rmsd_list(i) = rmsd
            base = (i-1) * natoms_local
            do j = 1, natoms_local
               c_atomperm_list(base+j) = atomperm1(j) - 1  ! 0-based for C
            end do
         end do
         c_natoms = natoms_local

      else

         call assign_atomperm_conformer(adjcs1, adjcs2, atomtypes, coords1w, coords2w, atomperm1)

         coords2r = coords2
         rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))

         c_occ_records = 1
         c_rmsd_list(1) = rmsd
         c_natoms = size(atomperm1)
         do j = 1, size(atomperm1)
            c_atomperm_list(j) = atomperm1(j) - 1  ! 0-based for C
         end do

      end if

   else

      ! Maintain input atom order
      call init_array(atomperm1, size(atoms1), identity)

      ! Abort if atoms do not match
      if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
         c_error_code = 2
         return
      end if

      ! Abort if bonds do not match
      if (adjacencydiff(atomset1, atomperm1, adjcs1, adjcs2) > 0) then
         c_error_code = 4
         return
      end if

      if (align_flag) then
         rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
         coords2r = rotated_coords(coords2, rotquat, center1)
         call build_homogeneous_transform(rotquat, center1, center2, c_transform_list(1:16))
      else
         coords2r = coords2
      end if

      rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))

      c_occ_records = 1
      c_rmsd_list(1) = rmsd
      c_natoms = size(atomperm1)
      do j = 1, size(atomperm1)
         c_atomperm_list(j) = atomperm1(j) - 1  ! 0-based for C
      end do

   end if

end subroutine conformsd_calculate

end module c_binding_conformsd