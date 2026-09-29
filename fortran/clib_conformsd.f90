module clib_conformsd
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
   use cbind_utils
   use flags
   use error_codes
   implicit none

contains

subroutine conformsd_calculate(                                              &
      c_n_atoms1,  c_atom_data1,  c_coords1,                                    &
      c_n_bonds1,  c_bond_data1,                                                  &
      c_n_atoms2,  c_atom_data2,  c_coords2,                                    &
      c_n_bonds2,  c_bond_data2,                                                  &
      c_align_flag, c_remap_flag, c_heavy_flag, c_mass_flag,                &
      c_mirror_flag, c_atomlabel_flag, c_bondtol_flag, c_bond_tol,          &
      c_bondtype_flag,                                                      &
      c_print_stats, c_print_assigntree, c_random_flag,                      &
      c_conv_freq, c_max_trials,                                           &
      c_n_records,                                                          &
      c_rmsd_list, c_mapping_list,                                              &
      c_transform_list, c_occ_records, c_error_code)                         &
      bind(C, name="conformsd_calculate")
! C-callable symmetry-corrected RMSD between conformers (see the conformsd
! program).
!
! Inputs, as flat C arrays:
!   c_atom_data1/2 : packed atom data, length c_n_atoms*2:
!                    [elnum0, label0, elnum1, label1, ...]
!                    label = 0 means unlabelled.
!   c_coords1/2    : XYZ coordinates, row-major (c_n_atoms x 3), length c_n_atoms*3
!   c_n_bonds1/2   : number of bonds (may be 0 when c_bondtol_flag is true, i.e. derive from geometry)
!   c_bond_data1/2 : flat bond array, length c_n_bonds*3, layout: [atom1, atom2, type, ...]
!                    (1-based atom indices as in the original file)
!   c_bond_tol     : bond detection tolerance (only used, and required, when
!                    c_bondtol_flag is true; no default)
!   c_bondtype_flag: use the bond types in the HNA refinement. Types
!                    are compared, never interpreted, so both molecules
!                    must use the same convention (read from the same
!                    file format with the same parser). Ignored
!                    when c_bondtol_flag is true (perceived bonds are untyped).
!
! Multiple ranked candidate solutions:
!   c_n_records requests up to that many ranked candidate solutions. Records
!   beyond the first are only ever produced when both c_align_flag and
!   c_remap_flag are true (the optimize_mapping_conformer search); in every
!   other case exactly one record is written regardless of c_n_records.
!   c_occ_records reports how many were actually written; only the first
!   c_occ_records entries of c_rmsd_list, c_mapping_list, and
!   c_transform_list are meaningful. All output arrays are flattened and
!   must be allocated by the caller with at least c_n_records elements per
!   record (n_padding for c_mapping_list, 16 for c_transform_list).
!
! Transforms:
!   Each c_transform_list record is a row-major 4x4 homogeneous matrix that
!   maps the input coordinates of molecule 2 onto molecule 1. It includes the
!   reflection when c_mirror_flag is true, also without alignment (then it is
!   just the reflection; otherwise the identity).
!
! Atom permutations and padding:
!   n_padding = max(c_n_atoms1, c_n_atoms2). The smaller molecule is padded with
!   dummy atoms appended after its real atoms, so each record is a true
!   permutation of 0..n_padding-1. Entry j (0-based) is the atom of molecule 2
!   that goes on line j of molecule 1. Values >= c_n_atoms2 denote padding
!   atoms of molecule 2; entries j >= c_n_atoms1 are padding lines of
!   molecule 1 and carry the extra atoms of molecule 2. The sizes can only
!   differ when c_heavy_flag is true.
!   Atoms excluded from the comparison (hydrogens with c_heavy_flag) are
!   paired afterwards: bonded to the image of their heavy neighbour first,
!   then by distance, then with padding atoms.
!
! c_error_code values: see the error_codes module (error_codes.f90) and its
! C mirror, error_codes.h. This function can return MOLALIGN_SUCCESS,
! MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER,
! MOLALIGN_ERROR_NOT_ISOMERS, MOLALIGN_ERROR_MISSING_BONDS,
! MOLALIGN_ERROR_ATOM_TYPE_MISMATCH and MOLALIGN_ERROR_BOND_MISMATCH (both
! only when c_remap_flag is false), and MOLALIGN_ERROR_NOT_CONFORMERS
! (only when c_remap_flag is true; passed through from the conformer
! assignment routines).

   ! Molecule 1
   integer(ik), intent(in), value :: c_n_atoms1
   integer(ik), dimension(c_n_atoms1*2), intent(in) :: c_atom_data1
   real(rk),    dimension(c_n_atoms1*3), intent(in) :: c_coords1
   integer(ik), intent(in), value :: c_n_bonds1
   integer(ik), dimension(c_n_bonds1*3), intent(in) :: c_bond_data1

   ! Molecule 2
   integer(ik), intent(in), value :: c_n_atoms2
   integer(ik), dimension(c_n_atoms2*2), intent(in) :: c_atom_data2
   real(rk),    dimension(c_n_atoms2*3), intent(in) :: c_coords2
   integer(ik), intent(in), value :: c_n_bonds2
   integer(ik), dimension(c_n_bonds2*3), intent(in) :: c_bond_data2

   ! Flags
   logical(lk), intent(in), value :: c_align_flag, c_remap_flag, c_heavy_flag, c_mass_flag
   logical(lk), intent(in), value :: c_mirror_flag, c_atomlabel_flag, c_bondtol_flag
   real(rk),    intent(in), value :: c_bond_tol
   logical(lk), intent(in), value :: c_bondtype_flag
   logical(lk), intent(in), value :: c_print_stats, c_print_assigntree, c_random_flag
   integer(ik), intent(in), value :: c_conv_freq, c_max_trials

   ! Requested number of ranked records
   integer(ik), intent(in), value :: c_n_records

   ! Outputs
   real(rk),    dimension(*), intent(out) :: c_rmsd_list
   integer(ik), dimension(*), intent(out) :: c_mapping_list
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
   real(rk), dimension(:),   allocatable :: weights1, weights2, unit_weights
   real(rk) :: transmat2(3,3)
   real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
   real(rk), dimension(:,:), allocatable :: full_coords1, full_coords2, full_coords2r
   integer(ik), dimension(:), allocatable :: atomset1, atomset2
   integer(ik), dimension(:), allocatable :: bondtypes
   integer(ik), dimension(:), allocatable :: mapping1, full_atomperm1
   integer(ik) :: n_records, n_padding
   integer(ik) :: i, j, base

   c_error_code = MOLALIGN_SUCCESS
   c_occ_records = 0

   ! Requested record count (must be at least 1)
   n_records = max(1_ik, c_n_records)

   ! Initialise all requested transform slots to identity so that, even on
   ! an early error return, every slot the caller allocated is well-defined.
   do i = 1, n_records
      call set_identity_transform(c_transform_list((i-1)*16+1:i*16))
   end do

   ! Set options
   align_flag       = c_align_flag
   remap_flag       = c_remap_flag
   heavy_flag       = c_heavy_flag
   mass_flag        = c_mass_flag
   mirror_flag      = c_mirror_flag
   atomlabel_flag   = c_atomlabel_flag
   bondtol_flag     = c_bondtol_flag
   bond_tol         = c_bond_tol
   bondtype_flag    = c_bondtype_flag .and. .not. c_bondtol_flag
   random_flag      = c_random_flag
   print_stats      = c_print_stats
   print_assigntree = c_print_assigntree

   ! Build atom_t and bond_t arrays from flat C arrays
   ! (returns MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER for out-of-range elnums)
   call build_atoms(c_n_atoms1, c_atom_data1, c_coords1, atoms1, c_error_code)
   if (c_error_code /= MOLALIGN_SUCCESS) return
   call build_atoms(c_n_atoms2, c_atom_data2, c_coords2, atoms2, c_error_code)
   if (c_error_code /= MOLALIGN_SUCCESS) return
   call build_bonds(c_n_bonds1, c_bond_data1, bonds1)
   call build_bonds(c_n_bonds2, c_bond_data2, bonds2)

   ! Pad the smaller molecule with dummy atoms appended at the end, so that
   ! both molecules have n_padding atoms and every atom permutation can be a
   ! bijection. Real atoms keep their indices; from here on padding atoms
   ! are recognised by index (> c_n_atoms1 or > c_n_atoms2).
   n_padding = max(c_n_atoms1, c_n_atoms2)
   call pad_atoms(atoms1, n_padding)
   call pad_atoms(atoms2, n_padding)

   ! Atom sets only ever contain real atoms
   if (heavy_flag) then
      ! Include only heavy atoms
      call include_heavy_atoms(atoms1(1:c_n_atoms1), atomset1)
      call include_heavy_atoms(atoms2(1:c_n_atoms2), atomset2)
   else
      ! Include all atoms
      call include_all_atoms(atoms1(1:c_n_atoms1), atomset1)
      call include_all_atoms(atoms2(1:c_n_atoms2), atomset2)
   end if

   ! From here on the comparison works on the included atoms only, numbered
   ! 1..size(atomset1) in molecule 1 and 1..size(atomset2) in molecule 2.
   ! mapping1 maps included atoms to included atoms in that compact
   ! numbering; complete_mapping turns it into full_atomperm1 over all
   ! (padded) atoms, which is what the caller receives.

   ! Collect atom types in a partition
   call collect_atomtypes(atoms1(atomset1), atoms2(atomset2), atomtypes)

   ! Abort if molecules are not isomers. Without heavy_flag the atom sets
   ! hold all real atoms, so different hydrogen counts are caught here.
   if (any(atomtypes%parts%n_items1 /= atomtypes%parts%n_items2)) then
      c_error_code = MOLALIGN_ERROR_NOT_ISOMERS
      return
   end if

   ! Reset bonds (padding atoms never get bonds)
   if (bondtol_flag) then
      call bonds_from_atoms(atoms1(1:c_n_atoms1), bonds1)
      call bonds_from_atoms(atoms2(1:c_n_atoms2), bonds2)
   end if

   ! Abort if either molecule has no bonds
   if (size(bonds1) < 1 .or. size(bonds2) < 1) then
      c_error_code = MOLALIGN_ERROR_MISSING_BONDS
      return
   end if

   ! Set adjacency lists of the included atoms (bonds to excluded atoms are
   ! dropped). With bondtype_flag, the edges used by the refinement carry
   ! the bond types, compacted jointly for both molecules.
   if (bondtype_flag) then
      bondtypes = distinct_bondtypes(bonds1, bonds2)
      call adjacency_from_bonds(atoms1(atomset1), extract_bonds(atomset1, n_padding, bonds1), adjcs1, bondtypes)
      call adjacency_from_bonds(atoms2(atomset2), extract_bonds(atomset2, n_padding, bonds2), adjcs2, bondtypes)
   else
      call adjacency_from_bonds(atoms1(atomset1), extract_bonds(atomset1, n_padding, bonds1), adjcs1)
      call adjacency_from_bonds(atoms2(atomset2), extract_bonds(atomset2, n_padding, bonds2), adjcs2)
   end if

   ! Weights of the included atoms, normalised to sum 1 over them so both
   ! molecules are scaled by the same factor
   if (mass_flag) then
      weights1 = atomic_masses(atoms1(atomset1)%elnum)
      weights2 = atomic_masses(atoms2(atomset2)%elnum)
   else
      allocate (weights1(size(atomset1)), source=1.0_rk)
      allocate (weights2(size(atomset2)), source=1.0_rk)
   end if
   weights1 = weights1/sum(weights1)
   weights2 = weights2/sum(weights2)

   ! Linear transformation applied to molecule 2
   if (mirror_flag) then
      transmat2 = MIRROR_MATRIX
   else
      transmat2 = IDENTITY_MATRIX
   end if

   ! Coordinates of all (padded) atoms, molecule 2 transformed. They are only
   ! used to align the whole molecule 2 and to pair the excluded atoms.
   allocate (unit_weights(n_padding), source=1.0_rk)
   full_coords1 = get_coords(atoms1, unit_weights, ORIGIN, IDENTITY_MATRIX)
   full_coords2 = get_coords(atoms2, unit_weights, ORIGIN, transmat2)

   ! Coordinates of the included atoms, in the compact numbering of mapping1
   coords1 = full_coords1(:, atomset1)
   coords2 = full_coords2(:, atomset2)

   ! Centers of the included atoms. They are taken from the coordinates, so
   ! center2 is in the transformed (e.g. mirrored) frame of molecule 2.
   if (align_flag) then
      center1 = get_centroid(coords1, weights1)
      center2 = get_centroid(coords2, weights2)
   else
      center1 = ORIGIN
      center2 = ORIGIN
   end if

   ! Weighted coordinates of the included atoms.
   ! coords2 is already transformed, hence IDENTITY_MATRIX for both.
   coords1w = get_coords(coords1, weights1, center1, IDENTITY_MATRIX)
   coords2w = get_coords(coords2, weights2, center2, IDENTITY_MATRIX)

   ! Move molecule 2 onto the center of molecule 1
   if (align_flag) then
      call translate_coords(full_coords2, center1 - center2)
      call translate_coords(coords2, center1 - center2)
   end if

   if (remap_flag) then

      if (align_flag) then

         call allocate_registry(registry, n_records)
         call optimize_mapping_conformer(adjcs1, adjcs2, atomtypes, &
               coords1w, coords2w, c_conv_freq, c_max_trials, registry, c_error_code)
         if (c_error_code /= MOLALIGN_SUCCESS) return

         if (print_stats) call print_records(registry)

         c_occ_records = registry%occ_records
         do i = 1, registry%occ_records
            mapping1 = registry%records(i)%mapping1
            rotquat = least_rotquat(mapping1, coords1w, coords2w)
            full_coords2r = rotated_coords(full_coords2, rotquat, center1)
            coords2r = rotated_coords(coords2, rotquat, center1)
            rmsd = sqrt(sqdistmean(mapping1, weights1, coords1, coords2r))
            call build_homogeneous_transform(rotquat, transmat2, center1, center2, &
                  c_transform_list((i-1)*16+1:i*16))

            ! Pair the excluded atoms in the aligned frame
            call complete_mapping(atomset1, atomset2, mapping1, atoms1, atoms2, &
                  c_n_atoms1, c_n_atoms2, bonds1, bonds2, full_coords1, full_coords2r, full_atomperm1)

            c_rmsd_list(i) = rmsd
            base = (i-1) * n_padding
            do j = 1, n_padding
               c_mapping_list(base+j) = full_atomperm1(j) - 1  ! 0-based for C
            end do
         end do

      else

         call assign_mapping_conformer(adjcs1, adjcs2, atomtypes, coords1w, coords2w, mapping1, c_error_code)
         if (c_error_code /= MOLALIGN_SUCCESS) return

         full_coords2r = full_coords2
         coords2r = coords2
         ! Unaligned: molecule 2 is only transformed (mirrored or not)
         call build_homogeneous_transform(IDENTITY_QUATERNION, transmat2, ORIGIN, ORIGIN, &
               c_transform_list(1:16))
         rmsd = sqrt(sqdistmean(mapping1, weights1, coords1, coords2r))

         call complete_mapping(atomset1, atomset2, mapping1, atoms1, atoms2, &
               c_n_atoms1, c_n_atoms2, bonds1, bonds2, full_coords1, full_coords2r, full_atomperm1)

         c_occ_records = 1
         c_rmsd_list(1) = rmsd
         do j = 1, n_padding
            c_mapping_list(j) = full_atomperm1(j) - 1  ! 0-based for C
         end do

      end if

   else

      ! Abort if atoms do not match: the same atoms must be included in
      ! both molecules, with the same types (the isomer check above
      ! guarantees that both atom sets have the same size)
      if (any(atomset1 /= atomset2)) then
         c_error_code = MOLALIGN_ERROR_ATOM_TYPE_MISMATCH
         return
      end if
      if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
         c_error_code = MOLALIGN_ERROR_ATOM_TYPE_MISMATCH
         return
      end if

      ! Maintain input atom order, both among the included atoms and in the
      ! full (padded) permutation. Identity is meaningful here: it is the
      ! mapping being tested.
      call init_array(mapping1, size(atomset1), identity)
      call init_array(full_atomperm1, size(atoms1), identity)

      ! Abort if bonds do not match
      if (adjacencydiff(mapping1, adjcs1, adjcs2) > 0) then
         c_error_code = MOLALIGN_ERROR_BOND_MISMATCH
         return
      end if

      if (align_flag) then
         rotquat = least_rotquat(mapping1, coords1w, coords2w)
         full_coords2r = rotated_coords(full_coords2, rotquat, center1)
         coords2r = rotated_coords(coords2, rotquat, center1)
         call build_homogeneous_transform(rotquat, transmat2, center1, center2, c_transform_list(1:16))
      else
         full_coords2r = full_coords2
         coords2r = coords2
         ! Unaligned: molecule 2 is only transformed (mirrored or not)
         call build_homogeneous_transform(IDENTITY_QUATERNION, transmat2, ORIGIN, ORIGIN, &
               c_transform_list(1:16))
      end if

      rmsd = sqrt(sqdistmean(mapping1, weights1, coords1, coords2r))

      c_occ_records = 1
      c_rmsd_list(1) = rmsd
      do j = 1, n_padding
         c_mapping_list(j) = full_atomperm1(j) - 1  ! 0-based for C
      end do

   end if

end subroutine conformsd_calculate

end module clib_conformsd
