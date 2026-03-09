module c_binding_conformsd
   use parameters
   use str_utils
   use molecule
   use euclidean
   use chemdata
   use permutation
   use assorting
   use pruning_atoms
   use recording
   use assignment_conformer
   use alignment_conformer
   use c_binding_utils
   use options
   implicit none

contains

! ---------------------------------------------------------------------------
! C-callable wrapper for conformsd functionality.
!
! Caller supplies pre-read molecular data as flat C arrays:
!   c_atom_data1/2 : packed atom data, length n_atoms*2:
!                     [elnum0, label0, elnum1, label1, ...]
!                   label = 0 means unlabelled.
!   coords1/2     : XYZ coordinates, row-major (n_atoms x 3), length n_atoms*3
!   n_bonds1/2    : number of bonds (may be 0 when bond_flag=1, i.e. derive from geometry)
!   c_bond_data1/2    : flat bond array, length n_bonds*3, layout: [atom1, atom2, type, ...]
!                   (1-based atom indices as in the original file)
! ---------------------------------------------------------------------------
subroutine conformsd_calculate(                                              &
      n_atoms1,  c_atom_data1,  c_coords1,                                    &
      n_bonds1,  c_bond_data1,                                                  &
      n_atoms2,  c_atom_data2,  c_coords2,                                    &
      n_bonds2,  c_bond_data2,                                                  &
      c_align_flag, c_remap_flag, c_heavy_flag, c_mass_flag,                &
      c_mirror_flag, c_label_flag, c_bond_flag,                             &
      c_stats_flag, c_random_flag,                                          &
      c_confo_thres, c_max_trials,                                           &
      c_rmsd, c_natoms, c_atomperm, c_transform, c_error_code)              &
      bind(C, name="conformsd_calculate")

   ! --- molecule 1 ---
   integer(ik), intent(in), value :: n_atoms1
   integer(ik), dimension(n_atoms1*2), intent(in) :: c_atom_data1
   real(rk),    dimension(n_atoms1*3), intent(in) :: c_coords1
   integer(ik), intent(in), value :: n_bonds1
   integer(ik), dimension(n_bonds1*3), intent(in) :: c_bond_data1

   ! --- molecule 2 ---
   integer(ik), intent(in), value :: n_atoms2
   integer(ik), dimension(n_atoms2*2), intent(in) :: c_atom_data2
   real(rk),    dimension(n_atoms2*3), intent(in) :: c_coords2
   integer(ik), intent(in), value :: n_bonds2
   integer(ik), dimension(n_bonds2*3), intent(in) :: c_bond_data2

   ! --- flags ---
   logical(lk), intent(in), value :: c_align_flag, c_remap_flag, c_heavy_flag, c_mass_flag
   logical(lk), intent(in), value :: c_mirror_flag, c_label_flag, c_bond_flag
   logical(lk), intent(in), value :: c_stats_flag, c_random_flag
   integer(ik), intent(in), value :: c_confo_thres, c_max_trials

   ! --- outputs ---
   real(rk),               intent(out) :: c_rmsd
   integer(ik),            intent(out) :: c_natoms
   integer(ik), dimension(*), intent(out) :: c_atomperm
   real(rk),               intent(out) :: c_transform(16)
   integer(ik),            intent(out) :: c_error_code

   ! --- local ---
   type(atom_t), dimension(:), allocatable :: atoms1, atoms2
   type(bond_t), dimension(:), allocatable :: bonds1, bonds2
   type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
   type(partition_t) :: atomtypes
   type(registry_t)  :: registry
   real(rk) :: rmsd, center1(3), center2(3), rotquat(4)
   real(rk), dimension(:),   allocatable :: weights1, weights2
   real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
   integer(ik),  dimension(:),   allocatable :: atomset1, atomset2, atomperm1
   integer(ik) :: i

   c_error_code = 0;  c_rmsd = 0.0_rk;  c_natoms = 0
   call set_identity_transform(c_transform)

   ! --- set option flags ---
   align_flag  = c_align_flag;    remap_flag  = c_remap_flag
   heavy_flag  = c_heavy_flag;    mass_flag   = c_mass_flag
   mirror_flag = c_mirror_flag;   label_flag  = c_label_flag
   bond_flag   = c_bond_flag;     stats_flag  = c_stats_flag
   random_flag = c_random_flag
   stoch_flag   = .TRUE.;   adaptive_flag = .TRUE.
   print_tree_flag     = .FALSE.

   confo_thres = c_confo_thres
   max_trials  = c_max_trials
   num_records = 1

   ! --- build atom_t / bond_t arrays from flat C arrays ---
   call build_atoms(n_atoms1, c_atom_data1, c_coords1, atoms1)
   call build_atoms(n_atoms2, c_atom_data2, c_coords2, atoms2)
   call build_bonds(n_bonds1, c_bond_data1, bonds1)
   call build_bonds(n_bonds2, c_bond_data2, bonds2)

   ! --- atom sets ---
   if (heavy_flag) then
      call include_heavy_atoms(atoms1, atomset1)
      call include_heavy_atoms(atoms2, atomset2)
   else
      call include_all_atoms(atoms1, atomset1)
      call include_all_atoms(atoms2, atomset2)
   end if

   call collect_atomtypes(atomset1, atomset2, atoms1, atoms2, atomtypes)

   if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
      c_error_code = 1;  return
   end if

   ! --- connectivity ---
   if (bond_flag) then
      call bonds_from_atoms( atoms1, bonds1)
      call bonds_from_atoms( atoms2, bonds2)
   end if

   if (size(bonds1) < 1 .or. size(bonds2) < 1) then
      c_error_code = 2;  return
   end if

   call adjacency_from_bonds( atomset1, atoms1, bonds1, adjcs1)
   call adjacency_from_bonds( atomset2, atoms2, bonds2, adjcs2)

   ! --- weights ---
   if (mass_flag) then
      weights1 = atomic_masses(atoms1%elnum)
      weights2 = atomic_masses(atoms2%elnum)
   else
      allocate (weights1(size(atoms1)), source=1.0_rk)
      allocate (weights2(size(atoms2)), source=1.0_rk)
   end if

   ! --- coordinates (mirror if requested) ---
   coords1 = get_coords(atoms1)
   if (mirror_flag) then
      coords2 = get_mirrored_coords(atoms2)
   else
      coords2 = get_coords(atoms2)
   end if

   ! --- alignment / RMSD ---
   if (align_flag) then
      center1 = get_centroid(atomset1, atoms1, weights1)
      center2 = get_centroid(atomset2, atoms2, weights2)
      call translate_coords(coords2, center1 - center2)

      coords1w = get_weighted_coords(atoms1, weights1, center1)
      coords2w = get_weighted_coords(atoms2, weights2, center2)

      if (remap_flag) then
         call allocate_registry(registry, num_records)
         call optimize_atomperm_conformer(atomset1, atomset2, adjcs1, adjcs2, atomtypes, &
                                          coords1w, coords2w, registry)
         if (stats_flag) call print_records(registry)
         atomperm1 = registry%records(1)%atomperm1
      else
         if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
            c_error_code = 3;  return
         end if
         call init_identity_permutation(size(atoms1), atomperm1)
      end if

      rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
      coords2r = rotated_coords(coords2, rotquat, center1)
      rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))
      call build_homogeneous_transform(rotquat, center1, center2, c_transform)
   else
      coords1w = get_weighted_coords(atoms1, weights1)
      coords2w = get_weighted_coords(atoms2, weights2)

      if (remap_flag) then
         call assign_atomperm_conformer(adjcs1, adjcs2, atomtypes, coords1w, coords2w, atomperm1)
      else
         if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
            c_error_code = 3;  return
         end if
         call init_identity_permutation(size(atoms1), atomperm1)
      end if
      rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2))
   end if

   c_rmsd   = rmsd
   c_natoms = size(atomperm1)
   do i = 1, size(atomperm1)
      c_atomperm(i) = atomperm1(i) - 1  ! 0-based for C
   end do

end subroutine conformsd_calculate

end module c_binding_conformsd
