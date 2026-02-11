module conformsd_c_binding
use, intrinsic :: iso_c_binding
use parameters
use molecule
use euclidean
use utils
use chemistry
use permutation
use file_reading
use file_writing
use assorting
use pruning_atoms
use recording
use assignment_conformer
use alignment_conformer
use options
implicit none

contains

! C-callable wrapper for conformsd functionality
subroutine conformsd_calculate(c_file1, c_file2, &
                                c_align, c_remap, c_heavy, c_mass, &
                                c_mirror, c_label, c_bond, &
                                c_stochastic, c_exhaustive, c_stats, &
                                c_confo_thres, c_max_trials, c_num_records, &
                                c_rmsd, c_natoms, c_atomperm, c_error_code) &
                                bind(C, name="conformsd_calculate")
   ! Input parameters (strings as C char arrays)
   character(kind=c_char), dimension(*), intent(in) :: c_file1
   character(kind=c_char), dimension(*), intent(in) :: c_file2
   
   ! Boolean flags (C int: 0=false, 1=true)
   integer(c_int), intent(in), value :: c_align
   integer(c_int), intent(in), value :: c_remap
   integer(c_int), intent(in), value :: c_heavy
   integer(c_int), intent(in), value :: c_mass
   integer(c_int), intent(in), value :: c_mirror
   integer(c_int), intent(in), value :: c_label
   integer(c_int), intent(in), value :: c_bond
   integer(c_int), intent(in), value :: c_stochastic
   integer(c_int), intent(in), value :: c_exhaustive
   integer(c_int), intent(in), value :: c_stats
   
   ! Integer parameters
   integer(c_int), intent(in), value :: c_confo_thres
   integer(c_int), intent(in), value :: c_max_trials
   integer(c_int), intent(in), value :: c_num_records
   
   ! Output parameters
   real(c_double), intent(out) :: c_rmsd
   integer(c_int), intent(out) :: c_natoms
   integer(c_int), dimension(*), intent(out) :: c_atomperm
   integer(c_int), intent(out) :: c_error_code
   
   ! Local Fortran variables
   character(:), allocatable :: file1, file2
   character(:), allocatable :: title1, title2
   character(:), allocatable :: typein
   type(atom_t), dimension(:), allocatable :: atoms1, atoms2
   type(bond_t), dimension(:), allocatable :: bonds1, bonds2
   type(adjc_t), dimension(:), allocatable :: adjcs1, adjcs2
   type(partition_t) :: atomtypes
   type(registry_t) :: registry
   real(rk) :: rmsd
   real(rk) :: center1(3), center2(3), rotquat(4)
   real(rk), dimension(:), allocatable :: weights1, weights2
   real(rk), dimension(:,:), allocatable :: coords1, coords2, coords1w, coords2w, coords2r
   integer, dimension(:), allocatable :: atomset1, atomset2
   integer, dimension(:), allocatable :: atomperm1
   integer :: unitin
   integer :: i
   
   ! Initialize error code
   c_error_code = 0
   c_rmsd = 0.0_c_double
   c_natoms = 0
   
   ! Convert C strings to Fortran strings
   call c_to_f_string(c_file1, file1)
   call c_to_f_string(c_file2, file2)
   
   ! Set options from C parameters
   align_flag = (c_align /= 0)
   remap_flag = (c_remap /= 0)
   heavy_flag = (c_heavy /= 0)
   mass_flag = (c_mass /= 0)
   mirror_flag = (c_mirror /= 0)
   label_flag = (c_label /= 0)
   bond_flag = (c_bond /= 0)
   stats_flag = (c_stats /= 0)
   
   if (c_exhaustive /= 0) then
      stoch_flag = .FALSE.
      adaptive_flag = .FALSE.
   else if (c_stochastic /= 0) then
      stoch_flag = .TRUE.
      adaptive_flag = .FALSE.
   else
      stoch_flag = .TRUE.
      adaptive_flag = .TRUE.
   end if
   
   ! Set other default flags
   test_flag = .FALSE.
   aligned_flag = .FALSE.
   write_aligned = .FALSE.
   tree_flag = .FALSE.
   random_flag = .FALSE.
   mapping_flag = .FALSE.
   
   ! Set integer parameters
   confo_thres = c_confo_thres
   max_trials = c_max_trials
   num_records = c_num_records
   
   ! Read input files
   call open2read(file1, typein, unitin)
   call read_file(unitin, typein, title1, atoms1, bonds1)
   close(unitin)
   
   call open2read(file2, typein, unitin)
   call read_file(unitin, typein, title2, atoms2, bonds2)
   close(unitin)
   
   ! Select atoms (heavy or all)
   if (heavy_flag) then
      call include_heavy_atoms(atoms1, atomset1)
      call include_heavy_atoms(atoms2, atomset2)
   else
      call include_all_atoms(atoms1, atomset1)
      call include_all_atoms(atoms2, atomset2)
   end if
   
   ! Collect atom types
   call collect_atomtypes(atomset1, atomset2, atoms1, atoms2, atomtypes)
   
   ! Check if molecules are isomers
   if (any(atomtypes%parts%num_items1 /= atomtypes%parts%num_items2)) then
      c_error_code = 1  ! Not isomers
      return
   end if
   
   ! Set adjacency lists
   if (bond_flag) then
      call adjacency_from_atoms(atomset1, atoms1, adjcs1)
      call adjacency_from_atoms(atomset2, atoms2, adjcs2)
   else
      if (size(bonds1) < 1 .or. size(bonds2) < 1) then
         c_error_code = 2  ! Missing bonds
         return
      end if
      call adjacency_from_bonds(atomset1, bonds1, size(atoms1), adjcs1)
      call adjacency_from_bonds(atomset2, bonds2, size(atoms2), adjcs2)
   end if
   
   ! Set weights
   if (mass_flag) then
      weights1 = atomic_masses(atoms1%elnum)
      weights2 = atomic_masses(atoms2%elnum)
   else
      weights1 = uniform_weights(1._rk, size(atoms1))
      weights2 = uniform_weights(1._rk, size(atoms2))
   end if
   
   ! Get coordinates
   coords1 = get_coords(atoms1)
   if (mirror_flag) then
      coords2 = get_mirrored_coords(atoms2)
   else
      coords2 = get_coords(atoms2)
   end if
   
   ! Calculate RMSD based on alignment flag
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
         
         ! Print statistics if requested
         if (stats_flag) then
            call print_records(registry)
         end if
         
         ! Get first (best) result
         atomperm1 = registry%records(1)%atomperm1
         rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
         coords2r = rotated_coords(coords2, rotquat, center1)
         rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))
      else
         if (any(atomtypes%itemdir1 /= atomtypes%itemdir2)) then
            c_error_code = 3  ! Atom types do not match
            return
         end if
         
         allocate(atomperm1(size(coords1w, 2)))
         call init_identity_permutation(atomperm1)
         rotquat = least_rotquat(atomset1, atomperm1, coords1w, coords2w)
         coords2r = rotated_coords(coords2, rotquat, center1)
         rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2r))
      end if
   else
      coords1w = get_weighted_coords(atoms1, weights1)
      coords2w = get_weighted_coords(atoms2, weights2)
      
      if (remap_flag) then
         call assign_atomperm_conformer(adjcs1, adjcs2, atomtypes, coords1w, coords2w, atomperm1)
         rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2))
      else
         allocate(atomperm1(size(coords1w, 2)))
         call init_identity_permutation(atomperm1)
         rmsd = sqrt(sqdistmean(atomset1, atomperm1, weights1, coords1, coords2))
      end if
   end if
   
   ! Set output values
   c_rmsd = real(rmsd, c_double)
   c_natoms = size(atomperm1)
   
   ! Copy atom permutation to C array (C uses 0-based indexing, so subtract 1)
   do i = 1, size(atomperm1)
      c_atomperm(i) = atomperm1(i) - 1
   end do
   
end subroutine conformsd_calculate

! Helper to convert C string to Fortran allocatable string
subroutine c_to_f_string(c_string, f_string)
   character(kind=c_char), dimension(*), intent(in) :: c_string
   character(:), allocatable, intent(out) :: f_string
   integer :: i, str_len
   
   ! Find string length (up to null terminator)
   str_len = 0
   do i = 1, huge(i)
      if (c_string(i) == c_null_char) exit
      str_len = str_len + 1
   end do
   
   ! Allocate and copy
   allocate(character(str_len) :: f_string)
   do i = 1, str_len
      f_string(i:i) = c_string(i)
   end do
end subroutine c_to_f_string

end module conformsd_c_binding
