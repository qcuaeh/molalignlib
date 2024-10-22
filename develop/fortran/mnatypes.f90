! Level up MNA types
subroutine levelup_mnatypes(atoms, types, subtypes)
   type(atom_type), dimension(:), intent(in) :: atoms
   type(partition_type), intent(in) :: types
   type(partition_type), intent(out) :: subtypes
   ! Local variables
   integer :: h, i, index_i
   type(hashtable_type) :: subtypedict
   type(pointertopart_type), allocatable :: subtypelist(:)
   integer, allocatable :: neighborhood(:)

   num_atoms = size(atoms0) + size(atoms1)
   max_part_size = max(types0%max_part_size, types1%max_part_size)

   call subtypes%init(num_atoms)
   call subtypedict%init(max_part_size)
   allocate (subtypelist(subtypedict%size))

   do h = 1, types%size
      do i = 1, types%parts(h)%size
         atom = types%parts(h)%items(i)
         neighborhood = types%index_part_map(adjlists(atom%index, atom%molnum)%indices, atom%molnum)
         if (.not. subtypedict%has_index(neighborhood)) then
            subtypelist(subtypedict%get_new_index(neighborhood))%ptr => subtypes%get_new_part()
         end if
         call subtypelist(subtypedict%get_index(neighborhood))%ptr%add(atom)
      end do
      call subtypedict%reset()
   end do

end subroutine
