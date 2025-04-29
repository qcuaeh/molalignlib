subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_container), intent(in) :: mnatypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(bipartition_container) :: submnatypes

   submnatypes = mnatypes
   call assign_atoms_conf_rec(submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)

end subroutine

recursive subroutine assign_atoms_conf_rec( submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_container), intent(inout) :: submnatypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(metapartition_type) :: metatypes
   integer :: h, i

   call collect_mnatypes(mol1, submnatypes%partition1(), metatypes)
!   call metatypes%print_parts()
!   call submnatypes%print_parts()

   do i = 1, metatypes%num_parts
!      write (stderr, *) 'loop:', i
      h = random_element(metatypes%parts(i)%items)
      call solve_lap(submnatypes%parts(h), coords1, coords2, atomperm, dist)
      call split_crossmnatypes(h, atomperm, submnatypes)
      call compute_crossmnatypes(mol1, mol2, submnatypes)
      call assign_atoms_conf_rec(submnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   end do

end subroutine

subroutine assign_atoms_conf( mnatypes, mol1, mol2, coords1, coords2, atomperm, dist)
   type(mol_type), intent(in) :: mol1, mol2
   type(bipartition_container), intent(in) :: mnatypes
   real(rk), dimension(:,:), intent(in) :: coords1, coords2
   integer, dimension(:), intent(out) :: atomperm
   real(rk), intent(out) :: dist
   ! Local variables
   type(bipartition_container) :: submnatypes
   type(metapartition_type) :: metatypes
   integer :: h, i, j, k

   write (stderr, *) repeat('*', 80)
   call mnatypes%print_parts()

   call collect_mnatypes(mol1, mnatypes%partition1(), metatypes)
   call metatypes%print_parts()
   submnatypes = mnatypes

   do while (metatypes%num_parts > 0)
      write (stderr, *) repeat('+', 80)
      do i = 1, metatypes%num_parts
         h = random_element(metatypes%parts(i)%items)
         do j = 1, metatypes%num_parts
            do k = 1, metatypes%parts(j)%num_items
               if (metatypes%parts(j)%items(k) > h) then
                  metatypes%parts(j)%items(k) = metatypes%parts(j)%items(k) + submnatypes%parts(h)%num_items1 - 1
               end if
            end do
         end do
         call solve_lap(submnatypes%parts(h), coords1, coords2, atomperm, dist)
         call split_crossmnatypes(h, atomperm, submnatypes)
         write (stderr, *)
         write (stderr, *) repeat(str(h)//'   ', 8)
         call submnatypes%print_parts()
         call metatypes%print_parts()
      end do
      write (stderr, *) repeat('-', 80)
      call compute_crossmnatypes(mol1, mol2, submnatypes)
      call submnatypes%print_parts()
      call collect_mnatypes(mol1, submnatypes%partition1(), metatypes)
      call metatypes%print_parts()
   end do

   call assign_atoms(submnatypes, coords1, coords2, atomperm, dist)

end subroutine
