module metapartitioning
use molecule
use partition
use partitioning
use metapartition
use partitiondict
implicit none

contains

! Split MNA types
subroutine split_mnatypes(h, mnatypes)
   integer, intent(in) :: h
   type(partition_type), intent(inout) :: mnatypes
   ! Local variables
   integer :: k, i
   type(part_type), pointer :: newtype
   type(partition_type) :: subtypes

   call subtypes%initialize(mnatypes%num_items)

   do k = 1, h - 1
      call subtypes%add_part(mnatypes%parts(k))
   end do

   do i = 1, mnatypes%parts(h)%part_size
      newtype => subtypes%new_part(mnatypes%parts(h)%part_size)
      call newtype%add(mnatypes%parts(h)%items(i))
   end do

   do k = h + 1, mnatypes%num_parts
      call subtypes%add_part(mnatypes%parts(k))
   end do

   mnatypes = subtypes

end subroutine

! Split MNA types
subroutine collect_mnatypes(mol, mnatypes, metatypes)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: mnatypes
   type(metapartition_type), intent(out) :: metatypes
   ! Local variables
   integer :: h, k
   type(partition_type) :: subtypes
   type(partitiondict_type) :: partitiondict
   type(metapartition_type) :: mnametatypes
   type(metapartpointer_type), allocatable :: partlist(:)
   type(partition_type), allocatable :: subpartitionlist(:)
   logical, allocatable :: primality(:)

   call mnametatypes%initialize(mnatypes%num_parts)
   call partitiondict%initialize(mnatypes%num_parts)
   allocate (partlist(partitiondict%num_slots))
   allocate (subpartitionlist(mnatypes%num_parts))

   do h = 1, mnatypes%num_parts
      if (mnatypes%parts(h)%part_size > 1) then
         subtypes = mnatypes
         call split_mnatypes(h, subtypes)
         call compute_mnatypes(mol, subtypes)
         if (.not. (subtypes .in. partitiondict)) then
!            call subtypes%print_parts()
            partlist(partitiondict%new_index(subtypes))%ptr => &
               mnametatypes%new_part(mnatypes%num_parts)
            subpartitionlist(mnametatypes%num_parts) = subtypes
         end if
         call partlist(partitiondict%get_index(subtypes))%ptr%add(h)
      end if
   end do

   allocate (primality(mnametatypes%num_parts))
   primality(:) = .true.

   do h = 1, mnametatypes%num_parts
      if (primality(h)) then
         do k = h + 1, mnametatypes%num_parts
            if (primality(k)) then
               if (subpartitionlist(h) < subpartitionlist(k)) then
                  primality(k) = .false.
               else if (subpartitionlist(k) < subpartitionlist(h)) then
                  primality(h) = .false.
               end if
            end if
         end do
      end if
   end do

   call metatypes%initialize(count(primality))

   do h = 1, mnametatypes%num_parts
      if (primality(h)) then
         call metatypes%add_part(mnametatypes%parts(h))
      end if
   end do

end subroutine

end module
