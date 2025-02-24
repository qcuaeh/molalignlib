module partitiondict
use parameters
use partition
implicit none
private
public operator (.in.)

real, parameter :: MAX_LOAD_FACTOR = 0.5

type, public :: partitiondict_type
   type(partition_type), allocatable :: partitions(:)
   logical, allocatable :: occupied(:)
   integer :: num_slots
   integer :: num_occupied
contains
   procedure :: get_index
   procedure :: new_index
   procedure :: initialize
   procedure :: reset
end type

interface operator (.in.)
   module procedure in_dict
end interface

contains

subroutine initialize( this, dict_size)
   class(partitiondict_type), intent(inout) :: this
   integer, intent(in) :: dict_size
   
   this%num_slots = int(dict_size / MAX_LOAD_FACTOR)
   
   allocate (this%occupied(this%num_slots))
   allocate (this%partitions(this%num_slots))
   
   this%occupied = .false.
   this%num_occupied = 0

end subroutine

subroutine reset(this)
   class(partitiondict_type), intent(inout) :: this

   this%occupied = .false.
   this%num_occupied = 0

end subroutine

function new_index( this, partition) result(index)
   class(partitiondict_type), intent(inout) :: this
   type(partition_type), intent(in) :: partition
   integer :: index

   index = modulo(hash(partition), this%num_slots) + 1

   do while (this%occupied(index))
      if (this%partitions(index) == partition) then
         error stop "Key already in hash table"
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   if (this%num_occupied >= this%num_slots*MAX_LOAD_FACTOR) then
      error stop "Hash table too full"
   end if

   this%partitions(index) = partition
   this%occupied(index) = .true.
   this%num_occupied = this%num_occupied + 1

end function

function get_index( this, partition) result(index)
   class(partitiondict_type), intent(in) :: this
   type(partition_type), intent(in) :: partition
   integer :: index

   index = modulo(hash(partition), this%num_slots) + 1

   do while (this%occupied(index))
      if (this%partitions(index) == partition) then
         return
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   error stop "Key not found"

end function

function in_dict( partition, dict)
   class(partitiondict_type), intent(in) :: dict
   type(partition_type), intent(in) :: partition
   logical :: in_dict
   integer :: index

   index = modulo(hash(partition), dict%num_slots) + 1

   do while (dict%occupied(index))
      if (dict%partitions(index) == partition) then
         in_dict = .true.
         return
      end if
      index = modulo(index, dict%num_slots) + 1
   end do

   in_dict = .false.

end function

integer function hash( partition)
   type(partition_type), intent(in) :: partition
   integer :: h, i, part_hash

   ! DJB2 partition hashing
   hash = 5381
   do h = 1, partition%num_parts
      ! DJB2 part hashing
      part_hash = 5381
      do i = 1, partition%parts(h)%part_size
         ! Update part hash
         part_hash = iand(part_hash * 33 + partition%parts(h)%items(i), 2**16 - 1)
      end do
      ! Update partition hash
      hash = iand(hash * 33 + part_hash, 2**16 - 1)
   end do

end function

end module
