module tupledict
use parameters
implicit none
private

public tupledict_type
public operator (.in.)

real, parameter :: MAX_LOAD_FACTOR = 0.5

abstract interface
   integer function hash_proc(tuple)
      integer, intent(in) :: tuple(:)
   end function
end interface

abstract interface
   logical function equality_proc(tuple1, tuple2)
      integer, intent(in) :: tuple1(:)
      integer, intent(in) :: tuple2(:)
   end function
end interface

type, public :: tuple_type
   integer, allocatable :: items(:)
end type

type, public :: tupledict_type
   type(tuple_type), allocatable :: tuples(:)
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

procedure(hash_proc), pointer :: hash
procedure(equality_proc), pointer :: equality

contains

subroutine initialize( this, dict_size, dict_mode)
   class(tupledict_type), intent(inout) :: this
   integer, intent(in) :: dict_size
   character(*), intent(in) :: dict_mode
   
   this%num_slots = int(dict_size / MAX_LOAD_FACTOR)
   
   allocate (this%occupied(this%num_slots))
   allocate (this%tuples(this%num_slots))
   
   this%occupied = .false.
   this%num_occupied = 0

   if (dict_mode == 'ordered') then
      hash => ordered_hash
      equality => ordered_equality
   else if (dict_mode == 'unordered') then
      hash => unordered_hash
      equality => unordered_equality
   else
      error stop 'Unknown dict mode'
   end if

end subroutine

subroutine reset(this)
   class(tupledict_type), intent(inout) :: this

   this%occupied = .false.
   this%num_occupied = 0

end subroutine

function new_index( this, tuple) result(index)
   class(tupledict_type), intent(inout) :: this
   integer, intent(in) :: tuple(:)
   integer :: index

   index = modulo(hash(tuple), this%num_slots) + 1

   do while (this%occupied(index))
      if (equality(this%tuples(index)%items, tuple)) then
         error stop "Key already in hash table"
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   if (this%num_occupied >= this%num_slots) then
      error stop "Hash table is full"
   end if

   this%tuples(index)%items = tuple
   this%occupied(index) = .true.
   this%num_occupied = this%num_occupied + 1

end function

function get_index( this, tuple) result(index)
   class(tupledict_type), intent(in) :: this
   integer, intent(in) :: tuple(:)
   integer :: index

   index = modulo(hash(tuple), this%num_slots) + 1

   do while (this%occupied(index))
      if (equality(this%tuples(index)%items, tuple)) then
         return
      end if
      index = modulo(index, this%num_slots) + 1
   end do

   error stop "Key not found"

end function

function in_dict( tuple, dict)
   integer, intent(in) :: tuple(:)
   class(tupledict_type), intent(in) :: dict
   logical :: in_dict
   integer :: index

   index = modulo(hash(tuple), dict%num_slots) + 1

   do while (dict%occupied(index))
      if (equality(dict%tuples(index)%items, tuple)) then
         in_dict = .true.
         return
      end if
      index = modulo(index, dict%num_slots) + 1
   end do

   in_dict = .false.

end function

integer function ordered_hash( tuple) result(hash)
! DJB2 hashing
   integer, intent(in) :: tuple(:)
   integer, parameter :: BIT_MASK = 2**16 - 1
   integer, parameter :: DJB2_INIT_HASH = 5381
   integer, parameter :: DJB2_MULTIPLIER = 33
   integer :: i

   hash = DJB2_INIT_HASH
   do i = 1, size(tuple)
      hash = iand(hash * DJB2_MULTIPLIER + tuple(i), BIT_MASK)
   end do

end function

function ordered_equality( tuple1, tuple2) result(equality)
   integer, intent(in) :: tuple1(:), tuple2(:)
   logical :: equality
   
   ! Check if sizes are equal
   if (size(tuple1) /= size(tuple2)) then
      equality = .false.
      return
   end if

   ! Check if tuples are equal
   if (all(tuple1 == tuple2)) then
      equality = .true.
   else
      equality = .false.
   end if

end function

integer function unordered_hash( tuple) result(hash)
   integer, intent(in) :: tuple(:)
   integer, parameter :: BIT_MASK = 2**16 - 1
   integer, parameter :: HASH_CONSTANT = 5381
   integer :: i

   hash = 1
   do i = 1, size(tuple)
      hash = iand(hash * (HASH_CONSTANT + 2 * tuple(i)), BIT_MASK)
   end do
   hash = hash / 2

end function

function unordered_equality( tuple1, tuple2) result(equality)
   integer, intent(in) :: tuple1(:), tuple2(:)
   logical :: equality
   integer :: i, j, matches
   
   ! Check if sizes are equal
   if (size(tuple1) /= size(tuple2)) then
      equality = .false.
      return
   end if

   ! Quick check for exact match first
   if (all(tuple1 == tuple2)) then
       equality = .true.
       return
   end if
   
   ! Early exit: Check if any value appears more times in one array
   ! than it does in the other
   do i = 1, size(tuple1)
       matches = 0
       do j = 1, size(tuple1)
           if (tuple1(i) == tuple2(j)) matches = matches + 1
           if (tuple1(i) == tuple1(j)) matches = matches - 1
       end do
       if (matches /= 0) then
           equality = .false.
           return
       end if
   end do

   equality = .true.

end function

end module
