module hashtable
use stdio
use bounds

implicit none

integer, parameter :: HASH_CONSTANT = 1779033703

type, public :: hashtable_type
   integer :: size
   integer :: num_collision_levels
   integer, allocatable :: keys(:,:)
   integer, allocatable :: key_lengths(:)
   integer, allocatable :: collisions(:)
   integer, allocatable :: occupations(:)
   logical, allocatable :: occupied(:)
contains
   procedure :: init => hashtable_init
   procedure :: reset => hashtable_reset
   procedure :: has_index => hashtable_has_index
   procedure :: get_index => hashtable_get_index
   procedure :: get_new_index => hashtable_get_new_index
   final :: hashtable_finalize
end type

contains

subroutine hashtable_init(self, minimum_size)
   class(hashtable_type), intent(inout) :: self
   integer, intent(in) :: minimum_size
   integer :: i, collision_level_size

   collision_level_size = 1
   do while (collision_level_size < minimum_size)
      collision_level_size = collision_level_size * 2
   end do
   self%size = collision_level_size * 2
   self%num_collision_levels = 0
   collision_level_size = collision_level_size / 2
   do while (collision_level_size > 0)
      self%num_collision_levels = self%num_collision_levels + 1
      collision_level_size = collision_level_size / 2
   end do

   allocate(self%occupied(self%size))
   allocate(self%key_lengths(self%size))
   allocate(self%keys(maxcoord, self%size))
   allocate(self%collisions(self%num_collision_levels))
   allocate(self%occupations(0:self%num_collision_levels))

   self%occupied = .false.
   self%collisions = 0
   self%occupations = 0

end subroutine

function hashtable_get_new_index(self, key) result(index)
   class(hashtable_type), intent(inout) :: self
   integer, intent(in) :: key(:)
   integer :: hash, index
   integer :: collision_level, collision_level_index, collision_level_offset, collision_level_size

   hash = compute_hash(key)
   collision_level_index = modulo(hash, self%size/2)

   collision_level = 0
   collision_level_offset = 0
   collision_level_size = self%size / 2

   do
      index = collision_level_index + collision_level_offset + 1
      if (self%occupied(index)) then
         if (are_keys_equivalent(self%keys(:self%key_lengths(index), index), key)) then
            return
         else
            collision_level = collision_level + 1
            if (collision_level > self%num_collision_levels) then
               error stop "Dictionary is full"
            end if
            self%collisions(collision_level) = self%collisions(collision_level) + 1
            collision_level_offset = collision_level_offset + collision_level_size
            collision_level_size = collision_level_size / 2
            collision_level_index = collision_level_index / 2
         end if
      else
         self%keys(:size(key), index) = key
         self%key_lengths(index) = size(key)
         self%occupied(index) = .true.
         self%occupations(collision_level) = self%occupations(collision_level) + 1
         return
      end if
   end do

end function

function hashtable_get_index(self, key) result(index)
   class(hashtable_type), intent(in) :: self
   integer, intent(in) :: key(:)
   integer :: hash, index
   integer :: collision_level, collision_level_index, collision_level_offset, collision_level_size

   hash = compute_hash(key)
   collision_level_index = modulo(hash, self%size/2)

   collision_level = 0
   collision_level_offset = 0
   collision_level_size = self%size / 2

   do
      index = collision_level_index + collision_level_offset + 1
      if (self%occupied(index)) then
         if (are_keys_equivalent(self%keys(:self%key_lengths(index), index), key)) then
            return
         else
            collision_level = collision_level + 1
            if (collision_level > self%num_collision_levels) then
               error stop "Key not found"
            end if
            collision_level_offset = collision_level_offset + collision_level_size
            collision_level_size = collision_level_size / 2
            collision_level_index = collision_level_index / 2
         end if
      else
         error stop "Key not found"
      end if
   end do

end function

function hashtable_has_index(self, key) result(has_index)
   class(hashtable_type), intent(in) :: self
   integer, intent(in) :: key(:)
   logical :: has_index
   integer :: hash, index
   integer :: collision_level_index, collision_level, collision_level_offset, collision_level_size

   hash = compute_hash(key)
   collision_level_index = modulo(hash, self%size/2)

   collision_level = 0
   collision_level_offset = 0
   collision_level_size = self%size / 2

   do
      index = collision_level_index + collision_level_offset + 1
      if (self%occupied(index)) then
         if (are_keys_equivalent(self%keys(:self%key_lengths(index), index), key)) then
            has_index = .true.
            return
         else
            collision_level = collision_level + 1
            if (collision_level > self%num_collision_levels) then
               has_index = .false.
               return
            end if
            collision_level_offset = collision_level_offset + collision_level_size
            collision_level_size = collision_level_size / 2
            collision_level_index = collision_level_index / 2
         end if
      else
         has_index = .false.
         return
      end if
   end do

end function

subroutine hashtable_reset(self)
   class(hashtable_type), intent(inout) :: self

   self%occupied = .false.
   self%collisions = 0
   self%occupations = 0

end subroutine

subroutine hashtable_finalize(self)
   type(hashtable_type), intent(inout) :: self

   if (allocated(self%keys)) deallocate(self%keys)
   if (allocated(self%key_lengths)) deallocate(self%key_lengths)
   if (allocated(self%occupied)) deallocate(self%occupied)
   if (allocated(self%collisions)) deallocate(self%collisions)
   if (allocated(self%occupations)) deallocate(self%occupations)

end subroutine

function are_keys_equivalent(key1, key2) result(equivalent)
   integer, intent(in) :: key1(:), key2(:)
   logical :: equivalent
   integer :: i

   if (size(key1) /= size(key2)) then
      equivalent = .false.
      return
   end if

   equivalent = .true.
   do i = 1, size(key1)
      if (count(key1 == key1(i)) /= count(key2 == key1(i))) then
         equivalent = .false.
         return
      end if
   end do

end function

function compute_hash(key) result(hash)
   integer, intent(in) :: key(:)
   integer :: hash
   integer :: i

   hash = 1
   do i = 1, size(key)
      hash = hash * (HASH_CONSTANT + 2 * key(i))
   end do
   hash = hash / 2

end function

end module
