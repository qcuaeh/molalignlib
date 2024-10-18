module dict_mod
   use, intrinsic :: iso_fortran_env, only: int32, int64
   implicit none
   private

   public :: dict, create_dict

   integer, parameter :: MAX_KEY_SIZE = 16
   integer, parameter :: HASH_CONSTANT = 1779033703

   type dict_item
      integer, allocatable :: key(:)
      integer :: value
      type(dict_item), pointer :: next => null()
   end type dict_item

   type, public :: dict
      type(dict_item), pointer, public :: buckets(:) => null()
      integer, public :: actual_capacity
      integer, public :: num_collision_levels
      integer, allocatable, public :: collisions(:)
      integer, allocatable, public :: occupations(:)
   contains
      procedure :: put => dict_put
      procedure :: get => dict_get
      procedure :: exists => dict_exists
      procedure :: reset => dict_reset
      final :: dict_finalize
   end type dict

contains

   function create_dict(requested_capacity) result(d)
      integer, intent(in) :: requested_capacity
      type(dict) :: d
      integer :: i, level_capacity

      d%actual_capacity = 1
      do while (d%actual_capacity < requested_capacity)
         d%actual_capacity = d%actual_capacity * 2
      end do

      d%num_collision_levels = 0
      level_capacity = d%actual_capacity / 2
      do while (level_capacity > 0)
         d%num_collision_levels = d%num_collision_levels + 1
         level_capacity = level_capacity / 2
      end do

      allocate(d%buckets(d%actual_capacity))
      do i = 1, d%actual_capacity
         nullify(d%buckets(i)%next)
      end do

      allocate(d%collisions(d%num_collision_levels))
      allocate(d%occupations(0:d%num_collision_levels))

      d%collisions = 0
      d%occupations = 0
   end function create_dict

   subroutine dict_put(this, key, value)
      class(dict), intent(inout) :: this
      integer, intent(in) :: key(:)
      integer, intent(in) :: value
      integer :: hash
      integer :: bucket_index, level
      type(dict_item), pointer :: new_item, current

      hash = compute_hash(key)
      bucket_index = modulo(hash, this%actual_capacity) + 1

      level = 0
      current => this%buckets(bucket_index)
      do while (associated(current%next))
         if (are_keys_equivalent(current%next%key, key)) then
            current%next%value = value  ! Update existing value
            return
         end if
         current => current%next
         level = level + 1
         if (level <= this%num_collision_levels) then
            this%collisions(level) = this%collisions(level) + 1
         end if
      end do

      ! Insert new item
      allocate(new_item)
      allocate(new_item%key(size(key)))
      new_item%key = key
      new_item%value = value
      new_item%next => current%next
      current%next => new_item

      this%occupations(min(level, this%num_collision_levels)) = &
         this%occupations(min(level, this%num_collision_levels)) + 1
   end subroutine dict_put

   function dict_get(this, key) result(value)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      integer :: value
      integer :: hash
      integer :: bucket_index
      type(dict_item), pointer :: current

      hash = compute_hash(key)
      bucket_index = modulo(hash, this%actual_capacity) + 1

      current => this%buckets(bucket_index)%next
      do while (associated(current))
         if (are_keys_equivalent(current%key, key)) then
            value = current%value
            return
         end if
         current => current%next
      end do

      error stop "Key not found"
   end function dict_get

   function dict_exists(this, key) result(exists)
      class(dict), intent(in) :: this
      integer, intent(in) :: key(:)
      logical :: exists
      integer :: hash
      integer :: bucket_index
      type(dict_item), pointer :: current

      hash = compute_hash(key)
      bucket_index = modulo(hash, this%actual_capacity) + 1

      current => this%buckets(bucket_index)%next
      do while (associated(current))
         if (are_keys_equivalent(current%key, key)) then
            exists = .true.
            return
         end if
         current => current%next
      end do

      exists = .false.
   end function dict_exists

   subroutine dict_reset(this)
      class(dict), intent(inout) :: this
      integer :: i
      type(dict_item), pointer :: current, next

      do i = 1, this%actual_capacity
         current => this%buckets(i)%next
         do while (associated(current))
            next => current%next
            deallocate(current%key)
            deallocate(current)
            current => next
         end do
         nullify(this%buckets(i)%next)
      end do

      this%collisions = 0
      this%occupations = 0
   end subroutine dict_reset

   subroutine dict_finalize(this)
      type(dict), intent(inout) :: this
      call this%reset()
      deallocate(this%buckets)
      deallocate(this%collisions)
      deallocate(this%occupations)
   end subroutine dict_finalize

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
   end function are_keys_equivalent

   function compute_hash(key) result(hash)
      integer, intent(in) :: key(:)
      integer :: hash
      integer :: i

      hash = 1
      do i = 1, size(key)
         hash = hash * (HASH_CONSTANT + 2 * key(i))
      end do
      hash = hash / 2
   end function compute_hash

end module dict_mod