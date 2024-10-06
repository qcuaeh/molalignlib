module partition

implicit none
private

type, public :: List
   integer, allocatable :: atomidx(:)
end type

type, public :: Partition
   integer, private :: step = 0
   integer, private :: currpartidx = 0
   integer, private :: curratomidx = 0
   integer, private :: currnatomidx = 0
   type(List), allocatable :: partidx(:)
contains
   procedure :: init_partition_list   ! uses a list of atom types to initilize
   procedure :: init_partition_id   ! initializes identity part for natoms
   generic   :: init_partition => init_partition_list, init_partition_id  ! creates a partition structure
   procedure :: get_atomidx   ! returns the current atomidx while looping
   procedure :: get_partidx   ! returns the current partidx while looping
   procedure :: get_step   ! returns the loop step counter (1, 2, ..., natomidx)
   procedure :: get_npartitions   ! returns the number of partition lists stored
   procedure :: get_natomidx_part   ! returns the number of indices in a partit
   procedure :: get_natomidx_total   ! returns the total number of indices
   generic :: get_natomidx => get_natomidx_part, get_natomidx_total  ! returns the number of indices in one/all partitions
   procedure :: loop  ! loops over all atomidx: do while (self%loop()) ... end do
   procedure :: looppart   ! loops over all atomidx within a specified part
   procedure :: print => print_partition   ! prints the contents of the partition
   procedure :: del_partition   ! deallocates partitions stored
end type

contains

subroutine init_partition_list(self, typelist)
   class(Partition), intent(inout) :: self
   integer, dimension(:), intent(in) :: typelist
   integer :: ntype, i
   integer, dimension(size(typelist)) :: typelengths
   integer, dimension(size(typelist), size(typelist)) :: atomidxmat

   self%step = 0
   self%currpartidx = 0
   self%curratomidx = 0
   self%currnatomidx = 0

   ntype = 0
   typelengths(:) = 0
   atomidxmat(:,:) = 0

! populate the temporary matrix atomidxmat([atomidx], [type])
   do i = 1, size(typelist)
      typelengths(typelist(i)) = typelengths(typelist(i)) + 1
      atomidxmat(typelengths(typelist(i)),typelist(i)) = i
      if (typelist(i) > ntype) then
         ntype = typelist(i)
      end if
   end do

   allocate(self%partidx(ntype))
! build partition structure
   do i = 1, ntype
      allocate(self%partidx(i)%atomidx(typelengths(i)))
      self%partidx(i)%atomidx = atomidxmat(:typelengths(i),i)
   end do

end subroutine init_partition_list

subroutine init_partition_id(self, natom)
   class(Partition), intent(inout) :: self
   integer, intent(in) :: natom
   integer :: n

   self%step = 0
   self%currpartidx = 0
   self%curratomidx = 0
   self%currnatomidx = 0

   allocate(self%partidx(1))   ! single partition (identity)
   allocate(self%partidx(1)%atomidx(natom))
   self%partidx(1)%atomidx = [(n, n = 1, natom)]
   
end subroutine init_partition_id

function get_atomidx(self) result(atomidx)
   class(Partition), intent(in) :: self
   integer :: atomidx

   if (self%currpartidx == 0) then
      atomidx = 0
   else
      atomidx = self%partidx(self%currpartidx)%atomidx(self%curratomidx)
   end if

end function get_atomidx

function get_partidx(self) result(partidx)
   class(Partition), intent(in) :: self
   integer :: partidx

   partidx = self%currpartidx

end function get_partidx

function get_step(self) result(stepord)
   class(Partition), intent(in) :: self
   integer :: stepord

   stepord = self%step

end function get_step

function get_npartitions(self) result(npart)
   class(Partition), intent(in) :: self
   integer :: npart

   npart = size(self%partidx)

end function get_npartitions

function get_natomidx_part(self, partidx) result(natomidx)
   class(Partition), intent(in) :: self
   integer, intent(in) :: partidx
   integer :: natomidx

   natomidx = 0
   if (partidx > self%get_npartitions()) then
      write (stderr, '(a)') "partidx exceedes number of partitions."
      return
   end if
   natomidx = size(self%partidx(partidx)%atomidx)

end function get_natomidx_part

function get_natomidx_total(self) result(natomidx)
   class(Partition), intent(in) :: self
   integer :: natomidx
   integer :: i

   natomidx = 0
   do i = 1, size(self%partidx)
      natomidx = natomidx + size(self%partidx(i)%atomidx)
   end do

end function get_natomidx_total

function loop (self, atomidx, partidx) result(continues)
   class(Partition), intent(inout) :: self
   integer, intent(out), optional :: atomidx, partidx   ! to optionally output current atomidx and partidx
   logical :: continues

   if (self%currpartidx == 0) then   ! starts loop from first block and atomidx
      self%step = 1
      self%currpartidx = 1
      self%curratomidx = 1
      self%currnatomidx = size(self%partidx(1)%atomidx)
      continues = .true.
   else   ! continues loop
      if (self%curratomidx /= self%currnatomidx) then   ! intermediate atomidx
         self%step = self%step + 1
         self%curratomidx = self%curratomidx + 1
         continues = .true.
      else   ! last position (atomidx) in a block
         if (self%currpartidx /= size(self%partidx)) then  ! intermediate block
            self%step = self%step + 1
            self%currpartidx = self%currpartidx + 1
            self%curratomidx = 1
            self%currnatomidx = size(self%partidx(self%currpartidx)%atomidx)
            continues = .true.
         else   ! ends loop (last atomidx in the last block)
            self%step = 0
            self%currpartidx = 0
            self%curratomidx = 0
            self%currnatomidx = 0
            continues = .false.
         end if
      end if
   end if

   if (present(atomidx)) then
      atomidx = self%curratomidx
   end if
   if (present(partidx)) then
      partidx = self%currpartidx
   end if

end function loop

function looppart(self, partidx, atomidx) result(continues)
   class(Partition), intent(inout) :: self
   integer, intent(in) :: partidx
   integer, intent(out), optional :: atomidx   ! to optionally output current atomidx
   logical :: continues

   if (self%curratomidx == 0) then   ! starts loop at first atomidx of specified partidx
      self%step = 1
      self%currpartidx = partidx
      self%curratomidx = 1
      self%currnatomidx = size(self%partidx(partidx)%atomidx)
      continues = .true.
   else   ! continues loop
      if (self%curratomidx /= self%currnatomidx) then   ! intermediate atomidx
         self%step = self%step + 1
         self%curratomidx = self%curratomidx + 1
         continues = .true.
      else   ! ends loop (last atomidx in the specified block)
         self%step = 0
         self%currpartidx = 0
         self%curratomidx = 0
         self%currnatomidx = 0
         continues = .false.
      end if
   end if

   if (present(atomidx)) then
      atomidx = self%curratomidx
   end if

end function looppart

subroutine print_partition(self)
   class(Partition), intent(inout) :: self
   integer :: npart, natomidx, i
   character(3) :: num
   
   npart = self%get_npartitions()
   do i = 1, npart
      natomidx = size(self%partidx(i)%atomidx)
      write (num, '(i0)') natomidx
      write (stderr, '(i0,a,'//num//'(i0,X))') i, ": ", self%partidx(i)%atomidx(:natomidx)
   end do

end subroutine print_partition

subroutine del_partition(self)
   class(Partition), intent(inout) :: self
   integer :: i

   do i = 1, size(self%partidx)
      deallocate(self%partidx(i)%atomidx)
   end do
   deallocate(self%partidx)

end subroutine del_partition

end module
