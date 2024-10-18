module point_cloud_module
  implicit none
  private
  public :: find_points_in_radius, get_closest_points

  ! Constants
  real, parameter :: R = 1.0  ! Fixed radius
  real, parameter :: s = 0.1  ! Grid size (s << R)
  integer, parameter :: MAX_POINTS_PER_CELL = 100  ! Maximum points per cell

  ! Type to store grid cell data
  type grid_cell
    integer :: num_points
    integer :: point_indices(MAX_POINTS_PER_CELL)
  end type grid_cell

  ! Grid data
  type(grid_cell), allocatable :: grid(:,:,:)
  real :: grid_min(3), grid_max(3)
  integer :: grid_dims(3)

contains

  subroutine find_points_in_radius(cloud_points)
    real, intent(in) :: cloud_points(:,:)
    integer :: i, ix, iy, iz, N

    N = size(cloud_points, 2)

    ! Initialize grid
    call initialize_grid(cloud_points)

    ! Assign points to grid cells
    do i = 1, N
      ix = floor((cloud_points(1, i) - grid_min(1)) / s) + 1
      iy = floor((cloud_points(2, i) - grid_min(2)) / s) + 1
      iz = floor((cloud_points(3, i) - grid_min(3)) / s) + 1

      if (ix >= 1 .and. ix <= grid_dims(1) .and. &
          iy >= 1 .and. iy <= grid_dims(2) .and. &
          iz >= 1 .and. iz <= grid_dims(3)) then
        call add_point_to_cell(ix, iy, iz, i)
      end if
    end do
  end subroutine find_points_in_radius

  subroutine initialize_grid(cloud_points)
    real, intent(in) :: cloud_points(:,:)
    integer :: i

    ! Find grid boundaries
    grid_min = minval(cloud_points, dim=2)
    grid_max = maxval(cloud_points, dim=2)

    ! Add padding
    grid_min = grid_min - R
    grid_max = grid_max + R

    ! Calculate grid dimensions
    grid_dims = ceiling((grid_max - grid_min) / s)

    ! Allocate grid
    allocate(grid(grid_dims(1), grid_dims(2), grid_dims(3)))

    ! Initialize grid cells
    do i = 1, product(grid_dims)
      grid(i,1,1)%num_points = 0
    end do
  end subroutine initialize_grid

  subroutine add_point_to_cell(ix, iy, iz, point_index)
    integer, intent(in) :: ix, iy, iz, point_index
    integer :: num_points

    num_points = grid(ix, iy, iz)%num_points

    if (num_points < MAX_POINTS_PER_CELL) then
      num_points = num_points + 1
      grid(ix, iy, iz)%num_points = num_points
      grid(ix, iy, iz)%point_indices(num_points) = point_index
    else
      print *, "Warning: Maximum number of points per cell reached. Some points may be omitted."
    end if
  end subroutine add_point_to_cell

  subroutine get_closest_points(r, close_points)
    real, intent(in) :: r(3)
    integer, allocatable, intent(out) :: close_points(:)
    integer :: ix, iy, iz, num_close_points

    ! Find the closest grid cell
    ix = floor((r(1) - grid_min(1)) / s) + 1
    iy = floor((r(2) - grid_min(2)) / s) + 1
    iz = floor((r(3) - grid_min(3)) / s) + 1

    ! Bound check
    ix = max(1, min(ix, grid_dims(1)))
    iy = max(1, min(iy, grid_dims(2)))
    iz = max(1, min(iz, grid_dims(3)))

    ! Get the number of close points
    num_close_points = grid(ix, iy, iz)%num_points

    ! Allocate and fill the close_points array
    allocate(close_points(num_close_points))
    close_points = grid(ix, iy, iz)%point_indices(1:num_close_points)
  end subroutine get_closest_points

end module point_cloud_module
