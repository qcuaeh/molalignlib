! MolAlignLib
! Copyright (C) 2025 José M. Vásquez

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module options
implicit none

! Flags
logical :: random_flag
logical :: mirror_flag
logical :: label_flag
logical :: heavy_flag
logical :: mass_flag
logical :: align_flag
logical :: remap_flag
logical :: stoch_flag
logical :: adaptive_flag
logical :: aligned_flag
logical :: stats_flag
logical :: tree_flag
logical :: bond_flag
logical :: atomorder_flag

! Bounds
integer :: max_trials
integer :: num_records
integer :: ato_thres
integer :: iso_thres
integer :: confo_thres

end module
