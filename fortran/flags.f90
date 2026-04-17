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

module flags
use parameters
implicit none

! Flags
logical(lk) :: bond_flag
logical(lk) :: random_flag
logical(lk) :: mirror_flag
logical(lk) :: label_flag
logical(lk) :: heavy_flag
logical(lk) :: mass_flag
logical(lk) :: align_flag
logical(lk) :: remap_flag
logical(lk) :: stoch_flag
logical(lk) :: adaptive_flag
logical(lk) :: stats_flag
logical(lk) :: print_tree_flag

end module
