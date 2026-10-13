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

! Global run options, set by the programs and C bindings
logical(LK) :: rebond_flag
logical(LK) :: usebondtype_flag
logical(LK) :: random_flag
logical(LK) :: mirror_flag
logical(LK) :: useatomtype_flag
logical(LK) :: heavy_flag
logical(LK) :: massweight_flag
logical(LK) :: align_flag
logical(LK) :: remap_flag
logical(LK) :: printstats_flag
logical(LK) :: printassigntree_flag

end module
