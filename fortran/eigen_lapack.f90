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

module eigen
use parameters

implicit none

external ssyev
external dsyev

private
public leasteigval
public leasteigvec

interface leasteigval
   module procedure leasteigval_sp
   module procedure leasteigval_dp
end interface

interface leasteigvec
   module procedure leasteigvec_sp
   module procedure leasteigvec_dp
end interface

contains

function leasteigval_sp(a) result(leasteigval)
   real(real32), intent(inout) :: a(4, 4)
   integer(ik) :: info
   real(real32) :: leasteigval
   real(real32) :: w(4), work(25)
   call ssyev('N', 'U', 4, a, 4, w, work, 25, info)
   leasteigval = w(1)
end function

function leasteigval_dp(a) result(leasteigval)
   real(real64), intent(inout) :: a(4, 4)
   integer(ik) :: info
   real(real64) :: leasteigval
   real(real64) :: w(4), work(25)
   call dsyev('N', 'U', 4, a, 4, w, work, 25, info)
   leasteigval = w(1)
end function

function leasteigvec_sp(a) result(leasteigvec)
   real(real32), intent(inout) :: a(4, 4)
   integer(ik) :: info
   real(real32) :: w(4), work(25)
   real(real32) :: leasteigvec(4)
   call ssyev('V', 'U', 4, a, 4, w, work, 25, info)
   leasteigvec = a(:, 1)
end function

function leasteigvec_dp(a) result(leasteigvec)
   real(real64), intent(inout) :: a(4, 4)
   integer(ik) :: info
   real(real64) :: leasteigvec(4)
   real(real64) :: w(4), work(25)
   call dsyev('V', 'U', 4, a, 4, w, work, 25, info)
   leasteigvec = a(:, 1)
end function 

end module
