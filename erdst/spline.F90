! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2019 Nobuyuki Matubayasi
! Copyright (C) 2010-2019 Shun Sakuraba
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

module spline
  implicit none
  real, allocatable :: coeff(:)
  integer :: order
contains

  subroutine spline_init(spline_order)
    integer, intent(in) :: spline_order
    integer :: i, k
    real :: factor
    order = spline_order
    allocate( coeff(0:order) )
    do i = 0, order
       factor = 1.0
       do k = 1, i ! pass thru when i == 0
          factor = factor * real(k)
       end do
       do k = 1, order - i ! pass thru when i == order
          factor = factor * real(k)
       end do
       factor = order / factor
       if(mod(i,2) == 1) factor = -factor
       coeff(i) = factor
    end do
  end subroutine spline_init

  ! FIXME: speed it up
  real function spline_value(rst)
    real, intent(in) :: rst
    integer :: i, k
    real :: f
    f = 0.0
    if((rst > 0.0) .and. (rst < order)) then
       k = int(rst)
       do i = 0, k
          f = f + coeff(i) * ((rst-i)**(order-1))
       end do
    endif
    spline_value = f
  end function spline_value

  ! never called in usual case
  subroutine spline_cleanup()
    deallocate(coeff)
  end subroutine spline_cleanup
end module spline
