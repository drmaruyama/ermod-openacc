! -*- F90 -*-
! ERmod - Energy Representation Module
! Copyright (C) 2000- The ERmod authors
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
  use precision_kinds, only: wp
  implicit none
  real(wp), allocatable :: coeff(:)
  integer :: order
contains

  subroutine spline_init(spline_order)
    integer, intent(in) :: spline_order
    integer :: i, k
    real(wp) :: factor
    order = spline_order
    allocate( coeff(0:order) )
    do i = 0, order
       factor = 1.0_wp
       do k = 1, i ! pass thru when i == 0
          factor = factor * real(k, wp)
       end do
       do k = 1, order - i ! pass thru when i == order
          factor = factor * real(k, wp)
       end do
       factor = order / factor
       if (mod(i,2) == 1) factor = -factor
       coeff(i) = factor
    end do
  end subroutine spline_init

  ! FIXME: speed it up
  pure real(wp) function spline_value(rst)
    real(wp), intent(in) :: rst
    integer :: i, k
    real(wp) :: f
    f = 0.0_wp
    if ((rst > 0.0_wp) .and. (rst < order)) then
       k = int(rst)
       do i = 0, k
          f = f + coeff(i) * ((rst-i)**(order-1))
       end do
    endif
    spline_value = f
  end function spline_value

  ! Computes all "order" nonzero cardinal B-spline weights for a
  ! fractional offset u (0 <= u < 1) in a single O(order) pass, using
  ! the standard Cox-de Boor recursion (as in Essmann et al. 1995's
  ! SPME paper, and the "fill_bspline" routine used by most PME
  ! implementations). This replaces calling spline_value(u), spline_value(u+1),
  ! ..., spline_value(u+order-1) separately -- each of which does
  ! O(order) work on its own closed-form truncated-power-function sum,
  ! for O(order^2) total. values(spi) on return is exactly what
  ! spline_value(u + spi) computes, for spi = 0, ..., order-1
  ! (verified numerically against the closed-form version to within
  ! floating-point rounding, for orders 4, 5, 6 and 8).
  pure subroutine spline_values_all(u, values)
    real(wp), intent(in) :: u
    real(wp), intent(out) :: values(0:order-1)
    real(wp) :: data(0:order-1)
    real(wp) :: div
    integer :: k, l

    data(:) = 0.0_wp
    data(1) = u
    data(0) = 1.0_wp - u
    do k = 3, order
       div = 1.0_wp / real(k - 1, wp)
       data(k-1) = div * u * data(k-2)
       do l = 1, k - 2
          data(k-l-1) = div * ((u + real(l, wp)) * data(k-l-2) &
               + (real(k-l, wp) - u) * data(k-l-1))
       end do
       data(0) = div * (1.0_wp - u) * data(0)
    end do

    ! the recursion above fills data(:) "back to front" relative to
    ! the spi convention used by callers; un-reverse it here.
    do k = 0, order - 1
       values(k) = data(order - 1 - k)
    end do
  end subroutine spline_values_all

  ! never called in usual case
  subroutine spline_cleanup()
    deallocate(coeff)
  end subroutine spline_cleanup
end module spline
