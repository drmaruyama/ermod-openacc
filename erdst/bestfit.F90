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

!----------------
! best fit RMSD minimization module
! Based on reviews of:
!   Horn B. K. P. "Closed-form solution of 
!     absolute orientation using unit quaternions." 
!     Opt. Soc. Am. A, 4(4), 629 (1987).
! see also:
!   Theobald D. L. "Rapid calculation of RMSDs"
!     Acta Cryst., A61, 478 (2005).
!----------------

module quaternion
  implicit none
contains
  subroutine array_of_quaternion(q, a) 
    implicit none
    real, intent(in)  :: q(0:3)
    real, intent(out) :: a(3)
    a = q(1:3)
  end subroutine array_of_quaternion
  
  subroutine quaternion_of_array(a, q)
    implicit none
    real, intent(in)  :: a(3)
    real, intent(out) :: q(0:3)
    q(0) = 0.
    q(1:3) = a
  end subroutine quaternion_of_array

  subroutine prod(q1, q2, p)
    implicit none
    real, intent(in)  :: q1(0:3), q2(0:3)
    real, intent(out) :: p(0:3)
    real :: temp(3)
    p(0) = q1(0) * q2(0) - dot_product(q1(1:3), q2(1:3))
    call cross_product(q1(1:3), q2(1:3), temp)
    p(1:3) = temp + q1(0) * q2(1:3) + q2(0) * q1(1:3)

  contains
    subroutine cross_product(u, v, r)
      implicit none
      real, intent(in) :: u(3), v(3)
      real, intent(out) :: r(3)
      r = cshift(u, 1) * cshift(v, -1) - cshift(u, -1) * cshift(v, 1)
    end subroutine cross_product
  end subroutine prod
  
  subroutine conjugate(q, r)
    implicit none
    real, intent(in)  :: q(0:3)
    real, intent(out) :: r(0:3)
    r(0)   =  q(0)
    r(1:3) = -q(1:3)
  end subroutine conjugate

  subroutine rotate(q, r, res)
    implicit none
    real, intent(in)  :: q(1:3), r(0:3)
    real, intent(out) :: res(1:3)
    real :: t(0:3), cj(0:3)
    real :: q2(0:3), res2(0:3)
    q2(0) = 0.
    q2(1:3) = q(:)
    call conjugate(r, cj)
    call prod(r, q2,  t)
    call prod(t, cj, res2)
    res(:) = res2(1:3)
  end subroutine rotate

  subroutine rotate_inplace(n, coord, rot_quaternion)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: rot_quaternion(0:3)
    real, intent(inout) :: coord(3, n)
    integer :: i
    real :: temp(3)
    do i = 1, n
       call rotate(coord(:, i), rot_quaternion, temp)
       coord(:, i) = temp(:)
    end do
  end subroutine rotate_inplace
end module quaternion

module bestfit
  implicit none
contains  
  ! find a rotation that satisfies 
  subroutine find_rotation_quaternion(n, refPt, movedPt, masses, rotation)
    implicit none
    external dsyev
    integer :: n, info
    real, intent(in) :: refPt(3, n), movedPt(3, n), masses(n)
    real, intent(out) :: rotation(0:3)
    real :: inner_prod(3, 3)
    real :: matmax(4, 4)
    real :: eigenvalue(4)
    integer, parameter :: lwork = 256
    real :: work(lwork)
    integer :: i, j, k

    do i = 1, 3
       do j = 1, 3
          inner_prod(j, i) = sum(refPt(i, :) * movedPt(j, :) * masses(:))
       end do
    end do
    
    !-- Calculate a linear operator that describes sums of products.
    !-- Regrettably, there is no succinct expression to write; we just write it "as-is" in page 635.
    !-- We need upper half of the matrix only, therefore omitting lower.
    !-- (Be aware that LAPACK's matrix is
    !-- a(1, 1) a(1, 2) a(1, 3) ...
    !-- a(2, 1) a(2, 2) a(2, 3) ...) 
    matmax(1, 1) = + inner_prod(1, 1) + inner_prod(2, 2) + inner_prod(3, 3)
    matmax(2, 2) = + inner_prod(1, 1) - inner_prod(2, 2) - inner_prod(3, 3)
    matmax(3, 3) = - inner_prod(1, 1) + inner_prod(2, 2) - inner_prod(3, 3)
    matmax(4, 4) = - inner_prod(1, 1) - inner_prod(2, 2) + inner_prod(3, 3)

    matmax(1, 2) = + inner_prod(2, 3) - inner_prod(3, 2)
    matmax(1, 3) = + inner_prod(3, 1) - inner_prod(1, 3)
    matmax(1, 4) = + inner_prod(1, 2) - inner_prod(2, 1)

    matmax(2, 3) = + inner_prod(1, 2) + inner_prod(2, 1)
    matmax(3, 4) = + inner_prod(2, 3) + inner_prod(3, 2)
    matmax(2, 4) = + inner_prod(3, 1) + inner_prod(1, 3)

    !-- solve eigenvalue and eigenvector
    call dsyev('V', 'U', 4, matmax, 4, eigenvalue, work, lwork, info)

    if (info /= 0) then
       print *, "error on LAPACK/DSYEV"
       print *, "info =", info
       stop
    end if

    ! value with maximum eigenvalue is the required vector
    rotation(0:3) = matmax(1:4, 4)
  end subroutine find_rotation_quaternion

  subroutine center_of_mass(n, points, masses, center)
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: points(3, n), masses(n)
    real, intent(out) :: center(3)
    real :: sumOfMasses
    integer :: i
    
    sumOfMasses = sum(masses)
    if(sumOfMasses == 0) call halt_with_error('bst_zrw')
    do i = 1, 3
       center(i) = dot_product(points(i, :), masses(:)) / sumOfMasses
    end do
  end subroutine center_of_mass
  
  subroutine com_shift(n, points, masses, center)
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: points(3, n)
    real, intent(in) :: masses(n)
    real, intent(out) :: center(3)
    real :: com(3)
    integer :: i

    call center_of_mass(n, points, masses, com)

    do i = 1, n
       points(:, i) = points(:, i) - com(:)
    end do

    center(:) = com(:)
  end subroutine com_shift

  subroutine com_unshift(n, points, masses, center)
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: points(3, n)
    real, intent(in) :: masses(n)
    real, intent(in) :: center(3)
    integer :: i

    do i = 1, n
       points(:, i) = points(:, i) + center(:)
    end do
  end subroutine com_unshift

  ! standard fit-to-structure procedure
  subroutine fit(n, refcoord, coord, masses, outcoord)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: refcoord(3, n), coord(3, n), masses(n)
    real, intent(out) :: outcoord(3, n)
    call fit_a_rotate_b(n, refcoord, coord, masses, n, coord, outcoord)
  end subroutine fit

  ! fit a to refa, and get corresponding position for b
  subroutine fit_a_rotate_b(na, refa, a, massa, nb, b, bout)
    use quaternion, only: rotate
    implicit none
    integer, intent(in) :: na, nb
    real, intent(in) :: refa(3, na), a(3, na), massa(na)
    real, intent(in) :: b(3, nb)
    real, intent(out) :: bout(3, nb)
    
    real :: workref(3, na), work(3, na)
    real :: com_refa(3), com_a(3)
    real :: bcrd(3), bcrdr(3)
    real :: rotation(0:3)
    integer :: i
    
    ! copy coordinate & align com
    workref(:, :) = refa(:, :)
    call com_shift(na, workref, massa, com_refa)
    work(:, :) = a(:, :)
    call com_shift(na, work, massa, com_a)
    
    call find_rotation_quaternion(na, workref, work, massa, rotation)
    
    do i = 1, nb
       bcrd(:) = b(:, i) - com_a(:)
       call rotate(bcrd, rotation, bcrdr)
       bout(:, i) = bcrdr(:) + com_refa(:)
    end do
  end subroutine fit_a_rotate_b

  ! RMSD calculation with best-fit
  real function rmsd_bestfit(n, refcoord, coord, masses)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: refcoord(3, n), coord(3, n), masses(n)
    real :: fitted_coord(3, n)
    call fit(n, refcoord, coord, masses, fitted_coord)
    rmsd_bestfit = rmsd_nofit(n, refcoord, fitted_coord, masses)
  end function rmsd_bestfit

  ! RMSD calculation without any fitting procedure
  real function rmsd_nofit(n, crdA, crdB, masses)
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: crdA(3, n), crdB(3, n), masses(n)
    integer :: i
    real :: sumOfMasses, disp, dx(3)
    sumOfMasses = sum(masses)
    if(sumOfMasses == 0) call halt_with_error('bst_zrw')
    disp = 0.0
    do i = 1, n
       dx(1:3) = crdB(1:3, i) - crdA(1:3, i)
       disp = disp + masses(i) * sum( dx(1:3) ** 2 )
    enddo
    rmsd_nofit = sqrt( disp / sumOfMasses )
  end function rmsd_nofit
end module bestfit
