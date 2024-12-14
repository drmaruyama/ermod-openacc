! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2019 Nobuyuki Matubayasi
! Copyright (C) 2010-2019 Shun Sakuraba
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! As a special exception, you may use this file as part of a free software
! without restriction.  Specifically, if other files instantiate
! templates or use macros or inline functions from this file, or you compile
! this file and link it with other files to produce an executable, this
! file does not by itself cause the resulting executable to be covered by
! the GNU General Public License.  
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!

module utility
  implicit none

contains

#ifdef HAVE_TRANSFER
  integer(8) function hash(arr, size) 
    implicit none
    integer, intent(in) :: size
    real, intent(in) :: arr(size)
    integer(8) :: ret
    integer :: i
    ret = 0

    do i = 1, size
       ret = ishftc(ret, 7)
       ret = ieor(ret, transfer(arr(i), ret))
    end do

    hash = ret
  end function hash

#else
! (== not HAVE_TRANSFER)
  integer(8) function hash(arr, size) 
    implicit none
    integer, intent(in) :: size
    real, intent(in) :: arr(size)
    integer(8) :: ret
    external hash_double, hash_float
    select case(kind(arr))
    case(4)
       call hash_float(arr, size, ret)
    case(8)
       call hash_double(arr, size, ret)
    case default
       stop "Error: hash(): unknown real type"
    end select
    hash = ret
  end function hash
#endif

  ! The following function is a snippet from Fortran wiki and in public domain.
  ! 
  ! This is a simple function to search for an available unit.
  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
  ! The UNIT value is returned by the function, and also by the optional
  ! argument. This allows the function to be used directly in an OPEN
  ! statement, and optionally save the result in a local variable.
  ! If no units are available, -1 is returned.
  integer function newunit(unit)
    implicit none
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=100, LUN_MAX=110
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

  ! convert cell-length & (alpha, beta, gamma) to cell vectors
  subroutine angles_to_cell_vector(cell_len, angles, out_cell_vectors)
    use engmain, only: PI
    implicit none
    real, intent(in) :: cell_len(3)
    real, intent(in) :: angles(3)
    real, intent(out) :: out_cell_vectors(3, 3)
    real :: alpha, beta, gamma, x, y, u, v, w

    alpha = angles(1) * PI / 180.0 ! for b-c axes
    beta  = angles(2) * PI / 180.0 ! for a-c axes
    gamma = angles(3) * PI / 180.0 ! for a-b axes

    ! ~a = (1, 0, 0)
    ! ~b = (x, y, 0)
    ! ~c = (u, v, w)
    ! ~a.~b = x = cos gamma
    ! |~a*~b| = y = sin gamma
    ! ~a.~c = u = cos beta
    ! ~b.~c = xu + yv = cos alpha

    x = cos(gamma)
    y = sin(gamma)
    u = cos(beta)
    v = (cos(alpha) - x * u) / y ! FIXME: potential underflow risk
    w = sqrt(1 - u * u - v * v)  ! FIXME: same above
    
    out_cell_vectors(1, 1) = cell_len(1)
    out_cell_vectors(2, 1) = 0.0
    out_cell_vectors(3, 1) = 0.0

    out_cell_vectors(1, 2) = cell_len(2) * x
    out_cell_vectors(2, 2) = cell_len(2) * y
    out_cell_vectors(3, 2) = 0.0

    out_cell_vectors(1, 3) = cell_len(3) * u
    out_cell_vectors(2, 3) = cell_len(3) * v
    out_cell_vectors(3, 3) = cell_len(3) * w

  end subroutine angles_to_cell_vector

  character(len=16) function itoa(x)
    integer, intent(in) :: x
    character(len=16) :: buf
    write(buf,"(I16)") x
    itoa = buf
  end function itoa
end module utility
  
