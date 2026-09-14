! -*- F90 -*-
! ERmod - Energy Representation Module
! Copyright (C) 2000- The ERmod authors
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
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_ptr, c_loc
  use, intrinsic :: iso_fortran_env, only: real32, real64, int64
  use precision_kinds, only: wp
  implicit none

#ifndef HAVE_TRANSFER
  ! Explicit, checked interfaces to the C routines in hash_real.c.
  ! The symbol names below (with the trailing underscore) are exactly
  ! what is compiled into hash_real.c, so no change to that file is
  ! required; this only replaces the old untyped "external" declarations
  ! with a real interface.
  !
  ! NOTE: the dummy "v" is deliberately an untyped type(c_ptr), not
  ! real(c_double)/real(c_float). hash() below is called at several
  ! sites (see engproc.F90) by passing a 2D array section into a 1D
  ! explicit-shape dummy via sequence association; that means hash()
  ! itself has to stay a single, non-generic procedure taking a
  ! real(wp) array (wp being ERmod's single build-wide precision, see
  ! precision.F90). If the two C entry points were declared with
  ! explicit real(c_double)/real(c_float) dummies, the *textual* call
  ! to whichever one doesn't match wp would fail to compile, even
  ! though it is never executed. Routing through type(c_ptr) sidesteps
  ! that: both calls type-check unconditionally, and the correct one
  ! is chosen at run time exactly as before.
  interface
     subroutine hash_double_c(v, elms, hash_out) bind(c, name="hash_double_")
       import :: c_ptr, c_int, c_int64_t
       type(c_ptr), value :: v
       integer(c_int), intent(in) :: elms
       integer(c_int64_t), intent(out) :: hash_out
     end subroutine hash_double_c

     subroutine hash_float_c(v, elms, hash_out) bind(c, name="hash_float_")
       import :: c_ptr, c_int, c_int64_t
       type(c_ptr), value :: v
       integer(c_int), intent(in) :: elms
       integer(c_int64_t), intent(out) :: hash_out
     end subroutine hash_float_c
  end interface
#endif

contains

  integer(int64) function hash(arr, n) result(hash_out)
    integer, intent(in) :: n
    real(wp), intent(in) :: arr(n)
#ifdef HAVE_TRANSFER
    integer(int64) :: ret
    integer :: i
    ret = 0_int64
    do i = 1, n
       ret = ishftc(ret, 7)
       ret = ieor(ret, transfer(arr(i), ret))
    end do
    hash_out = ret
#else
    real(wp), target :: local_arr(n)
    integer(c_int64_t) :: ret

    local_arr = arr
    select case(kind(arr))
    case(real32)
       call hash_float_c(c_loc(local_arr), int(n, c_int), ret)
    case(real64)
       call hash_double_c(c_loc(local_arr), int(n, c_int), ret)
    case default
       stop "Error: hash(): unknown real type"
    end select
    hash_out = ret
#endif
  end function hash

  ! convert cell-length & (alpha, beta, gamma) to cell vectors
  subroutine angles_to_cell_vector(cell_len, angles, out_cell_vectors)
    use engmain, only: PI
    implicit none
    real(kind=8), intent(in) :: cell_len(3)
    real(kind=8), intent(in) :: angles(3)
    real(kind=8), intent(out) :: out_cell_vectors(3, 3)
    real(kind=8) :: alpha, beta, gamma, x, y, u, v, w

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

  pure character(len=16) function itoa(x)
    integer, intent(in) :: x
    character(len=16) :: buf
    write(buf,"(I16)") x
    itoa = buf
  end function itoa
end module utility
