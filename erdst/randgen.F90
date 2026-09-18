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

! Pure-Fortran port of xoshiro256** 1.0 (Blackman & Vigna, public domain),
! replacing the previous bind(C) wrapper around xoshiro256ss.c. This file
! reproduces xoshiro256ss.c bit-for-bit (verified against it directly: same
! seed, same warm-up, same next()/next_double()/next_int31()/jump()/
! long_jump() outputs) but needs no C compiler, no interoperability layer,
! and therefore no compiler-specific behavior to worry about.
!
! Fortran has no unsigned integer type, which is why the previous version
! of this file called out to C instead of implementing xoshiro directly.
! That turns out not to be a real obstacle: xoshiro256** only ever uses
! its 64-bit words as bit patterns -- XOR, AND, shifts, rotates,
! multiplication and addition -- and never compares or divides them as
! numbers. On every mainstream Fortran compiler, INTEGER(8) arithmetic
! silently wraps using two's-complement, which is bit-for-bit identical to
! unsigned 64-bit arithmetic modulo 2^64; Fortran's ISHFT is a *logical*
! (zero-filling) shift regardless of sign, matching C's shifts on an
! unsigned type; and ISHFTC is exactly the "rotl" used throughout. So a
! bit-exact port is possible using only standard, portable Fortran.
!
! As a side effect, this also fixes the old "can't copy/assign randstate
! variables" limitation noted in the previous version of this file: since
! the state is now a plain array of integers instead of a C pointer to
! heap-allocated memory, ordinary Fortran assignment (copy = s) works
! correctly with no aliasing or double-free concerns.
module randgen
  use, intrinsic :: iso_fortran_env, only: int64, real64
  implicit none

  type, public :: randstate
     integer(int64) :: s(4) = 0_int64
  contains
     procedure, public :: init, dtor
     procedure, public :: next_double, next_real, next_int31
     procedure, public :: jump, long_jump
  end type randstate

  private
  public :: randstate

contains

  ! splitmix64, used only to expand a single 64-bit seed into four
  ! 64-bit words of well-mixed xoshiro256** initial state (same
  ! approach as splitmix64_next() in xoshiro256ss.c).
  function splitmix64_next(x) result(z)
    integer(int64), intent(inout) :: x
    integer(int64) :: z

    x = x + int(z'9E3779B97F4A7C15', int64)
    z = x
    z = ieor(z, ishft(z, -30))
    z = z * int(z'BF58476D1CE4E5B9', int64)
    z = ieor(z, ishft(z, -27))
    z = z * int(z'94D049BB133111EB', int64)
    z = ieor(z, ishft(z, -31))
  end function splitmix64_next

  ! Core xoshiro256** generator step: returns the next 64-bit output
  ! and advances this%s in place. Mirrors xoshiro256ss_next() exactly,
  ! with this%s(1:4) corresponding to C's s[0..3].
  function xoshiro_next(this) result(res)
    class(randstate), intent(inout) :: this
    integer(int64) :: res
    integer(int64) :: t

    res = ishftc(this%s(2) * 5_int64, 7) * 9_int64

    t = ishft(this%s(2), 17)

    this%s(3) = ieor(this%s(3), this%s(1))
    this%s(4) = ieor(this%s(4), this%s(2))
    this%s(2) = ieor(this%s(2), this%s(3))
    this%s(1) = ieor(this%s(1), this%s(4))

    this%s(3) = ieor(this%s(3), t)

    this%s(4) = ishftc(this%s(4), 45)
  end function xoshiro_next

  subroutine init(this, seed)
    class(randstate), intent(inout) :: this
    integer(int64), intent(in) :: seed
    integer(int64) :: x, discard
    integer :: i

    x = seed
    this%s(1) = splitmix64_next(x)
    this%s(2) = splitmix64_next(x)
    this%s(3) = splitmix64_next(x)
    this%s(4) = splitmix64_next(x)

    ! warm-up, matching xoshiro256ss_init_state_with_seed()
    do i = 1, 4
       discard = xoshiro_next(this)
    end do
  end subroutine init

  ! No-op: kept only for interface compatibility with the previous,
  ! heap-backed implementation. There is no longer any resource to
  ! release (this%s is a plain, stack/derived-type-resident array).
  subroutine dtor(this)
    class(randstate), intent(inout) :: this
    this%s = 0_int64
  end subroutine dtor

  ! returns [0, 1)
  function next_double(this) result(r)
    class(randstate), intent(inout) :: this
    real(real64) :: r

    r = real(ishft(xoshiro_next(this), -11), real64) * 2.0_real64 ** (-53)
  end function next_double

  function next_real(this) result(r)
    class(randstate), intent(inout) :: this
    real :: r

    r = real(next_double(this))
  end function next_real

  ! returns 0 <= x < 2^31
  function next_int31(this) result(r)
    class(randstate), intent(inout) :: this
    integer :: r

    r = int(iand(xoshiro_next(this), int(z'7FFFFFFF', int64)), kind(r))
  end function next_int31

  ! Shared implementation of jump()/long_jump(): advance this%s as if
  ! next() had been called an astronomical (2^128 / 2^192) number of
  ! times, using the standard "polynomial in GF(2)[x]" jump-ahead
  ! technique -- see xoshiro256ss_jump()/xoshiro256ss_long_jump() for
  ! the reference this mirrors.
  subroutine jump_with_poly(this, poly)
    class(randstate), intent(inout) :: this
    integer(int64), intent(in) :: poly(4)
    integer(int64) :: s0, s1, s2, s3, discard
    integer :: i, b

    s0 = 0_int64 ; s1 = 0_int64 ; s2 = 0_int64 ; s3 = 0_int64
    do i = 1, 4
       do b = 0, 63
          if (btest(poly(i), b)) then
             s0 = ieor(s0, this%s(1))
             s1 = ieor(s1, this%s(2))
             s2 = ieor(s2, this%s(3))
             s3 = ieor(s3, this%s(4))
          end if
          discard = xoshiro_next(this)
       end do
    end do
    this%s(1) = s0
    this%s(2) = s1
    this%s(3) = s2
    this%s(4) = s3
  end subroutine jump_with_poly

  ! Equivalent to 2^128 calls to next(): generates 2^128 non-overlapping
  ! subsequences for parallel computations.
  subroutine jump(this)
    class(randstate), intent(inout) :: this
    integer(int64), parameter :: JUMP_POLY(4) = [ &
         int(z'180EC6D33CFD0ABA', int64), int(z'D5A61266F0C9392C', int64), &
         int(z'A9582618E03FC9AA', int64), int(z'39ABDC4529B1661C', int64) ]

    call jump_with_poly(this, JUMP_POLY)
  end subroutine jump

  ! Equivalent to 2^192 calls to next(): generates 2^64 starting points,
  ! from each of which jump() generates 2^64 non-overlapping
  ! subsequences for parallel distributed computations (e.g. one
  ! long_jump() per MPI rank, as in insertion.F90's urand_init()).
  subroutine long_jump(this)
    class(randstate), intent(inout) :: this
    integer(int64), parameter :: LONG_JUMP_POLY(4) = [ &
         int(z'76E15D3EFEFDCBBF', int64), int(z'C5004E441C522FB3', int64), &
         int(z'77710069854EE241', int64), int(z'39109BB02ACBE635', int64) ]

    call jump_with_poly(this, LONG_JUMP_POLY)
  end subroutine long_jump

end module randgen
