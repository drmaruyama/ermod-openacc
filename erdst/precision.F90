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

! Single point of control for ERmod's floating-point working precision.
!
! Every real/complex declaration in ERmod is meant to spell out its kind
! explicitly as "real(wp)" / "complex(wp)", using the parameter defined
! here, instead of a bare "real"/"complex" whose meaning used to depend
! on a compiler flag (e.g. nvfortran's -Mr8) silently promoting the
! default real kind to 8 bytes. That flag is no longer needed: the
! precision is switched by the -DDP preprocessor definition alone, which
! configure.ac already adds for --enable-double.
module precision_kinds
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none

#ifdef DP
  integer, parameter :: wp = real64
#else
  integer, parameter :: wp = real32
#endif

end module precision_kinds
