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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION ""
#endif

subroutine enganal_init()
  use engmain, only: insorigin, insstructure, INSORG_REFSTR, INSSTR_RMSD
  use ptinsrt, only: load_refstructure
  use setconf, only: setparam
  use engproc, only: enginit
  implicit none

  call setparam
  call enginit
  ! getting the reference structure from file
  if ((insorigin == INSORG_REFSTR) .or. (insstructure == INSSTR_RMSD)) then
     call load_refstructure
  endif
end subroutine enganal_init

! FIXME: recover routine which runs as "combined with MD program"
!  connection to the main routine of trajectory generation is done in
!  setparam for parameter setting and getconf for configuration reading
subroutine enganal(stnum, nread)
  use engproc, only: engconst
  use mpiproc, only: nactiveproc
  implicit none
  integer, intent(in) :: stnum, nread
  nactiveproc = nread
  call engconst(stnum)
end subroutine enganal
!     
program trjmain
  use engmain, only: maxcnf, skpcnf, engdiv
  use OUTname, only: opentrj, closetrj, OUTinitial, initconf, finiconf
  use setconf, only: getconf_parallel
  use trajectory, only: init_trajectory, finish_trajectory
  use engproc, only: engclear, engstore, engproc_cleanup
  use mpiproc               ! MPI
  implicit none
  integer :: stnum, idiv, frames_per_div, nread, iframe

  call mpi_setup('init')
  if (myrank == 0) then
     print *, "ERmod " // PACKAGE_VERSION // ", Copyright (C) 2000-2024 Nobuyuki Matubayasi"
     print *, "                           2010-2024 Shun Sakuraba"
     print *, "                           2024-2024 Hidekazu Kojima"
     print *, "                           2024-2024 Yutaka Maruyama"
     print *, "ERmod comes with ABSOLUTELY NO WARRANTY."
     print *, "This is free software, and you can redistribute it"
     print *, "and/or modify it under certain conditions."
     print *, "See LICENSE file for details."
  end if
  call initconf()

  if (myrank == 0) then
     call init_trajectory()
     call opentrj()
  end if

  ! initialize
  call enganal_init()

  stnum = 0
  frames_per_div = maxcnf / skpcnf / engdiv
  if (frames_per_div <= 0) call halt_with_error("eng_par")

  do idiv = 1, engdiv
     call engclear

     do iframe = 1, frames_per_div, nprocs
        call getconf_parallel(frames_per_div - iframe + 1, nread)
        call enganal(stnum + myrank + 1, nread)
        stnum = stnum + nread
     end do
     call engstore(stnum)
  end do
  call engproc_cleanup

  if (myrank == 0) then
     call closetrj
     call finish_trajectory()
  end if

  call finiconf
  call mpi_setup('stop')    ! MPI
! stop  ! commented out to suppress harmless IEEE errors on Cygwin and MacOS X
end program trjmain
