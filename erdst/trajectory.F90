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

! module that governs trajectory I/O

module trajectory
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_char, c_float, c_int
  implicit none

  type handle
     ! Opaque handle to the VMD plugin's internal file-reader state.
     ! type(c_ptr) is the portable, standards-guaranteed way to hold a
     ! C pointer (it replaces the old integer(8), which only worked
     ! because pointers happened to be 8 bytes on the target platforms).
     type(c_ptr) :: vmdhandle = c_null_ptr
  end type handle

  ! Explicit, checked interfaces to the C routines in vmdfio.c.
  ! Symbol names (with the trailing underscore) match exactly what is
  ! already compiled into vmdfio.c, so no change to that file -- which
  ! also deals with the VMD plugin ABI -- is required here.
  interface
     subroutine vmdfio_init_traj() bind(c, name="vmdfio_init_traj_")
     end subroutine vmdfio_init_traj

     subroutine vmdfio_fini_traj() bind(c, name="vmdfio_fini_traj_")
     end subroutine vmdfio_fini_traj

     subroutine vmdfio_open_traj(vmdhandle, fname, fnamelen, status) &
          bind(c, name="vmdfio_open_traj_")
       import :: c_ptr, c_char, c_int
       type(c_ptr), intent(inout) :: vmdhandle
       character(kind=c_char), intent(in) :: fname(*)
       integer(c_int), intent(in) :: fnamelen
       integer(c_int), intent(out) :: status
     end subroutine vmdfio_open_traj

     subroutine vmdfio_close_traj(vmdhandle) bind(c, name="vmdfio_close_traj_")
       import :: c_ptr
       type(c_ptr), intent(inout) :: vmdhandle
     end subroutine vmdfio_close_traj

     subroutine vmdfio_read_traj_step(vmdhandle, xout, box, natoms, status) &
          bind(c, name="vmdfio_read_traj_step_")
       import :: c_ptr, c_float, c_int
       type(c_ptr), intent(in) :: vmdhandle
       real(c_float), intent(out) :: xout(*)
       real(c_float), intent(out) :: box(*)
       integer(c_int), intent(in) :: natoms
       integer(c_int), intent(out) :: status
     end subroutine vmdfio_read_traj_step
  end interface

contains
  subroutine init_trajectory()
    implicit none
    call vmdfio_init_traj()
  end subroutine init_trajectory

  subroutine finish_trajectory()
    implicit none
    call vmdfio_fini_traj()
  end subroutine finish_trajectory

  ! Open trajectory and returns handle as htraj. 
  ! Should open fail, the program abends.
  subroutine open_trajectory(htraj, fname)
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    character(kind=c_char) :: c_fname(len_trim(fname))
    integer(c_int) :: status
    integer :: i

    ! character(len=*) is not itself interoperable; marshal into a
    ! plain array of C characters (no NUL terminator needed, since the
    ! C side is given the explicit length and uses strncpy).
    do i = 1, len_trim(fname)
       c_fname(i) = fname(i:i)
    end do

    call vmdfio_open_traj(htraj%vmdhandle, c_fname, int(len_trim(fname), c_int), status)
    if (status /= 0) then
       stop "vmdfio_open_traj: unable to open trajectory. HISTORY must be a symlink"
    endif
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    call vmdfio_close_traj(htraj%vmdhandle)
  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, cell, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real, intent(out) :: crd(3, natom)
    real, intent(out) :: cell(3, 3)
    integer, intent(out) :: status

    ! The VMD plugin ABI (vmdfio.c) always speaks single precision,
    ! regardless of whether ERmod itself is built in single or double
    ! precision. Marshalling unconditionally through a c_float buffer
    ! here replaces the old "#ifdef DP" branch that duplicated this
    ! logic for the double-precision build.
    real(c_float) :: crd_tmp(3, natom)
    real(c_float) :: cell_tmp(3, 3)
    integer(c_int) :: c_status

    call vmdfio_read_traj_step(htraj%vmdhandle, crd_tmp, cell_tmp, int(natom, c_int), c_status)
    crd = real(crd_tmp, kind(crd))
    cell = real(cell_tmp, kind(cell))
    status = c_status
  end subroutine read_trajectory

end module trajectory
