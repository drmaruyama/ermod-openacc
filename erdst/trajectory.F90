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

! module that governs trajectory I/O

module trajectory
  type handle
     integer(8) :: vmdhandle
  end type handle
 
contains
  subroutine init_trajectory()
    implicit none
    external vmdfio_init_traj
    call vmdfio_init_traj()
  end subroutine init_trajectory

  subroutine finish_trajectory()
    implicit none
    external vmdfio_fini_traj
    call vmdfio_fini_traj()
  end subroutine finish_trajectory

  ! Open trajectory and returns handle as htraj. 
  ! Should open fail, the program abends.
  subroutine open_trajectory(htraj, fname)
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    integer :: status
    external vmdfio_open_traj

    call vmdfio_open_traj(htraj%vmdhandle, fname, len_trim(fname), status)
    if(status /= 0) then
       stop "vmdfio_open_traj: unable to open trajectory. HISTORY must be a symlink"
    endif
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    external vmdfio_close_traj
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
   
    external vmdfio_read_traj_step

    if(kind(crd) /= 8) stop "vmdfio: write interfacing wrapper"

    call vmdfio_read_traj_step(htraj%vmdhandle, crd, cell, natom, status)
  end subroutine read_trajectory

end module trajectory
