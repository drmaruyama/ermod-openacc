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

module fft_iface
  use precision_kinds, only: wp
  use cufft
  implicit none

  integer :: fftsize(3)

  type fft_handle
     integer :: plan
  end type fft_handle

contains 

  subroutine fft_set_size(fftsize_in)
    integer, intent(in) :: fftsize_in(3)
    fftsize(:) = fftsize_in(:)
  end subroutine fft_set_size

  ! Check the return status of a cuFFT call and abort (on every rank, under
  ! MPI) with a diagnostic message if it did not succeed. Every cufftPlan3D/
  ! cufftExec*/cufftDestroy call below used to discard its status entirely,
  ! so a GPU/library-level failure (e.g. an unsupported transform size, or
  ! running out of device memory) would silently continue with whatever
  ! garbage was left in the output array, instead of stopping right away
  ! with a message that says what went wrong.
  !
  ! [location] should identify the failing call (e.g. "cufftPlan3D (fft_init_rtc)").
  ! Status codes are documented at:
  ! https://docs.nvidia.com/cuda/cufft/index.html#cufftresult
  ! They are spelled out here as literals rather than named constants,
  ! since not every "use cufft" module exposes all of them by name.
  subroutine check_cufft_status(stat, location)
    use engmain, only: stdout
    use mpiproc, only: mpi_abend
    implicit none
    integer, intent(in) :: stat
    character(len=*), intent(in) :: location
    character(len=64) :: msg

    if (stat == CUFFT_SUCCESS) return

    select case(stat)
    case(1);  msg = "invalid plan handle (CUFFT_INVALID_PLAN)"
    case(2);  msg = "failed to allocate GPU or CPU memory (CUFFT_ALLOC_FAILED)"
    case(3);  msg = "invalid transform type (CUFFT_INVALID_TYPE)"
    case(4);  msg = "invalid pointer or parameter (CUFFT_INVALID_VALUE)"
    case(5);  msg = "internal driver error (CUFFT_INTERNAL_ERROR)"
    case(6);  msg = "failed to execute the FFT on the GPU (CUFFT_EXEC_FAILED)"
    case(7);  msg = "the cuFFT library failed to initialize (CUFFT_SETUP_FAILED)"
    case(8);  msg = "invalid transform size (CUFFT_INVALID_SIZE)"
    case(9);  msg = "unaligned data (CUFFT_UNALIGNED_DATA)"
    case(10); msg = "missing parameters in call (CUFFT_INCOMPLETE_PARAMETER_LIST)"
    case(11); msg = "plan executed on a different GPU than it was created on (CUFFT_INVALID_DEVICE)"
    case(12); msg = "internal plan database error (CUFFT_PARSE_ERROR)"
    case(13); msg = "no workspace provided prior to plan execution (CUFFT_NO_WORKSPACE)"
    case(14); msg = "functionality not implemented for the given parameters (CUFFT_NOT_IMPLEMENTED)"
    case(15); msg = "license error (CUFFT_LICENSE_ERROR)"
    case(16); msg = "operation not supported for the given parameters (CUFFT_NOT_SUPPORTED)"
    case default; msg = "unrecognized cuFFT status code"
    end select

    write(stdout, "(A,A,A,I0,A,A,A)") " cuFFT error in ", trim(location), ": status = ", stat, " (", trim(msg), ")"
    call mpi_abend()
    stop "cuFFT call failed"
  end subroutine check_cufft_status

  ! 3D-FFT, cufft version
  
  subroutine fft_init_rtc(handle, in, out)
    type(fft_handle), intent(out) :: handle
    real(wp), intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex(wp), intent(out) :: out(fftsize(1)/2+1, fftsize(2), fftsize(3))
    integer :: stat
#ifdef DP
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_D2Z)
#else
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_R2C)
#endif
    call check_cufft_status(stat, "cufftPlan3D (fft_init_rtc)")
  end subroutine fft_init_rtc

  subroutine fft_init_ctr(handle, in, out)
    type(fft_handle), intent(out) :: handle
    complex(wp), intent(in) :: in(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real(wp), intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
    integer :: stat
#ifdef DP
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_Z2D)
#else
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_C2R)
#endif
    call check_cufft_status(stat, "cufftPlan3D (fft_init_ctr)")
  end subroutine fft_init_ctr

  subroutine fft_rtc(handle, in, out)
    use openacc
    use cufft
    type(fft_handle), intent(in) :: handle
    real(wp), intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex(wp), intent(out) :: out(fftsize(1)/2+1, fftsize(2), fftsize(3))
    integer :: stat
    !$acc data present(in, out)
    !$acc host_data use_device(in, out)
#ifdef DP
    stat = cufftExecD2Z(handle%plan, in, out)
#else
    stat = cufftExecR2C(handle%plan, in, out)
#endif
    !$acc end host_data
    !$acc end data
    call check_cufft_status(stat, "cufftExecR2C/D2Z (fft_rtc)")
  end subroutine fft_rtc

  subroutine fft_ctr(handle, in, out)
    type(fft_handle), intent(in) :: handle
    complex(wp), intent(in) :: in(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real(wp), intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
    integer :: stat
    !$acc data present(in, out)
    !$acc host_data use_device(in, out)
#ifdef DP
    stat = cufftExecZ2D(handle%plan, in, out)
#else
    stat = cufftExecC2R(handle%plan, in, out)
#endif
    !$acc end host_data
    !$acc end data
    call check_cufft_status(stat, "cufftExecC2R/Z2D (fft_ctr)")
  end subroutine fft_ctr

  subroutine fft_cleanup_rtc(handle)
    type(fft_handle), intent(in) :: handle
    integer :: stat
    stat = cufftDestroy(handle%plan)
    call check_cufft_status(stat, "cufftDestroy (fft_cleanup_rtc)")
  end subroutine fft_cleanup_rtc

  subroutine fft_cleanup_ctr(handle)
    type(fft_handle), intent(in) :: handle
    integer :: stat
    stat = cufftDestroy(handle%plan)
    call check_cufft_status(stat, "cufftDestroy (fft_cleanup_ctr)")
  end subroutine fft_cleanup_ctr

end module fft_iface
