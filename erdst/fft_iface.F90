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

  ! 3D-FFT, cufft version
  
  subroutine fft_init_rtc(handle, in, out)
    type(fft_handle), intent(out) :: handle
    real, intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex, intent(out) :: out(fftsize(1)/2+1, fftsize(2), fftsize(3))
    integer :: stat
#ifdef DP
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_D2Z)
#else
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_R2C)
#endif
  end subroutine fft_init_rtc

  subroutine fft_init_ctr(handle, in, out)
    type(fft_handle), intent(out) :: handle
    complex, intent(in) :: in(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real, intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
    integer :: stat
#ifdef DP
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_Z2D)
#else
    stat = cufftPlan3D(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         CUFFT_C2R)
#endif
  end subroutine fft_init_ctr

  subroutine fft_rtc(handle, in, out)
    use openacc
    use cufft
    type(fft_handle), intent(in) :: handle
    real, intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex, intent(out) :: out(fftsize(1)/2+1, fftsize(2), fftsize(3))
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
  end subroutine fft_rtc

  subroutine fft_ctr(handle, in, out)
    type(fft_handle), intent(in) :: handle
    complex, intent(in) :: in(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real, intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
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
  end subroutine fft_ctr

  subroutine fft_cleanup_rtc(handle)
    type(fft_handle), intent(in) :: handle
    integer :: stat
    stat = cufftDestroy(handle%plan)
  end subroutine fft_cleanup_rtc

  subroutine fft_cleanup_ctr(handle)
    type(fft_handle), intent(in) :: handle
    integer :: stat
    stat = cufftDestroy(handle%plan)
  end subroutine fft_cleanup_ctr

end module fft_iface
