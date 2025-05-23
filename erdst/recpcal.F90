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

module reciprocal
  use fft_iface, only: fft_handle
  implicit none
  integer :: rc1min, rc1max, rc2min, rc2max, rc3min, rc3max
  integer :: ccesize, ccemax
  integer, allocatable :: slvtag(:)
  real,    allocatable :: engfac(:,:,:)
  real,    allocatable :: gf_b(:)
  complex, allocatable :: rcpslt(:,:,:)
  real,    allocatable :: splslv(:,:,:)
  integer, allocatable :: grdslv(:,:)
  real,    allocatable :: cnvslt(:,:,:), cnvslt_r(:,:,:,:)
  real,    allocatable :: splfc1(:), splfc2(:), splfc3(:)
  complex, allocatable :: fft_buf(:, :, :)

  real,    allocatable :: solute_self_energy_refs(:)
  real :: solute_self_energy

  type(fft_handle) :: handle_c2r, handle_r2c

contains
  subroutine recpcal_init(slvmax, tagpt)
    use engmain, only:  nummol, numsite, splodr, ms1max, ms2max, ms3max, &
         maxins, slttype, SLT_SOLN
    use spline, only: spline_init
    use fft_iface, only: fft_init_ctr, fft_init_rtc, fft_set_size
    implicit none
    integer, intent(in) :: slvmax, tagpt(:)
    integer :: m, k
    integer :: gridsize(3), ptrnk

    allocate( slvtag(nummol) )
    !$acc enter data create(slvtag)
    slvtag(:) = -1
    ptrnk = 0
    do k = 1, slvmax
       m = tagpt(k)
       slvtag(m) = ptrnk + 1
       ptrnk = ptrnk + numsite(m)
    enddo
    !$acc update device(slvtag)

    allocate(gf_b(splodr)) ! for PPPM Green's function
    call calc_gfb_pppm()   ! calc gf_b (function of splodr)
    rc1min = 0 ; rc1max = ms1max - 1
    rc2min = 0 ; rc2max = ms2max - 1
    rc3min = 0 ; rc3max = ms3max - 1
    ccesize = ms1max / 2 + 1; ccemax = ccesize - 1
    call spline_init(splodr)
    allocate(splslv(0:splodr-1, 3, ptrnk), grdslv(3, ptrnk))
    !$acc enter data create(splslv, grdslv)
    if (slttype == SLT_SOLN) then
       allocate(cnvslt(rc1min:rc1max, rc2min:rc2max, rc3min:rc3max))
       !$acc enter data create(cnvslt)
    else
       allocate(solute_self_energy_refs(maxins))
       allocate(cnvslt_r(rc1min:rc1max, rc2min:rc2max, rc3min:rc3max, maxins))
       !$acc enter data create(cnvslt_r)
    end if
    ! initialize spline table for all axes
    allocate( splfc1(rc1min: rc1max) )
    allocate( splfc2(rc2min: rc2max) )
    allocate( splfc3(rc3min: rc3max) )
    call init_spline_axis(rc1min, rc1max, splfc1(rc1min:rc1max))
    call init_spline_axis(rc2min, rc2max, splfc2(rc2min:rc2max))
    call init_spline_axis(rc3min, rc3max, splfc3(rc3min:rc3max))
    gridsize(1) = ms1max
    gridsize(2) = ms2max
    gridsize(3) = ms3max
    call fft_set_size(gridsize)
    allocate( engfac(rc1min:ccemax, rc2min:rc2max, rc3min:rc3max) )
    allocate( rcpslt(rc1min:ccemax, rc2min:rc2max, rc3min:rc3max) )
    !$acc enter data create(engfac, rcpslt)
    ! init fft
    call fft_init_rtc(handle_r2c, cnvslt, rcpslt)
    call fft_init_ctr(handle_c2r, rcpslt, cnvslt)
  end subroutine recpcal_init

  subroutine init_spline_axis(imin, imax, splfc)
    use engmain, only: splodr, PI
    use spline, only: spline_value
    implicit none
    integer, intent(in) :: imin, imax
    real, intent(out) :: splfc(imin:imax)
    real :: chr, factor, rtp2
    real :: cosk, sink
    complex :: rcpi
    integer :: rci, spi
    do rci = imin, imax
       rcpi = (0.0, 0.0)
       do spi = 0, splodr - 2
          chr = spline_value(real(spi + 1))
          rtp2 = 2.0 * PI * real(spi * rci) / real(imax + 1)
          cosk = chr * cos(rtp2)
          sink = chr * sin(rtp2)
          rcpi = rcpi + cmplx(cosk, sink)
       end do
       factor = real(rcpi * conjg(rcpi))
       splfc(rci) = factor
    end do
  end subroutine init_spline_axis

  subroutine recpcal_spline_greenfunc()
    use engmain, only: invcl, ms1max, ms2max, ms3max, splodr, volume, screen, PI
    implicit none
    integer :: rc1, rc2, rc3, rci, m, rcimax
    real :: factor, rtp2, chr
    real :: inm(3), xst(3)
    do rc3 = rc3min, rc3max
       do rc2 = rc2min, rc2max
          do rc1 = rc1min, ccemax
             factor = 0.0
             if (rc1 == 0 .and. rc2 == 0 .and. rc3 == 0) cycle
             do m = 1, 3
                if (m == 1) rci = rc1
                if (m == 2) rci = rc2
                if (m == 3) rci = rc3

                if (m == 1) rcimax = ms1max
                if (m == 2) rcimax = ms2max
                if (m == 3) rcimax = ms3max

                if ((mod(splodr, 2) == 1) .and. (2*abs(rci) == rcimax)) goto 3219
                if (rci <= rcimax / 2) then
                   inm(m) = real(rci)
                else
                   inm(m) = real(rci - rcimax)
                endif
             end do
             do m = 1, 3
                xst(m) = dot_product(invcl(:, m), inm(:))
             end do
             rtp2 = sum(xst(1:3) ** 2)
             chr = (PI ** 2) * rtp2 / (screen ** 2 )
             factor = exp(-chr) / rtp2 / PI / volume
             rtp2 = splfc1(rc1) * splfc2(rc2) * splfc3(rc3)
             factor = factor / rtp2
3219         continue
             engfac(rc1, rc2, rc3) = factor
          end do
       end do
    end do
    engfac(0, 0, 0) = 0.0
    !$acc update device(engfac)
  end subroutine recpcal_spline_greenfunc

  ! implemented only for PPPM
  function factorial(n) result(fact)
    implicit none
    integer, intent(in) :: n
    integer(kind=8) :: fact, i
    if (n < 0) stop "Error"
    fact = 1
    if (n > 0) then
       do i = 1, n
          fact = fact * i
       end do
    end if
  end function factorial

  ! calc gf_b needed by denominator in Green's func for PPPM
  subroutine calc_gfb_pppm()
    use engmain, only: splodr
    implicit none
    integer :: l, m
    real :: gaminv

    gf_b(2:splodr) = 0.
    gf_b(1) = 1.
    do m = 1, splodr-1
       do l = m, 1, -1
          gf_b(l+1)=4.*(gf_b(l+1)*(l-m)*(l-m-0.5)-gf_b(l)*(l-m-1)*(l-m-1))
       end do
       gf_b(1)=4.*(gf_b(1)*(l-m)*(l-m-0.5))
    end do
    gaminv = 1./ factorial(2*splodr-1)
    do m = 1, splodr
       gf_b(m) = gf_b(m) * gaminv
    end do
  end subroutine calc_gfb_pppm

  ! calc denominator in Green's func for PPPM
  function gf_denom(x,y,z) result(denom)
    use engmain, only: splodr
    implicit none
    real, intent(in) :: x, y, z
    real :: s(3), denom
    integer :: l
    s(1:3) = 0.
    do l = splodr, 1, -1
       s(1) = gf_b(l) + s(1) * x
       s(2) = gf_b(l) + s(2) * y
       s(3) = gf_b(l) + s(3) * z
    end do
    denom = product(s)**2
  end function gf_denom
    
  ! calc green's function (optimal influence function) for PPPM
  subroutine recpcal_pppm_greenfunc()
    use engmain, only: &
         invcl, ms1max, ms2max, ms3max, splodr, volume, screen, PI, cell
    implicit none
    integer :: rc1, rc2, rc3, i
    integer :: mx, my, mz, mx_max, my_max, mz_max, mN(3)
    real :: factor
    real :: inm(3) ! folded rc1, rc2 or rc3
    real :: k(3), sin2kh2(3), k2, sum_dru2, tmp, kmk(3), km(3)
    real :: gam(3), m_max(3), ukm(3)
    real,parameter :: EPS_HOC = 0.0000001 ! this value is adopted in LAMMPS 10Feb15
    m_max(1) = screen/(PI*ms1max)*((-log(EPS_HOC))**0.25) ! LAMMPS' form
    m_max(2) = screen/(PI*ms2max)*((-log(EPS_HOC))**0.25)
    m_max(3) = screen/(PI*ms3max)*((-log(EPS_HOC))**0.25)
    mx_max = int(dot_product(cell(:, 1), m_max(:)))
    my_max = int(dot_product(cell(:, 2), m_max(:)))
    mz_max = int(dot_product(cell(:, 3), m_max(:)))

!    write (6,*) "recpcal_pppm_spline_greenfunc"
!    write (6,*) m_max
!    write (6,*) cell(:,1)
!    write (6,*) cell(:,2)
!    write (6,*) cell(:,3)
!    write (6,*) dot_product(cell(:, 1), m_max(:))
!    write (6,*) 'Calculating PPPM Green function...'
!    write (6,*) ' Range for m vector:'
!    write (6,*) '  [-', mx_max,':',mx_max,']'
!    write (6,*) '  [-', my_max,':',my_max,']'
!    write (6,*) '  [-', mz_max,':',mz_max,']'

    if (mx_max < 0 .or. my_max < 0 .or. mz_max < 0) &
         stop "cannot set range of m vectors for Green\'s function correctly"
    if (mod(splodr,2) == 1) &
         stop "interpolation order for PPPM must be even"

    do rc3 = rc3min, rc3max
       if (rc3 <= ms3max/2) then
          inm(3) = rc3
       else
          inm(3) = rc3 - ms3max
       end if
       sin2kh2(3) = (sin(PI*inm(3)/ms3max))**2 !! sin^2(kz*hz/2)

       do rc2 = rc2min, rc2max
          if (rc2 <= ms2max/2) then
             inm(2) = rc2
          else
             inm(2) = rc2 - ms2max
          end if
          sin2kh2(2) = (sin(PI*inm(2)/ms2max))**2 !! sin^2(ky*hy/2)

          do rc1 = rc1min, ccemax
             if (rc1 <= ms1max/2) then
                inm(1) = rc1
             else
                inm(1) = rc1 - ms1max
             end if
             sin2kh2(1) = (sin(PI*inm(1)/ms1max))**2 !! sin^2(kx*hx/2)

             factor = 0.0
             if (rc1 == 0 .and. rc2 == 0 .and. rc3 == 0) cycle
             do i = 1, 3
                k(i) = 2 * PI * dot_product(invcl(:, i), inm(:)) ! k = 2*pi*m
             end do
             k2 = dot_product(k,k) ! k^2

             if (k2 > 0.) then
                sum_dru2 = 0.

                do mx = -mx_max, mx_max
                   tmp = PI*(inm(1)/ms1max+mx)
                   if (tmp == 0.) then
                      ukm(1) = 1.
                   else
                      ukm(1) = (sin(tmp)/tmp)**(2*splodr) ! product(ukm) = U(km)^2
                   end if
                   mN(1) = mx*ms1max

                   do my = -my_max, my_max
                      tmp = PI*(inm(2)/ms2max+my)
                      if (tmp == 0.) then
                         ukm(2) = 1.
                      else
                         ukm(2) = (sin(tmp)/tmp)**(2*splodr)
                      end if
                      mN(2) = my*ms2max

                      do mz = -mz_max, mz_max
                         tmp = PI*(inm(3)/ms3max+mz)
                         if (tmp == 0.) then
                            ukm(3) = 1.
                         else
                            ukm(3) = (sin(tmp)/tmp)**(2*splodr)
                         end if
                         mN(3) = mz*ms2max

                         do i = 1, 3
                            kmk(i) = 2 * PI * dot_product(invcl(:, i), mN(:))
                            ! km - k
                         end do
                         km(1:3) = k(1:3) + kmk(1:3) ! km
                         gam(1:3) = exp(-0.25*(km(1:3)/screen)**2)
                         ! product(gam) = gamma(km)

                         sum_dru2 = sum_dru2 &
                              + dot_product(k,km)*product(gam)*product(ukm) &
                              /dot_product(km,km)
                      end do
                   end do
                end do
                factor = sum_dru2*4.*PI &
                     / (k2*gf_denom(sin2kh2(1),sin2kh2(2),sin2kh2(3))) &
                     / volume
             end if
             engfac(rc1, rc2, rc3) = factor
          end do
       end do
    end do
    engfac(0, 0, 0) = 0.0
    !$acc update device(engfac)
  end subroutine recpcal_pppm_greenfunc

  ! note: this routine is named as "solvent", but may include solute molecule, when mutiple solute is used.
  subroutine recpcal_prepare_solvent(tagpt, slvmax)
    use engmain, only: numsite
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: tagpt(:), slvmax
    integer :: i, k, svi, stmax

    do k = 1, slvmax
       i = tagpt(k)
       svi = slvtag(i)
       if (svi <= 0) call halt_with_error('rcp_cns')

       stmax = numsite(i)
       call calc_spline_molecule(i, stmax, splslv(:,:,svi:svi+stmax-1), &
            grdslv(:,svi:svi+stmax-1))
    end do
    !$acc update device(splslv, grdslv)
  end subroutine recpcal_prepare_solvent

  subroutine recpcal_prepare_solute(tagslt)
    use engmain, only: ms1max, ms2max, ms3max, sitepos, invcl, numsite, splodr, charge, mol_begin_index
    use fft_iface, only: fft_ctr, fft_rtc
    implicit none
    integer, intent(in) :: tagslt
    integer :: i, j, k
    integer :: rc1, rc2, rc3, sid, ati, cg1, cg2, cg3, stmax
    real :: factor, chr
    real, allocatable, save :: splval(:,:,:)
    integer, allocatable, save :: grdval(:,:)
    logical, save :: initialized = .false.

    stmax = numsite(tagslt)
    if(.not. initialized) then
       allocate( splval(0:splodr-1, 3, stmax), grdval(3, stmax) )
       initialized = .true.
    end if
    call calc_spline_molecule(tagslt, stmax, splval(:,:,1:stmax), grdval(:,1:stmax))
    !$acc parallel present(cnvslt)
    cnvslt = 0.0
    !$acc end parallel
    !$acc parallel loop present(mol_begin_index, charge, cnvslt, rcpslt)
    do sid = 1, stmax
!       ati = specatm(sid, tagslt)
       ati = mol_begin_index(tagslt) + (sid - 1)
       chr = charge(ati)
       do cg3 = 0, splodr - 1
          do cg2 = 0, splodr - 1
             do cg1 = 0, splodr - 1
                rc1 = modulo(grdval(1, sid) - cg1, ms1max)
                rc2 = modulo(grdval(2, sid) - cg2, ms2max)
                rc3 = modulo(grdval(3, sid) - cg3, ms3max)
                factor = chr * splval(cg1, 1, sid) * splval(cg2, 2, sid) &
                     * splval(cg3, 3, sid)
                !$acc atomic update
                cnvslt(rc1, rc2, rc3) = cnvslt(rc1, rc2, rc3) + factor
             end do
          end do
       end do
    end do
    !$acc end parallel

    call fft_rtc(handle_r2c, cnvslt, rcpslt)                         ! 3D-FFT
    !$acc update self(rcpslt)

    ! original form is:
    ! 0.5 * sum(engfac(:, :, :) * real(rcpslt_c(:, :, :)) * conjg(rcpslt_c(:, :, :)))
    ! where rcpslt_c(rc1, rc2, rc3) = conjg(rcpslt_buf(ms1max - rc1, ms2max - rc2, ms3max - rc3))
    ! Here we use symmetry of engfac to calculate efficiently
    if (mod(ms1max, 2) == 0) then
       solute_self_energy = &
            sum(engfac(1:(ccemax-1), :, :) * real(rcpslt(1:(ccemax-1), :, :) * conjg(rcpslt(1:(ccemax-1), :, :)))) + &
            0.5 * sum(engfac(0,      :, :) * real(rcpslt(0,      :, :) * conjg(rcpslt(0,      :, :)))) + &
            0.5 * sum(engfac(ccemax, :, :) * real(rcpslt(ccemax, :, :) * conjg(rcpslt(ccemax, :, :))))
    else
       solute_self_energy = &
            sum(engfac(1:ccemax, :, :) * real(rcpslt(1:ccemax, :, :) * conjg(rcpslt(1:ccemax, :, :)))) + &
            0.5 * sum(engfac(0,      :, :) * real(rcpslt(0,      :, :) * conjg(rcpslt(0,      :, :))))
    endif

    !$acc parallel loop present (engfac, rcpslt)
    do concurrent (i = rc1min:ccemax, j = rc2min:rc2max, k = rc3min:rc3max)
       rcpslt(i, j, k) = engfac(i, j, k) * rcpslt(i, j, k)
    end do
    !$acc end parallel

    call fft_ctr(handle_c2r, rcpslt, cnvslt)                    ! 3D-FFT

  end subroutine recpcal_prepare_solute

  subroutine recpcal_prepare_solute_refs(tagslt, maxdst)
    use engmain, only: ms1max, ms2max, ms3max, sitepos, invcl, numsite, splodr, charge, mol_begin_index
    use fft_iface, only: fft_ctr, fft_rtc
    implicit none
    integer, intent(in) :: tagslt, maxdst
    integer :: i, j, k, cnt
    integer :: rc1, rc2, rc3, sid, ati, cg1, cg2, cg3, stmax
    real :: factor, chr
    real, allocatable, save :: splval(:,:,:)
    integer, allocatable, save :: grdval(:,:)
    logical, save :: initialized = .false.

    stmax = numsite(tagslt)
    if(.not. initialized) then
       allocate( splval(0:splodr-1, 3, stmax), grdval(3, stmax) )
       initialized = .true.
    end if
    !$acc parallel present(cnvslt_r)
    cnvslt_r = 0.0
    !$acc end parallel
    do cnt = 1, maxdst

       call calc_spline_molecule_refs(tagslt, cnt, stmax, &
            splval(:,:,1:stmax), grdval(:,1:stmax))
       !$acc parallel loop present(mol_begin_index, charge, cnvslt_r, rcpslt)
       do sid = 1, stmax
          ! ati = specatm(sid, tagslt)
          ati = mol_begin_index(tagslt) + (sid - 1)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             do cg2 = 0, splodr - 1
                do cg1 = 0, splodr - 1
                   rc1 = modulo(grdval(1, sid) - cg1, ms1max)
                   rc2 = modulo(grdval(2, sid) - cg2, ms2max)
                   rc3 = modulo(grdval(3, sid) - cg3, ms3max)
                   factor = chr * splval(cg1, 1, sid) * splval(cg2, 2, sid) &
                        * splval(cg3, 3, sid)
                   !$acc atomic update
                   cnvslt_r(rc1, rc2, rc3, cnt) = &
                        cnvslt_r(rc1, rc2, rc3, cnt) + factor
                end do
             end do
          end do
       end do
       !$acc end parallel

       call fft_rtc(handle_r2c, cnvslt_r(:,:,:,cnt), rcpslt)
       !$acc update self(rcpslt)

       ! original form is:
       ! 0.5 * sum(engfac(:, :, :) * real(rcpslt_c(:, :, :)) * conjg(rcpslt_c(:, :, :)))
       ! where rcpslt_c(rc1, rc2, rc3) = conjg(rcpslt_buf(ms1max - rc1, ms2max - rc2, ms3max - rc3))
       ! Here we use symmetry of engfac to calculate efficiently
       if (mod(ms1max, 2) == 0) then
          solute_self_energy_refs(cnt) = &
               sum(engfac(1:(ccemax-1), :, :) * real(rcpslt(1:(ccemax-1), :, :) * conjg(rcpslt(1:(ccemax-1), :, :)))) + &
               0.5 * sum(engfac(0,      :, :) * real(rcpslt(0,      :, :) * conjg(rcpslt(0,      :, :)))) + &
               0.5 * sum(engfac(ccemax, :, :) * real(rcpslt(ccemax, :, :) * conjg(rcpslt(ccemax, :, :))))
       else
          solute_self_energy_refs(cnt) = &
               sum(engfac(1:ccemax, :, :) * real(rcpslt(1:ccemax, :, :) * conjg(rcpslt(1:ccemax, :, :)))) + &
               0.5 * sum(engfac(0,      :, :) * real(rcpslt(0,      :, :) * conjg(rcpslt(0,      :, :))))
       endif

       !$acc parallel loop present (engfac, rcpslt)
       do concurrent (i = rc1min:ccemax, j = rc2min:rc2max, k = rc3min:rc3max)
          rcpslt(i, j, k) = engfac(i, j, k) * rcpslt(i, j, k)
       end do
       !$acc end parallel

       call fft_ctr(handle_c2r, rcpslt, cnvslt_r(:,:,:,cnt))

    end do
  end subroutine recpcal_prepare_solute_refs

  subroutine calc_spline_molecule(imol, stmax, store_spline, store_grid)
    use engmain, only: ms1max, ms2max, ms3max, splodr, specatm, sitepos, invcl
    use spline, only: spline_value
    implicit none

    integer, intent(in) :: imol, stmax
    real, intent(out) :: store_spline(0:splodr-1, 3, 1:stmax)
    integer, intent(out) :: store_grid(3, 1:stmax)
    
    integer :: sid, ati, rcimax, m, k, rci, spi
    real :: xst(3), inm(3)
    real :: factor, rtp2

    do sid = 1, stmax
       ati = specatm(sid, imol)
       xst(:) = sitepos(:, ati)
       do k = 1, 3
          factor = dot_product(invcl(k,:), xst(:))
          factor = factor - floor(factor)
          inm(k) = factor
       end do
       do m = 1, 3
          if (m == 1) rcimax = ms1max
          if (m == 2) rcimax = ms2max
          if (m == 3) rcimax = ms3max
          factor = inm(m) * real(rcimax)
          rci = int(factor)
          do spi = 0, splodr - 1
             rtp2 = factor - real(rci - spi)
             store_spline(spi, m, sid) = spline_value(rtp2)
          end do
          store_grid(m, sid) = rci
       end do
    end do
  end subroutine calc_spline_molecule

  subroutine calc_spline_molecule_refs(imol, cntdst, stmax, store_spline, store_grid)
    use engmain, only: ms1max, ms2max, ms3max, splodr, specatm, sitepos, invcl
    use spline, only: spline_value
    implicit none

    integer, intent(in) :: imol, cntdst, stmax
    real, intent(out) :: store_spline(0:splodr-1, 3, 1:stmax)
    integer, intent(out) :: store_grid(3, 1:stmax)

    integer :: sid, ati, ati_ext, rcimax, m, k, rci, spi
    real :: xst(3), inm(3)
    real :: factor, rtp2

    do sid = 1, stmax
       ati = specatm(sid, imol)
       ati_ext = ati + (cntdst - 1) * stmax
       xst(:) = sitepos(:, ati_ext)
       do k = 1, 3
          factor = dot_product(invcl(k,:), xst(:))
          factor = factor - floor(factor)
          inm(k) = factor
       end do
       do m = 1, 3
          if (m == 1) rcimax = ms1max
          if (m == 2) rcimax = ms2max
          if (m == 3) rcimax = ms3max
          factor = inm(m) * real(rcimax)
          rci = int(factor)
          do spi = 0, splodr - 1
             rtp2 = factor - real(rci - spi)
             store_spline(spi, m, sid) = spline_value(rtp2)
          end do
          store_grid(m, sid) = rci
       end do
    end do
  end subroutine calc_spline_molecule_refs

  function recpcal_self_energy() result(pairep)
    implicit none
    real :: pairep

    pairep = solute_self_energy
  end function recpcal_self_energy

  function recpcal_self_energy_refs(cnt) result(pairep)
    implicit none
    integer, intent(in) :: cnt
    real :: pairep

    pairep = solute_self_energy_refs(cnt)
  end function recpcal_self_energy_refs

  subroutine recpcal_energy_soln(tagslt, tagpt, slvmax, uvengy, cnt)
    use engmain, only: ms1max, ms2max, ms3max, splodr, numsite, sluvid, charge, mol_begin_index
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: tagslt, tagpt(:), slvmax, cnt
    real, intent(inout) :: uvengy(:, :)

    real :: pairep
    integer :: cg1, cg2, cg3, i, k
    integer :: rc1, rc2, rc3, ptrnk, sid, ati, svi, stmax
    real :: fac1, fac2, fac3, chr
    integer :: grid1
    complex :: rcpt

    if (sluvid(tagslt) == 0) stop  ! call halt_with_error('rcp_fst')

    !$acc parallel loop present(uvengy, mol_begin_index, tagpt, charge, numsite, sluvid, slvtag, splslv, grdslv, cnvslt)
    do k = 1, slvmax
       i = tagpt(k)
       if (i == tagslt) cycle

       pairep = 0.0
       svi = slvtag(i)
       if (svi <= 0) stop  ! call halt_with_error('rcp_cns')
       stmax = numsite(i)
       do sid = 1, stmax
          ptrnk = svi + sid - 1
          ati = mol_begin_index(i) + (sid - 1) ! = specatm(sid, i)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             fac1 = chr * splslv(cg3, 3, ptrnk)
             rc3 = modulo(grdslv(3, ptrnk) - cg3, ms3max)
             do cg2 = 0, splodr - 1
                fac2 = fac1 * splslv(cg2, 2, ptrnk)
                rc2 = modulo(grdslv(2, ptrnk) - cg2, ms2max)
                grid1 = grdslv(1, ptrnk)
                if (grid1 >= splodr-1 .and. grid1 < ms1max) then
                   !$acc loop seq
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = grid1 - cg1
                      pairep = pairep + fac3 * real(cnvslt(rc1, rc2, rc3))
                   enddo
                else
                   !$acc loop seq
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = mod(grid1 + ms1max - cg1, ms1max) ! speedhack
                      pairep = pairep + fac3 * real(cnvslt(rc1, rc2, rc3))
                   end do
                endif
             end do
          end do
       end do
       uvengy(k, cnt) = uvengy(k, cnt) + pairep
    end do
    !$acc end parallel
  end subroutine recpcal_energy_soln

  subroutine recpcal_energy_refs(tagslt, maxdst, slvmax, uvengy)
    use engmain, only: ms1max, ms2max, ms3max, splodr, numsite, sluvid, charge, mol_begin_index
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: tagslt, maxdst, slvmax
    real, intent(inout) :: uvengy(:, :)

    real :: pairep
    integer :: cg1, cg2, cg3, i, k, cnt
    integer :: rc1, rc2, rc3, ptrnk, sid, ati, svi, stmax
    real :: fac1, fac2, fac3, chr
    integer :: grid1
    complex :: rcpt

    if (sluvid(tagslt) == 0) stop  ! call halt_with_error('rcp_fst')

    !$acc parallel loop collapse(2) present(uvengy, mol_begin_index, charge, numsite, sluvid, slvtag, splslv, grdslv, cnvslt)
    do cnt = 1, maxdst
    do i = 1, slvmax

       pairep = 0.0
       svi = slvtag(i)
       if (svi <= 0) stop  ! call halt_with_error('rcp_cns')
       stmax = numsite(i)
       do sid = 1, stmax
          ptrnk = svi + sid - 1
          ati = mol_begin_index(i) + (sid - 1) ! = specatm(sid, i)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             fac1 = chr * splslv(cg3, 3, ptrnk)
             rc3 = modulo(grdslv(3, ptrnk) - cg3, ms3max)
             do cg2 = 0, splodr - 1
                fac2 = fac1 * splslv(cg2, 2, ptrnk)
                rc2 = modulo(grdslv(2, ptrnk) - cg2, ms2max)
                grid1 = grdslv(1, ptrnk)
                if (grid1 >= splodr-1 .and. grid1 < ms1max) then
                   !$acc loop seq
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = grid1 - cg1
                      pairep = pairep &
                           + fac3 * real(cnvslt_r(rc1, rc2, rc3, cnt))
                   enddo
                else
                   !$acc loop seq
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = mod(grid1 + ms1max - cg1, ms1max) ! speedhack
                      pairep = pairep &
                           + fac3 * real(cnvslt_r(rc1, rc2, rc3, cnt))
                   end do
                endif
             end do
          end do
       end do
       uvengy(i, cnt) = uvengy(i, cnt) + pairep
    end do
    end do
    !$acc end parallel
  end subroutine recpcal_energy_refs

end module reciprocal
