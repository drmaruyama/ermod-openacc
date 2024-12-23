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

module realcal
  implicit none
  
  ! "straight" coordinate system
  real(kind=8), allocatable :: sitepos_normal(:, :)
  real(kind=8) :: cell_normal(3, 3), invcell_normal(3), cell_len_normal(3)
  logical :: is_cuboid
  real, parameter :: check_rotate = 1e-8, cuboid_thres = 1e-8
  real, parameter :: cutoff_thres = 1e-4

contains
  subroutine realcal_prepare
    use engmain, only: numatm, sitepos, boxshp, SYS_PERIODIC
    implicit none
    logical, save :: initialized = .false.

    if (.not. initialized) then
       allocate(sitepos_normal(3, numatm))
       !$acc enter data create(sitepos_normal)
       initialized = .true.
    end if
    sitepos_normal(:, :) = sitepos(:, :)

    ! "Straighten" box, and normalize coordinate system
    if(boxshp == SYS_PERIODIC) call normalize_periodic
    !$acc update device(sitepos_normal)
  end subroutine realcal_prepare

  ! Calculate i-j interaction energy on GPU
  subroutine realcal_acc(tagslt, tagpt, slvmax, uvengy)
    use engmain, only:  boxshp, numsite, &
         elecut, lwljcut, upljcut, cltype, screen, charge, mol_begin_index, &
         ljswitch, ljtype, ljtype_max, ljene_mat, ljlensq_mat, &
         SYS_NONPERIODIC, SYS_PERIODIC, EL_COULOMB, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX
    implicit none
    integer, intent(in) :: tagslt, tagpt(:), slvmax
    real, intent(inout) :: uvengy(0:slvmax)

    integer :: i, k, is, js, ismax, jsmax, ati, atj
    real :: reelcut, pairep, rst, dis2, invr2, invr3, invr6
    real :: eplj, epcl, xst(3), half_cell(3)
    real :: lwljcut2, upljcut2
    real, save :: lwljcut3, upljcut3, lwljcut6, upljcut6
    real :: ljeps, ljsgm2, ljsgm3, ljsgm6, vdwa, vdwb, swth, swfac
    real, save :: repA, repB, repC, attA, attB, attC
    integer :: ljtype_i, ljtype_j
    logical, save :: initialized = .false.
    real, parameter :: infty = huge(infty)      ! essentially equal to infinity
    !
    if(boxshp == SYS_NONPERIODIC) reelcut = infty
    if(boxshp == SYS_PERIODIC) then
       reelcut = elecut
       half_cell(:) = 0.5 * cell_len_normal(:)
    else
       ! suppress warnings
       half_cell(:) = 0.0
    endif

    if(.not. initialized) then
       if(ljswitch == LJSWT_FRC_CHM) then       ! force switch (CHARMM type)
          lwljcut3 = lwljcut ** 3
          upljcut3 = upljcut ** 3
          lwljcut6 = lwljcut3 * lwljcut3
          upljcut6 = upljcut3 * upljcut3
       endif
       if(ljswitch == LJSWT_FRC_GMX) then       ! force switch (GROMACS type)
          call calc_gmx_switching_force_params(12, lwljcut, upljcut, &
               repA, repB, repC)
          call calc_gmx_switching_force_params(6,  lwljcut, upljcut, &
               attA, attB, attC)
       endif
       initialized = .true.
    end if

    ! calculated only when PME or PPPM, non-self interaction
    ismax = numsite(tagslt)
    !$acc data pcreate(xst) copyin(half_cell, cell_normal, invcell_normal)
    !$acc parallel loop collapse(2) gang vector present(uvengy, mol_begin_index, tagpt, sitepos_normal, ljlensq_mat, ljene_mat, charge, ljtype, numsite)
    do k = 1, slvmax
       do is = 1, ismax
          i = tagpt(k)
          if(i == tagslt) cycle

          pairep = 0.0
          do js = 1, numsite(i)
!             ati = specatm(is, tagslt)
!             atj = specatm(js, i)
             ati = mol_begin_index(tagslt) + (is - 1)
             atj = mol_begin_index(i) + (js - 1)
             ljtype_i = ljtype(ati)
             ljtype_j = ljtype(atj)
             xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
             if(boxshp == SYS_PERIODIC) then    ! when the system is periodic
                if(is_cuboid) then
                   xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
                else
                   ! Note some ops can be skipped because cell_normal
                   ! is upper triangular
                   xst(:) = xst(:) - cell_normal(:, 3) &
                        * anint(xst(3) * invcell_normal(3))
                   xst(:) = xst(:) - cell_normal(:, 2) &
                        * anint(xst(2) * invcell_normal(2))
                   xst(:) = xst(:) - cell_normal(:, 1) &
                        * anint(xst(1) * invcell_normal(1))
                end if
             endif
             dis2 = sum(xst(1:3) ** 2)
             rst = sqrt(dis2)
             if(rst > upljcut) then
                eplj = 0.0
             else
                ljeps = ljene_mat(ljtype_i, ljtype_j)
                ljsgm2 = ljlensq_mat(ljtype_i, ljtype_j)

                invr2 = ljsgm2 / dis2
                invr6 = invr2 * invr2 * invr2
                select case(ljswitch)
                case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
                   eplj = 4.0 * ljeps * invr6 * (invr6 - 1.0)
                   if(rst > lwljcut) then
                      select case(ljswitch)
                      case(LJSWT_POT_CHM)                  ! CHARMM type
                         lwljcut2 = lwljcut ** 2
                         upljcut2 = upljcut ** 2
                         swth = (2.0 * dis2 + upljcut2 - 3.0 * lwljcut2)      &
                              * ((dis2 - upljcut2) ** 2)                      &
                              / ((upljcut2 - lwljcut2) ** 3)
                      case(LJSWT_POT_GMX)                  ! GROMACS type
                         swfac = (rst - lwljcut) / (upljcut - lwljcut)
                         swth = 1.0 - 10.0 * (swfac ** 3)                     &
                              + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5)
                      case default
                         stop "Unknown ljswitch"
                      end select
                      eplj = swth * eplj
                   endif
                case(LJSWT_FRC_CHM)               ! force switch (CHARMM type)
                   ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                   if(rst <= lwljcut) then
                      vdwa = invr6 * invr6 &
                           - ljsgm6 *ljsgm6 / (lwljcut6 * upljcut6)
                      vdwb = invr6 - ljsgm6 / (lwljcut3 * upljcut3)
                   else
                      invr3 = sqrt(invr6)
                      ljsgm3 = sqrt(ljsgm6)
                      vdwa = upljcut6 / (upljcut6 - lwljcut6)                 &
                           * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
                      vdwb = upljcut3 / (upljcut3 - lwljcut3)                 &
                           * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
                   endif
                   eplj = 4.0 * ljeps * (vdwa - vdwb)
                case(LJSWT_FRC_GMX)               ! force switch (GROMACS type)
                   ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                   if(rst <= lwljcut) then
                      vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * repC
                      vdwb = invr6 - ljsgm6 * attC
                   else
                      swfac = rst - lwljcut
                      vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 *                &
                           (repA * (swfac ** 3) + repB * (swfac ** 4) + repC)
                      vdwb = invr6 - ljsgm6 *                                 &
                           (attA * (swfac ** 3) + attB * (swfac ** 4) + attC)
                   endif
                   eplj = 4.0 * ljeps * (vdwa - vdwb)
                case default
                   stop "Unknown ljswitch"
                end select
             endif
             if(rst >= reelcut) then
                epcl = 0.0
             else
                epcl = charge(ati) * charge(atj) &
                     * (1.0 - erf(screen * rst)) / rst
             endif
             pairep = pairep + eplj + epcl
          end do
          !$acc atomic update
          uvengy(k) = uvengy(k) + pairep
       end do
    end do
    !$acc end parallel
    !$acc end data
  end subroutine realcal_acc

  ! Calculate i-j interaction energy in the bare 1/r form
  subroutine realcal_bare(tagslt, tagpt, slvmax, uvengy)
    use engmain, only:  boxshp, numsite, &
         elecut, lwljcut, upljcut, cltype, screen, charge, mol_begin_index, &
         ljswitch, ljtype, ljtype_max, ljene_mat, ljlensq_mat, &
         SYS_NONPERIODIC, SYS_PERIODIC, EL_COULOMB, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX
    implicit none
    integer, intent(in) :: tagslt, tagpt(:), slvmax
    real, intent(inout) :: uvengy(0:slvmax)

    integer :: i, k, is, js, ismax, jsmax, ati, atj
    real :: reelcut, pairep, rst, dis2, invr2, invr3, invr6
    real :: eplj, epcl, xst(3), half_cell(3)
    real :: lwljcut2, upljcut2, lwljcut3, upljcut3, lwljcut6, upljcut6
    real :: ljeps, ljsgm2, ljsgm3, ljsgm6, vdwa, vdwb, swth, swfac
    real :: repA, repB, repC, attA, attB, attC
    integer :: ljtype_i, ljtype_j
    real, parameter :: infty = huge(infty)      ! essentially equal to infinity
    !
    if(cltype /= EL_COULOMB) stop "cannot happen: realcal_bare is called only when cltype is 'bare coulomb'."

    if(boxshp == SYS_NONPERIODIC) reelcut=infty
    if(boxshp == SYS_PERIODIC) then
       reelcut = elecut
       half_cell(:) = 0.5 * cell_len_normal(:)
    else
       ! suppress warnings
       half_cell(:) = 0.0
    endif

    if(ljswitch == LJSWT_FRC_CHM) then       ! force switch (CHARMM type)
       lwljcut3 = lwljcut ** 3
       upljcut3 = upljcut ** 3
       lwljcut6 = lwljcut3 * lwljcut3
       upljcut6 = upljcut3 * upljcut3
    endif
    if(ljswitch == LJSWT_FRC_GMX) then       ! force switch (GROMACS type)
       call calc_gmx_switching_force_params(12, lwljcut, upljcut, repA, repB, repC)
       call calc_gmx_switching_force_params(6,  lwljcut, upljcut, attA, attB, attC)
    endif

    ! Bare coulomb solute-solvent interaction
    ismax = numsite(tagslt)
    !$acc data pcreate(xst)
    !$acc parallel loop collapse(2) gang vector present(uvengy, mol_begin_index, tagpt, sitepos_normal, charge, numsite)
    do k = 1, slvmax
       do is = 1, ismax
          i = tagpt(k)
          if(i == tagslt) cycle

          pairep = 0.0
          do js = 1, numsite(i)
!             ati = specatm(is, tagslt)
!             atj = specatm(js, i)
             ati = mol_begin_index(tagslt) + (is - 1)
             atj = mol_begin_index(i) + (js - 1)
             ljtype_i = ljtype(ati)
             ljtype_j = ljtype(atj)
             xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
             if(boxshp == SYS_PERIODIC) then    ! when the system is periodic
                if(is_cuboid) then
                   xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
                else
                   ! Note some ops can be skipped because cell_normal is upper triangular
                   xst(:) = xst(:) - cell_normal(:, 3) &
                        * anint(xst(3) * invcell_normal(3))
                   xst(:) = xst(:) - cell_normal(:, 2) &
                        * anint(xst(2) * invcell_normal(2))
                   xst(:) = xst(:) - cell_normal(:, 1) &
                        * anint(xst(1) * invcell_normal(1))
                end if
             endif
             dis2 = sum(xst(1:3) ** 2)
             rst = sqrt(dis2)
             if(rst > upljcut) then
                eplj = 0.0
             else
                ljeps = ljene_mat(ljtype_i, ljtype_j)
                ljsgm2 = ljlensq_mat(ljtype_i, ljtype_j)

                invr2 = ljsgm2 / dis2
                invr6 = invr2 * invr2 * invr2
                select case(ljswitch)
                case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
                   eplj = 4.0 * ljeps * invr6 * (invr6 - 1.0)
                   if(rst > lwljcut) then
                      select case(ljswitch)
                      case(LJSWT_POT_CHM)                  ! CHARMM type
                         lwljcut2 = lwljcut ** 2
                         upljcut2 = upljcut ** 2
                         swth = (2.0 * dis2 + upljcut2 - 3.0 * lwljcut2)      &
                              * ((dis2 - upljcut2) ** 2)                      &
                              / ((upljcut2 - lwljcut2) ** 3)
                      case(LJSWT_POT_GMX)                  ! GROMACS type
                         swfac = (rst - lwljcut) / (upljcut - lwljcut)
                         swth = 1.0 - 10.0 * (swfac ** 3)                     &
                              + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5)
                      case default
                         stop "Unknown ljswitch"
                      end select
                      eplj = swth * eplj
                   endif
                case(LJSWT_FRC_CHM)               ! force switch (CHARMM type)
                   ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                   if(rst <= lwljcut) then
                      vdwa = invr6 * invr6 &
                           - ljsgm6 *ljsgm6 / (lwljcut6 * upljcut6)
                      vdwb = invr6 - ljsgm6 / (lwljcut3 * upljcut3)
                   else
                      invr3 = sqrt(invr6)
                      ljsgm3 = sqrt(ljsgm6)
                      vdwa = upljcut6 / (upljcut6 - lwljcut6)                 &
                           * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
                      vdwb = upljcut3 / (upljcut3 - lwljcut3)                 &
                           * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
                   endif
                   eplj = 4.0 * ljeps * (vdwa - vdwb)
                case(LJSWT_FRC_GMX)               ! force switch (GROMACS type)
                   ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                   if(rst <= lwljcut) then
                      vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * repC
                      vdwb = invr6 - ljsgm6 * attC
                   else
                      swfac = rst - lwljcut
                      vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 *                &
                           (repA * (swfac ** 3) + repB * (swfac ** 4) + repC)
                      vdwb = invr6 - ljsgm6 *                                 &
                           (attA * (swfac ** 3) + attB * (swfac ** 4) + attC)
                   endif
                   eplj = 4.0 * ljeps * (vdwa - vdwb)
                case default
                   stop "Unknown ljswitch"
                end select
             endif
             if(rst >= reelcut) then
                epcl = 0.0
             else
                epcl = charge(ati) * charge(atj) / rst
             endif
             pairep = pairep + eplj + epcl
          end do
          !$acc atomic update
          uvengy(k) = uvengy(k) + pairep
       end do
    end do
    !$acc end parallel
    !$acc end data
  end subroutine realcal_bare

  ! self-energy part, no LJ calculation performed
  subroutine realcal_self(i, pairep)
    use engmain, only: numsite, screen, cltype, charge, specatm, &
                       EL_COULOMB, PI
    implicit none
    integer, intent(in) :: i
    real, intent(inout) :: pairep
    integer :: is, js, ismax, ati, atj
    real :: rst, dis2, epcl, xst(3), half_cell(3)

    pairep = 0.0
    if(cltype == EL_COULOMB) return

    half_cell(:) = 0.5 * cell_len_normal(:) 

    ismax=numsite(i)

    do is = 1, ismax
       ati = specatm(is, i)

       ! Atom residual
       ! self (the same ati arguments for two charge variables below)
       epcl = - charge(ati) * charge(ati) * screen / sqrt(PI)
       pairep = pairep + epcl

       do js = is + 1, ismax
          atj = specatm(js, i)
 
          xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
          if(is_cuboid) then
             xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
          else
             ! Note some ops can be skipped because cell_normal is upper triangular
             xst(:) = xst(:) - cell_normal(:, 3) &
                  * anint(xst(3) * invcell_normal(3))
             xst(:) = xst(:) - cell_normal(:, 2) &
                  * anint(xst(2) * invcell_normal(2))
             xst(:) = xst(:) - cell_normal(:, 1) &
                  * anint(xst(1) * invcell_normal(1))
          end if

          dis2 = sum(xst(1:3) ** 2)
          rst = sqrt(dis2)

          ! distinct (different ati and atj arguments for two charge variables)
          epcl = - charge(ati) * charge(atj) * erf(screen * rst) / rst
          pairep = pairep + epcl
       enddo
    enddo

  end subroutine realcal_self

  ! Rotate box and coordinate so that the cell(:,:) is upper triangular:
  ! cell(2, 1) = cell(3, 1) = cell(3, 2) = 0
  ! (1st axis to be aligned to x-axis, 2nd axis to be within xy-plane)
  subroutine make_cell_uppertriangular
    use engmain, only: cell, sitepos
    use mpiproc, only: warning
    implicit none

    real :: q(3,3), r(3, 3), temp(3, 3), axis_len
    integer :: i, j

    if (cell(2, 1) == 0.0 .and. cell(3, 1) == 0.0 .and. cell(3, 2) == 0.0) then
       ! already upper triangular
       ! As long as we use modern MD software the program shall go this path.
       cell_normal(:, :) = cell(:, :)
       sitepos_normal(:, :) = sitepos(:, :)
       return
    end if

    call warning("cell")

    ! modified Gram-Schmidt
    q(:, :) = 0.0
    r(:, :) = 0.0
    temp(:, :) = cell(:, :)
    do i = 1, 3
       axis_len = sqrt(dot_product(temp(:, 1), temp(:, 1)))
       q(:, i) = temp(:, i) / axis_len
       r(i, i) = axis_len
       do j = i + 1, 3
          ! remove inner products
          r(i, j) = dot_product(q(:, i), temp(:, j))
          temp(:, j) = temp(:, j) - r(i, j) * q(:, i)
       end do
    end do

    ! Now cell = QR
    ! cell_normal = R = Q^T Q R
    cell_normal(:, :) = r(:, :)

    ! sitepos_normal = Q^T X
    do i = 1, size(sitepos, 2)
       sitepos_normal(:, i) = matmul(transpose(q), sitepos(:, i))
    end do
  end subroutine make_cell_uppertriangular

  ! Cell vectors are normalized by either inverting the cell vector or by adding integer times other vectors
  subroutine normalize_cell_vector
    use engmain, only: elecut, upljcut
    use mpiproc, only: warning
    integer :: i
    real :: cutoff

    do i = 1, 3
       if (cell_normal(1, 1) < 0.0) then
          cell_normal(:, i) = -cell_normal(:, i)
       end if
    end do

    ! check cell size restrictions
    if (abs(cell_normal(1, 2)) > 0.5 * cell_normal(1, 1) + cutoff_thres) then
       call warning("cel2")
       cell_normal(:, 2) = cell_normal(:, 2) - &
            cell_normal(:, 1) * anint(cell_normal(1, 2) * invcell_normal(1))
    end if
    if (abs(cell_normal(1, 3)) > 0.5 * cell_normal(1, 1) + cutoff_thres) then
       call warning("cel2")
       cell_normal(:, 3) = cell_normal(:, 3) - &
            cell_normal(:, 1) * anint(cell_normal(1, 3) * invcell_normal(1))
    end if
    if (abs(cell_normal(2, 3)) > 0.5 * cell_normal(2, 2) + cutoff_thres) then
       call warning("cel2")
       cell_normal(:, 3) = cell_normal(:, 3) - &
            cell_normal(:, 2) * anint(cell_normal(2, 3) * invcell_normal(2))
    end if

    cutoff = min(elecut, upljcut)
    ! check cutoff restrictions
    do i = 1, 3
       if (cutoff > cell_normal(i, i) * 0.5) then
          stop "One of axis in periodic cell is too small. This is either the box is too small compared to the cell, " // &
               "or the periodic cell is too skewed."
       endif
    end do
  end subroutine normalize_cell_vector

  subroutine normalize_periodic
    use engmain, only: cell, sitepos
    implicit none
    integer :: n, i

    call make_cell_uppertriangular

    call normalize_cell_vector

    n = size(sitepos, 2)

    do i = 1, 3
       cell_len_normal(i) = abs(cell_normal(i, i))
    end do
    invcell_normal(:) = 1 / cell_len_normal(:)

    if(abs(cell(1, 2)) > cuboid_thres .or. &
         abs(cell(1, 3)) > cuboid_thres .or. &
         abs(cell(2, 3)) > cuboid_thres ) then
       is_cuboid = .false.
    else
       is_cuboid = .true.
    end if

    ! move all particles inside the cuboid spanned by (0 .. cell(1, 1)), (0 .. cell(2,2)), (0 .. cell(3,3)).
    ! Shift with the order Z -> Y -> X
    do i = 1, n
       ! shift Z
       sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(:, 3) * floor(invcell_normal(3) * sitepos_normal(3, i))
       ! shift Y
       sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(:, 2) * floor(invcell_normal(2) * sitepos_normal(2, i))
       ! shift X
       sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(:, 1) * floor(invcell_normal(1) * sitepos_normal(1, i))

       if(sitepos_normal(1, i) < 0 .or. sitepos_normal(1, i) > cell_len_normal(1) .or.&
          sitepos_normal(2, i) < 0 .or. sitepos_normal(2, i) > cell_len_normal(2) .or.&
          sitepos_normal(3, i) < 0 .or. sitepos_normal(3, i) > cell_len_normal(3)) then
          print *, sitepos_normal(:, i), cell_len_normal
          stop "INVALID sitepos"
       endif
    end do
  end subroutine normalize_periodic

  ! get the coefficients for gromacs force switching
  subroutine calc_gmx_switching_force_params(pow, lwljcut, upljcut, coeffA, coeffB, coeffC)
    implicit none
    integer, intent(in) :: pow
    real, intent(in) :: lwljcut, upljcut
    real, intent(out) :: coeffA, coeffB, coeffC
    real :: dfljcut

    dfljcut = upljcut - lwljcut
    coeffA = - real(pow) * (real(pow + 4) * upljcut                   &
                          - real(pow + 1) * lwljcut)                  &
           / ((upljcut ** (pow + 2)) * (dfljcut ** 2)) / 3.0
    coeffB =   real(pow) * (real(pow + 3) * upljcut                   &
                          - real(pow + 1) * lwljcut)                  &
           / ((upljcut ** (pow + 2)) * (dfljcut ** 3)) / 4.0
    coeffC = 1.0 / (upljcut ** pow) - coeffA * (dfljcut ** 3)         &
                                    - coeffB * (dfljcut ** 4)
  end subroutine calc_gmx_switching_force_params

end module realcal
