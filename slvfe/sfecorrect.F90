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

module uvcorrect
   implicit none
   integer, dimension(:), allocatable :: ptsite

   ! In the following, all the variables are copied from engmain.F90
   !                   and their names are all the same as in engmain.F90
   character(len=*), parameter :: ene_confname = 'parameters_er'

   character(len=*), parameter :: solute_file = 'SltInfo'  ! solute species
   character(len=*), parameter :: solvent_file = 'MolPrm'  ! solvent species
   character(len=*), parameter :: ljtable_file = 'LJTable' ! table for LJ
   integer, parameter :: mol_io = 79        ! IO for SltInfo and MolPrmX
   integer, parameter :: ljtable_io = 70    ! IO for LJ table

   integer, parameter :: LJFMT_EPS_cal_SGM_nm = 0, LJFMT_EPS_Rminh = 1, &
      LJFMT_EPS_J_SGM_A = 2, LJFMT_A_C = 3, &
      LJFMT_C12_C6 = 4, LJFMT_TABLE = 5
   integer, parameter :: LJSWT_POT_CHM = 0, LJSWT_POT_GMX = 1, &
      LJSWT_FRC_CHM = 2, LJSWT_FRC_GMX = 3
   integer, parameter :: LJCMB_ARITH = 0, LJCMB_GEOM = 1

   integer                              :: ljtype_max
   integer, dimension(:,:), allocatable :: ljtype
   real, dimension(:,:),    allocatable :: ljlensq_mat, ljene_mat

   ! Declarations of the variables listed in namelist
   ! Only ljformat, ljswitch, cmbrule, lwljcut, and upljcut used in this program
   integer :: ljformat, ljswitch, cmbrule
   real :: lwljcut, upljcut

   ! The other variables are declared below since they appear in the namelist
   !     and the namelist is copied from engmain.F90
   integer :: iseed, skpcnf, corrcal, selfcal
   integer :: slttype, wgtslf, wgtsys, wgtins, boxshp, estype
   integer, parameter :: max_for_spec = 100    ! maximum number of species
   integer :: sltspec, hostspec(max_for_spec), refspec(max_for_spec)
   real :: lwreg, upreg, lwstr, upstr
   integer :: insposition, insorient, insstructure
   integer :: sltpick, refpick, inscnd, inscfg
   ! ljformat and ljswitch are declared above
   real :: inptemp, temp
   integer :: engdiv, maxins
   ! cmbrule, lwljcut, upljcut are declared above
   integer :: intprm, cltype, splodr
   real :: elecut, screen, ewtoler
   integer :: ew1max, ew2max, ew3max, ms1max, ms2max, ms3max
   integer :: ermax_limit
   real :: block_threshold
   logical :: force_calculation

   namelist /ene_param/ iseed, &
      skpcnf, corrcal, selfcal, &
      slttype, wgtslf, wgtsys, wgtins, boxshp, estype, &
      sltspec, hostspec, refspec, lwreg, upreg, lwstr, upstr, &
      insposition, insorient, insstructure, &
      sltpick, refpick, inscnd, inscfg, &                  ! deprecated
      ljformat, ljswitch, &
      inptemp, temp, &
      engdiv, maxins, &
      intprm, elecut, lwljcut, upljcut, &
      cmbrule, cltype, screen, ewtoler, splodr, &
      ew1max, ew2max, ew3max, ms1max, ms2max, ms3max, &
      ermax_limit, block_threshold, force_calculation

contains
   subroutine ljcorrect(cntrun)
      use sysvars, only: uvread, clcond, numslv, aveuv, blockuv
      implicit none
      real, dimension(:), allocatable, save :: ljcorr
      integer, intent(in) :: cntrun
      logical, save :: first_time = .true.
      integer :: pti
      if(first_time) then
         call set_keyparam
         write(6, '(A, f7.1, A)') '  Be sure that ' // &
            'the solvent distribution is homogeneous ' // &
            '(radial distribution function is essentially unity) ' // &
            'when the solvent molecule is separated beyond distance of', &
            lwljcut, ' Angstrom in any direction ' // &
            'from any atom within the solute molecule'
         write(6, *)
         call get_ljtable
         allocate( ljcorr(numslv) )
         do pti = 1, numslv
            call calc_ljlrc(pti, ljcorr(pti))
         end do
         write(6, '(A, 9999f12.4)') '  LJ long-range correction    =   ', &
            ljcorr(1:numslv)
         write(6, *)
         first_time = .false.
      endif
      aveuv(1:numslv) = aveuv(1:numslv) + ljcorr(1:numslv)
      if((uvread /= 'not') .and. (clcond == 'merge')) then
         blockuv(1:numslv, cntrun) = blockuv(1:numslv, cntrun) + ljcorr(1:numslv)
         blockuv(0, cntrun) = blockuv(0, cntrun) + sum( ljcorr(1:numslv) )
      endif
      return
   end subroutine ljcorrect

   subroutine set_keyparam
      use sysvars, only: avevolume, refsdirec
      implicit none
      character(len=80) :: keyfile
      integer, parameter :: iounit = 555
      real, parameter :: volm_min = 1.40e4
      integer :: stat

      ! Initialization of the variables
      !    that are in the ene_param namelist but are not used in this program
      iseed = 0 ; skpcnf = 1 ; corrcal = 0 ; selfcal = 0
      slttype = 1 ; wgtslf = 0 ; wgtsys = 0 ; wgtins = 0
      boxshp = 1 ; estype = 1
      sltspec = 1 ; hostspec(:) = 0 ; refspec(:) = 0
      lwreg = 0 ; upreg = lwreg ; lwstr = 0 ; upstr = lwstr
      insposition = 0 ; insorient = 0 ; insstructure = 0
      sltpick = sltspec ; refpick = 0 ; inscnd = 0 ; inscfg = 0
      inptemp = 300 ; temp = 0.592 ; engdiv = 1 ; maxins = 0
      intprm = 1 ; elecut = 12.0
      cltype = 2 ; screen = 0 ; ewtoler = 1.0e-5 ; splodr = 6
      ew1max = 0 ; ew2max = 0 ; ew3max = 0 ; ms1max = 0 ; ms2max = 0 ; ms3max = 0
      ermax_limit = 15000 ; block_threshold = 4.0 ; force_calculation = .false.
      ! The above initialization is done
      !    to prevent possible errors due to memory allocation

      ljformat = LJFMT_EPS_Rminh                    ! default setting
      ljswitch = LJSWT_POT_CHM                      ! default setting
      cmbrule = LJCMB_ARITH                         ! default setting
      upljcut = 12.0                                ! default setting
      lwljcut = upljcut - 2.0                       ! default setting

      keyfile = trim(refsdirec)//'/'//ene_confname
      open(unit = iounit, file = keyfile, action = 'read', status = 'old', iostat = stat)
      if(stat == 0) then
         read(iounit, nml = ene_param)
      else
         stop "The parameters_er file is not found in the refs directory"
      endif
      close(iounit)

      if(avevolume <= 0.0) then   ! no input of avevolume from parameters_fe
         write(6, '(A)') "  What is the average volume of reference solvent? (in Angstrom^3)"
         read(5, *) avevolume
      endif
      if(avevolume < volm_min) then
         write(6, '(A)') "  Warning: your input volume seems too small"
         write(6, '(A, f8.1)') "           This warning appears when your input is less than ",volm_min
         write(6, '(A)') "  Re-type the volume in Angstrom^3 (NOT in nm^3)"
         read(5, *) avevolume
      endif
   end subroutine set_keyparam

   ! setup of the Lennard-Jones parameters
   ! the program is taken from the corresponding part in setconf.F90
   subroutine get_ljtable
      use sysvars, only: numslv, refsdirec
      implicit none
      real, parameter :: sgmcnv = 1.7817974362806784e0 ! from Rmin/2 to sigma, 2.0**(5.0/6.0)
      real, parameter :: lencnv = 1.0e1                ! from nm to Angstrom
      real, parameter :: engcnv = 1.0e0/4.184e0        ! from kJ/mol to kcal/mol
      integer :: pti, sid, stmax, maxsite, i, m
      real :: factor, xst(3), mass
      integer, allocatable :: ljtype_temp(:)
      real, dimension(:), allocatable :: ljlen_temp, ljene_temp
      real, dimension(:), allocatable :: ljlen_temp_table, ljene_temp_table
      integer :: ljtype_found
      logical :: lj_is_new
      character(len=12) :: atmtype
      character(len=8) :: atmname
      character(len=80) :: molfile
      character(len=9) :: numbers = '123456789'
      character(len=120) :: linebuf
      integer :: ierr

      allocate( ptsite(0:numslv) )
      do pti = 0, numslv
         if(pti == 0) then
            molfile = solute_file                          ! solute
         else
            molfile = solvent_file // numbers(pti:pti)     ! solvent
         endif
         molfile = trim(refsdirec) // '/' // molfile
         open(unit = mol_io, file = molfile, status='old')
         stmax = 0
         do
            read(mol_io, *, end = 99) m
            stmax = stmax + 1
         end do
99       close(mol_io)
         ptsite(pti) = stmax
      end do

      ! large enough LJ table size
      allocate( ljlen_temp_table(1:sum(ptsite(:))), &
         ljene_temp_table(1:sum(ptsite(:))) )
      ! temporary set of LJ
      maxsite = maxval( ptsite(0:numslv) )
      allocate( ljtype_temp(maxsite), ljlen_temp(maxsite), ljene_temp(maxsite) )

      allocate( ljtype(maxsite, 0:numslv) )
      ljtype(:, :) = 0

      do pti = 0, numslv
         if(pti == 0) then
            molfile = solute_file                          ! solute
         else
            molfile = solvent_file // numbers(pti:pti)     ! solvent
         endif
         molfile = trim(refsdirec) // '/' // molfile
         stmax = ptsite(pti)
         open(unit = mol_io, file = molfile, status = 'old')
         do sid = 1, stmax
            read(mol_io, '(A)') linebuf ! read entire line into buf
            ! new format, covering both rigid and flexible
            read(linebuf, *, iostat = ierr) m, mass, atmtype, atmname, xst(1:3)
            if(ierr/=0) then
               ! old format, covering both rigid and flexible
               read(linebuf, *) m, atmtype, xst(1:3)
            end if

            if(ljformat == LJFMT_EPS_Rminh) xst(3) = sgmcnv * xst(3)
            if((ljformat == LJFMT_A_C) .or. (ljformat == LJFMT_C12_C6)) then
               if(xst(3) /= 0.0) then
                  factor = (xst(2) / xst(3)) ** (1.0 / 6.0)
                  xst(2) = xst(3) / (4.0 * (factor ** 6))
                  xst(3) = factor
               else
                  xst(2) = 0.0
               endif
            endif
            if((ljformat == LJFMT_EPS_J_SGM_A) .or. (ljformat == LJFMT_C12_C6)) then
               xst(2) = engcnv * xst(2)
               xst(3) = lencnv * xst(3)
            endif
            ljene_temp(sid) = xst(2)
            ljlen_temp(sid) = xst(3)
         end do
         close(mol_io)

         if(ljformat == LJFMT_TABLE) then
            ljtype_temp(1:stmax) = ljene_temp(1:stmax)
         else
            do sid = 1, stmax
               lj_is_new = .true.
               do i = 1, ljtype_max
                  ! linear search LJ table
                  if((ljlen_temp_table(i) == ljlen_temp(sid)) .and. &
                     (ljene_temp_table(i) == ljene_temp(sid))) then
                     ljtype_found = i
                     lj_is_new = .false.
                     exit
                  endif
               end do
               if(lj_is_new) then
                  ! new LJ type
                  ljtype_max = ljtype_max + 1
                  ljlen_temp_table(ljtype_max) = ljlen_temp(sid)
                  ljene_temp_table(ljtype_max) = ljene_temp(sid)
                  ljtype_found = ljtype_max
               endif
               ljtype_temp(sid) = ljtype_found
            end do
         endif

         ljtype(1:stmax, pti) = ljtype_temp(1:stmax)
      end do
      deallocate( ljlen_temp, ljene_temp, ljtype_temp )

      ! Fill LJ table
      if(ljformat == LJFMT_TABLE) then
         ! From table (directly)
         open(unit = ljtable_io, file = trim(refsdirec) // '/' // ljtable_file, status = 'old', action = 'read')
         read(ljtable_io, *) ljtype_max
         allocate( ljlensq_mat(ljtype_max, ljtype_max), &
            ljene_mat(ljtype_max, ljtype_max) )
         do i = 1, ljtype_max
            read (ljtable_io, *) ljlensq_mat(i, 1:ljtype_max)
            ljlensq_mat(i, 1:ljtype_max) = ljlensq_mat(i, 1:ljtype_max) ** 2
         end do
         do i = 1, ljtype_max
            read (ljtable_io, *) ljene_mat(i, 1:ljtype_max)
         end do
         close(ljtable_io)
      else
         ! From LJ data
         allocate( ljlensq_mat(ljtype_max, ljtype_max), &
            ljene_mat(ljtype_max, ljtype_max) )
         do i = 1, ljtype_max
            select case(cmbrule)
             case(LJCMB_ARITH)    ! arithmetic mean
               ljlensq_mat(1:ljtype_max, i) = (( ljlen_temp_table(1:ljtype_max) &
                  + ljlen_temp_table(i) ) / 2.0) ** 2
             case(LJCMB_GEOM)     ! geometric mean
               ljlensq_mat(1:ljtype_max, i) = ljlen_temp_table(1:ljtype_max) &
                  * ljlen_temp_table(i)
             case default
               stop "Incorrect cmbrule"
            end select
            ljene_mat(1:ljtype_max, i) = sqrt( ljene_temp_table(1:ljtype_max)  &
               * ljene_temp_table(i) )
         end do
      endif
      deallocate(ljlen_temp_table, ljene_temp_table)

      return
   end subroutine get_ljtable

   subroutine calc_ljlrc(pti, correction)
      use sysvars, only: avevolume, nummol
      implicit none
      integer, intent(in) :: pti
      real, intent(out) :: correction
      real :: dens, ljeps, ljsgm2
      integer :: ui, vi
      dens = nummol(pti) / avevolume
      correction = 0.0
      do ui = 1, ptsite(0)           ! sum over solute sites
         do vi = 1, ptsite(pti)      ! sum over solvent sites
            ljeps = ljene_mat( ljtype(ui, 0), ljtype(vi, pti) )
            ljsgm2 = ljlensq_mat( ljtype(ui, 0), ljtype(vi, pti) )
            correction = correction + dens * enelj(ljeps, ljsgm2)
         end do
      end do
      return
   end subroutine calc_ljlrc

   ! calculation of the Lennard-Jones long-range correction energy
   ! the program corresponds closely to the relevant parts
   !     the get_pair_energy_block subroutine in realcal.F90
   real function enelj(ljeps, ljsgm2)
      implicit none
      real, save :: rbin = 1.0e-3
      integer, save :: numbin
      real, intent(in) :: ljeps, ljsgm2
      real, parameter :: PI = 3.1415926535897932
      real :: ljint, edev, dist, r
      real :: invr2, invr3, invr6, ljsgm3, ljsgm6, vdwa, vdwb, swth, swfac
      real, save :: lwljcut2, upljcut2, lwljcut3, upljcut3, lwljcut6, upljcut6
      real, save :: repA, repB, repC, attA, attB, attC
      logical, save :: do_swth, first_time = .true.
      integer :: i

      if(first_time) then
         if(lwljcut > upljcut) then
            stop "Incorrect setting of lwljcut and upljcut (lwljcut > upljcut)"
         else
            numbin = nint((upljcut - lwljcut) / rbin)
            if(numbin >= 1) then
               do_swth = .true.
               rbin = (upljcut - lwljcut) / real(numbin)

               lwljcut2 = lwljcut ** 2
               upljcut2 = upljcut ** 2
               if(ljswitch == LJSWT_FRC_CHM) then    ! force switch (CHARMM type)
                  lwljcut3 = lwljcut ** 3
                  upljcut3 = upljcut ** 3
                  lwljcut6 = lwljcut3 * lwljcut3
                  upljcut6 = upljcut3 * upljcut3
               endif
               if(ljswitch == LJSWT_FRC_GMX) then    ! force switch (GROMACS type)
                  call calc_gmx_switching_force_params(12, lwljcut, upljcut, repA, repB, repC)
                  call calc_gmx_switching_force_params(6,  lwljcut, upljcut, attA, attB, attC)
               endif
            else
               do_swth = .false.
            endif
         endif

         first_time = .false.
      endif

      ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
      ljsgm3 = sqrt(ljsgm6)

      ljint = 0.0

      ! r < lwljcut
      select case(ljswitch)
       case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
         ! do nothing
       case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
         vdwa = ljsgm6 * ljsgm6 / (lwljcut6 * upljcut6)
         vdwb = ljsgm6 / (lwljcut3 * upljcut3)
         edev = 4.0 * ljeps * (vdwa - vdwb)
         ljint = ljint + (4.0 * PI /3.0) * lwljcut3 * edev
       case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
         vdwa = ljsgm6 * ljsgm6 * repC
         vdwb = ljsgm6 * attC
         edev = 4.0 * ljeps * (vdwa - vdwb)
         ljint = ljint + (4.0 * PI /3.0) * lwljcut3 * edev
       case default
         stop "Unknown ljswitch"
      end select

      ! lwljcut < r < upljcut
      if(do_swth) then
         do i = 1, numbin
            r = lwljcut + (real(i) - 0.5) * rbin
            dist = r * r
            invr2 = ljsgm2 / dist
            invr6 = invr2 * invr2 * invr2
            select case(ljswitch)
             case(LJSWT_POT_CHM)             ! potential switch (CHRAMM type)
               swth = (2.0 * dist + upljcut2 - 3.0 * lwljcut2)             &
                  * ((dist - upljcut2) ** 2) / ((upljcut2 - lwljcut2) ** 3)
               edev = 4.0 * ljeps * invr6 * (invr6 - 1.0) * (1.0 - swth)
             case(LJSWT_POT_GMX)             ! potential switch (GROMACS type)
               swfac = (r - lwljcut) / (upljcut - lwljcut)
               swth = 1.0 - 10.0 * (swfac ** 3)                            &
                  + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5)
               edev = 4.0 * ljeps * invr6 * (invr6 - 1.0) * (1.0 - swth)
             case(LJSWT_FRC_CHM)             ! force switch (CHARMM type)
               invr3 = sqrt(invr6)
               vdwa = upljcut6 / (upljcut6 - lwljcut6)                     &
                  * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
               vdwb = upljcut3 / (upljcut3 - lwljcut3)                     &
                  * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
               edev = 4.0 * ljeps * ( invr6 * (invr6 - 1.0) - (vdwa - vdwb) )
             case(LJSWT_FRC_GMX)             ! force switch (GROMACS type)
               swfac = r - lwljcut
               vdwa = ljsgm6 * ljsgm6 *                                    &
                  (repA * (swfac ** 3) + repB * (swfac ** 4) + repC)
               vdwb = ljsgm6 *                                             &
                  (attA * (swfac ** 3) + attB * (swfac ** 4) + attC)
               edev = 4.0 * ljeps * (vdwa - vdwb)
            end select
            ljint = ljint + 4.0 * PI * r * r * rbin * edev
         end do
      endif

      ! r > upljcut
      invr3 = ljsgm3 / (upljcut ** 3)
      ljint = ljint + ljeps * ljsgm3 &
         * (16.0 * PI / 3.0) * ( (invr3 ** 3 / 3.0) - invr3 )

      enelj = ljint
      return
   end function enelj

   ! get the coefficients for gromacs force switching
   ! the program is taken from realcal.F90
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
end module uvcorrect
