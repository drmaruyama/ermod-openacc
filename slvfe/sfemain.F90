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

module sysvars
  implicit none

  character(len=5) :: clcond = 'merge'
  character(len=3) :: uvread = 'yes',    slfslt = 'yes',   ljlrc = 'not'
  character(len=3) :: infchk = 'not',    meshread = 'not', cumuint = 'not'
  character(len=3) :: write_mesherror = 'cnd'
  character(len=4) :: zerosft = 'orig',  wgtfnform = 'harm'
  character(len=3) :: refmerge = 'yes',  extsln = 'lin'
  character(len=3) :: wgtf2smpl = 'yes', slncor = 'not'
  character(len=3) :: normalize = 'yes', showdst= 'not'
  character(len=3) :: wrtzrsft = 'not',  readwgtfl = 'yes'

  integer :: numprm = 0                           ! initialized to 0
  integer :: numprm_def_inf_yes = 11   ! default numprm at infchk = 'yes'
  integer :: numprm_def_inf_not = 5    ! default numprm at infchk = 'not'

  integer :: numsln = 0, numref = 0, numdiv = 0   ! initialized to 0
  integer :: maxsln, maxref, numrun, prmmax
  integer :: numslv, ermax

  real :: inptemp = 300.0              ! temperature in Kelvin, initialized
  real :: temp, kT, slfeng
  real :: avevolume = 0.0              ! average volume of system, initialized

  integer :: pickgr = 3
  integer :: msemin = 1, msemax = 5
  real :: mesherr = 0.1                ! allowed mesh error in kcal/mol

  integer :: extthres_soln = 1, extthres_refs = 1
  integer :: minthres_soln = 0, minthres_refs = 0
  real, parameter :: zero = 0.0
  real :: error = 1.0e-8, tiny = 1.0e-8
  integer :: ermax_limit = 15000
  integer :: large = 500000, itrmax = 100
  
  character(len=1024) :: solndirec = 'soln'
  character(len=1024) :: refsdirec = 'refs'
  character(len=1024) :: wgtslnfl  = 'weight_soln'
  character(len=1024) :: wgtreffl  = 'weight_refs'
  character(len=1024) :: slndnspf  = 'engsln'
  character(len=1024) :: slncorpf  = 'corsln'
  character(len=1024) :: refdnspf  = 'engref'
  character(len=1024) :: refcorpf  = 'corref'
  character(len=1024) :: aveuvfile = 'aveuv.tt'
  character(len=1024) :: engmeshfile = 'EngMesh'
  character(len=1024) :: cumuintfl = 'cumsfe'
  character(len=10), parameter :: numbers='0123456789'
  
  real, dimension(:),     allocatable :: nummol
  integer, dimension(:),  allocatable :: rduvmax, rduvcore
  real, dimension(:),     allocatable :: rdcrd, rddst, rddns
  real, dimension(:,:),   allocatable :: rdslc, rdcor
  integer, dimension(:),  allocatable :: rdspec
  real, dimension(:,:,:), allocatable :: chmpt
  real, dimension(:),     allocatable :: aveuv
  real, dimension(:,:),   allocatable :: uvene, blockuv
  integer, dimension(:),  allocatable :: svgrp, svinf
  real, dimension(:),     allocatable :: wgtsln, wgtref
  
  namelist /fevars/ clcond, numprm, numsln, numref, numdiv, &
       uvread, slfslt, infchk, meshread, zerosft, wgtfnform, &
       refmerge, extsln, extthres_soln, extthres_refs, &
       minthres_soln, minthres_refs, &
       wgtf2smpl, slncor, normalize, showdst, wrtzrsft, readwgtfl, &
       inptemp, pickgr, write_mesherror, msemin, msemax, mesherr, &
       ljlrc, avevolume, &
       solndirec, refsdirec, wgtslnfl, wgtreffl, &
       slndnspf, slncorpf, refdnspf, refcorpf, &
       aveuvfile, engmeshfile, cumuint, cumuintfl, &
       ermax_limit, large, itrmax, error, tiny

contains

  subroutine init_sysvars
    implicit none
    character(len=*), parameter :: parmfname = 'parameters_fe'
    character(len=3) :: file_suf
    character(len=1024) :: opnfile
    integer, parameter :: iounit = 191, sufmax = 99
    integer :: ioerr, count_suf, i, j, srcnt, count_soln, count_refs
    logical :: file_exist
    
    open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
    if(ioerr /= 0) goto 99
    read(iounit, nml = fevars)
    close(iounit)
99  continue

    if(clcond == 'merge') then
       do srcnt = 1, 2
          do count_suf = 1, sufmax
             i = count_suf / 10
             j = mod(count_suf, 10)
             file_suf = '.' // numbers(i+1:i+1) // numbers(j+1:j+1)
             select case(srcnt)
             case(1)
                opnfile = trim(solndirec) // '/' // trim(slndnspf) // file_suf
             case(2)
                opnfile = trim(refsdirec) // '/' // trim(refdnspf) // file_suf
             end select
             inquire(file = opnfile, exist = file_exist)
             if( file_exist ) then
                if(srcnt == 1) count_soln = count_suf
                if(srcnt == 2) count_refs = count_suf
             else
                exit
             endif
          enddo
       enddo

       open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
       if(ioerr /= 0) goto 91
       read(iounit, nml = fevars)
       close(iounit)
91     continue

       if((numsln <= 0) .or. (numsln > count_soln)) numsln = count_soln
       if((numref <= 0) .or. (numref > count_refs)) numref = count_refs

       if((numdiv <= 0) .or. (numdiv >= numsln)) numdiv = numsln
       if(mod(numsln, numdiv) /= 0) then
          do i = numdiv + 1, numsln      ! find the larger and closest divisor
             if(mod(numsln, i) == 0) exit
          enddo
          numdiv = i
       endif
       if(refmerge == 'not') then        ! see subroutine datread for refmerge
          if(numdiv > numref) stop " With refmerge = 'not', numdiv needs to be not larger than numref"
          if(mod(numref, numdiv) /= 0) then
             write(6, "(A,i2,A,i2,A)") " Note: only ", numdiv * (numref / numdiv), &
   &" files out of ", numref, " engref and corref files prepared"
             numref = numdiv * (numref / numdiv)
          endif
       endif
    endif

    if(numprm <= 0) then                 ! default setting
       if(infchk == 'yes') then
          numprm = numprm_def_inf_yes    ! default numprm at infchk = 'yes'
       else
          numprm = numprm_def_inf_not    ! default numprm at infchk = 'not'
       endif
    endif

    if(pickgr < msemin) stop " Incorrect setting: pickgr < msemin not allowed"
    if(pickgr > msemax) stop " Incorrect setting: pickgr > msemax not allowed"
    if(pickgr > numprm) stop " Incorrect setting: pickgr > numprm not allowed"

  end subroutine init_sysvars
end module sysvars
