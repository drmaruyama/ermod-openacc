! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2024 Nobuyuki Matubayasi
! Copyright (C) 2010-2024 Shun Sakuraba
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

module sysread
  use sysvars, only: clcond, uvread, slfslt, slncor, &
       maxsln, maxref, numrun, &
       numslv, ermax, slfeng, nummol, &
       rdcrd, rddst, rddns, rdslc, rdcor, rdspec, &
       aveuv, uvene, blockuv, wgtsln, wgtref
  implicit none
  character(len=1024) :: engfile(5)
contains

  subroutine defcond
    use sysvars, only: infchk, meshread, readwgtfl, &
         solndirec, refsdirec, wgtslnfl, wgtreffl, &
         slndnspf, aveuvfile, refdnspf, engmeshfile, &
         numprm, prmmax, numsln, numref, numdiv, &
         inptemp, temp, kT, &
         ermax_limit, tiny, large, &
         rduvmax, rduvcore, &
         chmpt, svgrp, svinf, &
         get_suffix, suffix_of_engsln_is_tt
    implicit none
    logical :: file_exist
    integer :: group, inft, prmcnt, iduv, i, k, pti, ioerr
    real :: factor, crdnow, crdprev, crddif_now, crddif_prev
    character(len=1024) :: opnfile, lineread

    select case(clcond)
    case('basic', 'range', 'merge')
    case default
       stop ' The clcond parameter is incorrect'
    end select
    !
    select case(clcond)
    case('basic', 'range')
       write(6, "(A)") " What is the energy distribution in solution?"
       read(5, *) engfile(1)
       if(slncor == 'yes') then
          write(6, "(A)") " What is the energy correlation in solution?"
          read(5, *) engfile(2)
       endif
       write(6, "(A)") " What is the energy density for insertion?"
       read(5, *) engfile(3)
       write(6, "(A)") " What is the energy correlation for insertion?"
       read(5, *) engfile(4)
       if((infchk == 'yes') .and. (meshread == 'yes')) then   ! EngMesh file
          write(6, "(A)") " Which file has the meshes for energy coordinates?"
          read(5, *) engfile(5)
       endif
       maxsln = 1
       maxref = 1
       numrun = 1
       if(clcond == 'basic') prmmax = 1
       if(clcond == 'range') prmmax = numprm
    case('merge')
       maxsln = numsln
       maxref = numref
       numrun = numdiv
       prmmax = numprm
    end select
    !
    if(clcond == 'merge') then
       opnfile = trim(solndirec) // '/' // trim(slndnspf) // '.' // get_suffix(1, suffix_of_engsln_is_tt)
    else
       opnfile = engfile(1)
    endif
    open(unit = 71, file = opnfile, status = 'old')
    read(71, *) ! skip header line
    ermax = 0
    numslv = 0
    k = 0
    do iduv = 1, large
       read(71, *, end = 781) factor, factor, pti ! only read species
       if(pti /= k) then
          numslv = numslv + 1
          k = pti
       endif
       ermax = ermax + 1
    end do
781 continue
    close(71)
    !
    allocate( nummol(numslv) )
    allocate( rduvmax(numslv), rduvcore(numslv) )
    allocate( rdcrd(ermax), rddst(ermax), rddns(ermax) )
    if(slncor == 'yes') allocate( rdslc(ermax, ermax) )
    allocate( rdcor(ermax,ermax) )
    allocate( rdspec(ermax) )
    allocate( chmpt(0:numslv, prmmax, numrun), aveuv(numslv) )
    if((uvread /= 'not') .and. (clcond == 'merge')) then
       allocate( uvene(numslv, maxsln), blockuv(0:numslv, numrun) )
    endif
    allocate( svgrp(prmmax), svinf(prmmax) )
    allocate( wgtsln(maxsln), wgtref(maxref) )
    !
    ! opnfile is still engfile(1) or soln/engsln.01
    open(unit = 71, file = opnfile, status = 'old')
    read(71, *) ! skip header
    k = 0
    do iduv = 1, ermax
       read(71, *) factor, crdnow, pti, factor
       if(pti /= k) then
          rduvmax(pti) = 1
          rduvcore(pti) = 0
          k = pti
       else
          if(rduvmax(pti) < 1) stop "Bug in counting rduvmax"
          rduvmax(pti) = rduvmax(pti) + 1
          crddif_now = crdnow - crdprev
          if(rduvmax(pti) > 2) then
             if(abs(crddif_now - crddif_prev) < tiny) then
                rduvcore(pti) = 0                      ! linear mesh
             else
                rduvcore(pti) = rduvcore(pti) + 1      ! logarithmic mesh
             endif
          endif
          crddif_prev = crddif_now
       endif
       crdprev = crdnow
    end do
    close(71)

    if((infchk == 'yes') .and. (meshread == 'yes')) then      ! EngMesh file
       if(clcond == 'merge') then
          opnfile = trim(solndirec) // '/' // trim(engmeshfile)
       else
          opnfile = engfile(5)
       endif
       do pti = 1, numslv
          open(unit = 72, file = opnfile, status = 'old')
          do i = 1, large
             read(72, *, end = 729) k
             if(k == pti) then
                backspace(72)
                read(72, *, end = 729) k, k, rduvcore(pti)
                exit
             endif
          enddo
729       close(72)
       enddo
    endif

    if(sum( rduvmax(1:numslv) ) /= ermax) stop ' The file format is incorrect'
    if(ermax > ermax_limit) stop ' The number of energy bins is too large'

    select case(clcond)
    case('basic')
       write(6, "(A)") " How many data are grouped into one?"
       read(5, *) group
       inft = 0
       if(infchk == 'yes') then
          write(6, "(A)") " How many large-energy meshes are merged ? (in %)"
          read(5, *) inft
       endif
       svgrp(1) = group
       svinf(1) = inft
       write(6, "(A)") " What is the temperature in Kelvin?"
       read(5, *) temp
    case('range', 'merge')
       do prmcnt = 1, prmmax                 ! inft in % (percentage)
          if(infchk == 'yes') then
             if(prmcnt == 1) then  ; group = 1  ; inft = 0   ; endif
             if(prmcnt == 2) then  ; group = 1  ; inft = 60  ; endif
             if(prmcnt == 3) then  ; group = 1  ; inft = 80  ; endif
             if(prmcnt == 4) then  ; group = 1  ; inft = 100 ; endif
             if(prmcnt == 5) then  ; group = 2  ; inft = 0   ; endif
             if(prmcnt == 6) then  ; group = 3  ; inft = 0   ; endif
             if(prmcnt == 7) then  ; group = 4  ; inft = 0   ; endif
             if(prmcnt == 8) then  ; group = 5  ; inft = 0   ; endif
             if(prmcnt == 9) then  ; group = 5  ; inft = 60  ; endif
             if(prmcnt == 10) then ; group = 5  ; inft = 80  ; endif
             if(prmcnt == 11) then ; group = 5  ; inft = 100 ; endif
             if(prmcnt == 12) then ; group = 8  ; inft = 0   ; endif
             if(prmcnt >= 13) then
                group = 10 + (prmcnt - 13) * 5  ; inft = 0   ; endif
          else
             inft = 0
             if(prmcnt <= 10) then
                group = prmcnt
             else
                group = 10 + (prmcnt - 10) * 2
             endif
          endif
          svgrp(prmcnt) = group
          svinf(prmcnt) = inft
       end do
       temp = inptemp
    end select
    kT = temp * 8.314510e-3 / 4.184              ! kcal/mol

    if(uvread /= 'not') then
       select case(clcond)
       case('basic', 'range')
          write(6, "(A)") " What is average solute-solvent energy in solution?"
          read(5, *) aveuv(1:numslv)
       case('merge')
          opnfile = trim(solndirec) // '/' // trim(aveuvfile)
          inquire(file = opnfile, exist = file_exist)
          if( file_exist ) then
             open(unit = 82, file = opnfile, status = 'old')
             do i = 1, maxsln
                read(82, *) k, uvene(1:numslv, i)
             end do
             close(82)
          else
             uvread = 'not'
             write(6, '(A,A,A)') " Warning: Although the uvread parameter " // &
                  "was set to 'yes', it is changed into 'not' " // &
                  "since the ", trim(opnfile), " file was not found"
          endif
       end select
    endif

    ! solution system
    wgtsln(1:maxsln) = 1.0
    if((clcond == 'merge') .and. (readwgtfl == 'yes')) then
       opnfile = trim(solndirec) // '/' // trim(wgtslnfl)
       open(unit = 81, file = opnfile, status = 'old')
       do i = 1, maxsln
          read(81, *) k, wgtsln(i)
       end do
       close(81)
    endif
    factor = sum( wgtsln(1:maxsln) )
    wgtsln(1:maxsln) = wgtsln(1:maxsln) / factor

    ! reference solvent system
    wgtref(1:maxref) = 1.0
    if((clcond == 'merge') .and. (readwgtfl == 'yes')) then
       opnfile = trim(refsdirec) // '/' // trim(wgtreffl)
       open(unit = 81, file = opnfile, status = 'old')
       do i = 1, maxref
          read(81, *) k, wgtref(i)
       end do
       close(81)
    endif
    factor = sum( wgtref(1:maxref) )
    wgtref(1:maxref) = wgtref(1:maxref) / factor

    if(slfslt == 'yes') then
       select case(clcond)
       case('basic', 'range')
          write(6, "(A)") " What is the solute self-energy?"
          read(5, *) slfeng
       case('merge')
          if(readwgtfl == 'not') then
             stop "readwgtfl needs to be yes when slfslt is yes"
          endif
          slfeng = 0.0
          opnfile = trim(refsdirec) // '/' // trim(wgtreffl)
          inquire(file = opnfile, exist = file_exist)
          if( .not. file_exist ) stop " weight_refs is absent although slfslt is set to yes"
          open(unit = 81, file = opnfile, status = 'old')
          do i = 1, maxref
             read(81, '(a)') lineread
             read(lineread, * , iostat = ioerr) k, factor, factor
             if(ioerr /= 0) then
                slfslt = 'not'
                write(6, '(A)') " Warning: Although the slfslt parameter " // &
                   "was set to 'yes', it is changed into 'not'. Maybe " // &
                   "the MD was done without periodic boundary condition, " // &
                   "with PME employed for the isolated solute, " // &
                   "or with Coulombic interaction in its bare form."
                exit
             endif
             slfeng = slfeng + wgtref(i) * factor
          end do
          close(81)
       end select
    endif

    return
  end subroutine

  subroutine datread(cntrun)
    use sysvars, only: refmerge, tiny, numbers, &
         solndirec, refsdirec, slndnspf, slncorpf, refdnspf, refcorpf, &
         get_suffix, suffix_of_engsln_is_tt, suffix_of_engref_is_tt
    implicit none
    integer, intent(in) :: cntrun
    integer :: slnini, slnfin, refini, reffin, ecmin, ecmax
    integer :: iduv, iduvp, i, k, m, pti, cnt
    real :: factor, ampl, leftbin
    real, allocatable :: bin_consistency_check(:)
    logical :: num_different, suffix_is_tt
    real, allocatable :: cormat_temp(:, :)
    character(len=1024) :: opnfile
    character(len=3) :: suffnum

    select case(clcond)
    case('basic', 'range')
       slnini = 1
       slnfin = 1
       refini = 1
       reffin = 1
    case('merge')
       if(maxsln >= numrun) then
          k = maxsln / numrun
          slnini = (cntrun - 1) * k + 1
          slnfin = cntrun * k
       else
          slnini = mod(cntrun - 1, maxsln) + 1
          slnfin = slnini
       endif
       if(refmerge == 'not') then
          if(maxref >= numrun) then
             m = maxref / numrun
             refini = (cntrun - 1) * m + 1
             reffin = cntrun * m
          else
             refini = mod(cntrun - 1, maxref) + 1
             reffin = refini
          endif
       else
          refini = 1
          reffin = maxref
       endif
    end select

    rddst(:) = 0.0
    if(slncor == 'yes') rdslc(:,:) = 0.0
    if((cntrun == 1) .or. (refmerge == 'not')) then
       rddns(:) = 0.0
       rdcor(:,:) = 0.0
    endif

    allocate(bin_consistency_check(ermax))
    bin_consistency_check(:) = 0.0

    ! FIXME: this part is kinda spaghetti and should be rewritten WITHOUT looping by cnt!
    do cnt = 1, 4
       if((cnt == 2) .and. (slncor /= 'yes')) cycle
       if((cnt >= 3) .and. (cntrun > 1) .and. (refmerge == 'yes')) cycle

       if(cnt <= 2) then                         ! solution
          ecmin = slnini
          ecmax = slnfin
          factor = sum( wgtsln(ecmin:ecmax) )
          wgtsln(ecmin:ecmax) = wgtsln(ecmin:ecmax) / factor
          suffix_is_tt = suffix_of_engsln_is_tt
       endif
       if(cnt >= 3) then                         ! reference solvent
          ecmin = refini
          ecmax = reffin
          factor = sum( wgtref(ecmin:ecmax) )
          wgtref(ecmin:ecmax) = wgtref(ecmin:ecmax) / factor
          suffix_is_tt = suffix_of_engref_is_tt
       endif

       do i = ecmin, ecmax
          select case(clcond)
          case('basic', 'range')
             opnfile = engfile(cnt)
          case('merge')
             suffnum = '.' // get_suffix(i, suffix_is_tt)
             if(cnt == 1) opnfile = trim(solndirec) // '/' // trim(slndnspf)
             if(cnt == 2) opnfile = trim(solndirec) // '/' // trim(slncorpf)
             if(cnt == 3) opnfile = trim(refsdirec) // '/' // trim(refdnspf)
             if(cnt == 4) opnfile = trim(refsdirec) // '/' // trim(refcorpf)
             opnfile = trim(opnfile) // suffnum
          end select
          if((cnt == 1) .or. (cnt == 3)) then    ! 1-D distribution
             open(unit = 71, file = opnfile, status = 'old')
             read(71, *) ! skip header
             k = 0
             m = 0
             do iduv = 1, ermax
                read(71, *) leftbin, rdcrd(iduv), pti, factor
                if(cnt == 1) then
                   bin_consistency_check(iduv) = leftbin
                elseif(cnt == 3) then
                   if(bin_consistency_check(iduv) /= leftbin) then
                      stop "Solution and reference system energy coordinates are inconsitent"
                   endif
                end if
                if(pti /= k) then
                   k = pti
                   m = m + 1
                endif
                if(cnt == 1) rddst(iduv) =rddst(iduv) + wgtsln(i) * factor
                if(cnt == 3) rddns(iduv) =rddns(iduv) + wgtref(i) * factor
                rdspec(iduv) = m
             enddo
             close(71)
          endif
          if((cnt == 2) .or. (cnt == 4)) then    ! 2-D correlation matrix
             allocate(cormat_temp(ermax, ermax))
             open(unit = 72, file = opnfile, status = 'old', form = "UNFORMATTED")
             read(72) cormat_temp
             close(72)
             if(cnt == 2) rdslc(:, :) = rdslc(:, :) + wgtsln(i) * cormat_temp(:, :)
             if(cnt == 4) rdcor(:, :) = rdcor(:, :) + wgtref(i) * cormat_temp(:, :)
             deallocate(cormat_temp)
          endif
       end do
    end do

    if(cntrun == 1) write(6, *)
    do pti = 1, numslv
       factor = sum( rddst, mask = (rdspec == pti) )
       ampl   = sum( rddns, mask = (rdspec == pti) )
       num_different = .false.
       if(abs(factor - ampl) > tiny) num_different = .true.
       if(cntrun > 1) then
          if(nint(factor) /= nint(nummol(pti))) num_different = .true.
          if(nint(ampl) /= nint(nummol(pti))) num_different = .true.
       endif
       if(num_different) then
          write(6, '(A,i4)') '  Incorrect normalization at ', pti
          stop
       endif
       if(cntrun == 1) then
          nummol(pti) = real(nint(factor))
          write(6, '(A,i3,A,i12)') '  Number of the ', pti, '-th solvent  = ', nint(nummol(pti))
       endif
    end do
    if(cntrun == 1) write(6, *)
!
    if((uvread /= 'not') .and. (clcond == 'merge')) then
       do pti = 1, numslv
          aveuv(pti) = sum( wgtsln(slnini:slnfin) * uvene(pti, slnini:slnfin) )
       end do
       blockuv(1:numslv, cntrun) = aveuv(1:numslv)
       blockuv(0, cntrun) = sum( blockuv(1:numslv, cntrun) )
       if(slfslt == 'yes') blockuv(0, cntrun) = blockuv(0, cntrun) + slfeng
    endif

    return
  end subroutine datread
end module sysread


module sfecalc
  use sysvars, only: invmtrx, zerosft, wgtfnform, slncor, &
       numslv, ermax, nummol, kT, itrmax, zero, error, tiny, &
       rduvmax, rduvcore, &
       rdcrd, rddst, rddns, rdslc, rdcor, rdspec
  implicit none
  integer, dimension(:), allocatable :: idrduv, uvmax
  real, dimension(:),    allocatable :: uvcrd, edist, edens
  real, dimension(:,:),  allocatable :: edscr, ecorr
  integer, dimension(:), allocatable :: uvspec
  real, dimension(:),    allocatable :: slncv, inscv, sdrcv
  real, dimension(:),    allocatable :: zrsln, zrref, zrsdr
  integer gemax
contains
  subroutine posv_wrap(n, mat, vec, info)
    use iso_c_binding
    use cuBlas_v2
    use cuSolverDn
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: mat(n, n)
    real, intent(inout) :: vec(n)
    integer, intent(out) :: info
    type(cuSolverDnHandle) :: h
    integer :: istat, lwork
    real, allocatable :: input_mat(:, :)
    real, allocatable :: input_vec(:), residual(:)
    real, allocatable :: work(:)

    istat = cuSolverDnCreate(h)

    lwork = -1
    !$acc data copyin(mat)
    !$acc host_data use_device(mat)
#ifdef DP
    istat = cusolverDnDpotrf_bufferSize(h, CUBLAS_FILL_MODE_UPPER, &
         n, mat, n, lwork)
#else
    istat = cusolverDnDpotrf_bufferSize(h, CUBLAS_FILL_MODE_UPPER, &
         n, mat, n, lwork)
#endif
    !$acc end host_data
    !$acc end data

    if (istat /= 0) then
       return
    endif

    allocate( input_mat(n, n), input_vec(n), residual(n) )
    input_mat(:, :) = mat(:, :)
    input_vec(:) = vec(:)
    allocate(work(lwork))
    !$acc enter data create(work)
    !$acc data copy(mat, vec) copyout(info)
    !$acc host_data use_device(mat, vec, work, info)
#ifdef DP
    istat = cusolverDnDpotrf(h, CUBLAS_FILL_MODE_UPPER, &
         n, mat, n, work, lwork, info)
    istat = cusolverDnDpotrs(h, CUBLAS_FILL_MODE_UPPER, &
         n, 1, mat, n, vec, n, info)
#else
    istat = cusolverDnSpotrf(h, CUBLAS_FILL_MODE_UPPER, &
         n, mat, n, work, lwork, info)
    istat = cusolverDnSpotrs(h, CUBLAS_FILL_MODE_UPPER, &
         n, 1, mat, n, vec, n, info)
#endif
    !$acc end host_data
    !$acc end data

    if (info == 0) then
       residual(:) = matmul( input_mat(:, :), vec(:) )
       residual(:) = abs( residual(:) - input_vec(:) )
       if(maxval( residual(:) ) > error) info = 1    ! invalid solution
    endif

    deallocate( input_mat, input_vec, residual )
  end subroutine posv_wrap

  subroutine syevr_wrap(n, mat, eigval, info)
    use cuBlas_v2
    use cuSolverDn
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: mat(n, n)
    real, intent(out) :: eigval(n)
    integer, intent(out) :: info
    type(cuSolverDnHandle) :: h
    integer :: istat, lwork
    real(kind=8), allocatable :: work(:)

    istat = cuSolverDnCreate(h)

    lwork = -1
    !$acc data copyin(mat, eigval)
    !$acc host_data use_device(mat, eigval)
#ifdef DP
    istat = cusolverDnDsyevd_bufferSize(h, CUSOLVER_EIG_MODE_VECTOR, &
         CUBLAS_FILL_MODE_UPPER, n, mat, n, eigval, lwork)
#else
    istat = cusolverDnSsyevd_bufferSize(h, CUSOLVER_EIG_MODE_VECTOR, &
         CUBLAS_FILL_MODE_UPPER, n, mat, n, eigval, lwork)
#endif
    !$acc end host_data
    !$acc end data

    if (istat /= 0) then
       return
    endif

    allocate(work(lwork))
    !$acc enter data create(work)
    !$acc data copy(mat) copyout(eigval, info)
    !$acc host_data use_device(mat, eigval, work, info)
#ifdef DP
    istat = cusolverDnDsyevd(h, CUSOLVER_EIG_MODE_VECTOR, &
         CUBLAS_FILL_MODE_UPPER, n, mat, n, &
         eigval, work, lwork, info)
#else
    istat = cusolverDnSsyevd(h, CUSOLVER_EIG_MODE_VECTOR, &
         CUBLAS_FILL_MODE_UPPER, n, mat, n, &
         eigval, work, lwork, info)
#endif
    !$acc end host_data
    !$acc end data
    !$acc exit data delete(work)
    deallocate(work)
  end subroutine syevr_wrap

  subroutine chmpot(prmcnt, cntrun)
    use sysvars, only: uvread, slfslt, ljlrc, normalize, showdst, wrtzrsft, &
         slfeng, chmpt, aveuv, svgrp, svinf, &
         numrun, pickgr, &
         minthres_soln, minthres_refs, &
         cumuint, cumuintfl, &
         get_suffix, &
         numbers
    use uvcorrect, only: ljcorrect
    implicit none
    integer, intent(in) :: prmcnt, cntrun
    integer :: group, inft
    integer :: iduv, iduvp, pti, cnt, j, k, m, cntdiv, ge_perslv
    real :: factor, ampl, slvfe, uvpot, lcent, lcsln, lcref
    real :: soln_zero, refs_zero
    integer, dimension(:), allocatable :: gpnum
    real, dimension(:,:), allocatable, save :: cumsfe
    real, dimension(:), allocatable :: cumu_coord, cumu_write
    logical :: cumu_process, cumu_homoform
    integer, parameter :: cumu_io = 51
    character(len=1024) :: opnfile

    group = svgrp(prmcnt)
    inft = svinf(prmcnt)

    allocate( idrduv(ermax), uvmax(numslv) )

    gemax = 0
    do pti = 1, numslv
       k = rduvcore(pti) * inft / 100
       m = (rduvmax(pti) - k) / group
       uvmax(pti) = m
       gemax = gemax + m
       if(pti == 1) then
          cnt = 1
          j = 0
       else
          cnt = 1 + sum( rduvmax(1:(pti - 1)) )
          j = sum( uvmax(1:(pti - 1)) )
       endif
       do iduv = cnt, cnt + rduvmax(pti) - 1
          k = (iduv - cnt) / group + 1
          if(k >= m) k = m
          idrduv(iduv) = k + j
       end do
       if((pti == numslv) .and. (cnt + rduvmax(pti) - 1 /= ermax)) then
          print *, "Error: The total no. of meshes does not match with input"
          print *, "(Sum should be", ermax, " but was", cnt + rduvmax(pti) - 1, ")"
          stop
       endif
    end do

    allocate( uvcrd(gemax), edist(gemax), edens(gemax), uvspec(gemax) )
    uvcrd(:) = 0.0
    edist(:) = 0.0
    edens(:) = 0.0

    if(slncor == 'yes') then
       allocate( edscr(gemax, gemax) )
       edscr(:,:) = 0.0
    endif
    allocate( ecorr(gemax, gemax) )
    ecorr(:,:) = 0.0

    allocate( slncv(gemax), inscv(gemax), zrsln(numslv), zrref(numslv) )
    slncv(:) = 0.0
    inscv(:) = 0.0
    zrsln(:) = 0.0
    zrref(:) = 0.0
    if(slncor == 'yes') then
       allocate( sdrcv(gemax), zrsdr(numslv) )
       sdrcv(:) = 0.0
       zrsdr(:) = 0.0
    endif

    allocate( gpnum(gemax) )
    gpnum(:) = 0
    do iduv = 1, ermax
       k = idrduv(iduv)
       uvcrd(k) = uvcrd(k) + rdcrd(iduv)
       gpnum(k) = gpnum(k) + 1
       uvspec(k) = rdspec(iduv)
    end do
    do k = 1, gemax
       if(gpnum(k) > 0) uvcrd(k) = uvcrd(k) / real(gpnum(k))
    end do
    cnt = rduvmax(1)
    k = uvmax(1)
    do pti = 1, numslv
       if(pti > 1) then
          cnt = cnt + rduvmax(pti)
          k = k + uvmax(pti)
       endif
       uvcrd(k) = rdcrd(cnt)
    end do
    deallocate( gpnum )

    do cnt = 1, 2
       do iduv = 1, ermax
          k = idrduv(iduv)
          if(cnt == 1) edist(k) = edist(k) + rddst(iduv)
          if(cnt == 2) edens(k) = edens(k) + rddns(iduv)
       end do
       if((cnt == 1) .and. (slncor /= 'yes')) goto 1115
       do iduv = 1, ermax
          do iduvp = 1, ermax
             k = idrduv(iduv)
             m = idrduv(iduvp)
             if(cnt == 1) edscr(m,k) = edscr(m,k) + rdslc(iduvp, iduv)
             if(cnt == 2) ecorr(m,k) = ecorr(m,k) + rdcor(iduvp, iduv)
          end do
       end do
1115   continue
    end do

    if(normalize == 'yes') call distnorm
    if(showdst == 'yes')   call distshow

    call getslncv
    call getinscv

    ! when uvread = 'not', aveuv is constructed for each set of group and inft
    !                      from the distribution function (edist) in solution
    ! see the datread subroutine for treatment of aveuv when uvread = 'yes'
    !    in this case, aveuv is common among all the sets of group and inft
    ! aveuv at uvread = 'not' is constructed only once at the finest mesh
    !    when the next if(prmcnt == 1) and corresponding endif are made active
    !    mesh error is then the same as the one at uvread = 'yes' (default)
    if(uvread == 'not') then
       !  if(prmcnt == 1) then       ! at the first call of the chmpot subroutine
       do pti = 1, numslv
          aveuv(pti) = sum( uvcrd * edist, mask = (uvspec == pti) )
       end do
       ! LJ long-range correction at each construction of aveuv
       if(ljlrc == 'yes') call ljcorrect(cntrun)
       !  endif
    else  ! default setting with aveuv read in the datread subroutine
       ! LJ long-range correction at the first call of the chmpot subroutine
       if(prmcnt == 1) then       ! at the first call of the chmpot subroutine
          if(ljlrc == 'yes') call ljcorrect(cntrun)
       endif
    endif

    cumu_process = .false.
    if((cumuint == 'yes') .and. (group == pickgr) .and. (inft == 0)) then
       cumu_process = .true.                  ! cumulative integral treated
       if(cntrun == 1) allocate( cumsfe(gemax, 0:numrun) )
    endif

    soln_zero = minthres_soln * minval( rddst, mask = (rddst > zero) )
    refs_zero = minthres_refs * minval( rddns, mask = (rddns > zero) )

    do pti = 1, numslv
       uvpot = 0.0
       slvfe = 0.0
       do iduv = 1, gemax
          if(uvspec(iduv) == pti) then
             if((edist(iduv) <= soln_zero) .and. &
                  (edens(iduv) <= refs_zero)) goto 5009
             uvpot = uvpot + uvcrd(iduv) * edist(iduv)
             slvfe = slvfe - kT * (edist(iduv) - edens(iduv))

             ! kT * log(edist/edens)
             lcent = - (slncv(iduv) + zrsln(pti) + uvcrd(iduv))
             if((slncor == 'yes') .and. (edist(iduv) > soln_zero) &
                  .and. (edens(iduv) <= refs_zero)) then
                ! special case to be examined and fixed later
                ampl = lcent * edens(iduv) / edist(iduv)
                lcent = ampl - (zrsln(pti) + uvcrd(iduv)) &
                     * (1.0 - edens(iduv) / edist(iduv))
             endif
             slvfe = slvfe + lcent * edist(iduv)

             lcsln = pyhnc(slncv(iduv), 1)    ! solution
             lcref = pyhnc(inscv(iduv), 2)    ! reference solvent
             if((slncor == 'yes') .and. (edist(iduv) > soln_zero) &
                  .and. (edens(iduv) <= refs_zero)) then
                ! special case to be examined and fixed later
                lcsln = pyhnc(sdrcv(iduv) + zrsdr(pti), 3)
             endif
             ampl = sfewgt(edist(iduv), edens(iduv))
             factor = ampl * lcsln + (1.0 - ampl) * lcref
             slvfe = slvfe + kT * factor * (edist(iduv) - edens(iduv))
5009         continue

             ! cumulative integral
             if(cumu_process) cumsfe(iduv, cntrun) = uvpot + slvfe
          endif
       end do
       chmpt(pti, prmcnt, cntrun) = slvfe + aveuv(pti)
    end do

    ! cumulative integral stored
    if(cumu_process .and. (cntrun == numrun)) then
       do iduv = 1, gemax                     ! averaged cumsfe
          cumsfe(iduv, 0) = sum( cumsfe(iduv, 1:numrun) ) / real(numrun)
       enddo

       do cntdiv = 0, numrun
          ! only cntdiv = 0 when numrun = 1, and cntdiv >= 0 when numrun > 1
          if((numrun == 1) .and. (cntdiv /= 0)) cycle
          if(cntdiv == 0) then                ! averaged cumsfe
             opnfile = trim(cumuintfl)
          else                                ! cumsfe in each block
             opnfile = trim(cumuintfl) // get_suffix(cntdiv)
          endif
          open(unit = cumu_io, file = opnfile, status = 'replace')
          if(numslv == 1) then
             do iduv = 1, gemax
                write(cumu_io, '(g15.5, f12.5)') &
                     uvcrd(iduv), cumsfe(iduv, cntdiv)
             enddo
          else
             ge_perslv = gemax / numslv
             if(all(uvmax(1:numslv) == ge_perslv)) then
                cumu_homoform = .true.
                allocate( cumu_coord(numslv) )
                do iduv = 1, ge_perslv
                   do j = 1, numslv
                      cumu_coord(j) = uvcrd(iduv + (j - 1) * ge_perslv)
                   enddo
                   if(all(cumu_coord(2:numslv) == cumu_coord(1))) then
                      ! do nothing
                   else             ! unless all the coordinates are identical
                      cumu_homoform = .false.
                   endif
                enddo
                deallocate( cumu_coord )
             else
                cumu_homoform = .false.
             endif
             if(cumu_homoform) then ! same coordinates for all the solvent
                allocate( cumu_write(numslv) )
                do iduv = 1, ge_perslv
                   do j = 1, numslv
                      cumu_write(j) = cumsfe(iduv + (j - 1) * ge_perslv, cntdiv)
                   enddo
                   factor = sum( cumu_write(1:numslv) )
                   write(cumu_io, '(g15.5, 999f12.5)') &
                        uvcrd(iduv), factor, cumu_write(1:numslv)
                enddo
                deallocate( cumu_write )
             else
                do iduv = 1, gemax
                   write(cumu_io, '(g15.5, i5, f12.5)') &
                        uvcrd(iduv), uvspec(iduv), cumsfe(iduv, cntdiv)
                enddo
             endif
          endif
          endfile(cumu_io)
          close(cumu_io)
       enddo
       deallocate( cumsfe )
    endif

    chmpt(0, prmcnt, cntrun) = sum( chmpt(1:numslv, prmcnt, cntrun) )
    if(slfslt == 'yes') chmpt(0, prmcnt, cntrun) = chmpt(0, prmcnt, cntrun) + slfeng

    if(wrtzrsft == 'yes') then
       write(6, 381) zrsln(1:numslv)
       write(6, 382) zrref(1:numslv)
       if(slncor == 'yes') write(6, 383) zrsdr(1:numslv)
381    format('  Zero shift for solution             = ', 9999f12.4)
382    format('  Zero shift for reference solvent    = ', 9999f12.4)
383    format('  Zero shift for solution correlation = ', 9999f12.4)
    endif

    deallocate( slncv, inscv, zrsln, zrref )
    if(slncor == 'yes') deallocate( sdrcv, zrsdr, edscr )
    deallocate( uvcrd, edist, edens, ecorr, uvspec, idrduv, uvmax )

    return
  end subroutine chmpot

  subroutine getslncv
    use sysvars, only: extsln, extthres_soln, extthres_refs
    implicit none
    integer :: iduv, iduvp, pti, j, k, m
    real :: cvzero, factor, mat11, mat22, mat12, mat21, min_rddst, min_rddns
    real, dimension(:), allocatable :: work
    integer, parameter :: ofdmp = 10 ! factor to suppress the integer overflow
    logical, dimension(:), allocatable :: ext_target

    min_rddst = minval( rddst, mask = (rddst > zero) )
    min_rddns = minval( rddns, mask = (rddns > zero) )
    allocate( ext_target(gemax) )
    ext_target(:) = .true.

    do iduv = 1, gemax
       factor = edist(iduv) / min_rddst
       m = ofdmp * extthres_soln
       if(factor > real(m)) then
          j = m
       else
          j = nint(factor)
       endif
       factor = edens(iduv) / min_rddns
       m = ofdmp * extthres_refs
       if(factor > real(m)) then
          k = m
       else
          k = nint(factor)
       endif
       if((j < extthres_soln) .or. (k < extthres_refs)) ext_target(iduv) = .false.
    enddo

    do iduv = 1, gemax
       if(ext_target(iduv)) then
          slncv(iduv) = - kT * log(edist(iduv) / edens(iduv)) - uvcrd(iduv)
       endif
    end do

    do iduv = 1, gemax
       if(.not. ext_target(iduv)) then
          if(edist(iduv) <= zero) then
             slncv(iduv) = 0.0
             cycle
          endif
          pti = uvspec(iduv)
          m = 1
          k = gemax
          do iduvp = 1, iduv - 1
             if((uvspec(iduvp) == pti) .and. (ext_target(iduvp)) .and. &
                  (m < iduvp)) m = iduvp
          enddo
          do iduvp = gemax, iduv + 1, -1
             if((uvspec(iduvp) == pti) .and. (ext_target(iduvp)) .and. &
                  (k > iduvp)) k = iduvp
          enddo
          !
          if(extsln == 'sim') then
             if(abs(m - iduv) <  abs(k - iduv)) factor = slncv(m)
             if(abs(m - iduv) >  abs(k - iduv)) factor = slncv(k)
             if(abs(m - iduv) == abs(k - iduv)) then
                factor = (slncv(m) + slncv(k)) / 2.0
             endif
          else
             if(ext_target(m) .and. ext_target(k)) then
                j = k
                if(abs(m - iduv) < abs(k - iduv)) j = m
             elseif(ext_target(m) .and. (.not. ext_target(k))) then
                j = m
             elseif((.not. ext_target(m)) .and. ext_target(k)) then
                j = k
             else  ! (.not. ext_target(m)) .and. (.not. ext_target(k))
                write(6, *) "Extrapolation is not possible at ", uvcrd(iduv)
                stop
             endif
             allocate( work(gemax) )
             work(:) = 0.0
             do iduvp = 1, gemax
                if((uvspec(iduvp) == pti) .and. (ext_target(iduvp))) then
                   if(iduvp == iduv) stop ' A bug in program or data'
                   factor = uvcrd(iduvp) - uvcrd(j)
                   if(iduvp < iduv) then
                      factor = - factor - 2.0 * (uvcrd(j) - uvcrd(iduv))
                   endif
                   work(iduvp) = exp(- factor / kT) &
                        * wgtdst(iduvp, 1, 'extsl', wgtfnform)
                endif
             enddo
             factor = sum( work, mask = (work > zero) )
             do iduvp = 1, gemax
                work(iduvp) = work(iduvp) / factor
             enddo
             mat11 = sum( work * uvcrd, mask = (work > zero) )
             mat22 = sum( work * uvcrd * uvcrd, mask = (work > zero) )
             mat12 = sum( work * slncv, mask = (work > zero) )
             mat21 = sum( work * uvcrd * slncv, mask = (work > zero) )
             work(1) = (mat22 * mat12 - mat11 * mat21) / (mat22 - mat11 ** 2)
             work(2) = (mat21 - mat11 * mat12) / (mat22 - mat11 ** 2)
             factor = work(1) + work(2) * uvcrd(iduv)
             deallocate( work )
          endif
          slncv(iduv) = factor
       endif
    enddo

    deallocate( ext_target )

    do pti = 1, numslv
       select case(zerosft)
       case('eczr', 'orig')
          cvzero = 0.0
       case('mxco')
          factor = wgtmxco(pti)
          cvzero = factor* cvfcen(pti, 1, 'slncv', wgtfnform, 'not') &
               - (1.0 - factor) * cvfcen(pti, 2, 'uvcrd', 'smpl', 'yes')
       case('zero')
          cvzero = cvfcen(pti, 1, 'slncv', wgtfnform, 'yes')
       case('cntr')
          cvzero = cvfcen(pti, 1, 'slncv', wgtfnform, 'not')
       case default
          stop ' zerosft not properly set '
       end select
       do iduv = 1, gemax
          if(uvspec(iduv) == pti) slncv(iduv) = slncv(iduv) - cvzero
       end do
       zrsln(pti) = cvzero
    end do

    return
  end subroutine getslncv

  subroutine getinscv
    implicit none
    integer :: iduv, iduvp, pti, cnt, k, inv_info
    real :: dns, dnsp, dmcr, ddiff, cvzero, factor
    real, dimension(:),   allocatable :: edvec, work, zerouv, egnval
    real, dimension(:,:), allocatable :: edmcr, edmcr_nonzeroec
    character(len=5) :: invmtrx_cnt
    logical, save :: first_time = .true.
    !
    if(first_time) then
       first_time = .false.
       if(invmtrx /= 'syevr') invmtrx = 'posv'    ! default is 'posv'
    endif

    do cnt = 1, 2     ! cnt = 1: solution   cnt = 2: reference solvent
       if((cnt == 1) .and. (slncor /= 'yes')) cycle
       invmtrx_cnt = invmtrx
       allocate( edvec(gemax), work(gemax), zerouv(numslv) )
       allocate( edmcr(gemax, gemax) )

       ! construction of the distribution and the correlation matrix
       select case(cnt)
       case(1)       ! solution
          edvec(:) = edist(:)
          edmcr(:,:) = edscr(:,:)
       case(2)       ! reference solvent
          edvec(:) = edens(:)
          edmcr(:,:) = ecorr(:,:)
       case default
          stop "Unknown type of the system"
       end select
       do iduv = 1, gemax
          do iduvp = 1, gemax
             dns  = edvec(iduv)
             dnsp = edvec(iduvp)
             if( (dnsp > zero) .and. (dns > zero) ) then
                dmcr = edmcr(iduvp, iduv) - dnsp * dns
             else
                ! edmcr is a unit matrix for the part with edvec = 0
                if(iduv == iduvp) then
                   dmcr = 1.0
                else
                   dmcr = 0.0
                endif
             endif
             edmcr(iduvp, iduv) = dmcr
          end do
       end do

       if(invmtrx_cnt == 'posv') then
          ! correlation matrix only for non-zero solute-solvent energy part
          ! unit matrix and null vector for the zero energy part
          ! matrix inversion over the non-zero energy part
          allocate( edmcr_nonzeroec(gemax, gemax) )
          edmcr_nonzeroec(:,:) = edmcr(:,:)
          work(:) = - kT * ( edist(:) - edens(:) )
          where(edvec(:) <= zero) work(:) = 0.0
          do pti = 1, numslv
             k = zeroec(pti, cnt)        ! coordinate at which energy = 0
             edmcr_nonzeroec(k, :) = 0.0
             edmcr_nonzeroec(:, k) = 0.0
             edmcr_nonzeroec(k, k) = 1.0
             work(k) = 0.0
          enddo
          call posv_wrap(gemax, edmcr_nonzeroec, work, inv_info)
          if(inv_info == 0) then
             if(cnt == 1) then           ! solution
                sdrcv(:) = work(:)
             else                        ! reference solvent
                inscv(:) = work(:)
             endif
          else
             ! write(6, *) "Failed in solving the linear equations"
             invmtrx_cnt = 'syevr'       ! alternative solution with syevr
          endif
          deallocate( edmcr_nonzeroec )
       endif

       if(invmtrx_cnt == 'syevr') then
          ! pseudoinversion to solve the linear equations
          allocate( egnval(gemax) )
          call syevr_wrap(gemax, edmcr, egnval, inv_info)
          if (inv_info /= 0) stop "Failed inversion of correlation matrix"
          ! - kT * edmcr^(-1) * (edist - edens)
          ! first numslv eigenvalues skipped in the pseudoinversion
          pti = numslv + 1
          do iduv = pti, gemax
             ! inner product between (edist - edens) and eigenvector
             factor = 0.0
             do iduvp = 1, gemax
                if(edvec(iduvp) > zero) factor = factor &
                     + (edist(iduvp) - edens(iduvp)) * edmcr(iduvp, iduv)
             end do
             ! multiplication by -kT and division by non-zero eigenvalue
             work(iduv) = - kT * factor / egnval(iduv)
          end do
          do iduv = 1, gemax
             ! sum over eigenvectors
             factor = dot_product( edmcr(iduv, pti:gemax), work(pti:gemax) )
             if(cnt == 1) then           ! solution
                sdrcv(iduv) = factor
             else                        ! reference solvent
                inscv(iduv) = factor
             endif
          end do
          deallocate( egnval )
       endif

       ! reference solute-solvent energy to fix the additive constant below
       do pti = 1, numslv
          select case(zerosft)
          case('eczr', 'orig')
             k = zeroec(pti, cnt)        ! coordinate at which energy = 0
             if(cnt == 1) cvzero = sdrcv(k)
             if(cnt == 2) cvzero = inscv(k)
          case('mxco')
             cvzero = cvfcen(pti, cnt, 'inscv', 'smpl', 'yes')
          case default
             cvzero = 0.0
          end select
          zerouv(pti) = cvzero
       end do

       ! conversion from solute-solvent potential to indirect PMF
       do iduv = 1, gemax
          ! - kT (edist - edens) / edvec - (solute-solvent potential)
          if(edvec(iduv) > zero) then
             ddiff = - kT * (edist(iduv) - edens(iduv)) / edvec(iduv)
             if(cnt == 1) factor = ddiff - sdrcv(iduv)
             if(cnt == 2) factor = ddiff - inscv(iduv)
          else
             factor = 0.0
          endif
          if(cnt == 1) then              ! solution
             sdrcv(iduv) = factor
          else                           ! reference solvent
             inscv(iduv) = factor
          endif
       end do

       ! fixing the additive constant for indirect PMF
       do pti = 1, numslv
          select case(zerosft)
          case('eczr', 'orig')
             ! negative of solute-solvent energy at epsilon = 0
             cvzero = - zerouv(pti)
          case('mxco')
             factor = wgtmxco(pti)
             cvzero = factor * cvfcen(pti, cnt, 'inscv', wgtfnform, 'not') &
                  - (1.0 - factor) * zerouv(pti)
          case('zero')
             cvzero = cvfcen(pti, cnt, 'inscv', wgtfnform, 'yes')
          case('cntr')
             cvzero = cvfcen(pti, cnt, 'inscv', wgtfnform, 'not')
          case default
             stop ' zerosft not properly set '
          end select
          do iduv = 1, gemax
             if(uvspec(iduv) == pti) then
                if(cnt == 1) then        ! solution
                   sdrcv(iduv) = sdrcv(iduv) - cvzero
                else                     ! reference solvent
                   inscv(iduv) = inscv(iduv) - cvzero
                endif
             endif
          end do
          if(cnt == 1) then              ! solution
             zrsdr(pti) = cvzero
          else                           ! reference solvent
             zrref(pti) = cvzero
          endif
       end do

       deallocate( edvec, work, zerouv, edmcr )
    end do

    return
  end subroutine getinscv

  real function wgtmxco(pti)
    implicit none
    integer, intent(in) :: pti
    real :: numpt
    numpt = nummol(pti)
    wgtmxco = 1.0 / numpt
    return
  end function wgtmxco

  real function cvfcen(pti, cnt, systype, wgttype, engtype)
    implicit none
    ! pti : identifier of solvent species
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: pti, cnt
    character(len=5), intent(in) :: systype
    character(len=4), intent(in) :: wgttype
    character(len=3), intent(in) :: engtype
    integer :: iduv
    logical :: errtag
    real :: factor, cvfnc
    real, dimension(:), allocatable :: weight
    allocate( weight(gemax) )
    call getwght(weight, pti, cnt, systype, wgttype, engtype)
    factor = 0.0
    errtag = .false.
    do iduv = 1, gemax
       if(uvspec(iduv) == pti) then
          select case(systype)
          case('slncv')
             if(cnt == 1) cvfnc = slncv(iduv)
             if(cnt == 2) errtag = .true.
          case('inscv')
             if(cnt == 1) cvfnc = sdrcv(iduv)
             if(cnt == 2) cvfnc = inscv(iduv)
          case('uvcrd')
             cvfnc = uvcrd(iduv)
          case default
             errtag = .true.
          end select
          if(errtag) stop ' Bug in the program'
          factor = factor + cvfnc * weight(iduv)
       endif
    end do
    cvfcen = factor
    deallocate( weight )
    return
  end function cvfcen

  subroutine getwght(weight, pti, cnt, systype, wgttype, engtype)
    implicit none
    real, intent(out) :: weight(gemax)
    ! pti : identifier of solvent species
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: pti, cnt
    character(len=5), intent(in) :: systype
    character(len=4), intent(in) :: wgttype
    character(len=3), intent(in) :: engtype
    integer :: iduv
    real :: minuv, ampl
    weight(:) = 0.0
    do iduv = 1, gemax
       if(uvspec(iduv) == pti) then
          weight(iduv) = wgtdst(iduv, cnt, systype, wgttype)
       endif
    end do
    if(engtype == 'yes') then
       minuv = minval( abs(uvcrd(:)), &
            mask = ((uvspec == pti) .and. (weight > zero)) )
       do iduv = 1, gemax
          if(uvspec(iduv) == pti) then
             !  ampl = abs(uvcrd(iduv)) - minuv
             ampl = nummol(pti) * (abs(uvcrd(iduv)) - minuv)
             weight(iduv) = exp(- ampl / kT) * weight(iduv)
          endif
       end do
    endif
    ampl = sum( weight, mask = (uvspec == pti) )
    if(ampl > zero) then
       weight(:) = weight(:) / ampl
    else
       write(6, *) ' Zero weight at ', pti
       stop
    endif
  end subroutine getwght

  integer function zeroec(pti, cnt)
    implicit none
    ! pti : identifier of solvent species
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: pti, cnt
    integer :: iduv, k
    real :: factor, ampl, lcsln, lcref
    do iduv = 1, gemax - 1
       if(uvspec(iduv) == pti) then
          if((uvcrd(iduv) <= 0.0) .and. (uvcrd(iduv + 1) >= 0.0)) then
             factor = abs(uvcrd(iduv))
             ampl = uvcrd(iduv + 1)
             if(ampl > factor + tiny) then
                k = iduv
             elseif(ampl < factor - tiny) then
                k = iduv + 1
             elseif(abs(ampl - factor) <= tiny) then
                if(cnt == 1) then
                   lcsln = edist(iduv)
                   lcref = edist(iduv + 1)
                endif
                if(cnt == 2) then
                   lcsln = edens(iduv)
                   lcref = edens(iduv + 1)
                endif
                if(lcsln >= lcref) then
                   k = iduv
                else
                   k = iduv + 1
                endif
             else
                stop 'Bug in the program'
             endif
          endif
       endif
    end do
    zeroec = k
    return
  end function zeroec


  real function wgtdst(iduv, cnt, systype, wgttype)
    use sysvars, only: wgtf2smpl
    implicit none
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: iduv, cnt
    character(len=5), intent(in) :: systype
    character(len=4), intent(in) :: wgttype
    real :: fsln, fref, wght, factor
    logical :: errtag
    integer :: jdg
    fsln = edist(iduv)
    fref = edens(iduv)
    errtag = .false.
    jdg = 1
    if(cnt == 2) then
       if(systype == 'slncv')  jdg = 9
       if(systype == 'extsl')  jdg = 9
       if(wgttype == 'smpl')   jdg = 2
       if(wgtf2smpl == 'yes')  jdg = 2
    endif
    if((jdg /= 1) .and. (jdg /= 2)) errtag = .true.
    if(jdg == 2) wght = fref
    if(jdg == 1) then
       wght = 0.0
       select case(wgttype)
       case('smpl')
          wght = fsln
          if(cnt == 2) errtag = .true.
       case('geom')
          factor = fsln * fref
          if(factor > zero) wght = sqrt(factor)
       case default              ! corresponding to wgttype = 'harm'
          factor = fsln + fref
          if(factor > zero) wght = fsln * fref / factor
       end select
    endif
    if(errtag) stop ' Bug in the program'
    wgtdst = wght
    return
  end function wgtdst


  real function sfewgt(fsln, fref)
    implicit none
    real, intent(in) :: fsln, fref
    real :: wght, factor
    if(fsln >= fref) then
       wght = 1.0
    else
       factor = (fsln - fref) / (fsln + fref)
       wght = 1.0 - factor ** 2
    endif
    sfewgt = wght
    return
  end function sfewgt


  real function pyhnc(indpmf, cnt)
    implicit none
    real, intent(in) :: indpmf
    integer, intent(in) :: cnt
    real :: intg, factor
    factor = indpmf / kT
    select case(cnt)
    case(1, 2)      ! usual case   1 : solution  2 : referecen solvent
       if(factor < - zero) then    ! PY
          if(cnt == 1) intg = factor + factor / (exp(- factor) - 1.0)
          if(cnt == 2) intg = log(1.0 - factor) * (1.0 / factor - 1.0)
          intg = intg + 1.0
       else                        ! HNC
          intg = factor / 2.0
       endif
    case(3)         ! special case to be examined carefully
       if(factor > zero) then
          intg = 1.0 - log(1.0 + factor) * (1.0 / factor + 1.0)
       else
          intg = - factor / 2.0
       endif
    case default
       stop "Incorrct cnt argument in pyhnc"
    end select
    pyhnc = intg
    return
  end function pyhnc

  subroutine distnorm
    implicit none
    integer :: iduv, iduvp, pti, cnt, itrcnt
    real :: factor, ampl, lcsln, lcref, errtmp
    real, dimension(:), allocatable :: correc, edhst
    real, dimension(:,:), allocatable :: edmcr
    allocate( correc(gemax), edhst(gemax), edmcr(gemax, gemax) )
    do cnt = 1, 2     ! cnt = 1: solution   cnt = 2: reference solvent
       if(cnt == 1) then
          edhst(:) = edist(:)
          if(slncor == 'yes') edmcr(:,:) = edscr(:,:)
       endif
       if(cnt == 2) then
          edhst(:) = edens(:)
          edmcr(:,:) = ecorr(:,:)
       endif

       do pti = 1, numslv
          factor = sum( edhst, mask = (uvspec == pti) )
          if(factor > zero) then
             factor = nummol(pti) / factor
          else
             factor = 0.0
          endif
          do iduv = 1, gemax
             if(uvspec(iduv) == pti) then
                edhst(iduv) = factor * edhst(iduv)
             endif
          end do
       end do

       if((cnt == 1) .and. (slncor /= 'yes')) goto 5555

       errtmp = error + 1.0
       itrcnt = 0
       correc(:) = 1.0
       do while((errtmp > error) .and. (itrcnt <= itrmax))
          do iduv = 1, gemax
             lcsln = 0.0
             do pti = 1, numslv
                ampl = sum( correc(:) * edmcr(:, iduv), mask = (uvspec == pti) )
                if(ampl > zero) lcsln = lcsln + nummol(pti) / ampl
             end do
             lcsln = lcsln / real(numslv)
             correc(iduv) = lcsln * edhst(iduv)
          end do
          do iduv = 1, gemax
             do iduvp = 1, gemax
                ampl = correc(iduv) * correc(iduvp)
                edmcr(iduvp,iduv) = ampl * edmcr(iduvp,iduv)
             end do
          end do
          errtmp = maxval( abs(correc(:) - 1.0), mask = (edhst(:) > zero) )
          itrcnt = itrcnt + 1
          if(itrcnt >= itrmax) then
             write(6, *) ' The optimization of the correlation matrix'
             write(6, *) '  did not converge with an error of ', errtmp
             stop
          endif
       end do

5555   continue
       if(cnt == 1) then
          edist(:) = edhst(:)
          if(slncor == 'yes') edscr(:,:) = edmcr(:,:)
       endif
       if(cnt == 2) then
          edens(:) = edhst(:)
          ecorr(:,:) = edmcr(:,:)
       endif
    end do

    deallocate( correc, edhst, edmcr )

    return
  end subroutine distnorm

  subroutine distshow
    implicit none
    integer :: iduv, pti, cnt, ecmin, ecmax, k, ilist(gemax)
    real :: factor, ratio
    real, dimension(:), allocatable :: edhst

    allocate( edhst(gemax) )
    do iduv = 1, gemax
       ilist(iduv) = iduv
    end do
    write(6, *)

    do pti = 1, numslv
       do cnt = 1, 2     ! cnt = 1: solution   cnt = 2: reference solvent
          if(cnt == 1) then
             edhst(:) = edist(:)
             if(numslv == 1) then
                write(6, "(A)") " SOLUTION"
             else
                write(6, "(A,i4,A)") " SOLUTION for", pti, "-th species"
             endif
          endif
          if(cnt == 2) then
             edhst(:) = edens(:)
             if(numslv == 1) then
                write(6, "(A)") " INSERTION"
             else
                write(6, "(A,i4,A)") " INSERTION for", pti, "-th species"
             endif
          endif

          ecmin = minval( ilist, &
               mask = ((uvspec == pti) .and. (edhst > zero)) )
          ecmax = maxval( ilist, &
               mask = ((uvspec == pti) .and. (edhst > zero)) )
          k = count( mask = (edhst(ecmin:ecmax) > zero) )
          ratio = real(k) / real(ecmax - ecmin + 1)
          factor = sum( edhst(ecmin:ecmax), mask = (edhst(ecmin:ecmax) > zero) )

          do iduv = ecmin, ecmax
             if((cnt == 2) .and. (edhst(iduv) <= zero)) then
                write(6, "(A,i5,A,g14.6)") "     No sampling at ", iduv, &
                     " with energy ", uvcrd(iduv)
             endif
          end do

          if(cnt == 1) write(6, "(A,g12.4)") "     Nonzero component ratio in solution  = ", ratio
          if(cnt == 2) write(6, "(A,g12.4)") "     Nonzero component ratio at insertion = ", ratio
          write(6, "(A,g12.4)") "          Number of interacting molecules = ", factor

          if(cnt == 1) then
             write(6, "(A,g15.7)") "     Minimum energy in solution   =", uvcrd(ecmin)
             write(6, "(A,g15.7)") "     Maximum energy in solution   =", uvcrd(ecmax)
          endif
          if(cnt == 2) then
             write(6, "(A,g15.7)") "     Minimum energy at insertion  =", uvcrd(ecmin)
             write(6, "(A,g15.7)") "     Maximum energy at insertion  =", uvcrd(ecmax)
          endif
       end do
    end do

    deallocate( edhst )

    return
  end subroutine distshow
end module sfecalc

module opwrite
  use sysvars, only: clcond, uvread, slfslt, infchk, &
       prmmax, numrun, numslv, zero, &
       pickgr, write_mesherror, msemin, msemax, mesherr, &
       slfeng, chmpt, aveuv, blockuv, svgrp, svinf
  implicit none
  integer :: grref
  real :: fe_stat_error     ! 95% error of the solvation free energy
  real, dimension(:), allocatable :: mshdif
contains

  subroutine wrtresl
    implicit none
    integer :: prmcnt, pti, k, group, inft
    real :: totuv, differ, mesh_error, valcp
    !
    if(slfslt == 'yes') write(6, "(A,f12.4,A)") "  Self-energy of the solute   =   ", slfeng, "  kcal/mol"

    if((clcond == 'basic') .or. (clcond == 'range')) then
       write(6,*)
       if((numslv > 1) .and. (uvread == 'not')) write(6, 331) aveuv(1:numslv)
331    format('  Solute-solvent energy       =   ', 9999f12.4)
       totuv = sum( aveuv(1:numslv) )
       if(slfslt == 'yes') totuv = totuv + slfeng
       write(6, "(A,f12.4,A)") "  Total solvation energy      =   ", totuv, "  kcal/mol"
    endif

    if(clcond == 'basic') then
       if(numslv > 1) write(6, 351) chmpt(1:numslv, 1, 1)
351    format('  Solvation free energy       =   ', 9999f12.4)
       write(6, "(A,f12.4,A)") "  Total solvation free energy =   ", chmpt(0,1,1), "  kcal/mol"
    endif

    if((clcond == 'range') .or. (clcond == 'merge')) then
       do prmcnt = 1, prmmax
          group = svgrp(prmcnt)
          if(group == pickgr) then
             grref = prmcnt
             exit
          endif
       end do
       allocate( mshdif(msemin:msemax) )
       mshdif(:) = -1.0
    endif

    if(clcond == 'range') then
       if(numslv == 1) then
          k = 0
       else
          k = numslv
       endif
       do pti = 0, k
          write(6, *)
          if(infchk == 'yes') then
             if(pti == 0) then
                write(6, "(A)") " group  inft  solvation free energy difference"
             else
                write(6, "(A,i3,A)") "               ", pti, "-th component difference"
             endif
          else
             if(pti == 0) then
                write(6, "(A)") " group    solvation free energy   difference"
             else
                write(6, "(A,i3,A)") "           ", pti, "-th component       difference"
             endif
          endif
          do prmcnt = 1, prmmax
             group = svgrp(prmcnt)
             inft = svinf(prmcnt)
             valcp = chmpt(pti, prmcnt, 1)
             differ = valcp - chmpt(pti, grref, 1)
             if(infchk == 'yes') then
                write(6, "(i4,i7,f17.5,f18.5)") group, inft, valcp, differ
             else
                write(6, "(i4,f20.5,f18.5)") group, valcp, differ
             endif
             if((pti == 0) .and. (inft == 0) .and. &
                  (group >= msemin) .and. (group <= msemax)) then
                mshdif(group) = abs(differ)
             endif
          end do
       end do
    endif

    if(clcond == 'merge') call wrtmerge

    if((clcond == 'range') .or. (clcond == 'merge')) then
       mesh_error = maxval( mshdif(msemin:msemax), &
            mask = (mshdif(msemin:msemax) > zero) )
       ! mesh error not written in any case
       if(write_mesherror == 'not') then
          ! do nothing
          ! mesh error written in any case
       elseif(write_mesherror == 'yes') then
          write(6, *)
          write(6, "(A,f8.3,A)") " Mesh error is ", mesh_error, " kcal/mol"
          ! write_mesherror is other than 'not' or 'yes'
          ! no output for mesh error when mesh_error < mesherr
       elseif(mesh_error < mesherr) then
          ! do nothing
          ! clcond = 'merge', mesh_error > fe_stat_error
       elseif((clcond == 'merge') .and. (mesh_error > fe_stat_error)) then
          write(6, *)
          write(6, "(A,f8.3,A,f8.3,A)") &
               " Warning: mesh error is ", mesh_error, &
               " kcal/mol and is larger than the 95% error of ", &
               fe_stat_error, " kcal/mol"
          ! clcond = 'range' or mesh_error < fe_stat_error,
          ! mesherr is set to be positive, mesh_error > mesherr
       elseif((mesherr > zero) .and. (mesh_error > mesherr)) then
          write(6, *)
          write(6, "(A,f8.3,A,f8.3,A)") &
               " Warning: mesh error is ", mesh_error, &
               " kcal/mol and is larger than the threshold value of ", &
               mesherr, " kcal/mol"
       endif
       deallocate( mshdif )
    endif

    return
  end subroutine wrtresl

  subroutine wrtmerge
    implicit none
    integer :: prmcnt, cntrun, group, inft, pti, i, j, k, m
    real :: avecp, stdcp, avcp0, recnt, slvfe
    real, dimension(:),   allocatable :: showcp
    real, dimension(:,:), allocatable :: wrtdata

    allocate( showcp(numrun), wrtdata(0:numslv, numrun) )
    if(uvread /= 'not') then
       wrtdata(0:numslv, 1:numrun) = blockuv(0:numslv, 1:numrun)
       write(6, *)
       write(6, *)
       write(6, "(A)") " cumulative average & 95% error for solvation energy"
       call wrtcumu(wrtdata)
       write(6, *)
    endif

    do pti = 0, numslv
       if((numslv == 1) .and. (pti /= 0)) cycle
       avcp0 = sum( chmpt(pti, grref, 1:numrun) ) / real(numrun)
       do prmcnt = 1, prmmax
          group = svgrp(prmcnt)
          inft = svinf(prmcnt)
          avecp = 0.0
          stdcp = 0.0
          do cntrun = 1, numrun
             slvfe = chmpt(pti, prmcnt, cntrun)
             avecp = avecp + slvfe
             stdcp = stdcp + slvfe ** 2
          end do
          recnt = real(numrun)
          avecp = avecp / recnt
          if(numrun > 1) then
             stdcp = sqrt(recnt / (recnt - 1.0)) &
                  * sqrt(stdcp / recnt - avecp ** 2)    ! standard deviation
             stdcp = 2.0 * stdcp / sqrt(recnt)           ! 95% error
          endif
          if(prmcnt == 1) then
             write(6, *)
             if(pti == 0) then
                if(infchk == 'yes') then
                   if(numrun == 1) then
                      write(6, "(A)") " group  inft  solvation free energy     difference"
                   else
                      write(6, "(A)") " group  inft  solvation free energy     error          difference"
                   endif
                else
                   if(numrun == 1) then
                      write(6, "(A)") " group    solvation free energy     difference"
                   else
                      write(6, "(A)") " group    solvation free energy     error          difference"
                   endif
                endif
             endif
             if(numslv > 1) then
                if(pti == 0) then
                   write(6, "(A)") "  total solvation free energy"
                else
                   write(6, "(A,i2,A)") "  contribution from ", pti, "-th solvent component"
                endif
             endif
          endif
          if(infchk == 'yes') then
             if(numrun == 1) then
                write(6, "(i4,i7,f17.5,f21.5)") group, inft, avecp, (avecp - avcp0)
             else
                write(6, "(i4,i7,f17.5,2f18.5)") group, inft, avecp, stdcp, (avecp - avcp0)
             endif
          else
             if(numrun == 1) then
                write(6, "(i4,f20.5,f21.5)") group, avecp, (avecp - avcp0)
             else
                write(6, "(i4,f20.5,2f18.5)") group, avecp, stdcp, (avecp - avcp0)
             endif
          endif

          if((pti == 0) .and. (inft == 0) .and. &
               (group >= msemin) .and. (group <= msemax)) then
             mshdif(group) = abs(avecp - avcp0)
          endif
       end do
    end do

    if(numrun > 1) then
       write(6, *)
       write(6, *)
       do pti = 0, numslv
          if((numslv == 1) .and. (pti /= 0)) cycle
          do prmcnt = 1, prmmax
             group = svgrp(prmcnt)
             inft = svinf(prmcnt)
             showcp(1:numrun) = chmpt(pti, prmcnt, 1:numrun)
             if(prmcnt == 1) then
                if(infchk == 'yes') then
                   if(numslv == 1) then
                      write(6, "(A)") " group  inft   Estimated free energy (kcal/mol)"
                   else
                      if(pti == 0) then
                         write(6, "(A)") " group  inft   Estimated free energy: total (kcal/mol)"
                      else
                         write(6, *)
                         write(6, "(A,i2,A)") " group  inft   Estimated free energy:", pti, "-th solvent contribution (kcal/mol)"
                      endif
                   endif
                else
                   if(numslv == 1) then
                      write(6, "(A)") " group   Estimated free energy (kcal/mol)"
                   else
                      if(pti == 0) then
                         write(6, "(A)") " group   Estimated free energy: total (kcal/mol)"
                      else
                         write(6, *)
                         write(6, "(A,i2,A)") " group   Estimated free energy:", pti, "-th solvent contribution (kcal/mol)"
                      endif
                   endif
                endif
             endif

             k = (numrun - 1) / 5
             if(infchk == 'yes') then
                if(k == 0) then
                   write(6, "(i4,i7,5f13.4)") group, inft, showcp(1:numrun)
                else
                   write(6, "(i4,i7,5f13.4)") group, inft, showcp(1:5)
                   if(k > 1) then
                      do i = 1, k - 1
                         write(6, "('           ',5f13.4)") &
                              (showcp(5 * i + m), m = 1, 5)
                      end do
                   endif
                   if(5 * k < numrun) then
                      write(6, "('           ',5f13.4)") &
                           (showcp(m), m = 5 * k + 1, numrun)
                   endif
                endif
             else
                if(k == 0) then
                   write(6, "(i4,'  ',5f13.4)") group, showcp(1:numrun)
                else
                   write(6, "(i4,'  ',5f13.4)") group, showcp(1:5)
                   if(k > 1) then
                      do i = 1, k - 1
                         write(6, "('      ',5f13.4)") &
                              (showcp(5 * i + m), m = 1, 5)
                      end do
                   endif
                   if(5 * k < numrun) then
                      write(6, "('      ',5f13.4)") &
                           (showcp(m), m = 5 * k + 1, numrun)
                   endif
                endif
             endif
          end do
       end do
    endif

    wrtdata(0:numslv, 1:numrun) = chmpt(0:numslv, grref, 1:numrun)
    write(6, *)
    write(6, *)
    write(6, "(A)") " cumulative average & 95% error for solvation free energy"
    call wrtcumu(wrtdata, fe_stat_error)
    deallocate( showcp, wrtdata )

    return
  end subroutine wrtmerge

  subroutine wrtcumu(wrtdata, stat_error)
    implicit none
    real, intent(in) :: wrtdata(0:numslv, numrun)
    real, intent(out), optional :: stat_error
    integer :: cntrun, pti
    real :: avecp, factor, slvfe, recnt
    real, dimension(:), allocatable :: runcp, runer, wrtcp

    allocate( runcp(0:numslv), runer(0:numslv), wrtcp(2 * numslv + 2) )
    runcp(:) = 0.0
    runer(:) = 0.0

    if(numslv >= 2) then               ! mixed-solvent system
       write(6, "(A)", advance='no') "              total             1st component         2nd component"
       if(numslv >= 3) write(6, "(A)", advance='no') "         3rd component"
       if(numslv >= 4) then
          do pti = 4, numslv
             select case(pti / 10)
             case(0)                   ! 4 <= pti < 10
                write(6, "(A,i1,A)", advance='no') "         ", pti, "th component"
             case(1)                   ! 10 <= pti < 100
                write(6, "(A,i2,A)", advance='no') "        ", pti, "th component"
             case default              ! 100 <= pti
                write(6, "(A,i4,A)", advance='no') "        ", pti, "th component"
             end select
          enddo
       endif
       write(6, '(A)') ""
    endif

    do cntrun = 1, numrun
       recnt = real(cntrun)
       do pti = 0, numslv
          slvfe = wrtdata(pti, cntrun)
          runcp(pti) = runcp(pti) + slvfe
          runer(pti) = runer(pti) + slvfe ** 2
          avecp = runcp(pti) / recnt
          wrtcp(2 * pti + 1) = avecp
          if(cntrun > 1) then
             factor = runer(pti) / recnt - avecp ** 2
             if(factor <= zero) then
                wrtcp(2 * pti + 2) = 0.0
             else
                wrtcp(2 * pti + 2) = (2.0 / sqrt(recnt) ) &
                     *sqrt(recnt / (recnt - 1.0)) * sqrt(factor)
             endif
          endif
       end do
       if(cntrun == 1) then
          do pti = 0, numslv
             wrtcp(pti + 1) = wrtcp(2 * pti + 1)
          end do
          if(numslv == 1) then
             write(6, "(i3,f11.4)") cntrun, wrtcp(1)
          else
             write(6, 771) cntrun, wrtcp(1), wrtcp(2:numslv+1)
771          format(i3,f11.4,9999f22.4)
          endif
       else
          if(numslv == 1) then
             write(6, "(i3,2f11.4)") cntrun, wrtcp(1:2)
          else
             write(6, 772) cntrun, wrtcp(1:2*numslv+2)
772          format(i3,9999f11.4)
          endif
          ! 95% error of the solvation free energy
          ! used only when the mesh error for free energy is assessed
          if(present(stat_error) .and. (cntrun == numrun)) stat_error = wrtcp(2)
       endif
    end do

    deallocate( runcp, runer, wrtcp )

    return
  end subroutine wrtcumu
end module opwrite


program sfemain
  use sysvars, only: numrun, prmmax, init_sysvars
  use sysread, only: defcond, datread
  use sfecalc, only: chmpot
  use opwrite, only: wrtresl
  implicit none
  integer :: cntrun, prmcnt
  call init_sysvars
  call defcond
  do cntrun = 1, numrun
     call datread(cntrun)
     do prmcnt = 1, prmmax
        call chmpot(prmcnt, cntrun)
     end do
  end do
  call wrtresl
  ! stop  ! commented out to suppress harmless IEEE errors on Cygwin and MacOS X
end program sfemain
