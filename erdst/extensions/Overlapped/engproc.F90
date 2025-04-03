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

module engproc
   implicit none
   character(len=10), parameter :: numbers='0123456789'
   real, parameter :: tiny = 1.0e-30
   integer :: cntdst, slvmax
   integer :: maxdst
   integer :: tagslt
   integer, allocatable :: tagpt(:)

   ! flceng needs to be output in-order
   logical, allocatable :: flceng_stored(:)
   real, allocatable :: flceng(:, :)

   ! start of the extension for the method with the overlapped state
   integer :: do_overlapped
   character(len=1024) :: refs_filename = '../refs/engref'
   logical, dimension(:), allocatable :: zerodns_crd
   real, dimension(:),    allocatable :: edovlp, uvovlp
   real, dimension(:,:),  allocatable :: avovlp
   real :: soln_total, soln_ovlp, norm_ovlp
   ! end of the extension for the method with the overlapped state

contains
   !
   !  procedure for constructing energy distribution functions
   !
   subroutine enginit
      use engmain, only: numtype, nummol, engdiv, corrcal, selfcal, slttype, &
         moltype, sluvid, &
         ermax, numslv, uvmax, uvsoft, esmax, uvspec, &
         uvcrd, edens, ecorr, escrd, eself, &
         voffset, &
         aveuv, slnuv, avediv, minuv, maxuv, numslt, sltlist, &
         ene_confname, io_paramfile, io_flcuv, &
         ecdinfo_file, ecdinfo_io, ecdmesh_file, ecdmesh_io, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, PT_SOLVENT, NO, YES, &
         ermax_limit
      use mpiproc, only: halt_with_error, warning, myrank
      implicit none
      real(kind=8) :: ecdmin, ecfmns, ecmns0, ecdcen, ecpls0, ecfpls, eccore, ecdmax
      real(kind=8) :: eclbin, ecfbin, ec0bin, finfac, ectmvl
      integer :: pecore
      ! pemax :  number of discretization of the solute-solvent energy
      ! pesoft : number of discretization in the soft interaction region
      !    they are constructed from other parameters (pemax = pesoft + pecore)
      integer pemax, pesoft
      integer :: ecprread, meshread, peread
      !
      real, parameter :: infty = huge(infty)    ! essentially equal to infinity
      integer, parameter :: rglmax = 5, large = 10000
      real(kind=8) :: factor, incre, cdrgvl(0:rglmax+1), ecpmrd(large)
      integer :: solute_moltype
      integer :: iduv, i, q, pti, regn, minrg, maxrg, uprgcd(0:rglmax+1), dummy
      integer, dimension(:), allocatable :: tplst
      real(kind=8), dimension(:,:), allocatable  :: ercrd
      !
      integer :: param_err
      logical :: check_ok, start_line
      ! start of the extension for the method with the overlapped state
      logical :: edens_exist
      logical, dimension(:), allocatable :: exist_total
      character(len=1024) :: opnfile
      integer, parameter :: sufmax = 99, refs_io = 81
      character(len=3), parameter :: suffix_tt = '.tt'
      character(len=3) :: count_suf
      namelist /overlapped/ do_overlapped, refs_filename
      ! end of the extension for the method with the overlapped state
      namelist /hist/ ecdmin, ecfmns, ecdcen, eccore, ecdmax, &
         eclbin, ecfbin, ec0bin, finfac, pecore, &
         ecprread, meshread, peread
      !
      allocate( tplst(nummol) )
      solute_moltype = 0 ! suppress warning
      numslt = 0
      do i = 1, nummol
         if(sluvid(i) > 0) then
            numslt = numslt + 1
            tplst(numslt) = i
            solute_moltype = moltype(i)
         endif
      enddo
      ! solute must have moltype value equal to solute_moltype
      if(any(sluvid(:) /= PT_SOLVENT .and. moltype(:) /= solute_moltype)) call halt_with_error('eng_typ')
      ! solvent must have moltype value not equal to solute_moltype
      if(any(sluvid(:) == PT_SOLVENT .and. moltype(:) == solute_moltype)) call halt_with_error('eng_typ')
      !
      ! consistency check between slttype and numslt (number of solute molecules)
      check_ok = .true.
      if(numslt <= 0) check_ok = .false.
      if((slttype == SLT_REFS_RIGID) .or. (slttype == SLT_REFS_FLEX)) then
         if(numslt /= 1) check_ok = .false.
      endif
      if(.not. check_ok) call halt_with_error('eng_num')
      !
      allocate( sltlist(numslt) )
      sltlist(1:numslt) = tplst(1:numslt)   ! list of solute molecules
      deallocate( tplst )
      !
      ! solute needs to be the last particle in reference system
      if((slttype == SLT_REFS_RIGID) .or. (slttype == SLT_REFS_FLEX)) then
         if(sltlist(1) /= nummol) call halt_with_error('eng_ins')
      endif
      !
      ! number of solvent species
      if(numslt == 1) then
         numslv = numtype - 1
      elseif(numslt > 1) then     ! solute can also be a solvent species
         numslv = numtype
      else
         call halt_with_error('eng_bug')
      endif
      !
      allocate( uvspec(nummol) )
      uvspec(1:nummol) = moltype(1:nummol)
      if(numslt == 1) then        ! solute is not present in reference solvent
         where(sluvid(:) /= PT_SOLVENT) uvspec(:) = 0 ! solute
         ! after the solute, slide the value of molecule type
         where(sluvid(:) == PT_SOLVENT .and. moltype(:) > solute_moltype) uvspec(:) = uvspec(:) - 1
      endif
      if(numslv /= maxval(uvspec(:))) call halt_with_error('eng_bug')
      !
      allocate( uvmax(numslv), uvsoft(numslv), ercrd(large, 0:numslv) )
      !
      peread = NO                           ! deprecated
      ecprread = NO
      meshread = NO
      ermax = 0
      do pti = 0, numslv
         open(unit = io_paramfile, file = ene_confname, action = "read", iostat = param_err)
         if(param_err == 0) then
            read(io_paramfile, nml = hist)
            close(io_paramfile)
         else
            call halt_with_error('eng_per')
         end if
         if(peread == YES) ecprread = YES   ! deprecated

         ! read coordinate parameters from EcdInfo
         if(ecprread == YES) then
            open(unit = ecdinfo_io, file = ecdinfo_file, action = 'read', iostat = param_err)
            if(param_err /= 0) call halt_with_error('eng_eci')
            read(ecdinfo_io, *)             ! comment line
            do i = 1, large
               read(ecdinfo_io, *, end = 3109) q
               if(q == pti) then
                  backspace(ecdinfo_io)
                  if(pti == 0) then         ! solute self-energy
                     read(ecdinfo_io, *) iduv, ecpmrd(1:8)
                     pecore = 0             ! no core region for self-energy
                  else                      ! solute-solvent interaction energy
                     read(ecdinfo_io, *) iduv, ecpmrd(1:9), pecore
                     ecdmax = ecpmrd(9)
                     if(pecore <= 1) call halt_with_error('eng_pcr')
                  endif
                  eclbin = ecpmrd(1)
                  ecfbin = ecpmrd(2)
                  ec0bin = ecpmrd(3)
                  finfac = ecpmrd(4)
                  ecdmin = ecpmrd(5)
                  ecfmns = ecpmrd(6)
                  ecdcen = ecpmrd(7)
                  eccore = ecpmrd(8)
                  if(eccore < tiny) pecore=0
                  if(pecore == 1) call halt_with_error('eng_pcr')
                  exit
               endif
            enddo
3109        continue
            close(ecdinfo_io)
         end if

         ectmvl = finfac * ecfbin
         ecdmin = ecdmin - ectmvl
         ecfmns = ecfmns - ectmvl
         ecmns0 = ecdcen - ectmvl
         ecpls0 = 2.0 * ecdcen - ecmns0
         ecfpls = 2.0 * ecdcen - ecfmns
         eccore = eccore + ecfpls - ecdcen

         cdrgvl(0) = ecdmin
         cdrgvl(1) = ecfmns
         cdrgvl(2) = ecmns0
         cdrgvl(3) = ecpls0
         cdrgvl(4) = ecfpls
         cdrgvl(5) = eccore
         cdrgvl(6) = ecdmax

         uprgcd(0) = 1
         do regn = 1, rglmax
            if((regn == 1) .or. (regn == 5)) factor = eclbin
            if((regn == 2) .or. (regn == 4)) factor = ecfbin
            if(regn == 3) factor = ec0bin
            incre = cdrgvl(regn) - cdrgvl(regn - 1)
            if(incre < factor) call halt_with_error('eng_ecd')
            iduv = nint(incre / factor)
            uprgcd(regn) = uprgcd(regn - 1) + iduv
         enddo

         pesoft = uprgcd(rglmax) - uprgcd(0)
         pemax = pesoft + pecore
         uprgcd(rglmax + 1) = pemax
         ercrd(uprgcd(0:(rglmax + 1)), pti) = cdrgvl(0:(rglmax + 1))

         if(pemax > large) call halt_with_error('eng_siz')

         if(pecore == 0) i = rglmax         ! no core region
         if(pecore >  0) i = rglmax + 1     ! explicit treatment of core region
         do regn = 1, i
            minrg = uprgcd(regn - 1)
            maxrg = uprgcd(regn)
            if(regn <= rglmax) then
               incre = ercrd(maxrg, pti) - ercrd(minrg, pti)
            endif
            if(regn == (rglmax + 1)) then   ! effective only when pecore > 0
               incre = log(ercrd(maxrg, pti) / ercrd(minrg, pti))
            endif
            factor = incre / real(maxrg - minrg)
            do iduv = minrg, (maxrg - 1)
               incre = factor * real(iduv - minrg)
               if(regn <= rglmax) then
                  ercrd(iduv, pti) = ercrd(minrg, pti) + incre
               endif
               if(regn == (rglmax + 1)) then
                  ercrd(iduv, pti) = ercrd(minrg, pti) * exp(incre)
               endif
            enddo
         enddo

         ! read coordinate meshes from EcdMesh
         if(meshread == YES) then
            open(unit = ecdmesh_io, file = ecdmesh_file, action = 'read', iostat = param_err)
            if(param_err /= 0) call halt_with_error('eng_ecm')
            start_line = .true.
            do i = 1, large
               read(ecdmesh_io, *, end = 3209) q
               if(q == pti) then
                  backspace(ecdmesh_io)
                  if( start_line ) then     ! 0-th line of the q-th species
                     start_line = .false.
                     iduv = 0
                     if(pti == 0) then      ! solute self-energy
                        read(ecdmesh_io, *) dummy, pemax
                        pecore = 0          ! no core region for self-energy
                     else                   ! solute-solvent interaction energy
                        read(ecdmesh_io, *) dummy, pemax, pecore
                     endif
                     pesoft = pemax - pecore
                  else
                     iduv = iduv + 1        ! iduv-th line of the q-th species
                     read(ecdmesh_io, *) dummy, dummy, ercrd(iduv, pti)
                  endif
               endif
            enddo
3209        continue
            close(ecdmesh_io)
            check_ok = .true.
            if(iduv /= pemax) check_ok = .false.
            if(pemax <= pecore) check_ok = .false.
            do iduv = 1, pemax - 1
               if(ercrd(iduv, pti) >= ercrd(iduv + 1, pti)) check_ok = .false.
            enddo
            if(.not. check_ok) call halt_with_error('eng_emf')
         endif

         if(pti == 0) then                  ! solute self-energy
            esmax = pesoft
         else                               ! solute-solvent interaction energy
            uvmax(pti) = pemax
            uvsoft(pti) = pesoft
            ermax = ermax + pemax
         endif
      enddo

      allocate( uvcrd(ermax), edens(ermax) )
      i = 0
      do pti = 1, numslv
         pemax = uvmax(pti)
         do iduv = 1, pemax
            i = i + 1
            uvcrd(i) = ercrd(iduv, pti)
         enddo
      enddo

      if(corrcal == YES) then
         if(ermax > ermax_limit) call warning('emax')
         allocate( ecorr(ermax, ermax) )
      endif

      if(selfcal == YES)  then
         allocate( escrd(esmax), eself(esmax) )
         escrd(1:esmax) = ercrd(1:esmax,0)
      endif
      deallocate( ercrd )

      if(slttype == SLT_SOLN) allocate( aveuv(engdiv, numslv), slnuv(numslv) )
      allocate( avediv(engdiv, 2) )
      allocate( minuv(0:numslv), maxuv(0:numslv) )

      do pti = 0, numslv
         minuv(pti) = infty
         maxuv(pti) = -infty
      enddo
      voffset = -infty

      ! start of the extension for the method with the overlapped state
      do_overlapped = NO                 ! default = don't do it
      if(slttype == SLT_SOLN) then
         open(unit = io_paramfile, file = ene_confname, action = "read", iostat = param_err)
         if(param_err == 0) then
            read(io_paramfile, nml = overlapped, iostat = param_err)
            close(io_paramfile)
         endif
         if(do_overlapped == YES) then
            allocate( zerodns_crd(ermax), edovlp(ermax) )
            allocate( avovlp(engdiv, numslv), uvovlp(numslv) )
            zerodns_crd(:) = .true.
            soln_total = 0.0
            soln_ovlp = 0.0
            opnfile = trim(refs_filename) // suffix_tt
            inquire(file = opnfile, exist = edens_exist)
            if( edens_exist ) then
               regn = 1
            else
               regn = sufmax
            endif
            allocate( exist_total(regn) )
            exist_total(:) = .false.
            do pti = 1, regn
               if(regn == sufmax) then
                  i = pti / 10
                  q = mod(pti, 10)
                  count_suf = '.' // numbers(i+1:i+1) // numbers(q+1:q+1)
                  opnfile = trim(refs_filename) // count_suf
                  inquire(file = opnfile, exist = edens_exist)
               endif
               if( edens_exist ) then
                  open(unit = refs_io, file = opnfile, action = "read")
                  read(refs_io, *)    ! skip the comment line
                  do iduv = 1, ermax
                     read(refs_io, *) factor, factor, i, factor
                     if(factor > tiny) zerodns_crd(iduv) = .false.
                  enddo
                  close(refs_io)
                  exist_total(pti) = .true.
               endif
            enddo
            if( all(exist_total(:) == .false.) ) stop 'engref files did not exist'
            deallocate( exist_total )
         endif
      endif
      ! end of the extension for the method with the overlapped state

      call engclear

      ! Output for energy fluctuation
      if(myrank == 0) then
         if(slttype == SLT_SOLN) then
            ! open flcuv file
            open(unit = io_flcuv, file = 'flcuv.tt', status = 'new', action = 'write')
         else
            ! open progress file
            open(unit = io_flcuv, file = 'progress.tt', status = 'new', action = 'write')
         endif
      endif

      return
   end subroutine enginit

   subroutine engclear
      use engmain, only: corrcal, selfcal, slttype, SLT_SOLN, YES, &
         edens, ecorr, eself, slnuv, avslf, engnorm, engsmpl
      implicit none
      edens(:) = 0.0
      if(corrcal == YES) ecorr(:,:) = 0.0
      if(selfcal == YES) eself(:) = 0.0
      if(slttype == SLT_SOLN) slnuv(:) = 0.0
      avslf = 0.0
      engnorm = 0.0
      engsmpl = 0.0
      ! start of the extension for the method with the overlapped state
      if(do_overlapped == YES) then
         edovlp(:) = 0.0
         uvovlp(:) = 0.0
         norm_ovlp = 0.0
      endif
      ! end of the extension for the method with the overlapped state
      return
   end subroutine engclear

   subroutine engproc_cleanup
      use engmain, only: io_flcuv
      use mpiproc, only: myrank
      implicit none
      if(myrank == 0) then
         endfile(io_flcuv)
         close(io_flcuv)
      endif
   end subroutine engproc_cleanup
   !
   !
   subroutine engconst(stnum)
      use engmain, only: nummol, skpcnf, slttype, sluvid, &
         maxins, numslv, numslt, cltype, cell, &
         io_flcuv, &
         SYS_NONPERIODIC, SYS_PERIODIC, &
         EL_COULOMB, EL_PME, EL_PPPM, ES_NVT, ES_NPT, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
         PT_SOLVENT, PT_SOLUTE
      use reciprocal, only: recpcal_init, recpcal_spline_greenfunc, &
         recpcal_prepare_solvent, recpcal_pppm_greenfunc
      use mpiproc                                                      ! MPI
      implicit none
      integer, intent(in) :: stnum
      integer :: i, k, irank
      real :: stat_weight_solute
      integer, dimension(:), allocatable :: tplst
      real, dimension(:),    allocatable :: uvengy
      logical, allocatable :: flceng_stored_g(:,:)
      real, allocatable :: flceng_g(:,:,:)
      integer :: flceng_mpikind
      real, save :: prevcl(3,3)
      logical :: skipcond
      logical, save :: pme_initialized = .false. ! NOTE: this variable is also used when PPPM is selected.

      call mpi_init_active_group(nactiveproc)
      call sanity_check_sluvid()

      ! for soln: maxdst is number of solute molecules
      ! for refs: maxdst is number of insertions
      select case(slttype)
       case(SLT_SOLN)
         maxdst = numslt
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
         maxdst = maxins
      end select

      allocate( tplst(nummol) )
      slvmax = 0
      do i = 1, nummol
         ! particle exists in trajectory (not a test particle)
         if((sluvid(i) == PT_SOLVENT) .or. (sluvid(i) == PT_SOLUTE)) then
            slvmax = slvmax + 1
            tplst(slvmax) = i
         end if
      enddo
      allocate( tagpt(slvmax) )
      allocate( uvengy(0:slvmax) )
      tagpt(1:slvmax) = tplst(1:slvmax)  ! and copied from tplst
      deallocate( tplst )

      allocate( flceng_stored(maxdst) )
      allocate( flceng(numslv, maxdst) )
      flceng_stored(:) = .false.
      flceng(:,:) = 0

      if(myrank < nactiveproc) then
         ! Initialize reciprocal space - grid and charges
         call get_inverted_cell
         if(cltype == EL_PME .or. cltype == EL_PPPM) then
            if(.not. pme_initialized) call recpcal_init(slvmax, tagpt)

            ! check whether cell size changes
            ! recpcal is called only when cell size differ
            call perf_time("kgrn")
            if((.not. pme_initialized) .or. &
               (any(prevcl(:,:) /= cell(:,:)))) then
               if (cltype == EL_PME) then
                  call recpcal_spline_greenfunc()
               elseif (cltype == EL_PPPM) then
                  call recpcal_pppm_greenfunc()
               endif
            end if
            call perf_time()
            prevcl(:,:) = cell(:,:)

            pme_initialized = .true.

            call perf_time("kpre")
            !$omp parallel do schedule(dynamic) private(k, i)
            do k = 1, slvmax
               i = tagpt(k)
               call recpcal_prepare_solvent(i)
            enddo
            call perf_time()
         endif

         ! cntdst is the pick-up no. of solute molecule from plural solutes (soln)
         ! cntdst is the iteration no. of insertion (refs)
         do cntdst = 1, maxdst
            call get_uv_energy(stnum, stat_weight_solute, uvengy(0:slvmax), skipcond)
            if(skipcond) cycle
            call update_histogram(stat_weight_solute, uvengy(0:slvmax))
         enddo

         select case(slttype)
          case(SLT_SOLN)                           ! for soln: output flceng
            allocate( flceng_stored_g(maxdst, nactiveproc) )
            allocate( flceng_g(numslv, maxdst, nactiveproc) )

#ifdef MPI
            ! gather flceng values to rank 0
            if(kind(flceng) /= kind(flceng_g)) then
               stop "bug in the kind(real) for flceng or flceng_g"
            else
               call get_mympi_realkind(kind(flceng), flceng_mpikind)
            endif
            call mpi_gather(flceng_stored, maxdst, mpi_logical, &
               flceng_stored_g, maxdst, mpi_logical, &
               0, mpi_comm_activeprocs, ierror)
            call mpi_gather(flceng, numslv * maxdst, flceng_mpikind, &
               flceng_g, numslv * maxdst, flceng_mpikind, &
               0, mpi_comm_activeprocs, ierror)
#else
            flceng_stored_g(:,1) = flceng_stored(:)
            flceng_g(:,:,1) = flceng(:,:)
#endif

            if(myrank == 0) then
               do irank = 1, nactiveproc
                  do cntdst = 1, maxdst
                     if(flceng_stored_g(cntdst, irank)) then
                        if(maxdst == 1) then
                           write(io_flcuv, 911) &
                              (stnum + irank - 1) * skpcnf, &
                              flceng_g(1:numslv, cntdst, irank)
                        else
                           write(io_flcuv, 912) cntdst, &
                              (stnum + irank - 1) * skpcnf, &
                              flceng_g(1:numslv, cntdst, irank)
                        endif
911                     format(i9, 9999f15.5)
912                     format(2i9, 9999f15.5)
                     endif
                  enddo
               enddo
            endif
            deallocate( flceng_g, flceng_stored_g )

          case(SLT_REFS_RIGID, SLT_REFS_FLEX)      ! for refs: output progress
            if(myrank == 0) write(io_flcuv, *) ( &
               (stnum + irank - 1) * skpcnf, &
               irank = 1, nactiveproc)
         end select
#ifdef HAVE_FLUSH
         if (myrank == 0) flush(io_flcuv)
#endif
      endif

      deallocate( tagpt, uvengy )
      deallocate( flceng, flceng_stored )

      call mpi_finish_active_group()
      return
   end subroutine engconst
   !
   !
   subroutine engstore(stnum)
      !
      use engmain, only: maxcnf, skpcnf, engdiv, corrcal, selfcal, &
         slttype, wgtslf, &
         ermax, numslv, esmax, temp, &
         edens, ecorr, eself, &
         aveuv, slnuv, avediv, avslf, minuv, maxuv, &
         engnorm, engsmpl, voffset, voffset_initialized, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, NO, YES
      use mpiproc                                                      ! MPI
      implicit none
      integer :: stnum, pti, j, k, iduv, division, reduce_mpikind
      character(len=9) :: engfile
      character(len=3) :: suffeng
      integer, parameter :: eng_io = 51, cor_io = 52, slf_io = 53
      integer, parameter :: ave_io = 54, wgt_io = 55, uvr_io = 56
      real(kind=8) :: voffset_local, voffset_scale
      real :: factor
      real(kind=8) :: invwt, leftbin, middlebin
      ! start of the extension for the method with the overlapped state
      integer, parameter :: prb_io = 81
      ! end of the extension for the method with the overlapped state
      call mpi_rank_size_info                                          ! MPI

      ! synchronize voffset
      if(wgtslf == YES) then
         voffset_local = voffset
#ifdef MPI
         if(kind(voffset) /= kind(voffset_local)) then
            stop "bug in the kind(real) for voffset or voffset_local"
         else
            call get_mympi_realkind(kind(voffset), reduce_mpikind)
            call mpi_allreduce(voffset_local, voffset, 1, &
               reduce_mpikind, mpi_max, mpi_comm_world, ierror)
         endif

         ! if uninitialized, use voffset so as not to pollute results with NaN
         if(.not. voffset_initialized) voffset_local = voffset

         ! scale histograms accoording to the maximum voffset
         select case(slttype)
          case(SLT_SOLN)
            voffset_scale = exp((voffset_local - voffset)/temp)
          case(SLT_REFS_RIGID, SLT_REFS_FLEX)
            voffset_scale = exp(-(voffset_local - voffset)/temp)
         end select

         engnorm = engnorm * voffset_scale
         if(selfcal == YES) eself(:) = eself(:) * voffset_scale
         if(slttype == SLT_SOLN) slnuv(:) = slnuv(:) * voffset_scale
         edens(:) = edens(:) * voffset_scale
         if(corrcal == YES) ecorr(:, :) = ecorr(:, :) * voffset_scale
#endif
      endif

#ifdef MPI
      ! Gather all information to master node
      ! Only the variables which are initialized in the subsequent engclear
      !    can be subject to the mympi_reduce_real_* procedures
      call mympi_reduce_real8_scalar(avslf, mpi_sum, 0)
      call mympi_reduce_real8_scalar(engnorm, mpi_sum, 0)
      call mympi_reduce_real8_scalar(engsmpl, mpi_sum, 0)
      if(selfcal == YES) call mympi_reduce_real8_array(eself, esmax, mpi_sum, 0)
      if(slttype == SLT_SOLN) call mympi_reduce_real_array(slnuv, numslv, mpi_sum, 0)
      call mympi_reduce_real8_array(edens, ermax, mpi_sum, 0)
      if(corrcal == YES) call mympi_reduce_real8_array(ecorr, (ermax * ermax), mpi_sum, 0)
#endif


      ! index for the division of the trajectory
      division = stnum / (maxcnf / skpcnf / engdiv)

      if(division == engdiv) then
#ifdef MPI
         call mympi_reduce_real_array(minuv, (numslv + 1), mpi_min, 0)
         call mympi_reduce_real_array(maxuv, (numslv + 1), mpi_max, 0)
#endif
      endif

      ! start of the extension for the method with the overlapped state
      if(do_overlapped == YES) then
#ifdef MPI
         if(wgtslf == YES) then
            soln_total = soln_total * voffset_scale
            soln_ovlp = soln_ovlp * voffset_scale
            norm_ovlp = norm_ovlp * voffset_scale
            edovlp(:) = edovlp(:) * voffset_scale
            uvovlp(:) = uvovlp(:) * voffset_scale
         endif

         call mympi_reduce_real_scalar(norm_ovlp, mpi_sum, 0)
         call mympi_reduce_real_array(uvovlp, numslv, mpi_sum, 0)
         call mympi_reduce_real_array(edovlp, ermax, mpi_sum, 0)
         if(division == engdiv) then
            call mympi_reduce_real_scalar(soln_total, mpi_sum, 0)
            call mympi_reduce_real_scalar(soln_ovlp, mpi_sum, 0)
         endif
#endif
         if(myrank == 0) then
            if(engdiv == 1) then
               suffeng = '.tt'
            else
               j = division / 10
               k = mod(division, 10)
               suffeng = '.' // numbers(j+1:j+1) // numbers(k+1:k+1)
            endif

            edovlp(:) = edovlp(:) / norm_ovlp

            engfile = 'engovl' // suffeng
            open(unit = eng_io, file = engfile, form = "FORMATTED", action = 'write')
            write(eng_io, '(A2,A13,A15,A5,A25,A25)') "# ", "bin_left", "bin_middle", "spec", "ratio", "dens"
            do iduv = 1, ermax
               call representative_bin_info(iduv, leftbin, middlebin, invwt, pti)
               write(eng_io, '(g15.7,g15.7,i5,g25.15,g25.15)') leftbin, middlebin, pti, edovlp(iduv), edovlp(iduv) * invwt
            enddo
            endfile(eng_io)
            close(eng_io)

            avovlp(division, :) = uvovlp(:) / norm_ovlp

            if(division == engdiv) then
               factor = soln_ovlp / soln_total
               open(unit = prb_io, file = 'Prob_Result', action = 'write')
               write(prb_io, '(a,f12.6)') ' Probability to find the system in the overlapped state = ', factor
               write(prb_io, '(a)') ' Free-energy change from the overlapped state to the solution state'
               write(prb_io, '(a,f12.6,a)') '                                                = ', temp * log(factor), ' kcal/mol'
               endfile(prb_io)
               close(prb_io)

               open(unit = ave_io, file = 'avovl.tt', action = 'write')
               do k = 1, engdiv
                  write(ave_io, 751) k, avovlp(k, 1:numslv)
               enddo
               endfile(ave_io)
               close(ave_io)
            endif
         endif
      endif
      ! end of the extension for the method with the overlapped state

      if(myrank /= 0) return                                            ! MPI
      ! data to be stored; only the master node matters

      edens(:) = edens(:) / engnorm
      if(corrcal == YES) ecorr(:,:) = ecorr(:,:) / engnorm
      if(selfcal == YES) eself(:) = eself(:) / engnorm

      avslf = avslf / engnorm
      avediv(division, 1) = engnorm / engsmpl

      select case(slttype)
       case(SLT_SOLN)
         aveuv(division, 1:numslv) = slnuv(1:numslv) / engnorm
         avediv(division, 2) = voffset - temp * log(avslf)
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
         avediv(division, 2) = voffset + temp * log(avslf)
      end select
      !
      if(division == engdiv) then
         if(slttype == SLT_SOLN) then
            open(unit = ave_io, file = 'aveuv.tt', action = 'write')
            do k = 1, engdiv
               write(ave_io, 751) k, aveuv(k, 1:numslv)
751            format(i5, 9999f15.5)
            enddo
            endfile(ave_io)
            close(ave_io)
         endif

         select case(slttype)
          case(SLT_SOLN)
            open(unit = wgt_io, file = 'weight_soln', action = 'write')
          case(SLT_REFS_RIGID, SLT_REFS_FLEX)
            open(unit = wgt_io, file = 'weight_refs', action = 'write')
         end select
         do k = 1, engdiv
            select case(wgtslf)
             case(NO)
               write(wgt_io, '(i5,e20.8)') k, avediv(k,1)
             case(YES)
               write(wgt_io, '(i5,e20.8,e20.7)') k, avediv(k,1), avediv(k,2)
            end select
         enddo
         endfile(wgt_io)
         close(wgt_io)

         open(uvr_io, file = 'uvrange.tt', action = 'write')
         write(uvr_io,'(A)') ' species     minimum        maximum'
         do pti = 0, numslv
            if(maxuv(pti) < 1.0e5) then
               write(uvr_io, '(i5,2f15.5)') pti, minuv(pti), maxuv(pti)
            else
               write(uvr_io, '(i5,f15.5,g18.5)') pti, minuv(pti), maxuv(pti)
            endif
         enddo
         endfile(uvr_io)
         close(uvr_io)
      endif

      if(engdiv == 1) then
         suffeng = '.tt'
      else
         j = division / 10
         k = mod(division, 10)
         suffeng = '.' // numbers(j+1:j+1) // numbers(k+1:k+1)
      endif
      !
      ! storing the solute-solvent pair-energy distribution function
      select case(slttype)
       case(SLT_SOLN)
         engfile = 'engsln' // suffeng
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
         engfile = 'engref' // suffeng
      end select
      open(unit = eng_io, file = engfile, form = "FORMATTED", action = 'write')
      write(eng_io, '(A2,A13,A15,A5,A25,A25)') "# ", "bin_left", "bin_middle", "spec", "ratio", "dens"
      do iduv = 1, ermax
         call representative_bin_info(iduv, leftbin, middlebin, invwt, pti)
         write(eng_io, '(g15.7,g15.7,i5,g25.15,g25.15)') leftbin, middlebin, pti, edens(iduv), edens(iduv) * invwt
      enddo
      endfile(eng_io)
      close(eng_io)
      !
      ! storing the solvent-solvent correlation matrix in energy representation
      if(corrcal == YES) then
         select case(slttype)
          case(SLT_SOLN)
            engfile = 'corsln' // suffeng
          case(SLT_REFS_RIGID, SLT_REFS_FLEX)
            engfile = 'corref' // suffeng
         end select
         open(unit = cor_io, file = engfile, form = "UNFORMATTED", action = 'write')
         write(cor_io) ecorr
         endfile(cor_io)
         close(cor_io)
      endif
      !
      ! storing the distribution function of the self-energy of solute
      if(selfcal == YES) then
         engfile = 'slfeng' // suffeng
         open(unit = slf_io, file = engfile, form = "FORMATTED", action = 'write')
         do iduv = 1, esmax
            factor = representative_bin_selfenergy(iduv)
            write(slf_io, '(g15.7,g25.15)') factor, eself(iduv)
         enddo
         endfile(slf_io)
         close(slf_io)
      endif
      !
      return
   end subroutine engstore

   ! Calculate interaction energy between solute and solvent
   subroutine get_uv_energy(stnum, stat_weight_solute, uvengy, out_of_range)
      use engmain, only: maxcnf, skpcnf, slttype, sltlist, cltype, &
         EL_COULOMB, EL_PME, EL_PPPM, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX
      use ptinsrt, only: instslt
      use realcal, only: realcal_prepare, realcal_proc, realcal_self, &
         realcal_bare, realcal_cleanup
      use reciprocal, only: recpcal_prepare_solute, recpcal_self_energy, &
         recpcal_energy
      use mpiproc                                                      ! MPI
      implicit none
      integer, intent(in) :: stnum
      real, intent(out) :: uvengy(0:slvmax), stat_weight_solute
      logical, intent(out) :: out_of_range

      integer :: i, k
      integer(8) :: current_solute_hash
      integer(8) :: solute_hash = 0
      real :: pairep, residual, factor, uvrecp, pair_real, pair_recp
      real, save :: usreal
      logical, save :: initialized = .false.

      out_of_range = .false.
      ! determine / pick solute structure
      select case(slttype)
       case(SLT_SOLN)
         tagslt = sltlist(cntdst)
         stat_weight_solute = 1.0
         call check_mol_configuration(out_of_range)
         if(out_of_range) return
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
         tagslt = sltlist(1)
         if(.not. initialized) call instslt('init')
         initialized = .true.
         call instslt('proc', cntdst, stat_weight_solute)
         if((stnum == (maxcnf / skpcnf)) .and. (cntdst == maxdst)) call instslt('last')
      end select

      ! At this moment all coordinate in the system is determined
      call realcal_prepare

      uvengy(:) = 0
      uvrecp = 0.
      ! Calculate system-wide values
      if(cltype == EL_PME .or. cltype == EL_PPPM) then
         call recpcal_prepare_solute(tagslt)
         call perf_time("rblk")
         call realcal_proc(tagslt, tagpt, slvmax, uvengy)
         call perf_time()
         call perf_time("kslf")
         uvrecp = recpcal_self_energy()
         call perf_time()
      endif

      ! solute-solute self energy
      pairep = 0.0
      residual = 0.0
      current_solute_hash = get_solute_hash() ! FIXME: if this tuns into a bottleneck, add conditionals
      if(current_solute_hash == solute_hash .or. &
         (slttype == SLT_REFS_RIGID .and. solute_hash /= 0)) then
         ! For refs calculation, the configuration of solute may change with
         ! random translation and/or rotation upon insertion,
         ! though the self energy will not change.
         ! usreal is kept over subroutine calls during "save"
      else
         call perf_time("rslf")
         usreal = realcal_self(tagslt) ! calculate self-interaction
         call perf_time()
      endif
      solute_hash = current_solute_hash
      residual = residual_ene(tagslt, tagslt)
      uvengy(0) = uvrecp + usreal + residual

      ! solute-solvent pair
      !$omp parallel do schedule(dynamic) private(i, k, pairep, pair_real, pair_recp)
      do k = 1, slvmax
         i = tagpt(k)
         if(i == tagslt) cycle

         pair_real = 0
         pair_recp = 0
         if(cltype == EL_PME .or. cltype == EL_PPPM) then
            ! called only when PME or PPPM, non-self interaction
            pair_real = residual_ene(tagslt, i)
            pair_recp = recpcal_energy(tagslt, i)
            pairep = pair_real + pair_recp
         else                      ! Bare coulomb solute-solvent interaction
            pair_real = realcal_bare(tagslt, i)
            pairep = pair_real
         endif
         !$omp atomic
         uvengy(k) = uvengy(k) + pairep
      enddo

      call realcal_cleanup
   end subroutine get_uv_energy

   subroutine update_histogram(stat_weight_solute, uvengy)
      use engmain, only: wgtslf, estype, slttype, corrcal, selfcal, ermax, &
         volume, temp, uvspec, &
         slnuv, avslf, minuv, maxuv, &
         edens, ecorr, eself, &
         stat_weight_system, engnorm, engsmpl, &
         voffset, voffset_initialized, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
         ES_NVT, ES_NPT, NO, YES
      use mpiproc
      implicit none
      real, intent(in) :: uvengy(0:slvmax), stat_weight_solute
      integer, allocatable :: insdst(:)
      integer :: i, k, q, iduv, iduvp, pti
      real :: factor, pairep, total_weight
      real(kind=8) :: engnmfc

      allocate( insdst(ermax) )

      engnmfc = 0.0

      select case(wgtslf)
       case(NO)
         engnmfc = 1.0
       case(YES)
         if(.not. voffset_initialized) then
            voffset = uvengy(0)
            voffset_initialized = .true.
         endif
         factor = uvengy(0) - voffset  ! shifted by offset
         select case(slttype)
          case(SLT_SOLN)
            engnmfc = exp(factor / temp)
          case(SLT_REFS_RIGID, SLT_REFS_FLEX)
            engnmfc = exp(- factor / temp)
         end select
       case default
         stop "Unknown wgtslf"
      end select
      total_weight = stat_weight_system * stat_weight_solute
      if(estype == ES_NPT) call volcorrect(total_weight)
      engnmfc = engnmfc * total_weight
      !
      engnorm = engnorm + engnmfc      ! normalization factor
      engsmpl = engsmpl + 1.0          ! number of sampling
      avslf = avslf + total_weight     ! normalization without solute self-energy

      ! self energy histogram
      if(selfcal == YES) then
         call getiduv(0, uvengy(0), iduv)
         eself(iduv) = eself(iduv) + engnmfc
      endif
      minuv(0) = min(minuv(0), uvengy(0))
      maxuv(0) = max(maxuv(0), uvengy(0))

      insdst(:) = 0                              ! instantaneous histogram
      flceng(:, cntdst) = 0.0                    ! sum of solute-solvent energy
      flceng_stored(cntdst) = .true.
      do k = 1, slvmax
         i = tagpt(k)
         if(i == tagslt) cycle
         pti = uvspec(i)
         if(pti <= 0) call halt_with_error('eng_cns')

         pairep = uvengy(k)
         call getiduv(pti, pairep, iduv)

         ! instantaneous histogram of solute-solvent energy
         insdst(iduv) = insdst(iduv) + 1
         ! sum of solute-solvent energy
         flceng(pti, cntdst) = flceng(pti, cntdst) + pairep

         ! minimum and maxmium of solute-solvent energy
         minuv(pti) = min(minuv(pti), pairep)
         maxuv(pti) = max(maxuv(pti), pairep)
      enddo

      if(slttype == SLT_SOLN) then
         slnuv(:) = slnuv(:) + flceng(:, cntdst) * engnmfc
      endif

      do iduv = 1, ermax
         k = insdst(iduv)
         if(k == 0) cycle
         edens(iduv) = edens(iduv) + engnmfc * real(k)
      enddo
      if(corrcal == YES) then
         do iduv = 1, ermax
            k = insdst(iduv)
            if(k == 0) cycle
            do iduvp = 1, ermax
               q = insdst(iduvp)
               if(q == 0) cycle
               ecorr(iduvp,iduv) = ecorr(iduvp,iduv) + engnmfc * real(k) * real(q)
            enddo
         enddo
      endif

      ! start of the extension for the method with the overlapped state
      if(do_overlapped == YES) then
         soln_total = soln_total + engnmfc
         if(any( (insdst(:) > 0) .and. zerodns_crd(:) )) then
         else
            soln_ovlp = soln_ovlp + engnmfc
            norm_ovlp = norm_ovlp + engnmfc
            do iduv = 1, ermax
               k = insdst(iduv)
               if(k > 0) edovlp(iduv) = edovlp(iduv) + engnmfc * real(k)
            enddo
            uvovlp(:) = uvovlp(:) + flceng(:, cntdst) * engnmfc
         endif
      endif
      ! end of the extension for the method with the overlapped state

      deallocate( insdst )
   end subroutine update_histogram

   !
   function residual_ene(i, j) result(pairep)
      use engmain, only: screen, volume, mol_charge, cltype, EL_COULOMB, PI
      implicit none
      integer, intent(in) :: i, j
      real :: pairep
      real :: epcl
      if(cltype == EL_COULOMB) return
      epcl = PI * mol_charge(i) * mol_charge(j) / screen / screen / volume
      if(i == j) epcl = epcl / 2.0   ! self-interaction
      pairep = - epcl
   end function residual_ene
   !
   subroutine volcorrect(weight)
      use engmain, only:  sluvid, cltype, screen, mol_charge, volume, temp, &
         EL_EWALD, EL_PME, EL_PPPM, PT_SOLVENT, PT_SOLUTE, PI
      implicit none
      real, intent(inout) :: weight
      real :: total_charge, factor
      weight = weight * volume
      if((cltype == EL_EWALD) .or. (cltype == EL_PME) &
         .or. (cltype == EL_PPPM)) then  ! Ewald, PME and PPPM
         ! sum of charges over physical particles
         total_charge = sum( mol_charge, mask = ((sluvid == PT_SOLVENT) &
            .or. (sluvid == PT_SOLUTE)) )
         factor = PI * (total_charge ** 2)  / (screen ** 2 ) / volume / 2.0
         weight = weight * exp(factor / temp)
      endif
      return
   end subroutine volcorrect
   !
   subroutine update_cell_info()
      use engmain, only: cell, celllen
      implicit none
      integer :: i
      do i = 1, 3
         celllen(i) = sqrt(sum(cell(:, i) ** 2))
      enddo
   end subroutine update_cell_info

   subroutine get_inverted_cell
      use engmain, only:  cell, invcl, volume
      implicit none
      volume = cell(1,1) * cell(2,2) * cell(3,3) &
         + cell(1,2) * cell(2,3) * cell(3,1) &
         + cell(1,3) * cell(2,1) * cell(3,2) &
         - cell(1,3) * cell(2,2) * cell(3,1) &
         - cell(1,2) * cell(2,1) * cell(3,3) &
         - cell(1,1) * cell(2,3) * cell(3,2)
      invcl(1,1) = cell(2,2) * cell(3,3) - cell(2,3) * cell(3,2)
      invcl(1,2) = cell(1,3) * cell(3,2) - cell(1,2) * cell(3,3)
      invcl(1,3) = cell(1,2) * cell(2,3) - cell(1,3) * cell(2,2)
      invcl(2,1) = cell(2,3) * cell(3,1) - cell(2,1) * cell(3,3)
      invcl(2,2) = cell(1,1) * cell(3,3) - cell(1,3) * cell(3,1)
      invcl(2,3) = cell(1,3) * cell(2,1) - cell(1,1) * cell(2,3)
      invcl(3,1) = cell(2,1) * cell(3,2) - cell(2,2) * cell(3,1)
      invcl(3,2) = cell(1,2) * cell(3,1) - cell(1,1) * cell(3,2)
      invcl(3,3) = cell(1,1) * cell(2,2) - cell(1,2) * cell(2,1)

      invcl(1:3, 1:3) = invcl(1:3, 1:3) / volume
      call update_cell_info
   end subroutine get_inverted_cell


   ! binsearch returns the smallest index (ret) which satisfies
   ! coord(ret) >= v
   ! Special cases are as follows:
   ! v < coord(1)  ==>  ret = 0 (out of range)
   ! v > coord(n)  ==>  ret = n (infinite)
   subroutine binsearch(coord, n, v, ret)
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: coord(n)
      real, intent(in) :: v
      integer, intent(out) :: ret
      integer :: rmin, rmax, rmid
      if(v < coord(1)) then
         ret = 0
         return
      endif
      if(v > coord(n)) then
         ret = n
         return
      endif

      rmin = 1
      rmax = n + 1
      do
         if(rmax - rmin <= 1) then
            exit
         endif
         rmid = ishft((rmin + rmax - 1), -1) ! equal to (rmin + rmax - 1) / 2
         if(v > coord(rmid)) then
            rmin = rmid + 1
         else
            rmax = rmid + 1
         endif
      enddo
      ret = rmin - 1
   end subroutine binsearch

   ! returns the position of the bin corresponding to energy coordinate value
   subroutine getiduv(pti, engval, iduv)
      use engmain, only: slttype, uvmax, uvsoft, uvcrd, esmax, escrd, &
         stdout, SLT_SOLN
      use mpiproc, only: halt_with_error, warning
      implicit none
      integer, intent(in) :: pti
      real, intent(in) :: engval
      integer, intent(out) :: iduv
      integer :: idmin, idnum, idpti
      real, parameter :: warn_threshold = 1.0e3
      logical, save :: warn_bin_firsttime = .true.

      idmin = 0 ! clear to supress warning
      idpti = 0
      idnum = 0

      if(pti >  0) then            ! solute-solvent interaction
         if(pti == 1) then         ! first solvent species
            idmin = 0
         elseif(pti > 1) then      ! second and later solvent species
            idmin = sum( uvmax(1:(pti - 1)) )
         else                      ! bug --- pti must not be smaller than 1
            call halt_with_error('eng_bug')
         endif
         idnum = uvmax(pti)
         call binsearch(uvcrd((idmin + 1):(idmin + idnum)), idnum, engval, idpti)
      elseif(pti == 0) then        ! solute self-energy
         idmin = 0
         idnum = esmax
         call binsearch(escrd(1:idnum), idnum, engval, idpti)
      else                         ! bug --- pti must be non-negative
         call halt_with_error('eng_bug')
      endif

      ! inappropriate setting of energy coordinate
      ! smaller than the minimum energy mesh
      if(idpti <= 0) then
         write(stdout, '(A,g12.4,A,i3,A)') '  energy of ', engval, ' for ', pti, '-th species'
         call halt_with_error('eng_min')
      endif
      ! larger than the maximum energy mesh
      if(idpti > idnum) call halt_with_error('eng_bug')
      ! only for solute-solvent interaction in solution system
      ! larger than the maximum of soft part (linear-graduation part)
      if((slttype == SLT_SOLN) .and. (pti > 0)) then
         if(idpti > uvsoft(pti)) then
            write(stdout, '(A,g12.4,A,i3,A)') '  energy of ', engval, ' for ', pti, '-th species'
            call halt_with_error('eng_sft')
         endif
      endif

      ! Warning if the energy exceeds the maximum binning region and pecore = 0
      ! Since it is hard to see the pecore value at this point,
      ! the energy value is simply shown as a warning
      if((idpti == idnum) .and. (engval < warn_threshold) &
         .and. (warn_bin_firsttime)) then
         write(stdout, '(A,g12.4,A,i3,A)') '  energy of ', engval, ' for ', pti, '-th species'
         call warning('mbin')
         warn_bin_firsttime = .false.
      endif

      iduv = idmin + idpti

      return
   end subroutine getiduv

   ! Check system consistency: either test particle or solvent must exist
   subroutine sanity_check_sluvid()
      use engmain, only: slttype, nummol, sluvid, SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX
      use mpiproc, only: halt_with_error
      implicit none

      ! sanity check
      if(any(sluvid(:) < 0) .or. any(sluvid(:) > 3)) call halt_with_error('eng_bug')
      select case(slttype)
       case(SLT_SOLN)
         ! sluvid should be 0 (solvent) or 1 (solute)
         if(any(sluvid(:) >= 2)) call halt_with_error('eng_par')
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
         ! sluvid should be 0 (solvent), 2, 3 (test particles)
         if(any(sluvid(:) == 1)) call halt_with_error('eng_par')
         ! solvent must exist
         if(all(sluvid(:) /= 0)) call halt_with_error('eng_par')
      end select

      ! solute / test particle must exist
      if(all(sluvid(:) == 0)) call halt_with_error('eng_par')
   end subroutine sanity_check_sluvid

   ! Check whether molecule is within specified region
   subroutine check_mol_configuration(out_of_range)
      use engmain, only: insposition, insstructure, lwreg, upreg, lwstr, upstr, &
         numsite, mol_begin_index, mol_end_index, &
         sitepos, boxshp, invcl, celllen, &
         SYS_NONPERIODIC, &
         INSPOS_RANDOM, INSPOS_NOCHANGE, &
         INSPOS_SPHERE, &
         INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
         INSPOS_RMSD, INSPOS_GAUSS, &
         INSSTR_NOREJECT, INSSTR_RMSD
      use ptinsrt, only: refhost_natom, refhost_specatm, &
         refhost_crd, refhost_weight, &
         refslt_natom, refslt_crd, refslt_weight, &
         insscheme
      use bestfit, only: fit_a_rotate_b, rmsd_nofit, rmsd_bestfit
      use mpiproc, only: halt_with_error
      implicit none
      logical, intent(out) :: out_of_range
      real :: dx(3), distance
      integer :: i, ptb, pte
      real, dimension(:,:), allocatable :: hostcrd, refslt_bestfit

      out_of_range = .false.

      select case(insposition)
       case(INSPOS_RANDOM)                               ! random
         ! do nothing
       case(INSPOS_NOCHANGE)                             ! fixed configuration
         ! do nothing
       case(INSPOS_SPHERE)                               ! sphere geometry
         ! The following has the same structure as the corresponding part
         ! in the set_shift_com subroutine within insertion.F90
         call relative_com(tagslt, dx)
         distance = sqrt(dot_product(dx, dx))
         if((lwreg > distance) .or. (distance > upreg)) then
            out_of_range = .true.     ! configuration is rejected
            return
         endif
       case(INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC)  ! slab configuration
         ! constrained to z-axis
         ! The following has the same structure as the corresponding part
         ! in the set_shift_com subroutine within insertion.F90
         if(boxshp == SYS_NONPERIODIC) call halt_with_error('eng_slb')
         call relative_com(tagslt, dx)
         distance = dot_product(invcl(3,:), dx(:)) * celllen(3)
         if(insposition == INSPOS_SLAB_SYMMETRIC) then  ! symmetric bilayer
            distance = abs(distance)
         endif
         if((lwreg > distance) .or. (distance > upreg)) then
            out_of_range = .true.     ! configuration is rejected
            return
         endif
       case(INSPOS_RMSD)                                 ! comparison to reference
         ptb = mol_begin_index(tagslt)
         pte = mol_end_index(tagslt)
         if(numsite(tagslt) /= refslt_natom) call halt_with_error('eng_bug')
         allocate( hostcrd(3, refhost_natom), refslt_bestfit(3, refslt_natom) )
         ! The following has the same structure as the corresponding part
         ! in the reffit subroutine within insertion.F90
         do i = 1, refhost_natom
            hostcrd(1:3, i) = sitepos(1:3, refhost_specatm(i))
         enddo
         call fit_a_rotate_b(refhost_natom, hostcrd, &
            refhost_crd, refhost_weight, &
            refslt_natom, refslt_crd, refslt_bestfit)
         ! The following has the same structure as the corresponding part
         ! in the check_solute_configuration subroutine within insertion.F90
         distance = rmsd_nofit(refslt_natom, refslt_bestfit, &
            sitepos(1:3, ptb:pte), refslt_weight)
         deallocate( hostcrd, refslt_bestfit )
         if((lwreg > distance) .or. (distance > upreg)) then
            out_of_range = .true.     ! configuration is rejected
            return
         endif
       case(INSPOS_GAUSS)                                ! comparison to reference
         ! Exeperimental and under construction...  What is to be written?
       case default
         stop "Unknown insposition in check_mol_configuration"
      end select

      select case(insstructure)
       case(INSSTR_NOREJECT)    ! no rejection of solute structure
         ! do nothing
       case(INSSTR_RMSD)        ! solute structure rejection with RMSD
         ! The following has the same structure as the corresponding part
         ! in the getsolute subroutine within insertion.F90
         ptb = mol_begin_index(tagslt)
         pte = mol_end_index(tagslt)
         if(numsite(tagslt) /= refslt_natom) call halt_with_error('eng_bug')
         distance = rmsd_bestfit(refslt_natom, refslt_crd, &
            sitepos(1:3, ptb:pte), refslt_weight)
         if((lwstr > distance) .or. (distance > upstr)) then
            out_of_range = .true.     ! structure is rejected
            return
         endif
       case default
         stop "Unknown insstructure in check_mol_configuration"
      end select

      call insscheme(tagslt, out_of_range)
      return

   contains
      subroutine relative_com(tagpt, dx)
         use engmain, only: numsite, mol_begin_index, mol_end_index, &
            sitemass, sitepos
         use ptinsrt, only: com_aggregate
         use bestfit, only: center_of_mass
         implicit none
         integer, intent(in) :: tagpt
         real, intent(out) :: dx(3)
         integer ptb, pte, stmax
         real :: solute_com(3), aggregate_com(3)
         real, dimension(:), allocatable   :: ptmass
         real, dimension(:,:), allocatable :: ptsite
         stmax = numsite(tagpt)
         ptb = mol_begin_index(tagpt)
         pte = mol_end_index(tagpt)
         allocate( ptmass(stmax), ptsite(3, stmax) )
         ptmass(1:stmax) = sitemass(ptb:pte)
         ptsite(1:3, 1:stmax) = sitepos(1:3, ptb:pte)
         call center_of_mass(stmax, ptsite, ptmass, solute_com)  ! solute COM
         call com_aggregate(aggregate_com)                       ! aggregate COM
         dx(1:3) = solute_com(1:3) - aggregate_com(1:3)
         deallocate( ptmass, ptsite )
         return
      end subroutine relative_com
   end subroutine check_mol_configuration

   ! get the hashed function of solute coordinate
   integer(8) function get_solute_hash()
      use utility, only: hash
      use engmain, only: sitepos, mol_begin_index, numsite
      implicit none

      get_solute_hash = hash(sitepos(1:3, mol_begin_index(tagslt):(mol_begin_index(tagslt+1) - 1)), numsite(tagslt) * 3)
   end function get_solute_hash


   function representative_bin_selfenergy(iduv) result(engcoord)
      use engmain, only: esmax, escrd
      use mpiproc, only: halt_with_error
      implicit none
      integer, intent(in) :: iduv
      real :: engcoord
      if(iduv <  esmax) engcoord = (escrd(iduv) + escrd(iduv+1)) / 2.0
      if(iduv == esmax) engcoord = escrd(esmax)
      if(iduv >  esmax) call halt_with_error('eng_ecd')
   end function representative_bin_selfenergy

   subroutine representative_bin_info(iduv, engleft, engmiddle, enginvwidth, solvent_spec)
      use engmain, only: numslv, uvmax, uvsoft, uvcrd, esmax, escrd
      use mpiproc, only: halt_with_error
      implicit none
      integer, intent(in) :: iduv
      real(kind=8), intent(out) :: engleft, engmiddle, enginvwidth
      integer, intent(out) :: solvent_spec
      integer :: idpt, cnt, idmin, idmax, idsoft
      idmin = 0
      do cnt = 1, numslv
         idpt = idmin + uvmax(cnt)
         if(iduv <= idpt) exit
         idmin = idpt
      enddo
      solvent_spec = cnt
      idsoft = uvsoft(solvent_spec)
      idmax = uvmax(solvent_spec)
      idpt = iduv - idmin
      if((idpt < 1) .or. (idpt > idmax)) call halt_with_error('eng_ecd')
      engleft = uvcrd(iduv)
      if(idpt <= idsoft) then   ! linear graduation
         engmiddle = (uvcrd(iduv) + uvcrd(iduv+1)) / 2.0
         enginvwidth = 1. / (uvcrd(iduv + 1) - uvcrd(iduv))
      else                      ! logarithmic graduation
         if(idpt <  idmax) then 
            engmiddle = sqrt(uvcrd(iduv) * uvcrd(iduv+1))
            enginvwidth = 1. / (uvcrd(iduv + 1) - uvcrd(iduv))
         elseif(idpt == idmax) then
            engmiddle = uvcrd(iduv)
            enginvwidth = 0.
         end if
      endif
   end subroutine representative_bin_info
end module engproc
