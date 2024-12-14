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
  integer :: nsolu_atom, nsolv_atom
  integer, allocatable :: block_solu(:, :), block_solv(:, :)
  integer, allocatable :: belong_solu(:), belong_solv(:)
  integer, allocatable :: atomno_solu(:), atomno_solv(:)
  integer, allocatable :: counts_solu(:, :, :), counts_solv(:, :, :)
  integer, allocatable :: psum_solu(:), psum_solv(:)
  real, allocatable :: sitepos_solu(:, :), sitepos_solv(:, :)
  integer, allocatable :: ljtype_solu(:, :), ljtype_solv(:, :)

  integer :: max_solu_block, max_solv_block

  real, allocatable :: ljeps_lowlj(:, :), ljsgm2_lowlj(:, :), dist_lowlj(:, :)
  integer, allocatable :: belong_lowlj(:, :)
  real, allocatable :: ljeps_switch(:, :), ljsgm2_switch(:, :), dist_switch(:, :)
  integer, allocatable :: belong_switch(:, :)
  real, allocatable :: charge_el(:, :), dist_el(:, :)
  integer, allocatable :: belong_el(:, :)
  real, allocatable :: e_t(:, :)

  ! subcell_neighbour only stores the neighbour list on y-z plane. x direction is stored on subcell_xlen.
  ! (looks like sub"pillar" rather than cell)
  ! each subcell spans x=(-subcell_xlen(ix) .. 1+subcell_xlen(ix)),  y=subcell_neighbour(2, ix), z=subcell_neighbour(3, ix)
  integer, allocatable :: subcell_neighbour(:, :) ! of (2:3, subcell_num_neighbour)
  integer, allocatable :: subcell_xlen(:)
  integer :: subcell_num_neighbour

  integer :: block_size(3)
  
  ! "straight" coordinate system
  real, allocatable :: sitepos_normal(:, :)
  real :: cell_normal(3, 3), invcell_normal(3), cell_len_normal(3)
  real :: invcell(3, 3)
  logical :: is_cuboid
  real, parameter :: check_rotate = 1e-8, cuboid_thres = 1e-8

contains
  subroutine realcal_proc(target_solu, tagpt, slvmax, uvengy)
    use engmain, only: numsite
    !$ use omp_lib, only: omp_get_num_procs
    implicit none
    integer, intent(in) :: target_solu, tagpt(:), slvmax
    real, intent(out) :: uvengy(0:slvmax)
    real, allocatable :: eng(:, :)
    integer :: lsize, i, j
    integer :: npar

    ! print *, "DEBUG: relcal_proc called"
    ! FIXME: fix calling convention & upstream call tree
    ! to calculate several solutes at once
    nsolu_atom = numsite(target_solu)
    nsolv_atom = count_solv(target_solu, tagpt, slvmax)

    call set_block_info()

    allocate(block_solu(3, nsolu_atom), block_solv(3, nsolv_atom))
    allocate(belong_solu(nsolu_atom), belong_solv(nsolv_atom))
    allocate(atomno_solu(nsolu_atom), atomno_solv(nsolv_atom))
    allocate(counts_solu(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))
    allocate(counts_solv(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))
    allocate(psum_solu(0:block_size(1) * block_size(2) *  block_size(3)))
    allocate(psum_solv(0:block_size(1) * block_size(2) *  block_size(3)))

    call set_solv_atoms(target_solu, tagpt, slvmax)
    call set_solu_atoms(target_solu)

    ! assertion
    ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged"

    call blockify(nsolu_atom, atomno_solu, block_solu)
    call blockify(nsolv_atom, atomno_solv, block_solv)

    call sort_block(block_solu, nsolu_atom, belong_solu, atomno_solu, counts_solu, psum_solu)
    call sort_block(block_solv, nsolv_atom, belong_solv, atomno_solv, counts_solv, psum_solv)

    allocate(sitepos_solv(3, nsolv_atom + 1))
    sitepos_solv(1:3, 1:nsolv_atom) = sitepos_normal(1:3, atomno_solv(1:nsolv_atom))
    sitepos_solv(1:3, nsolv_atom+1) = 0

    max_solu_block = maxval(counts_solu)
    max_solv_block = 0
    do i = 0, block_size(2) - 1
       do j = 0, block_size(3) - 1
          max_solv_block = max(max_solv_block, sum(counts_solv(:, i, j)))
       end do
    end do
    lsize = max_solu_block * max_solv_block

    npar = 1
    !$ npar = omp_get_num_procs()
    
    allocate(ljeps_lowlj(lsize, npar), ljsgm2_lowlj(lsize, npar), dist_lowlj(lsize, npar), belong_lowlj(lsize, npar))
    allocate(ljeps_switch(lsize, npar), ljsgm2_switch(lsize, npar), dist_switch(lsize, npar), belong_switch(lsize, npar))
    allocate(charge_el(lsize, npar), dist_el(lsize, npar), belong_el(lsize, npar))
    allocate(e_t(lsize, npar))
    ! assertion
    ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged after sorting"

    allocate(eng(1:slvmax, npar))
    eng(:, :) = 0.0
    call get_pair_energy(eng)

    uvengy(1:slvmax) = sum(eng(1:slvmax, 1:npar), 2)

    deallocate(ljeps_lowlj, ljsgm2_lowlj, dist_lowlj, belong_lowlj)
    deallocate(ljeps_switch, ljsgm2_switch, dist_switch, belong_switch)
    deallocate(charge_el, dist_el, belong_el)
    deallocate(e_t)

    deallocate(eng)
    deallocate(block_solu, belong_solu, atomno_solu, counts_solu, psum_solu)
    deallocate(block_solv, belong_solv, atomno_solv, counts_solv, psum_solv)
    deallocate(sitepos_solv)
    deallocate(subcell_neighbour, subcell_xlen)
  end subroutine realcal_proc

  subroutine realcal_prepare
    use engmain, only: numatm, sitepos, boxshp, SYS_PERIODIC
    implicit none
    
    allocate(sitepos_normal(3, numatm))
    sitepos_normal(:, :) = sitepos(:, :)

    ! "Straighten" box, and normalize coordinate system
    if(boxshp == SYS_PERIODIC) call normalize_periodic
  end subroutine realcal_prepare

  subroutine realcal_cleanup
    implicit none
    deallocate(sitepos_normal)
  end subroutine realcal_cleanup


  ! Calculate i-j interaction energy in the bare 1/r form
  subroutine realcal_bare(i, j, pairep)
    use engmain, only:  boxshp, numsite, &
         elecut, lwljcut, upljcut, cltype, screen, charge, specatm, &
         ljswitch, ljtype, ljtype_max, ljene_mat, ljlensq_mat, &
         SYS_NONPERIODIC, SYS_PERIODIC, EL_COULOMB, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX
    implicit none
    integer :: i, j, is, js, ismax, jsmax, ati, atj
    real :: reelcut, pairep, rst, dis2, invr2, invr3, invr6
    real :: eplj, epcl, xst(3), half_cell(3)
    real :: lwljcut2, upljcut2, lwljcut3, upljcut3, lwljcut6, upljcut6
    real :: ljeps, ljsgm2, ljsgm3, ljsgm6, vdwa, vdwb, swth, swfac
    real :: repA, repB, repC, attA, attB, attC
    integer :: ljtype_i, ljtype_j
    real, parameter :: infty = 1.0e50      ! essentially equal to infinity
    !
    if(i == j) stop "cannot happen: two particle arguments should not be the same"
    if(cltype /= EL_COULOMB) stop "cannot happen: realcal_bare is called only when cltype is 'bare coulomb'."

    if(boxshp == SYS_NONPERIODIC) reelcut=infty
    if(boxshp == SYS_PERIODIC) then
       reelcut = elecut
       half_cell(:) = 0.5 * cell_len_normal(:)
    endif

    pairep = 0.0
    ismax = numsite(i)
    jsmax = numsite(j)

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

    do is = 1, ismax
       do js = 1, jsmax
          ati = specatm(is,i)
          atj = specatm(js,j)
          ljtype_i = ljtype(ati)
          ljtype_j = ljtype(atj)
          xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
          if(boxshp == SYS_PERIODIC) then    ! when the system is periodic
             if(is_cuboid) then
                xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
             else
                xst(:) = xst(:) - matmul(cell_normal, anint(matmul(invcell, xst)))
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
             case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
                ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                if(rst <= lwljcut) then
                   vdwa = invr6 * invr6 - ljsgm6 *ljsgm6 / (lwljcut6 * upljcut6)
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
             case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
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
    end do
    !
    return
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

    !$omp parallel do private(is,js,ati,atj,epcl,rst,dis2,xst) reduction(+:pairep)
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
             xst(:) = xst(:) - matmul(cell_normal, anint(matmul(invcell, xst)))
          end if

          dis2 = sum(xst(1:3) ** 2)
          rst = sqrt(dis2)

          ! distinct (different ati and atj arguments for two charge variables)
          epcl = - charge(ati) * charge(atj) * erf(screen * rst) / rst
          pairep = pairep + epcl
       enddo
    enddo

  end subroutine realcal_self


  integer function count_solv(solu, tagpt, slvmax)
    use engmain, only: numsite
    implicit none
    integer, intent(in) :: solu, tagpt(:), slvmax
    integer :: i, j, cnt
    cnt = 0
    do i = 1, slvmax
       j = tagpt(i)
       if(j == solu) cycle
       cnt = cnt + numsite(j)
    end do
    count_solv = cnt
  end function count_solv

  subroutine set_solu_atoms(solu)
    use engmain, only: numsite, mol_begin_index
    implicit none
    integer, intent(in) :: solu
    integer :: i
    do i = 1, numsite(solu)
       atomno_solu(i) = mol_begin_index(solu) + (i - 1)
       belong_solu(i) = solu
    end do
  end subroutine set_solu_atoms

  subroutine set_solv_atoms(solu, tagpt, slvmax)
    use engmain, only: numsite, mol_begin_index
    implicit none
    integer, intent(in) :: solu, tagpt(:), slvmax
    integer :: i, j, k, cnt
    cnt = 1
    do i = 1, slvmax
       j = tagpt(i)
       if(j == solu) cycle
       do k = 1, numsite(j)
          atomno_solv(cnt) = mol_begin_index(j) + (k - 1)
          belong_solv(cnt) = i
          cnt = cnt + 1
       end do
    end do
  end subroutine set_solv_atoms

  subroutine set_block_info()
    use engmain, only: block_threshold, upljcut, elecut
    implicit none
    real :: unit_axes(3), cut2, l
    integer :: i, j, k, bmax, ix, xlen
    real, allocatable :: grid_dist(:, :)
    real, allocatable :: box_dist(:, :)
    
    ! get the length of axes
    ! assumes cell's 1st axis being x-axis, 2nd axis on x-y plane

    ! set block size
    block_size(:) = ceiling(cell_len_normal(:) / block_threshold)

    ! pre-calculate grid-grid distance and box-box distance
    bmax = maxval(block_size(:))
    allocate(grid_dist(0:bmax-1, 3))
    allocate(box_dist(0:bmax-1, 3))
    unit_axes(:) = cell_len_normal(:) / block_size(:)

    ! only have to calculate axis-wise distance (because of orthogonality)
    do j = 1, 3
       do i = 0, block_size(j) - 1
          grid_dist(i, j) = min(i, modulo(-i, block_size(j))) * unit_axes(j)
       end do
    end do
    do j = 1, 3
       do i = 0, block_size(j) - 1
          box_dist(i, j) = min(grid_dist(i, j),&
               grid_dist(modulo(i - 1, block_size(j)), j), &
               grid_dist(modulo(i + 1, block_size(j)), j))
       end do
    end do

    deallocate(grid_dist)
    
    ! make a subcell list
    cut2 = max(upljcut, elecut) ** 2
    ! count...
    subcell_num_neighbour = 0
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          l = box_dist(j, 2) ** 2 + box_dist(k, 3) ** 2
          if(l < cut2) subcell_num_neighbour = subcell_num_neighbour + 1
       end do
    end do
    allocate(subcell_neighbour(2:3, subcell_num_neighbour))
    allocate(subcell_xlen(subcell_num_neighbour))

    ! then generate the list
    ix = 1
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          l = box_dist(j, 2) ** 2 + box_dist(k, 3) ** 2
          if(l < cut2) then
             xlen = ceiling(sqrt(cut2 - l) / unit_axes(1))
             subcell_neighbour(2, ix) = j
             subcell_neighbour(3, ix) = k
             subcell_xlen(ix) = xlen ! this subpillar contains subcell [-xlen .. xlen]; note the end of region is inclusive
             ix = ix + 1
          endif
       end do
    end do
    deallocate(box_dist)
  end subroutine set_block_info

  subroutine blockify(natom, atomlist, blk)
    implicit none
    integer, intent(in) :: natom, atomlist(:)
    integer, intent(out) :: blk(:, :)
    integer :: i, j, a, blktmp(3)

    do i = 1, natom
       a = atomlist(i)
       blktmp(:) = int(floor(sitepos_normal(:, a) * invcell_normal(:) * block_size(:)))
       do j = 1, 3
          blk(j, i) = modulo(blktmp(j), block_size(j))
          if(blk(j,i) < 0) then
             print *, cell_len_normal(:), blktmp(:), block_size(:)
             STOP "INVLBLK"
          endif
       end do
    end do
  end subroutine blockify

  subroutine sort_block(blk, nmol, belong, atomno, counts, psum)
    implicit none
    integer, intent(inout) :: blk(:, :)
    integer, intent(in) :: nmol
    integer, intent(inout) :: belong(:)
    integer, intent(inout) :: atomno(:)
    integer, intent(inout) :: counts(0:block_size(1) - 1, 0:block_size(2) - 1, 0:block_size(3) - 1)
    integer, intent(out) :: psum(0:block_size(1) * block_size(2) * block_size(3))
    integer, allocatable :: buffer(:, :) ! FIXME: ugly!
    integer, allocatable :: pnum(:, :, :)
    integer :: a, b, c
    integer :: i, j, k, partialsum, pos
    
    counts(:, :, :) = 0

    do i = 1, nmol
       a = blk(1, i)
       b = blk(2, i)
       c = blk(3, i)
       if(a < 0 .or. b < 0 .or. c < 0) STOP "INVL"
       counts(a, b, c) = counts(a, b, c) + 1
    end do

    allocate(pnum(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))

    partialsum = 0
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          do i = 0, block_size(1) - 1
             pnum(i, j, k) = partialsum
             partialsum = partialsum + counts(i, j, k)
          end do
       end do
    end do

    allocate(buffer(5, nmol))
    do i = 1, nmol
       a = blk(1, i)
       b = blk(2, i)
       c = blk(3, i)
       pos = pnum(a, b, c) + 1  ! buffer (and output) index starts from 1
       pnum(a, b, c) = pos
       buffer(1:3, pos) = blk(1:3, i)
       buffer(4, pos) = belong(i)
       buffer(5, pos) = atomno(i)
    end do
    deallocate(pnum)

    blk(1:3, :) = buffer(1:3, :)
    belong(:) = buffer(4, :)
    atomno(:) = buffer(5, :)

    deallocate(buffer)

    partialsum = 0
    pos = 0
    do k = 0, block_size(3) - 1
       do j = 0, block_size(2) - 1
          do i = 0, block_size(1) - 1
             psum(pos) = partialsum + 1 ! output index starts from 1
             pos = pos + 1
             partialsum = partialsum + counts(i, j, k)
          end do
       end do
    end do
    psum(pos) = partialsum + 1
    ! print *, psum
  end subroutine sort_block

  subroutine get_pair_energy(energy_vec)
    ! calculate for each subcell
    ! cut-off by subcell distance
    implicit none
    real, intent(inout) :: energy_vec(:, :)
    integer :: u1, u2, u3
    integer :: vbs(3)
    integer :: i, upos, vpos_base, vpos_line_end, vpos_begin, vpos_end
    integer :: xlen

    !$omp parallel &
    !$omp   private(u1, u2, u3, upos, i, vbs, vpos_begin, vpos_base, vpos_end, vpos_line_end, xlen) &
    !$omp   shared(energy_vec)
    !$omp single
    do u3 = 0, block_size(3) - 1
       do u2 = 0, block_size(2) - 1
          do u1 = 0, block_size(1) - 1
             upos = u1 + block_size(1) * (u2 + block_size(2) * u3)
             if(psum_solu(upos + 1) /= psum_solu(upos)) then ! if solute have atoms in the block
                !$omp task
                do i = 1, subcell_num_neighbour
                   vbs(2) = mod(u2 + subcell_neighbour(2, i) , block_size(2))
                   vbs(3) = mod(u3 + subcell_neighbour(3, i) , block_size(3))
                   vpos_base = block_size(1) * (vbs(2) + block_size(2) * vbs(3))
                   
                   xlen = subcell_xlen(i)
                   if(2 * xlen + 1 >= block_size(1)) then
                      ! spans all over x-axis
                      call get_pair_energy_block(upos, vpos_base, vpos_base + block_size(1), energy_vec)
                   else
                      vpos_begin = vpos_base + u1 - xlen
                      vpos_end = vpos_base + u1 + xlen + 1
                      vpos_line_end = vpos_base + block_size(1)

                      if(vpos_begin < vpos_base) then
                         ! spans over periodic boundary, case 1
                         call get_pair_energy_block(upos, vpos_begin + block_size(1), vpos_line_end, energy_vec)
                         call get_pair_energy_block(upos, vpos_base, vpos_end, energy_vec)
                      elseif(vpos_end > vpos_line_end) then
                         ! spans over periodic boundary, case 2
                         call get_pair_energy_block(upos, vpos_begin, vpos_line_end, energy_vec)
                         call get_pair_energy_block(upos, vpos_base, vpos_end - block_size(1), energy_vec)
                      else
                         ! standard case
                         call get_pair_energy_block(upos, vpos_begin, vpos_end, energy_vec)
                      endif
                   endif
                end do
                !$omp end task
             end if
          end do
       end do
    end do
    !$omp end single
    !$omp end parallel
  end subroutine get_pair_energy

  ! Computational kernel to calculate distance between particles
  subroutine get_pair_energy_block(upos, vpos_b, vpos_e, energy_vec)
    use engmain, only: cltype, boxshp, &
         upljcut, lwljcut, elecut, screen, charge,&
         ljswitch, ljtype, ljtype_max, ljene_mat, ljlensq_mat, &
         SYS_NONPERIODIC, EL_COULOMB, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX
    !$ use omp_lib, only: omp_get_thread_num
    implicit none
    integer, intent(in) :: upos, vpos_b, vpos_e
    real, intent(inout) :: energy_vec(:, :)
    integer :: ui, vi, ua, va, i, curp
    integer :: n_lowlj, n_switch, n_el
    integer :: belong_u, belong_v, ljtype_u, ljtype_v
    real :: crdu(3), crdv(3), d(3), dist, r, dist_next, invr2, invr3, invr6
    real :: lwljcut2, upljcut2, lwljcut3, upljcut3, lwljcut6, upljcut6
    real :: ljeps, ljsgm2, ljsgm3, ljsgm6, vdwa, vdwb, swfac
    real :: repA, repB, repC, attA, attB, attC
    real :: elecut2, half_cell(3)
    
    if(cltype == EL_COULOMB) stop "realcal%get_pair_energy_block: cltype assertion failure"
    if(boxshp == SYS_NONPERIODIC) stop "realcal%get_pair_energy_block: boxshp assertion failure"

    curp = 1
    !$ curp = omp_get_thread_num() + 1

    half_cell(:) = 0.5 * cell_len_normal(:)

    n_lowlj = 0
    n_switch = 0
    n_el = 0

    lwljcut2 = lwljcut ** 2
    upljcut2 = upljcut ** 2
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

    elecut2 = elecut ** 2

    ! TODO optimize:
    ! if you sort / reorder atomno, ljtype etc.
    ! this loop can be vectorized
    ! TODO optimize:
    ! software pipelining
    do ui = psum_solu(upos), psum_solu(upos + 1) - 1
       ua = atomno_solu(ui)
       belong_u = belong_solu(ui) ! FIXME: not used in later calculation
       ljtype_u = ljtype(ua)
       crdu(:) = sitepos_normal(:, ua)

       ! hide latency by calculating distance of next coordinate set
       d(:) = crdu(:) - sitepos_solv(:, psum_solv(vpos_b))
       if(is_cuboid) then
          d(:) = half_cell(:) - abs(half_cell(:) - abs(d(:)))
       else
          d(:) = d(:) - matmul(cell_normal, anint(matmul(invcell, d)))
       end if

       dist_next = sum(d(:) ** 2)
       
       do vi = psum_solv(vpos_b), psum_solv(vpos_e) - 1
          va = atomno_solv(vi)
          belong_v = belong_solv(vi)
          ljtype_v = ljtype(va)
          
          crdv(:) = sitepos_solv(:, vi + 1)

          d(:) = crdv(:) - crdu(:)
          if(is_cuboid) then
             d(:) = half_cell(:) - abs(half_cell(:) - abs(d(:))) ! get nearest image
          else
             d(:) = d(:) - matmul(cell_normal, anint(matmul(invcell, d)))
          end if

          ! assumes that only a single image matters for both electrostatic and LJ.
          ! if the box is very small and strongly anisotropic,
          ! there is a risk that second nearest image still being inside the cutoff length.
          ! But it's not considered in this case ...
          
          dist = dist_next
          dist_next = sum(d(:) ** 2) ! CHECK: any sane compiler will expand and unroll

          ! lines up all variables, to enable vectorization in 2nd phase
          ljeps = ljene_mat(ljtype_v, ljtype_u)
          if(ljeps > 0) then
             if(dist <= lwljcut2) then
                n_lowlj = n_lowlj + 1
                ljeps_lowlj (n_lowlj, curp) = ljeps
                ljsgm2_lowlj(n_lowlj, curp) = ljlensq_mat(ljtype_v, ljtype_u)
                dist_lowlj(n_lowlj, curp) = dist
                belong_lowlj(n_lowlj, curp) = belong_v
             elseif(dist <= upljcut2) then
                n_switch = n_switch + 1
                ljeps_switch (n_switch, curp) = ljeps
                ljsgm2_switch(n_switch, curp) = ljlensq_mat(ljtype_v, ljtype_u)
                dist_switch(n_switch, curp) = dist
                belong_switch(n_switch, curp) = belong_v
             end if
          end if
          
          if(dist <= elecut2) then
             n_el = n_el + 1
             charge_el(n_el, curp) = charge(ua) * charge(va)
             dist_el(n_el, curp) = dist
             belong_el(n_el, curp) = belong_v
          end if
       end do
    end do

    ! 2nd phase: calculate actual values with vectorized loop

    ! TODO optimize:
    ! explicit vectorization may increase performance

    ! LJ inside low cutoff
    do i = 1, n_lowlj
       ljeps = ljeps_lowlj(i, curp)
       ljsgm2 = ljsgm2_lowlj(i, curp)
       dist = dist_lowlj(i, curp)
       invr2 = ljsgm2 / dist
       invr6 = invr2 * invr2 * invr2
       select case(ljswitch)
       case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
          e_t(i, curp) = 4.0 * ljeps * invr6 * (invr6 - 1.0)
       case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
          ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
          vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 / (lwljcut6 * upljcut6)
          vdwb = invr6 - ljsgm6 / (lwljcut3 * upljcut3)
          e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
       case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
          ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
          vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * repC
          vdwb = invr6 - ljsgm6 * attC
          e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
       case default
          stop "Unknown ljswitch"
       end select
    end do
    do i = 1, n_lowlj
       energy_vec(belong_lowlj(i, curp), curp) = energy_vec(belong_lowlj(i, curp), curp) + e_t(i, curp)
    end do

    ! LJ switching region
    do i = 1, n_switch
       ljeps = ljeps_switch(i, curp)
       ljsgm2 = ljsgm2_switch(i, curp)
       dist = dist_switch(i, curp)
       invr2 = ljsgm2 / dist
       invr6 = invr2 * invr2 * invr2
       select case(ljswitch)
       case(LJSWT_POT_CHM)                   ! potential switch (CHRAMM type)
          e_t(i, curp) = 4.0 * ljeps * invr6 * (invr6 - 1.0)             &
                       * (2.0 * dist + upljcut2 - 3.0 * lwljcut2)        &
                       * ((dist - upljcut2) ** 2) / ((upljcut2 - lwljcut2) ** 3)
       case(LJSWT_POT_GMX)                   ! potential switch (GROMACS type)
          swfac = (sqrt(dist) - lwljcut) / (upljcut - lwljcut)
          e_t(i, curp) = 4.0 * ljeps * invr6 * (invr6 - 1.0)             &
                       * (1.0 - 10.0 * (swfac ** 3)                      &
                              + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5) )
       case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
          invr3 = sqrt(invr6)
          ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
          ljsgm3 = sqrt(ljsgm6)
          vdwa = upljcut6 / (upljcut6 - lwljcut6)          &
               * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
          vdwb = upljcut3 / (upljcut3 - lwljcut3)          &
               * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
          e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
       case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
          ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
          swfac = sqrt(dist) - lwljcut
          vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * &
                         (repA * (swfac ** 3) + repB * (swfac ** 4) + repC)
          vdwb = invr6 - ljsgm6 * &
                         (attA * (swfac ** 3) + attB * (swfac ** 4) + attC)
          e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
       case default
          stop "Unknown ljswitch"
       end select
    end do
    do i = 1, n_switch
       energy_vec(belong_switch(i, curp), curp) = energy_vec(belong_switch(i, curp), curp) + e_t(i, curp)
    end do

    ! ewald electrostatic
    do i = 1, n_el
       r = sqrt(dist_el(i, curp))
       e_t(i, curp) = charge_el(i, curp) * (1.0 - erf(screen * r)) / r
    end do
    do i = 1, n_el
       energy_vec(belong_el(i, curp), curp) = energy_vec(belong_el(i, curp), curp) + e_t(i, curp)
    end do
  end subroutine get_pair_energy_block

  ! Rotate box and coordinate so that the cell(:,:) is upper triangular:
  ! cell(2, 1) = cell(3, 1) = cell(3, 2) = 0
  ! (1st axis to be aligned to x-axis, 2nd axis to be within xy-plane)
  subroutine normalize_periodic
    use engmain, only: cell, sitepos
    implicit none
    integer :: n, i, perm(3), info, lwork
    real :: qr(3,3), scale(3), newcell(3, 3)
    real, allocatable :: work(:)
    real :: dummy

    n = size(sitepos, 2)
    lwork = max(3 * 3 + 1, n)
    allocate(work(lwork))

    ! QR-factorize box vector
    qr(:, :) = cell
    perm(:) = 0

    select case(kind(dummy))
    case(8)
       call dgeqp3(3, 3, qr, 3, perm, scale, work, lwork, info)
    case(4)
       call sgeqp3(3, 3, qr, 3, perm, scale, work, lwork, info)
    case default
       stop "The libraries are used only at real or double precision"
    end select
    if(info /= 0) stop "normalize_periodic: failed to factorize box vector"

    ! reorganize R
    ! AP = QR    <=>  Q^T AP = R
    newcell(:, :) = cell(:, perm(:))

    select case(kind(dummy))
    case(8)
       call dormqr('L', 'T', 3, 3, 3, qr, 3, scale, newcell, 3, work, lwork, info)
    case(4)
       call sormqr('L', 'T', 3, 3, 3, qr, 3, scale, newcell, 3, work, lwork, info)
    end select
    if(info /= 0) stop "normalize_periodic: failed to rotate cell"

    cell_normal(:, :) = newcell(:, :)
    if(abs(newcell(2, 1)) > check_rotate .or. &
       abs(newcell(3, 1)) > check_rotate .or. &
       abs(newcell(3, 2)) > check_rotate ) then
       print *, newcell
       stop "normalize_periodic: assertion failed, box rotation is bugged"
    endif

    ! rotate coordinates
    ! Note: sitepos does not need permutation. (the order of cell axis is not important)
    select case(kind(dummy))
    case(8)
       call dormqr('L', 'T', 3, n, 3, qr, 3, scale, sitepos_normal, 3, work, lwork, info)
    case(4)
       call sormqr('L', 'T', 3, n, 3, qr, 3, scale, sitepos_normal, 3, work, lwork, info)
    end select
    if(info /= 0) stop "normalize_periodic: failed to rotate coordinate"

    deallocate(work)

    do i = 1, 3
       cell_len_normal(i) = abs(cell_normal(i, i))
    end do
    invcell_normal(:) = 1 / cell_len_normal(:)

    info = 0
    invcell(:, :) = cell_normal(:, :)
    select case(kind(dummy))
    case(8)
       call dtrtri('U', 'N', 3, invcell, 3, info)
    case(4)
       call strtri('U', 'N', 3, invcell, 3, info)
    end select
    if(info /= 0) stop "normalize_periodic: failed to invert box vector"

    if(abs(cell(1, 2)) > cuboid_thres .or. &
       abs(cell(1, 3)) > cuboid_thres .or. &
       abs(cell(2, 3)) > cuboid_thres ) then
       is_cuboid = .false.
       stop "A non-cuboidal cell is not supported. Wait for ver 0.4"
    else
       is_cuboid = .true.
    end if

    ! normalize coordinates into a periodic parallelpiped cell
    do i = 1, n
       sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            matmul(cell_normal, floor(matmul(invcell, sitepos_normal(1:3, i))))
    end do

    ! move all particles inside the cuboid spanned by (0 .. cell(1, 1)), (0 .. cell(2,2)), (0 .. cell(3,3)).
    ! Z values are already within the range (only 1 axis exists within parallelpiped cell)
    ! Y values are bit tricky, shifted by the second vector
    ! X values are simply shifted along the first vector
    do i = 1, n
       sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(1:3, 2) * floor(dot_product(invcell(1:3, 2), sitepos_normal(1:3, i)))
       sitepos_normal(1, i) = sitepos_normal(1, i) - &
            cell_normal(1, 1) * floor(invcell(1, 1) * sitepos_normal(1, i))
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
