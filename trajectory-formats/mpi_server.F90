! -*- F90 -*-
! Trajectory module based on 
module trajectory

  type handle
     integer :: comm
     integer :: clientrank
  end type handle

  character(len=80) :: service_name
  character(len=*), parameter :: mpi_server_conffile = "MpiConn"
  integer :: total_num_open = 0

contains

  subroutine init_trajectory()
    implicit none
    integer :: ioconf, ioerr
    
    ! default service name
    service_name = "ermod-trajio"
    open(newunit = ioconf, file = mpi_server_conffile, action = "READ", iostat = ioerr)
    if(ioerr == 0) then ! if failed to open, use default
       do
          read(ioconf, iostat = ioerr) service_name
          if(ioerr /= 0) exit
          exit
       end do
       close(ioconf)
    end if

  end subroutine init_trajectory

  subroutine finish_trajectory()
    implicit none
  end subroutine finish_trajectory

  ! Open trajectory and returns handle as htraj. 
  ! Should open fail, the program abends.
  ! MPI vcersion can only be called once. (Hence cannot be used for refs)
  subroutine open_trajectory(htraj, fname)
    use mpi
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname
    character(len=MPI_MAX_PORT_NAME) :: port_name
    integer :: ierr
    integer :: csize, myrank
    
    total_num_open = total_num_open + 1
    if(total_num_open > 1) stop "Error: MPI trajctory module, 'open_trajectory' is called more than once" 

    ! Open MPI port
    call MPI_Open_port(MPI_INFO_NULL, port_name, ierr)
    if(ierr /= 0) stop "Error: failed to open MPI port"
    ! Publish ports with a service name
    call MPI_Publish_name(service_name, MPI_INFO_NULL, port_name, ierr)
    if(ierr /= 0) stop "Error: failed to publish MPI port"

    ! Accept connections
    call MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, htraj%comm, ierr)
    if(ierr /= 0) stop "Error: failed to accept communicator"

    ! Close after accepting port.
    call MPI_Close_port(port_name, ierr)

    ! ... and unpublish.
    call MPI_Unpublish_name(service_name, MPI_INFO_NULL, port_name, ierr)
    
    ! Assumes that client connects to server with MPI_COMM_SELF. Check consistency and emits warning..
    call MPI_Comm_Size(htraj%comm, csize, ierr)
    if(csize /= 2) then
       write(6, *) "Warning: MPI server trajectory: communicator size is expected to be 2, but was ", csize
    end if

    call MPI_Comm_Rank(htraj%comm, myrank, ierr)
    if(myrank /= 0) then
       write(6, *) "Warning: MPI server trajectory: rank is expected to be 0, but was ", myrank
    endif
    htraj%clientrank = 1 - myrank
    
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    use mpi
    implicit none
    type(handle), intent(inout) :: htraj
    integer :: ierr

    call MPI_Comm_disconnect(htraj%comm, ierr)
    call MPI_Comm_free(htraj%comm, ierr)

  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, cell, status)
    use mpi
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real(8), intent(out) :: crd(3,natom)
    real(8), intent(out) :: cell(3,3)
    integer, intent(out) :: status
    integer :: ierr

    do
       if(is_periodic) then
          call MPI_Recv(cell, 3 * 3, MPI_DOUBLE_PRECISION, &
               htraj%clientrank, MPI_ANY_TAG, &
               htraj%comm, ierr)
          if(ierr /= 0) exit
       end if
          
       ! Note for future implementers: you may rewrite here by MPI_FLOAT to reduce bandwidth.
       call MPI_Recv(crd, 3 * natom, MPI_DOUBLE_PRECISION, &
            htraj%clientrank, MPI_ANY_TAG, &
            htraj%comm, ierr)
       if(ierr /= 0) exit

       status = 0
       return
    end do

    status = 1
    return
  end subroutine read_trajectory

end module trajectory
