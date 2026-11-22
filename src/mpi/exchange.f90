subroutine exchange(phi) 
!
!***********************************************************************
!
!   Exchanges field values between connected processes from their
!   respective buffers using NON-BLOCKING COMMUNICATION.
!
!***********************************************************************
!

  use types
  use parameters
  use geometry
  use my_mpi_module
  use mpi

  implicit none

  !
  ! > Arguments 
  !
  real(dp), intent(inout) :: phi(numTotal)

  !
  ! > Locals
  !
  integer :: i, ijp, ijn, ib, iface
  integer :: iConnection, iDFriend, iStart, iEnd
  integer :: sendtag, rectag, length
  integer :: ierr

  integer :: numRequests
  integer, allocatable :: request(:)
  integer, allocatable :: status(:,:)
  integer, allocatable :: recv_requests(:)
  integer, allocatable :: send_requests(:)

  !
  ! > Fill the buffers with new values
  !

  ipro = 0
  do ib=1,numBoundaries
    if (bctype(ib) == 'process') then ! Faces on process boundary
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ipro = ipro + 1
        buffer(ipro) = phi(ijp)
      enddo
    endif 
  enddo

  !
  ! > Exchange the values using non-blocking communication
  !
  allocate(request(2 * numConnections))
  allocate(status(MPI_STATUS_SIZE, 2 * numConnections))
  allocate(recv_requests(numConnections))
  allocate(send_requests(numConnections))

  numRequests = 0

  ! Post non-blocking receives.
  ! Go trough all the connections this domains has - there is numConnections of them.
  do iConnection = 1, numConnections 
    iDFriend = neighbProcNo(iConnection)
    rectag   = 200 + MYID + iDFriend  ! Unique receive tag

    iStart = neighbProcOffset(iConnection)   
    iEnd   = neighbProcOffset(iConnection + 1) - 1
    length = iEnd - iStart + 1

    numRequests = numRequests + 1
    call MPI_IRECV(buffer(iStart), length, MPI_DOUBLE_PRECISION, &
                   iDFriend, rectag, MPI_COMM_WORLD, recv_requests(numRequests), ierr)
  end do

  ! Post non-blocking sends.
  ! Go trough all the connections this domains has - there is numConnections of them.
  do iConnection = 1, numConnections
    iDFriend = neighbProcNo(iConnection)
    sendtag  = 100 + MYID + iDFriend  ! Unique send tag

    iStart = neighbProcOffset(iConnection)
    iEnd   = neighbProcOffset(iConnection + 1) - 1
    length = iEnd - iStart + 1

    numRequests = numRequests + 1
    call MPI_ISEND(buffer(iStart), length, MPI_DOUBLE_PRECISION, &
                   iDFriend, sendtag, MPI_COMM_WORLD, send_requests(numRequests), ierr)
  end do

  ! Wait for all sends and receives to complete
  call MPI_WAITALL(2 * numConnections, request, status, ierr)

  !
  ! Move data from the buffer to variable array where it belongs.
  !
  ipro = 0
  do ib = 1, numBoundaries
    if (bctype(ib) == 'process') then
      do i = 1, nfaces(ib)
        ijn = iBndValueStart(ib) + i
        ipro = ipro + 1
        phi(ijn) = buffer(ipro)
      enddo
    endif 
  enddo

  ! Clean up
  deallocate(request, status, recv_requests, send_requests)

end subroutine exchange
