!***********************************************************************
!
subroutine exchange(phi) 
!
!***********************************************************************
!
!   Exchanges field values between connected processes from their
!   respective buffers. 
!   Executes a blocking send and receive. The same buffer is used both 
!   for the send and for the receive, so that the message sent is 
!   replaced by the message received.                   
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
  integer :: i,ijp,ijn,ib,ipro,iface
  integer :: iConnection,iDFriend,iStart,iEnd
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)

 
  !
  ! > Fill the buffers with new values
  !

  ! ! Moze i da se pojednostavi punjenje buffera jer je sve vec podeseno
  ! do i=1,lenbuf
  !   buffer(i) = phi( bufind(i) )
  ! enddo 

  ipro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'process' ) then

      ! Faces on process boundary
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ipro = ipro+1

        buffer( ipro ) = phi( ijp )

      enddo
    endif 

  enddo


  !  syncronization before communication
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !
  ! > Exchange the values
  !

  !
  ! > Idi po svim domenima sa kojima je ovaj konektovan.
  !
  do iConnection = 1,numConnections

    iDFriend = neighbProcNo(iConnection)

    sendtag = 123 + MYID + iDFriend   ! tag for sending
    rectag  = sendtag                 ! tag for receiving

    iStart = neighbProcOffset(iConnection)
    iEnd   = neighbProcOffset(iConnection+1)-1

    length = iEnd-iStart+1

    call MPI_SENDRECV_REPLACE & 
     (buffer(iStart),   &     ! buffer
      length,           &     ! length   
      MPI_DOUBLE_PRECISION, & ! datatype  
      iDFriend,          &    ! dest,      
      sendtag,          &     ! sendtag,    
      iDFriend,         &     ! source,      
      rectag,           &     ! recvtag,      
      MPI_COMM_WORLD,   &     ! communicator
      status,           &     ! status
      ierr)                   ! error

  end do


  ! 
  ! Move data from the buffer to variable array where it belongs.
  !

  ipro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'process' ) then

      ! Faces on process boundary
      do i=1,nfaces(ib)

        ijn = iBndValueStart(ib) + i
        ipro = ipro+1

        phi( ijn ) = buffer( ipro )

      enddo
    endif 

  enddo
 
end subroutine exchange

