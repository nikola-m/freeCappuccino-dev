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


  integer :: i,ijn,ib,ipro
  integer :: iDomain,iDFriend,iStart,iEnd
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)

  real(dp), intent(inout) :: phi(numTotal)
 
  !
  ! > Fill the buffers with new values
  !

  ! Moze i da se pojednostavi punjenje buffera jer je sve vec podeseno
  do i=1,lenbuf
    buffer(i) = phi( bufind(i) )
  enddo 


! >> Exchange the values

  ! Idi po svim domenima sa kojima je ovaj konektovan
  do iDomain = 1,numConnections

    iDFriend = neighbProcNo(iDomain)

    sendtag = 123 + this + iDFriend ! tag for sending
    rectag  = sendtag               ! tag for receiving

    iStart = neighbProcOffset(iDomain)
    iEnd   = neighbProcOffset(iDomain+1)-1

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

  ! Prebaci iz buffera u phi polje na odgovarajuce mesto 
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


!***********************************************************************
!
subroutine exchange__(phi) 
!
!***********************************************************************
!
!   Exchanges field values between connected processes from their
!   respective buffers. 
!   Executes a blocking send and receive. The same buffer is used both 
!   for the send and for the receive, so that the message sent is 
!   replaced by the message received.
!   !!!
!   This version of the exchange routine exchanges arrays of size (numPcells).                   
!   !!!
!***********************************************************************
!

  use types
  use parameters
  use geometry
  use my_mpi_module
  use mpi

  implicit none



  integer :: i
  integer :: iDomain,iDFriend,iStart,iEnd
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)

  real(dp), intent(inout) :: phi(numPCells)
 
  !
  ! >> Fill the buffers with new values
  !

  ! Moze i da se pojednostavi punjenje buffera jer je sve vec podeseno
  do i=1,lenbuf
    buffer( i ) = phi( bufind(i) )
  enddo 


! >> Exchange the values

  ! Idi po svim domenima sa kojima je ovaj konektovan
  do iDomain = 1,numConnections

    iDFriend = neighbProcNo(iDomain)

    sendtag = 123 + this + iDFriend ! tag for sending
    rectag  = sendtag ! tag for receiving

    iStart = neighbProcOffset(iDomain)
    iEnd   = neighbProcOffset(iDomain+1)-1

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

  ! Prebaci iz buffera u phi polje na odgovarajuce mesto
  ! Paznja: ovde je kraci niz pa vrednosti idu odmah nakon numCells mesta u nizu
  do i=1,lenbuf
    phi( numCells + i ) = buffer(i)
  enddo
 
end subroutine exchange__