subroutine recirculate_flow

  use types
  use parameters
  use geometry
  use variables
  use my_mpi_module
  use mpi

  implicit none

  integer :: i,ib,ini,iface
  integer :: asize
  integer :: iDFriend
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)
  real(dp) :: buf(5*lenbuf)


  asize = 4089

  ! Exchange
  if (myid.eq.0)  iDFriend = 3
  if (myid.eq.3)  iDFriend = 0

  sendtag = 123 + MYID + iDFriend   ! tag for sending
  rectag  = sendtag                 ! tag for receiving

  length = 5*asize

  !  syncronization before communication
  call MPI_Barrier(MPI_COMM_WORLD,ierr)


  if ( myid.eq.3 ) then  
    do ib=1,numBoundaries
      if ( bctype(ib) == 'outlet' ) then
        asize = nfaces(ib)
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          ini = iBndValueStart(ib) + i
          buf(i)         = u(ini)
          buf(  asize+i) = v(ini)
          buf(2*asize+i) = w(ini)
          buf(3*asize+i) = te(ini)
          buf(4*asize+i) = vis(ini) 
        end do
      endif 
    enddo

    call MPI_SEND(          &
      buf,                  &
      length,               &
      MPI_DOUBLE_PRECISION, &
      iDFriend,             &
      sendtag,              &
      MPI_COMM_WORLD,       &
      ierr)

  endif


  if ( myid.eq.0 ) then

    call MPI_RECV(      & 
      buf,              &     ! buffer
      length,           &     ! length   
      MPI_DOUBLE_PRECISION, & ! datatype     
      iDFriend,         &     ! source,      
      rectag,           &     ! recvtag,      
      MPI_COMM_WORLD,   &     ! communicator
      status,           &     ! status
      ierr)                   ! error

    do ib=1,numBoundaries    
      if ( bctype(ib) == 'inlet' ) then
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          ini = iBndValueStart(ib) + i
          u(ini) = buf(i) 
          v(ini) = buf(  asize+i) 
          w(ini) = buf(2*asize+i)
          te(ini) = buf(3*asize+i)
          vis(ini) = buf(4*asize+i) 
          flmass(iface) = den(ini)*(arx(iface)*u(ini)+ary(iface)*v(ini)+arz(iface)*w(ini))
        end do
      endif 
    enddo
  endif  

  if (myid == 0) write(*,'(a)') '  **Recycled flow on periodic boundaries.'


end subroutine