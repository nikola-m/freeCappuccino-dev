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
  real(dp) :: buf(7*lenbuf)
  integer :: indx, itemp
  real(dp) :: ytemp, ztemp, mindist, dist

  asize = 4089

  ! Exchange
  if (myid.eq.0)  iDFriend = 3
  if (myid.eq.3)  iDFriend = 0

  sendtag = 123 + MYID + iDFriend   ! tag for sending
  rectag  = sendtag                 ! tag for receiving

  length = 7*asize

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
          buf(5*asize+i) = yf(iface) 
          buf(6*asize+i) = zf(iface)            
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

!
! OK sad smo kod tog nekog fejsa, treba da vidimo koji mu je poslati podatak 'predodredjen'
! u tom smilsu da trazimo neku distancu po xf i yf i kada nadjemo poklapanje taj je indeks 
!
          mindist = 1e+30
          indx = 0
          do itemp=1,asize
            ytemp = buf(5*asize+itemp)
            ztemp = buf(6*asize+itemp)       
            dist = abs(yf(iface)-ytemp) + abs(zf(iface)-ztemp)
            if ( dist < mindist ) then
             mindist = dist
             indx  = itemp 
            endif
          enddo

          ini = iBndValueStart(ib) + i
          u(ini) = buf(indx) 
          v(ini) = buf(  asize+indx) 
          w(ini) = buf(2*asize+indx)
          te(ini) = buf(3*asize+indx)
          vis(ini) = buf(4*asize+indx) 
          flmass(iface) = den(ini)*(arx(iface)*u(ini)+ary(iface)*v(ini)+arz(iface)*w(ini))
        end do
      endif 
    enddo
  endif  

  if (myid == 0) write(*,'(a)') '  **Recycled flow on periodic boundaries.'


end subroutine