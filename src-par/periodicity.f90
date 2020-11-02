module periodicity
!
! Written for periodic flows. We map flow from one end of the domain to the other.
! These boundary regions are usually not on same process in parallel run.
! We need to arrange pairs of face indices for sender and receiver first,
! which enables mapped data to get to the right place.
! We are dealing with unstructured meshes so there is no certainty how are faces
! ordered on both ends.
!
  use types
  use parameters
  use geometry
  use variables
  use my_mpi_module
  use mpi

  implicit none

  integer, parameter :: procSender   = 3 ! Sending process
  integer, parameter :: procReceiver = 0 ! Receiving process
  integer, parameter :: asize = 27825 !14905 !21825  ! Number of faces in section, change for every case.
  integer, dimension(:), allocatable :: mappedIndx

  private

  public :: face_mapping, recirculate_flow


contains


subroutine face_mapping

  implicit none


  integer :: i,ib,iface
  integer :: iDFriend
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)
  real(dp) :: buf(2*asize)

  integer :: indx, itemp
  real(dp) :: ytemp, ztemp, mindist, dist

  ! Allocate array for mapped indices
  allocate( mappedIndx(asize) )

  ! Exchange
  if (myid.eq.procReceiver)  iDFriend = procSender
  if (myid.eq.procSender)    iDFriend = procReceiver

  sendtag = 123 + procSender + procReceiver   ! tag for sending
  rectag  = sendtag                           ! tag for receiving

  length = 2*asize

  !  syncronization before communication
  call MPI_Barrier(MPI_COMM_WORLD,ierr)


  if ( myid.eq.procSender ) then  
    do ib=1,numBoundaries
      if ( bctype(ib) == 'outlet' ) then ! NOTE, we use 'bcname' here.
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          buf(        i) = yf(iface) 
          buf(  asize+i) = zf(iface)                    
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


  if ( myid.eq.procReceiver ) then

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
            ytemp = buf(      itemp)
            ztemp = buf(asize+itemp)       
            dist = abs(yf(iface)-ytemp) + abs(zf(iface)-ztemp)
            if ( dist < mindist ) then
             mindist = dist
             indx  = itemp 
            endif
          enddo

          mappedIndx(i) = indx

        end do
      endif 
    enddo
  endif  

  if (myid == 0) write(*,'(a)') '  **Mapped faces for flow recirculation.'


end subroutine


subroutine recirculate_flow

  implicit none

  integer :: i,ib,ini,iface
  integer :: iDFriend
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)
  real(dp) :: buf(5*asize)

  ! Exchange
  if (myid.eq.procReceiver)  iDFriend = procSender
  if (myid.eq.procSender)    iDFriend = procReceiver

  sendtag = 123 + procSender + procReceiver   ! tag for sending
  rectag  = sendtag                           ! tag for receiving

  length = 5*asize

  !  syncronization before communication
  call MPI_Barrier(MPI_COMM_WORLD,ierr)


  if ( myid.eq.procSender ) then  
    do ib=1,numBoundaries
      if ( bctype(ib) == 'outlet' ) then  ! NOTE, we use 'bcname' here.
        ! asize = nfaces(ib)
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


  if ( myid.eq.procReceiver ) then

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
          u(ini)   = buf(         mappedIndx(i) ) 
          v(ini)   = buf(   asize+mappedIndx(i) ) 
          w(ini)   = buf( 2*asize+mappedIndx(i) )
          te(ini)  = buf( 3*asize+mappedIndx(i) )
          vis(ini) = buf( 4*asize+mappedIndx(i) ) 
          flmass(iface) = den(ini)*(arx(iface)*u(ini)+ary(iface)*v(ini)+arz(iface)*w(ini))

        end do
      endif 
    enddo
  endif  

  if (myid == 0) write(*,'(a)') '  **Recycled flow on periodic boundaries.'


end subroutine


end module