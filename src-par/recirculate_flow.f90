subroutine recirculate_flow

  use types
  use parameters
  use geometry
  use variables
  use my_mpi_module
  use mpi

  implicit none

  integer :: i,ib,ini,iface
  integer :: IDInlet, IDOutlet,idtmp
  integer :: asize
  integer :: iDFriend
  integer :: rectag,sendtag
  integer :: length
  integer :: status(mpi_status_size)
  real(dp) :: buf(7*lenbuf)
  logical :: IhaveInlet = .false.
  logical :: IhaveOutlet = .false.

  asize = 0
  IDOutlet= -1
  IDInlet = -1

  !  syncronization before communication
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Idi po granicama, ako imas OUTLET onda si ti taj koji treba da posalje
  ! obavesti ostale ko si preko broadcasta.
  do ib=1,numBoundaries
    if ( bctype(ib) == 'outlet' ) then

      IhaveOutlet = .true.  
      IDOutlet = myid

      asize = nfaces(ib)
      ! Javio si svoj identitet sad napuni buffer i pripremi se za komunikaciju
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ini = iBndValueStart(ib) + i
        buf(i)         = u(ini)
        buf(  asize+i) = v(ini)
        buf(2*asize+i) = w(ini)
        buf(3*asize+i) = te(ini)
        buf(4*asize+i) = ed(ini) 
        buf(5*asize+i) = vis(ini)
        buf(6*asize+i) = flmass(iface)
      end do

    endif 
  enddo

  call mpi_allreduce      &               
   (IDOutlet,             & ! send buffer
    idtmp,                & ! recv buffer 
    1,                    & ! length     
    mpi_integer,          & ! datatype  
    mpi_max,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  IDOutlet = idtmp

! print*, 'We know who is outlet:', IDOutlet

  ! Idi po granicama, ako imas INLET onda si ti taj koji treba da posalje
  ! obavesti ostale ko si preko broadcasta.
  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then

    IhaveInlet = .true.
    IDInlet = myid

    asize = nfaces(ib)
    ! Javio si svoj identitet sad napuni buffer i pripremi se za komunikaciju
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ini = iBndValueStart(ib) + i
      buf(i)         = u(ini)
      buf(  asize+i) = v(ini)
      buf(2*asize+i) = w(ini)
      buf(3*asize+i) = te(ini)
      buf(4*asize+i) = ed(ini) 
      buf(5*asize+i) = vis(ini)
      buf(6*asize+i) = flmass(iface)
      end do
    endif 
  enddo

  ! Communicate to others
  call mpi_allreduce      &               
   (IDInlet,             & ! send buffer
    idtmp,                & ! recv buffer 
    1,                    & ! length     
    mpi_integer,          & ! datatype  
    mpi_max,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  IDInlet = idtmp
  ! print*, 'We know who is inlet:', IDInlet

  ! Ucestvujes u exchange-u ako si bilo koji od ova dva
  if (IhaveInlet .or. IhaveOutlet) then

    if (IhaveInlet)  iDFriend = IDOutlet
    if (IhaveOutlet) iDFriend = IDInlet

    sendtag = 123 + MYID + iDFriend   ! tag for sending
    rectag  = sendtag                 ! tag for receiving

    length = 7*asize
    ! print*, 'packed now sending!'

    call MPI_SENDRECV_REPLACE & 
     (buf,              &     ! buffer
      length,           &     ! length   
      MPI_DOUBLE_PRECISION, & ! datatype  
      iDFriend,          &    ! dest,      
      sendtag,          &     ! sendtag,    
      iDFriend,         &     ! source,      
      rectag,           &     ! recvtag,      
      MPI_COMM_WORLD,   &     ! communicator
      status,           &     ! status
      ierr)                   ! error

  endif

  ! print*, 'Done exchange.'

  if ( IhaveOutlet ) then
    do ib=1,numBoundaries    
      if ( bctype(ib) == 'outlet' ) then
        ! print*, 'Inside Outlet.'
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          ini = iBndValueStart(ib) + i
          u(ini) = buf(i) 
          v(ini) = buf(  asize+i) 
          w(ini) = buf(2*asize+i)
          te(ini) = buf(3*asize+i)
          ed(ini)  = buf(4*asize+i) 
          vis(ini)  = buf(5*asize+i)
          flmass(iface) = buf(6*asize+i)
        end do
      endif 
    enddo
  endif  

  if ( IhaveInlet ) then
    do ib=1,numBoundaries    
      if ( bctype(ib) == 'inlet' ) then
        ! print*, 'Inside inlet.'
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          ini = iBndValueStart(ib) + i
          u(ini) = buf(i) 
          v(ini) = buf(  asize+i) 
          w(ini) = buf(2*asize+i)
          te(ini) = buf(3*asize+i)
          ed(ini)  = buf(4*asize+i) 
          vis(ini)  = buf(5*asize+i)
          ! flmass(iface) = buf(6*asize+i)
          flmass(iface) = den(ini)*(arx(iface)*u(ini)+ary(iface)*v(ini)+arz(iface)*w(ini))
        end do
      endif 
    enddo
  endif  

  if (myid == 0) write(*,'(a)') '  **Recycled flow on periodic boundaries.'

end subroutine