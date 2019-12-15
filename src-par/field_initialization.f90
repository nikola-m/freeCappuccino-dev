module field_initialization
!
! Initializes fields by reading appropriate file from 0/ folder in the case directory.
! Initialization amounts to initializing inner field and boundary field.
! Inner field may be initialized by uniform values or nonuniform list of values.
! Boundary field is initialized following the rules defined by the BC type:
! Dirichlet (fixed value), Neumann (fixed gradient, including zero gradient Neumann), Robin (mixed) and Cauchy. 
!

use types
use geometry
use utils, only: get_unit
use my_mpi_module

implicit none

!   interface initialize
!     module procedure initialize_vector_field
!     module procedure initialize_scalar_field
!   end interface

! private

! public :: initialize

contains


subroutine initialize_vector_field(u,v,w, dUdxi, field_name)
!
! Reads 'field_name' file form 0/ folder
!
implicit none

  real(dp), dimension(numTotal) :: u,v,w
  real(dp), dimension(3,numTotal) :: dUdxi
  character(len=*) :: field_name

  !
  ! Locals
  !
  integer :: i,ib,inp,ijp,ijb,iface
  integer :: input_status, input_unit
  character( len=80 ) :: key,field_type
  character( len=15 ) :: name,type,distribution
  character( len=5 ) :: nproc_char
  real(dp) :: u0,v0,w0

  if(myid .eq. 0) then

  write(*,'(a)') ' ' 
  write(*,'(a)') '  Initializing internal field and boundaries (reading 0/'//trim(field_name)//' ):'
  write(*,'(a)') ' '

  endif

  ! 
  ! > Initialize vector field, e.g. Velocity

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

  call get_unit ( input_unit )
  ! open ( unit = input_unit, file = 'processor'//adjustl(nproc_char)//'/0/U')
  open ( unit = input_unit, file = 'processor'//trim(nproc_char)//'/0/'//trim(field_name))  
  if(myid .eq. 0) write(*,'(a)') '  0/'//trim(field_name)
  rewind input_unit

  do
  read(input_unit,'(a)',iostat = input_status) key
  if(input_status /= 0) exit

    if( key=='internalField' ) then

      if(myid .eq. 0) write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if( field_type=='uniform' ) then

          read(input_unit,*) u0,v0,w0

          if(myid .eq. 0) write(*,'(4x,a,3e15.7)') 'uniform',u0,v0,w0

          do inp = 1,numCells
            u(inp) = u0
            v(inp) = v0
            w(inp) = w0
          enddo

      elseif( field_type=='nonuniform' ) then

          if(myid .eq. 0) write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) u(inp),v(inp),w(inp)
          enddo

      endif

    elseif( key =='boundaryField') then

      if(myid .eq. 0) write(*,'(2x,a)') 'boundaryField'      

      !
      !*** Initial values for boundary faces and prescription of the BC type 
      !***   which determines how will boundary face values be treated during the code run!
      !

      do ib=1,numBoundaries - size(neighbProcNo) ! values at process boundary are not set here

      read(input_unit,*) name   ! repeat bcname(ib)
        if(myid .eq. 0) write(*,'(2x,a)') name
        read(input_unit,*) type ! eg. Dirichlet, Neumann, Robin, Cauchy, empty
          if(myid .eq. 0) write(*,'(4x,a)') type

          if ( type == 'Dirichlet' ) then
            read(input_unit,*) distribution
              if ( distribution == 'uniform' ) then
                read(input_unit,*) u0,v0,w0
                if(myid .eq. 0) write(*,'(8x,a,3e15.7)') 'uniform',u0,v0,w0  
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  u(ijb) = u0
                  v(ijb) = v0
                  w(ijb) = w0
                end do        
              else ! nonuniform
                if(myid .eq. 0) write(*,'(8x,a)') 'nonuniform'
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  read(input_unit,*) u(ijb),v(ijb),w(ijb)
                  ! write(*,'(8x,3e15.7)') u(ijb),v(ijb),w(ijb)
                end do
              endif
          endif

          if ( type == 'Neumann' ) then
            read(input_unit,*) distribution
              if(myid .eq. 0) write(*,'(6x,a)') distribution
              if ( distribution == 'uniform' ) then
                read(input_unit,*) u0,v0,w0
                if(myid .eq. 0) write(*,'(8x,3e15.7)') u0,v0,w0  
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  dUdxi(1,ijb) = u0
                  dUdxi(2,ijb) = v0
                  dUdxi(3,ijb) = w0
                end do    
              elseif ( distribution == 'zeroGradient' ) then
                do i=1,nfaces(ib)
                  iface = startFace(ib) + i
                  ijp = owner(iface)
                  ijb = iBndValueStart(ib) + i
                  u(ijb) = u(ijp)
                  v(ijb) = v(ijp)
                  w(ijb) = w(ijp)
                end do       
              else ! nonuniform
                if(myid .eq. 0) write(*,'(8x,a,3e15.7)') 'nonuniform'
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  read(input_unit,*) dUdxi(1,ijb),dUdxi(2,ijb),dUdxi(3,ijb)
                end do
              endif
          endif

        enddo

    endif  

  enddo
  close ( input_unit )

end subroutine


subroutine initialize_scalar_field(T, dTdxi, field_name)
!
! Reads 'field_name' file form 0/ folder
!
implicit none

  real(dp), dimension(numTotal) :: T
  real(dp), dimension(numTotal) :: dTdxi
  character(len=*) :: field_name

  !
  ! Locals
  !
  integer :: i,ib,inp,ijp,ijb,iface
  integer :: input_status, input_unit
  character( len=80 ) :: key,field_type
  character( len=15 ) :: name,type,distribution
  character( len=5 ) :: nproc_char
  real(dp) :: t0
 
  if(myid .eq. 0) then

  write(*,'(a)') ' '
  write(*,'(a)') '  Initializing internal field and boundaries (reading 0/'//trim(field_name)//' ):'
  write(*,'(a)') ' '

  endif

  ! 
  ! > Initialize scalar field


  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

  call get_unit ( input_unit )
  open ( unit = input_unit, file = 'processor'//adjustl(nproc_char)//'/0/'//trim(field_name))
  if(myid .eq. 0) write(*,'(a)') '  0/'//trim(field_name)
  rewind input_unit

  do
  read(input_unit,'(a)',iostat = input_status) key
  if(input_status /= 0) exit

    if( key=='internalField' ) then

      if(myid .eq. 0) write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if( field_type=='uniform' ) then

          read(input_unit,*) t0

          if(myid .eq. 0) write(*,'(4x,a,3e15.7)') 'uniform',t0

          do inp = 1,numCells
            T(inp) = t0
          enddo

      elseif( field_type=='nonuniform' ) then

          if(myid .eq. 0) write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) T(inp)
          enddo

      endif

    elseif( key =='boundaryField') then

      if(myid .eq. 0) write(*,'(2x,a)') 'boundaryField'      

      !
      !*** Initial values for bounary faces and prescription of the BC type 
      !***   which determines how will boundary face values be treated during the code run!
      !

      do ib=1,numBoundaries-size(neighbProcNo) ! values at process boundary are not set here

      read(input_unit,*) name   ! repeat bcname(ib)
        if(myid .eq. 0) write(*,'(2x,a)') name
        read(input_unit,*) type ! eg. Dirichlet, Neumann, Robin, Cauchy, empty
          if(myid .eq. 0) write(*,'(4x,a)') type

          if ( type == 'Dirichlet' ) then
            read(input_unit,*) distribution
              if(myid .eq. 0) write(*,'(6x,a)') distribution
              if ( distribution == 'uniform' ) then
                read(input_unit,*) t0
                if(myid .eq. 0) write(*,'(8x,a,3e15.7)') 'uniform',t0  
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  T(ijb) = t0
                end do        
              else ! nonuniform
                if(myid .eq. 0) write(*,'(8x,a,3e15.7)') 'nonuniform'
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  read(input_unit,*) T(ijb)
                  ! write(*,'(8x,e15.7)') T(ijb)
                end do
              endif
          endif

          if ( type == 'Neumann' ) then
            read(input_unit,*) distribution
              if(myid .eq. 0) write(*,'(6x,a)') distribution
              if ( distribution == 'uniform' ) then
                read(input_unit,*) t0
                if(myid .eq. 0) write(*,'(8x,a,3e15.7)') 'uniform',t0
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  dTdxi(ijb) = t0
                end do    
              elseif ( distribution == 'zeroGradient' ) then
                do i=1,nfaces(ib)
                  iface = startFace(ib) + i
                  ijp = owner(iface)
                  ijb = iBndValueStart(ib) + i
                  T(ijb) = T(ijp)
                end do       
              else ! nonuniform
                if(myid .eq. 0) write(*,'(8x,a,3e15.7)') 'nonuniform'
                do i=1,nfaces(ib)
                  ijb = iBndValueStart(ib) + i
                  read(input_unit,*) dTdxi(ijb)
                end do
              endif
          endif

        enddo

    endif  

  enddo
  close ( input_unit )

end subroutine



end module
