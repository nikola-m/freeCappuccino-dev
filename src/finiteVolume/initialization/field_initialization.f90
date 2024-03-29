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
use utils,    only: get_unit

implicit none

!   interface initialize
!     module procedure initialize_vector_field
!     module procedure initialize_scalar_field
!   end interface

! private

! public :: initialize

contains


subroutine initialize_vector_field(u,v,w, dUdxi, iphi, field_name)
!
! Reads 'field_name' file form 0/ folder
!

implicit none

  real(dp), dimension(numTotal) :: u,v,w
  real(dp), dimension(3,numTotal) :: dUdxi
  integer, intent(in) :: iphi
  character(len=*) :: field_name

  !
  ! Locals
  !
  integer :: i,ib,inp,ijp,ijb,iface
  integer :: input_status, input_unit
  character(len=80) :: key,field_type
  character(len=15) :: name,type,distribution
  real(dp) :: u0,v0,w0

  integer :: numVals,k
  real(dp) :: ic
  real(dp), dimension(:), allocatable :: yi, ui, vi, wi 

  write(*,'(a)') ' ' 
  write(*,'(a)') '  Initializing internal field and boundaries (reading 0/'//trim(field_name)//' ):'
  write(*,'(a)') ' '

  ! 
  ! > Initialize vector field, e.g. Velocity


  call get_unit ( input_unit )
  open ( unit = input_unit, file = '0/'//trim(field_name))
  write(*,'(a)') '  0/'//trim(field_name)
  rewind input_unit

  do
  read(input_unit,'(a)',iostat = input_status) key
  if(input_status /= 0) exit

    if( key=='internalField' ) then

      write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if( field_type=='uniform' ) then

          read(input_unit,*) u0,v0,w0

          write(*,'(4x,a,3e15.7)') 'uniform',u0,v0,w0

          do inp = 1,numCells
            u(inp) = u0
            v(inp) = v0
            w(inp) = w0
          enddo

      elseif( field_type=='nonuniform' ) then

          write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) u(inp),v(inp),w(inp)
          enddo

      endif

    elseif( key =='boundaryField') then

      write(*,'(2x,a)') 'boundaryField'      

      !
      !*** Initial values for bounary faces and prescription of the BC type 
      !***   which determines how will boundary face values be treated during the code run!
      !

      do ib=1,numBoundaries

      read(input_unit,*) name   ! repeat bcname(ib)
        write(*,'(2x,a)') name
        read(input_unit,*) type ! eg. Dirichlet, Neumann, Robin, Cauchy, empty
          write(*,'(4x,a)') type

          if ( type == 'Dirichlet' ) then

            bcDFraction(iphi,ib) = 1.0_dp !!! Set Dirichlet fraction for value update to one.

            read(input_unit,*) distribution
            if ( distribution == 'uniform' ) then
              read(input_unit,*) u0,v0,w0
              write(*,'(8x,a,3e15.7)') 'uniform',u0,v0,w0  
              do i=1,nfaces(ib)
                ijb = iBndValueStart(ib) + i
                u(ijb) = u0
                v(ijb) = v0
                w(ijb) = w0
              end do   
 
            ! When you have a sample of points (for now along y-axis which is assumed to span the inlet in 2D and 3D)
            ! you can use this type of boundary condition, where values are linearly interpolated from sample distribution.
            elseif(distribution =='interpolated') then
              write(*,'(8x,a)') 'interpolated'
              read(input_unit, *) numVals ! Number of entries for values that need to be interpolated
              allocate(yi(0:numVals), ui(0:numVals), vi(0:numVals), wi(0:numVals))
              yi(0)=0.; ui(0)=0.; vi(0) = 0.; wi(0) = 0. ! Zero velocity at zero ordinate
              do i = 1, numVals
                read(input_unit, *) yi(i), ui(i), vi(i), wi(i)
              enddo
              face_loop: do i=1,nfaces(ib)
                ijb = iBndValueStart(ib) + i
                iface = startFace(ib) + i
                span_loop: do k=1,numVals ! Loop over interpolation points and find the interval where yf(iface) belongs
                  if(yf(iface) > yi(k-1) .and. yf(iface) < yi(k)) then
                    ic = (yf(iface) - yi(k-1))/(yi(k) - yi(k-1)) ! linear interpolation coefficient
                    u(ijb) = ui(k-1) + ic*(ui(k) - ui(k-1))
                    v(ijb) = vi(k-1) + ic*(vi(k) - vi(k-1))
                    w(ijb) = wi(k-1) + ic*(wi(k) - wi(k-1))
                    exit span_loop
                  endif
                enddo span_loop
                write(*,'(8x,3e15.7)') u(ijb),v(ijb),w(ijb)
              end do face_loop
              deallocate(yi,ui,vi,wi)

            else ! nonuniform
              write(*,'(8x,a)') 'nonuniform'
              do i=1,nfaces(ib)
                ijb = iBndValueStart(ib) + i
                read(input_unit,*) u(ijb),v(ijb),w(ijb)
                ! write(*,'(8x,3e15.7)') u(ijb),v(ijb),w(ijb)
              end do
            endif
          endif

          if ( type == 'Neumann' ) then
            read(input_unit,*) distribution
            write(*,'(6x,a)') distribution
            if ( distribution == 'uniform' ) then
              read(input_unit,*) u0,v0,w0
              write(*,'(8x,3e15.7)') u0,v0,w0  
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
              write(*,'(8x,a,3e15.7)') 'nonuniform'
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


subroutine initialize_scalar_field(T, dTdxi, iphi, field_name)
!
! Reads 'field_name' file form 0/ folder
!
use geometry, only: bcDFraction

implicit none

  real(dp), dimension(numTotal), intent(inout) :: T
  real(dp), dimension(numTotal), intent(inout):: dTdxi
  integer, intent(in) :: iphi
  character(len=*) :: field_name

  !
  ! Locals
  !
  integer :: i,ib,inp,ijp,ijb,iface
  integer :: input_status, input_unit
  character(len=80) :: key,field_type
  character(len=15) :: name,type,distribution
  real(dp) :: t0
 
  write(*,'(a)') ' '
  write(*,'(a)') '  Initializing internal field and boundaries (reading 0/'//trim(field_name)//' ):'
  write(*,'(a)') ' '

  ! 
  ! > Initialize scalar field


  call get_unit ( input_unit )
  open ( unit = input_unit, file = '0/'//trim(field_name))
  write(*,'(a)') '  0/'//trim(field_name)
  rewind input_unit

  do
  read(input_unit,'(a)',iostat = input_status) key
  if(input_status /= 0) exit

    if( key=='internalField' ) then

      write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if( field_type=='uniform' ) then

          read(input_unit,*) t0

          write(*,'(4x,a,3e15.7)') 'uniform',t0

          do inp = 1,numCells
            T(inp) = t0
          enddo

      elseif( field_type=='nonuniform' ) then

          write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) T(inp)
          enddo

      endif

    elseif( key =='boundaryField') then

      write(*,'(2x,a)') 'boundaryField'      

      !
      !*** Initial values for bounary faces and prescription of the BC type 
      !***   which determines how will boundary face values be treated during the code run!
      !

      do ib=1,numBoundaries

      read(input_unit,*) name   ! repeat bcname(ib)
        write(*,'(2x,a)') name
        read(input_unit,*) type ! eg. Dirichlet, Neumann, Robin, Cauchy, empty
          write(*,'(4x,a)') type

          if ( type == 'Dirichlet' ) then

            bcDFraction(iphi,ib) = 1.0_dp !!! Set Dirichlet fraction for value update to one.

            read(input_unit,*) distribution
            write(*,'(6x,a)') distribution

            if ( distribution == 'uniform' ) then
              read(input_unit,*) t0
              write(*,'(8x,a,3e15.7)') 'uniform',t0  
              do i=1,nfaces(ib)
                ijb = iBndValueStart(ib) + i
                T(ijb) = t0
              end do        
            else ! nonuniform
              write(*,'(8x,a,3e15.7)') 'nonuniform'
              do i=1,nfaces(ib)
                ijb = iBndValueStart(ib) + i
                read(input_unit,*) T(ijb)
                ! write(*,'(8x,e15.7)') T(ijb)
              end do
            endif
          endif

          if ( type == 'Neumann' ) then
            read(input_unit,*) distribution
            write(*,'(6x,a)') distribution
            if ( distribution == 'uniform' ) then
              read(input_unit,*) t0
              write(*,'(8x,a,3e15.7)') 'uniform',t0
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
              write(*,'(8x,a,3e15.7)') 'nonuniform'
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


subroutine pipe_disturbances(x,y,z,ux,uy,uz)
!
! Initial disturbances in pipe flow.
! Note: some values are specific for the case I was running - change it appropriately.
!
  use parameters, only: pi
  ! use utils, only: init_random_seed

  implicit none

  real(dp), intent(in) :: x,y,z
  real(dp), intent(inout) :: ux,uy,uz

  real(dp), parameter :: pipeRadius = 0.5
  real(dp), parameter :: pipeLength = 5.0
  real(dp) :: zr,yr,rr,th,xo,uxx
  real(dp) :: amp_x,freq_x,freq_t,amp_tht,amp_clip,blt,phase_x,arg_tht,amp_sin,rand,big

  zr = z/pipeRadius
  yr = y/pipeRadius
  rr = zr*zr + yr*yr
  if (rr.gt.0) rr=sqrt(rr)
  th = atan2(y,z)
  xo = 2*pi*x/pipeLength

  uxx = 6.*(1-rr**6)/5.

  ! Assign a wiggly shear layer near the wall
  amp_x    = 0.35  ! Fraction of 2pi for z-based phase modification
  freq_x   = 4     ! Number of wiggles in axial-direction
  freq_t   = 9     ! Frequency of wiggles in azimuthal-direction

  amp_tht  = 5     ! Amplification factor for clipped sine function
  amp_clip = 0.2   ! Clipped amplitude

  blt      = 0.07  ! Fraction of boundary layer with momentum deficit

  phase_x = amp_x*(2*pi)*sin(freq_x*xo)

  arg_tht = freq_t*th + phase_x
  amp_sin = 5*sin(arg_tht)
  if (amp_sin.gt. amp_clip) amp_sin =  amp_clip
  if (amp_sin.lt.-amp_clip) amp_sin = -amp_clip

  if (rr.gt.(1-blt)) uxx = uxx + amp_sin

  ! Quick P-independent randomizer
  big  = 1.e3*x + 1.e2*y + 1.e1*z
  rand = sin(big)

  ux   = ux * (uxx + .01*rand)
  uy   = .10*rand*rand*rand
  uz   = .05*rand*rand

end subroutine


subroutine channel_disturbances(x,y,z,ux,uy,uz)
!
! Initial disturbances in channel flow.
! Note: some values are specific for the case I was running - change it appropriately.
!
  use parameters, only: pi
  ! use utils, only: init_random_seed

  implicit none

  real(dp), intent(in) :: x,y,z
  real(dp), intent(inout) :: ux,uy,uz

  real(dp), parameter :: channelHalfWidth = 1.0
  real(dp), parameter :: channelLength = 6.28 ! 4.0
  real(dp) :: yr,rr,th,xo,uxx
  real(dp) :: amp_x,freq_x,freq_t,amp_tht,amp_clip,blt,phase_x,arg_tht,amp_sin,rand,big

  yr = y !yr = (y-1.)/channelHalfWidth ! NOTE: yr = y if y in [-1,1]; yr = y-1.0 if y in [0,2] in vertical direction.
  rr = yr*yr
  if (rr.gt.0) rr=sqrt(rr)
  th = atan2(y,z)
  xo = 2*pi*x/channelLength

  uxx = 6.*(1-rr**6)/5.

  ! Assign a wiggly shear layer near the wall
  amp_x    = 0.35  ! Fraction of 2pi for z-based phase modification
  freq_x   = 4     ! Number of wiggles in axial-direction
  freq_t   = 9     ! Frequency of wiggles in azimuthal-direction

  amp_tht  = 5     ! Amplification factor for clipped sine function
  amp_clip = 0.2   ! Clipped amplitude

  blt      = 0.07  ! Fraction of boundary layer with momentum deficit

  phase_x = amp_x*(2*pi)*sin(freq_x*xo)

  arg_tht = freq_t*th + phase_x
  amp_sin = 5*sin(arg_tht)
  if (amp_sin.gt. amp_clip) amp_sin =  amp_clip
  if (amp_sin.lt.-amp_clip) amp_sin = -amp_clip

  if (rr.gt.(1-blt)) uxx = uxx + amp_sin !<<<<

  ! Quick P-independent randomizer
  big  = 1.e3*x + 1.e2*y + 1.e1*z
  rand = sin(big)

  ux   = ux * (uxx + 0.01*rand) ! to scale it to the proper level of u
  uy   = 0.10*rand*rand*rand
  uz   = 0.10*rand*rand*rand !0.05*rand*rand

end subroutine



subroutine initialize_tgv(x,y,z,ux,uy,uz,p)
!
! Initialize field values for Taylor-Green vortex case
!
  use parameters, only: pi

  implicit none

  real(dp), intent(in) :: x,y,z
  real(dp), intent(inout) :: ux,uy,uz,p

  real(dp), parameter :: rho = 1.0_dp
  real(dp), parameter :: Vo = 1.0_dp
  real(dp), parameter :: po = 0.0_dp

  ux = Vo*sin(x)*cos(y)*cos(z)
  uy =-Vo*cos(x)*sin(y)*cos(z)
  uz = 0.0_dp

  p = po + rho*Vo**2/16.0 * ( cos(2*x)+cos(2*y) ) * ( cos(2*z) + 2.0_dp )

end subroutine



end module
