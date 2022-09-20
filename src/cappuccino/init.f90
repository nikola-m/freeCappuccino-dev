! ***********************************************************************
!
subroutine init
!
!***********************************************************************
!   Contents:
!
!   Various initialisations
!       Parameter Initialisation
!       Field Initialisation
!   Read Restart File And Set Field Values
!   Initial Gradient Calculation
!   Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use sparse_matrix
  use utils, only: get_unit
  use field_initialization
  use output
  use temperature, only: calcT
  use energy, only: calcEn
  use TurbModelData
  use mhd
  use wall_distance
  use node_interpolation
  ! use fieldManipulation, only: add_random_noise_to_field

  implicit none

  ! 
  ! Local variables 
  !

  integer :: i, ijp, ijn, ib, iface
  ! integer :: inp
  real(dp) :: fxp, fxn, ui, vi, wi
  real(dp) :: outare, uav, are

  ! real(dp) :: cosa,sina,velmag,distance,xvel,zvel
  ! real(dp) :: nx,ny,xvel,yvel
!
!***********************************************************************
!

!
! Various initialisations
!

  ! Initial time iz zero.
  if(.not.lread) time = 0.0_dp

  ! Set to zero cumulative error in continuity
  cumulativeContErr = 0.0_dp


! 1.2)  Field Initialisation

  call initialize_vector_field(u,v,w,dUdxi,1,'U')

    !
    ! Definition of initial velocity field directly - recompile fCp for specific case.
    ! Sometimes it is better just to print data obtained below and copy to 0/U file.
    ! Like this you can generate profiles for inlet boundaries, for inner field, etc.
    !

    ! ! > Rotating cylinders-Taylor-Couette; velocities at boundary faces for constant angular velocity.
    ! ! Loop over inlet boundaries
    ! do ib=1,numBoundaries   
    !   if ( bcname(ib) == 'wallInner' ) then
    !     do i=1,nfaces(ib)
    !       iface = startFace(ib) + i
    !       ijp = owner(ib)
    !       ijn = iBndValueStart(ib) + i
    !       are = sqrt(arx(iface)**2 + ary(iface)**2 + arz(iface)*2)
    !       nx = arx(iface)/are
    !       ny = ary(iface)/are
    !       xvel = -ny
    !       yvel = nx
    !       write(*,*) xvel, yvel, 0.0
    !     end do
    !   endif 
    ! enddo

    ! > Rotating lid - velocities at boundary faces for constant angular velocity.
    ! ! Loop over inlet boundaries
    ! do ib=1,numBoundaries  
    !   if ( bcname(ib) == 'lid' ) then
    !     do i=1,nfaces(ib)
    !       iface = startFace(ib) + i
    !       ijp = owner(ib)
    !       ijn = iBndValueStart(ib) + i
    !       distance = sqrt(xf(iface)**2 + zf(iface)**2)
    !       cosa = xf(iface)/(distance+1e-20)
    !       sina = zf(iface)/(distance+1e-20)
    !       velmag = distance
    !       xvel = velmag*sina
    !       zvel = -velmag*cosa
    !       write(*,*) xvel, 0.0, zvel
    !     end do
    !   endif 
    ! enddo

    ! Create initial disturbances for turbulent channel and pipe (below) flow
    ! do inp = 1,numCells
    !   call channel_disturbances(xc(inp),yc(inp),zc(inp),u(inp),v(inp),w(inp))
    ! enddo
    !** 
    ! do inp = 1,numCells
    !   call pipe_disturbances(xc(inp),yc(inp),zc(inp),u(inp),v(inp),w(inp))
    ! enddo
    !\***

    ! Curved tube case
    ! ! Loop over boundaries
    ! do ib=1,numBoundaries 
    !   if ( bctype(ib) == 'inlet' ) then
    !     do i=1,nfaces(ib)
    !       iface = startFace(ib) + i
    !       ijp = owner(ib)
    !       ijn = iBndValueStart(ib) + i
    !       U(ijn) = magUbar*1.5*( 1. - ( sqrt( yf(iface)**2+zf(iface)**2 )/0.004 )**2 )
    !     end do
    !   endif 
    ! enddo

    ! ! Initialize field values for Taylor-Green vortex
    ! do ijp=1,numCells
    !   call initialize_tgv( xc(ijp),yc(ijp),zc(ijp),u(ijp),v(ijp),w(ijp),p(ijp) )
    ! enddo



  ! > Pressure
  if(compressible) call initialize_scalar_field(p,dPdxi,4,'p')

  ! > TE Turbulent kinetic energy, 
  if(solveTKE) call initialize_scalar_field(te,dTEdxi,5,'k')

  ! > Nutilda - for Spallart-Allmaras,
  if (solveNuTilda) call initialize_scalar_field(te,dTEdxi,5,'nutilda')

  ! > ED Specific turbulent kinetic energy dissipation rate, also turbulence frequency - omega
  if (solveOmega) call initialize_scalar_field(ed,dEDdxi,6,'omega')
  if (solveEpsilon) call initialize_scalar_field(ed,dEDdxi,6,'epsilon')


  ! 
  ! > Temperature
  !
  if( calcT .or. compressible )   call initialize_scalar_field(t,dTdxi,7,'T')

  if( calcEn )   call initialize_scalar_field(En,dEndxi,7,'Energy')

  ! 
  ! > Magnetic field
  ! 
  if ( calcEpot ) call initialize_vector_field(bmagx,bmagy,bmagz,dEpotdxi,10,'B')

  ! Density
  den = densit

  ! Effective viscosity
  vis = viscos

  ! Set effective viscosity at inlet for appropriate turbulence model
  if(lturb) call modify_viscosity_inlet

  ! Wall effective viscosity
  visw = viscos  
  
!
!******* Initialization of mass fluxes at inner and boundary faces ******
!

  ! Initialize mass flow on inner faces
  do i=1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    fxn = facint(i)
    fxp = 1.0_dp-facint(i)

    ui = u(ijp)*fxp + u(ijn)*fxn
    vi = v(ijp)*fxp + v(ijn)*fxn
    wi = w(ijp)*fxp + w(ijn)*fxn

    flmass(i) = den(ijp)*(arx(i)*ui+ary(i)*vi+arz(i)*wi)

  enddo

  ! Initialize mass flow over inlet faces
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijn = iBndValueStart(ib) + i

        flmass(iface) = den(ijn)*(arx(iface)*u(ijn)+ary(iface)*v(ijn)+arz(iface)*w(ijn))

        ! Face normal vector is faced outwards, while velocity vector at inlet
        ! is faced inwards. That means their scalar product will be negative,
        ! so minus signs here is to turn net mass influx - flomas, into positive value.
        flomas = flomas - flmass(iface)         

      end do

    endif 

  enddo

  ! Outlet area
  outare = 0.0_dp

  ! Loop over outlet boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        outare = outare + sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

      end do

    endif 

  enddo

  ! Average velocity at outlet boundary
  uav = flomas/(densit*outare)

  ! Mass flow trough outlet faces using Uav velocity

  ! Loop over outlet boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijn = iBndValueStart(ib) + i

        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        u(ijn) = uav * arx(iface)/are
        v(ijn) = uav * ary(iface)/are
        w(ijn) = uav * arz(iface)/are

        flmass(iface)=den(ijn)*(arx(iface)*u(ijn)+ary(iface)*v(ijn)+arz(iface)*w(ijn))

      end do

    endif 

  enddo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialized mass flow.'
  write ( *, '(a)' ) ' '
  
  ! write ( *, '(a)' ) ' '
  ! write ( *, '(a)' ) '  Inlet boundary condition information:'
  ! write ( *, '(a,e12.6)' ) '  Mass inflow: ', flomas
  ! write ( *, '(a)' ) ' '

!
!******* END: Initialization of mass fluxes at inner and boundary faces ******
!

!
! Read Restart File And Set Field Values
!
  if(lread) call readfiles


!
! Set interpolation weights for cell-to-node interpolation
!
  call set_interpolation_weights

  
!
! Initial Gradient Calculation
!

  ! Set matrix with geometric coefficients for Least-squares gradient reconstruction.
  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_grad_matrix(U,dUdxi)
  endif

  ! Initial velocity gradient calculation with given velocity field.
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)


!
! Distance to the nearest wall (needed for some turbulence models) for all cells via Poisson equation.
!

  if (lturb) call wall_distance_poisson


end subroutine
