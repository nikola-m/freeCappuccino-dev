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
  use turbulence
  use mhd
  use wall_distance
  use fieldManipulation, only: add_random_noise_to_field

  implicit none

  ! 
  ! Local variables 
  !

  integer :: i, ijp, ijn
  real(dp) :: fxp, fxn, ui, vi, wi

  ! integer :: ib, iface
  ! real(dp) :: cosa,sina,velmag,distance,xvel,zvel

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


    ! > Rotating lid - velocitiesa at boundary faces for constant angular velocity.
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

  ! ! > Initialize inner cells for turbulent channel flow
  ! do inp=1,numCells
  ! U(inp) = magUbar*(1.2*(1-yc(inp)**6))  
  ! enddo
  ! call add_random_noise_to_field(U,20)
  !***
  !Create initial disturbances
  ! do inp = 1,numCells
  !   call channel_disturbances(xc(inp),yc(inp),zc(inp),u(inp),v(inp),w(inp))
  ! enddo  
  !\***

  ! > Pressure
  if(AllSpeedsSIMPLE) call initialize_scalar_field(p,dPdxi,4,'p')

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
  if( calcT )   call initialize_scalar_field(t,dTdxi,7,'T')

  if( calcEn )   call initialize_scalar_field(t,dTdxi,7,'Energy')

  ! 
  ! > Magnetic field
  ! 
  if ( calcEpot ) call initialize_vector_field(bmagx,bmagy,bmagz,dEpotdxi,10,'B')

  ! Density
  den = densit

  ! Effective viscosity
  vis = viscos
  ! if ( lturb ) call initialize_scalar_field(vis,8,'mueff') ! no gradient, it should be optional argument
  visw = viscos  
  

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

!
! Read Restart File And Set Field Values
!
  if(lread) call readfiles

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
