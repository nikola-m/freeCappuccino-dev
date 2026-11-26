module k_epsilon_zeta_f
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    k-epsilon-zeta-f turbulence model (Hanjalic, Popovac, Hadziabdic, 2004), written on top of 'turb_transport'.
!
!  Discussion:
!
!    Four scalar equations, TurbModel%Scalar(1:4) = k, epsilon, zeta, f:
!
!      Dk/Dt     = P - eps                              + div[ (nu + nut/sigma_k) grad k ]
!      Deps/Dt   = ( Ceps1 P - Ceps2 eps ) / T          + div[ (nu + nut/sigma_e) grad eps ]
!      Dzeta/Dt  = f - zeta P / k                       + div[ (nu + nut/sigma_z) grad zeta ]
!      L^2 lapl(f) - f = (1/T) ( C1 - 1 + C2' P/eps ) ( zeta - 2/3 )            (elliptic relaxation)
!
!      nut = Cmu zeta k T
!      T   = max( min( k/eps, 0.6/(sqrt(6) Cmu S zeta) ), CT sqrt(nu/eps) )
!      L   = CL max( min( k^1.5/eps, k^0.5/(sqrt(6) Cmu S zeta) ), Ceta (nu^3/eps)^0.25 )
!      Ceps1 = 1.4 (1 + 0.012/zeta),   S = sqrt(2 Sij Sij)
!
!    Wall: zeta_w = 0, k_w = 0 (zero flux), first-cell epsilon = 2 nu k / y^2 and f = -2 nu zeta / y^2 (fixed cell values).
!    If the first cell is in the log layer (y* > ctrans) the wall-function branch is used: epsilon = Cmu^0.75 k^1.5/(kappa y),
!    production from the wall shear stress, wall viscosity from the log law. This is a hard switch, NOT the blended
!    compound wall treatment of Popovac & Hanjalic (2007) - intended use is y+ ~ 1.
!
!    Changes with respect to the previous implementation of this module (patch 0033):
!      * assembly delegated to turb_transport (1392 -> ~450 lines),
!      * zeta is initialised (2/3 inside, 0 on walls); before it was allocated but never set,
!      * zeta wall BC: zeta_w = 0 Dirichlet; before the first wall cell was overwritten with 2 nu k / y^2,
!      * epsilon sink: Ceps2 rho eps / T (the factor 'eps' appeared twice in the implicit coefficient),
!      * division by zero protection of the time scale for zeta -> 0.
!
!  Reference:
!    [1] K. Hanjalic, M. Popovac, M. Hadziabdic, A Robust Near-Wall Elliptic-Relaxation Eddy-Viscosity Turbulence Model for CFD,
!        Int. J. Heat and Fluid Flow, Vol. 25, 2004, pp.1047-1051.
!    [2] M. Popovac, K. Hanjalic, Compound Wall Treatment for RANS Computation of Complex Turbulent Flows and Heat Transfer,
!        Flow, Turbulence and Combustion, Vol. 78, 2007, pp.177-202.
!    [3] Langley Research Center-Turbulence Modeling Resource, turbmodels.larc.nasa.gov
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Author:
!
!    Nikola Mirkov (original module, 17th February 2023); rewritten on turb_transport, November 2026.
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use TurbModelData, only: TurbModel
  use turb_transport

  implicit none

  ! Turbulence model constants
  real(dp), parameter :: CMU = 0.22_dp
  real(dp), parameter :: SIXCMU = sqrt(6.0_dp)*CMU
  real(dp), parameter :: C1 = 1.4_dp
  real(dp), parameter :: CEPS2 = 1.9_dp
  real(dp), parameter :: C2P = 0.65_dp
  real(dp), parameter :: CT = 6.0_dp
  real(dp), parameter :: CL = 0.36_dp
  real(dp), parameter :: CETHA = 85.0_dp

  real(dp), parameter :: sigma_k = 1.0_dp
  real(dp), parameter :: sigma_epsilon = 1.3_dp
  real(dp), parameter :: sigma_zeta = 1.2_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(cmu))
  real(dp), parameter :: CMU75 = cmu25**3

  ! Scalar slots in TurbModel%Scalar
  integer, parameter :: iK = 1, iEps = 2, iZeta = 3, iF = 4

  real(dp), dimension(:), allocatable :: frlx                  ! elliptic relaxation function f
  real(dp), dimension(:), allocatable :: zeta, zetao, zetaoo   ! zeta = v^2/k and its old time levels
  real(dp), dimension(:,:), allocatable :: dZetadxi, dFrlxdxi
  real(dp), dimension(:), allocatable :: Tt, Ll2               ! turbulent time scale T and squared length scale L^2
  real(dp), dimension(:), allocatable :: dif                   ! diffusivity work array
  real(dp), dimension(:), allocatable :: suv, spv              ! volume sources work arrays
  logical, dimension(:), allocatable :: fixed
  real(dp), dimension(:), allocatable :: fixval
  integer :: stamp = -1

  private 

  public :: modify_viscosity_k_epsilon_zeta_f
  public :: modify_viscosity_inlet_k_epsilon_zeta_f

contains


subroutine modify_viscosity_k_epsilon_zeta_f
!
! Purpose: 
!   Main module routine to solve turbulence model equations and update effective viscosity.
!
  implicit none

  if(.not.allocated(zeta)) call allocate_and_initialise

  call shift_time_levels(zeta,zetao,zetaoo,stamp)

  call time_and_length_scales

  call calcsc_tke     ! Assemble and solve turbulence kinetic energy eq.
  call calcsc_epsilon ! Assemble and solve dissipation rate of tke eq.
  call calcsc_frelax  ! Assemble and solve elliptic relaxation eqn.
  call calcsc_zeta    ! Assemble and solve zeta equation.
  call modify_mu_eff  ! Update efective viscosity

end subroutine


subroutine allocate_and_initialise
!
! Fields owned by the model. zeta starts from its isotropic value 2/3 (wall value 0), f from zero.
!
  implicit none
  integer :: ib, i, ijb

  allocate( frlx(numTotal), zeta(numTotal), zetao(numTotal), zetaoo(numTotal) )
  allocate( dZetadxi(3,numTotal), dFrlxdxi(3,numTotal) ) ! 'grad' fills boundary values too
  allocate( Tt(numCells), Ll2(numCells), dif(numTotal), suv(numCells), spv(numCells) )
  allocate( fixed(numCells), fixval(numCells) )

  frlx = 0.0_dp
  zeta = twothirds
  do ib=1,numBoundaries
    if ( bctype(ib) == 'wall' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        zeta(ijb) = 0.0_dp
      enddo
    endif
  enddo
  zetao = zeta
  zetaoo = zeta

  write(*,'(a)') ' '
  write(*,'(a)') '  **Allocated and initialised zeta and f fields for k-epsilon-zeta-f model.'
  write(*,'(a)') ' '

end subroutine


subroutine time_and_length_scales
!
! Time scale T and squared length scale L^2 of the model, in every cell.
!
  implicit none
  integer :: inp
  real(dp) :: nu, szeta, lenscale

  nu = viscos/densit

  do inp=1,numCells

    szeta = SIXCMU*magStrain(inp)*max(zeta(inp),small) + small ! sqrt(6) Cmu S zeta

    Tt(inp) = max( min( te(inp)/(ed(inp)+small), 0.6_dp/szeta ), CT*sqrt( nu/(ed(inp)+small) ) )

    lenscale = CL*max( min( te(inp)**1.5_dp/(ed(inp)+small), sqrt(te(inp))/szeta ), &
                       CETHA*sqrt(sqrt( nu**3/(ed(inp)+small) )) )
    Ll2(inp) = lenscale*lenscale

  enddo

end subroutine


subroutine calcsc_tke
!
! Purpose:
!  Assemble and solve turbulence kinetic energy equation.
!
  implicit none
  integer :: inp, ib, i, iface, ijp, ijb, iWall
  real(dp) :: genp, genn, prtr

  prtr = 1.0_dp/sigma_k

  ! Standard production
  gen(1:numCells) = abs(vis(1:numCells)-viscos)*magStrain(1:numCells)**2
  ! Production limiter P <= 10 rho eps (as in SST): the impulsive start from a uniform field gives P/eps ~ 1e6 at the wall
  gen(1:numCells) = min( gen(1:numCells), 10.0_dp*den(1:numCells)*ed(1:numCells) )

  ! Wall-adjacent cells in the log layer: production from the wall shear stress
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      if ( ypl(iWall) > ctrans ) gen(ijp) = abs(tau(iWall))*cmu25*sqrt(te(ijp))/(dnw(iWall)*cappa)
    enddo
  enddo

  do inp=1,numCells
    genp = max(gen(inp),zero)
    genn = min(gen(inp),zero)
    suv(inp) = genp*vol(inp)
    spv(inp) = ( den(inp)*ed(inp) - genn )*vol(inp)/(te(inp)+small)
  enddo

  dif = viscos + (vis-viscos)*prtr

  call transport_scalar( te, dTEdxi, teo, teoo, iK, 'k', suv, spv, dif )

end subroutine


subroutine calcsc_epsilon
!
! Purpose:
!  Assemble and solve dissipation rate equation. Wall-adjacent cell values are prescribed.
!
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall
  real(dp) :: genp, genn, prtr, ceps1

  prtr = 1.0_dp/sigma_epsilon

  do inp=1,numCells
    genp = max(gen(inp),zero)
    genn = min(gen(inp),zero)
    ceps1 = 1.4_dp*( 1.0_dp + 0.012_dp/(zeta(inp)+small) )

    suv(inp) = ceps1*genp*vol(inp)/Tt(inp)
    spv(inp) = ( CEPS2*den(inp) - ceps1*genn/(ed(inp)+small) )*vol(inp)/Tt(inp)
  enddo

  fixed = .false.
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      fixed(ijp) = .true.
      if ( ypl(iWall) > ctrans ) then
        fixval(ijp) = cmu75*te(ijp)**1.5_dp/(cappa*dnw(iWall))       ! log layer
      else
        fixval(ijp) = 2.0_dp*viscos/densit*te(ijp)/dnw(iWall)**2     ! viscous sublayer
      endif
    enddo
  enddo

  dif = viscos + (vis-viscos)*prtr

  call transport_scalar( ed, dEDdxi, edo, edoo, iEps, 'epsilon', suv, spv, dif, fixed=fixed, fixval=fixval )

end subroutine


subroutine calcsc_frelax
!
! Purpose: 
!   Assemble and solve elliptic relaxation equation   f - L^2 lapl(f) = - (1/T)(C1 - 1 + C2' P/eps)(zeta - 2/3).
!   No convection, no time derivative. Wall-adjacent cell: f = -2 nu zeta / y^2.
!
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall

  do inp=1,numCells
    suv(inp) = -( C1 - 1.0_dp + C2P*gen(inp)/(densit*ed(inp)+small) )*( zeta(inp)-twothirds )/Tt(inp)*vol(inp)
    spv(inp) = vol(inp)
  enddo

  fixed = .false.
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      fixed(ijp) = .true.
      fixval(ijp) = -2.0_dp*viscos/densit*zeta(ijp)/dnw(iWall)**2
    enddo
  enddo

  dif(1:numCells) = Ll2
  call copy_owner_to_boundary( dif )

  call transport_scalar( frlx, dFrlxdxi, frlx, frlx, iF, 'f', suv, spv, dif, &
                         lconv=.false., ltrans=.false., fixed=fixed, fixval=fixval, positive=.false. )

end subroutine


subroutine calcsc_zeta
!
! Purpose: 
!   Assemble and solve zeta = v^2/k equation. zeta_w = 0 on walls (Dirichlet, with diffusive flux).
!
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall
  real(dp) :: prtr

  prtr = 1.0_dp/sigma_zeta

  do inp=1,numCells
    suv(inp) = den(inp)*frlx(inp)*vol(inp)
    spv(inp) = gen(inp)*vol(inp)/(te(inp)+small)
  enddo

  ! Wall-adjacent cell: f = -2 nu zeta / y^2 is proportional to zeta, so it is taken implicitly (a sink
  ! coefficient) - treating it as a lagged explicit source makes the iteration oscillate and drives zeta negative.
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      suv(ijp) = 0.0_dp
      spv(ijp) = spv(ijp) + 2.0_dp*viscos/dnw(iWall)**2*vol(ijp)
    enddo
  enddo

  dif = viscos + (vis-viscos)*prtr

  call transport_scalar( zeta, dZetadxi, zetao, zetaoo, iZeta, 'zeta', suv, spv, dif, wallDirichlet=.true. )

  zeta(1:numCells) = min( zeta(1:numCells), 2.0_dp )   ! realizability: v'^2 <= 2k
  call updateBoundary( zeta )

end subroutine


subroutine modify_mu_eff
!
! Purpose: 
!   Update turbulent and effective viscosity, and wall viscosity, y+ and wall shear stress.
!
  implicit none

  integer :: inp
  real(dp) :: urf,visold

  urf = TurbModel%urfVis

  do inp=1,numCells
    visold = vis(inp)
    vis(inp) = viscos + den(inp)*cmu*zeta(inp)*te(inp)*Tt(inp)
    vis(inp) = urf*vis(inp)+(1.0_dp-urf)*visold
  enddo

  ! Every boundary type except 'wall', which is treated below.
  call updateBoundary( vis )

  call wall_viscosity_loglaw( cmu25 )

end subroutine modify_mu_eff


subroutine modify_viscosity_inlet_k_epsilon_zeta_f
!
! Update turbulent and effective viscosity at inlet.
!
  implicit none
  integer :: i,ib,ijb

  if(.not.allocated(zeta)) call allocate_and_initialise

  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        zeta(ijb) = twothirds
        vis(ijb) = viscos + den(ijb)*zeta(ijb)*cmu*te(ijb)**2/(ed(ijb)+small)
      end do
    endif
  enddo 

end subroutine

end module
