module earsm_gatski_speziale
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Gatski-Speziale explicit algebraic Reynolds stress model with the standard k-epsilon transport equations,
!    on top of 'turb_transport' and 'easm_closure' (HIGH-REYNOLDS version with wall functions).
!
!  Discussion:
!
!    Dk/Dt   = P - rho eps                                 + div[ (mu + mut/sigma_k) grad k ]
!    Deps/Dt = (eps/k) ( Ceps1 P - Ceps2 rho eps )          + div[ (mu + mut/sigma_e) grad eps ]
!
!    time scale tau = max( k/eps, 6 sqrt(nu/eps) );  mut = rho C*mu k tau;  P = T_ij dU_i/dx_j with the algebraic stress.
!    C*mu (variable, replaces the constant 0.09) and the excess anisotropy 'aij' come from 'easm_cell'; 'aij' goes into the
!    momentum equations as an explicit source (asm_stress_terms, activated through 'lasm').
!    Constants: Ceps1 = 1.44, Ceps2 = 1.83 (they enter the closure through gamma_1*), sigma_k = 1, sigma_e = 1.3.
!    Walls: standard wall functions - first cell in the log layer (y* > 11.63): epsilon = Cmu^3/4 k^3/2/(kappa y), production
!    from the wall shear stress; otherwise the viscous-sublayer formulas (eps = 2 nu k/y^2). Without near-wall damping
!    functions the model is NOT meant for y+ ~ 1 meshes - use 'easm_rumsey_gatski_komega' there. The wall function uses
!    the constant Cmu = 0.09, not the variable C*mu.
!
!  Reference:
!    T.B. Gatski, C.G. Speziale, On explicit algebraic stress models for complex turbulent flows, J. Fluid Mech. 254, 59-78, 1993.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use TurbModelData, only: TurbModel
  use turb_transport
  use easm_closure, only: easm_cell, CEPS1, CEPS2

  implicit none

  real(dp), parameter :: CMU_WALL = 0.09_dp
  real(dp), parameter :: CMU25 = sqrt(sqrt(CMU_WALL))
  real(dp), parameter :: CMU75 = CMU25**3
  real(dp), parameter :: SIGMA_K = 1.0_dp, SIGMA_E = 1.3_dp

  integer, parameter :: iK = 1, iEps = 2

  real(dp), dimension(:), allocatable :: cmustar, tscale, pk, dif, suv, spv
  logical, dimension(:), allocatable :: fixed
  real(dp), dimension(:), allocatable :: fixval

  private

  public :: modify_viscosity_earsm_gatski_speziale
  public :: modify_viscosity_inlet_earsm_gatski_speziale

contains


subroutine modify_viscosity_earsm_gatski_speziale
  implicit none

  if (.not.allocated(cmustar)) then
    allocate( cmustar(numCells), tscale(numCells), pk(numCells), suv(numCells), spv(numCells), dif(numTotal) )
    allocate( fixed(numCells), fixval(numCells) )
    cmustar = CMU_WALL
  endif

  call algebraic_stress
  call calcsc_tke
  call calcsc_epsilon
  call modify_mu_eff

end subroutine


subroutine algebraic_stress
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall
  real(dp), dimension(3,3) :: g
  real(dp) :: aG, eta2, mut, nu

  nu = viscos/densit

  do inp=1,numCells

    g(1,:) = dUdxi(:,inp)
    g(2,:) = dVdxi(:,inp)
    g(3,:) = dWdxi(:,inp)

    tscale(inp) = max( te(inp)/max(ed(inp),small), 6.0_dp*sqrt(nu/max(ed(inp),small)) )

    call easm_cell( g, tscale(inp), cmustar(inp), aij(:,inp) )

    eta2 = 0.25_dp*sum( (g + transpose(g))**2 )
    aG = aij(1,inp)*g(1,1) + aij(4,inp)*g(2,2) + aij(6,inp)*g(3,3) &
       + aij(2,inp)*(g(1,2)+g(2,1)) + aij(3,inp)*(g(1,3)+g(3,1)) + aij(5,inp)*(g(2,3)+g(3,2))
    mut = den(inp)*cmustar(inp)*te(inp)*tscale(inp)
    pk(inp) = 2.0_dp*mut*eta2 - den(inp)*te(inp)*aG
    pk(inp) = max( pk(inp), 0.0_dp )

  enddo

  ! Wall-adjacent cells in the log layer: production from the wall shear stress.
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      if ( ypl(iWall) > ctrans ) pk(ijp) = abs(tau(iWall))*cmu25*sqrt(te(ijp))/(dnw(iWall)*cappa)
    enddo
  enddo

  gen(1:numCells) = pk

end subroutine


subroutine calcsc_tke
  implicit none
  integer :: inp

  do inp=1,numCells
    suv(inp) = pk(inp)*vol(inp)
    spv(inp) = den(inp)*ed(inp)/(te(inp)+small)*vol(inp)
  enddo

  dif = viscos + (vis-viscos)/SIGMA_K

  call transport_scalar( te, dTEdxi, teo, teoo, iK, 'k', suv, spv, dif )

end subroutine


subroutine calcsc_epsilon
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall

  do inp=1,numCells
    suv(inp) = CEPS1*pk(inp)*ed(inp)/(te(inp)+small)*vol(inp)
    spv(inp) = CEPS2*den(inp)*ed(inp)/(te(inp)+small)*vol(inp)
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
        fixval(ijp) = cmu75*te(ijp)**1.5_dp/(cappa*dnw(iWall))
      else
        fixval(ijp) = 2.0_dp*viscos/densit*te(ijp)/dnw(iWall)**2
      endif
    enddo
  enddo

  dif = viscos + (vis-viscos)/SIGMA_E

  call transport_scalar( ed, dEDdxi, edo, edoo, iEps, 'epsilon', suv, spv, dif, fixed=fixed, fixval=fixval )

end subroutine


subroutine modify_mu_eff
  implicit none
  integer :: inp
  real(dp) :: urf, visold, nu

  urf = TurbModel%urfVis
  nu = viscos/densit

  do inp=1,numCells
    visold = vis(inp)
    tscale(inp) = max( te(inp)/max(ed(inp),small), 6.0_dp*sqrt(nu/max(ed(inp),small)) )
    vis(inp) = viscos + den(inp)*cmustar(inp)*te(inp)*tscale(inp)
    vis(inp) = urf*vis(inp) + (1.0_dp-urf)*visold
  enddo

  call updateBoundary( vis )
  call wall_viscosity_loglaw( cmu25 )

end subroutine


subroutine modify_viscosity_inlet_earsm_gatski_speziale
  implicit none
  integer :: i, ib, ijb

  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        vis(ijb) = viscos + den(ijb)*CMU_WALL*te(ijb)**2/(ed(ijb)+small)
      end do
    endif
  enddo

end subroutine

end module
