module easm_rumsey_gatski_komega
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Rumsey-Gatski explicit algebraic stress model with the k-omega transport equations (EASMko2003),
!    on top of 'turb_transport' and 'easm_closure'.
!
!  Discussion:
!
!    Dk/Dt     = P - f_beta* rho omega k                   + div[ (mu + sigma_k mut) grad k ]
!    Domega/Dt = gamma (omega/k) P - beta rho omega^2      + div[ (mu + sigma_w mut) grad omega ]
!
!    P = T_ij dU_i/dx_j is computed with the FULL algebraic stress (not mut S^2), limited to 20 f_beta* rho omega k,
!    mut = rho C*mu k/omega with C*mu from the cubic equation of the closure ('easm_cell', lowest real part root).
!    The excess anisotropy is stored in 'aij' and enters the momentum equations as an explicit source
!    (asm_stress_terms, activated through 'lasm' in turbulence_registry).
!
!    Constants (TMR): Cmu = 0.0895, sigma_k = 1, gamma = 0.53, beta = 0.83, kappa = 0.41,
!                     sigma_w = sqrt(Cmu)(beta - gamma)/kappa^2.
!    Wall: low-Re, k_w = 0 (Dirichlet), omega in the wall cell = 10 * 6 nu/(beta y^2); y+ ~ 1 mesh needed.
!
!    The cross-diffusion damping chi_k = Cmu^2 (1/omega^3) grad k . grad omega: omega here is eps/k = Cmu omega_Wilcox,
!    so chi_k of Wilcox (1/omega_w^3 grad k . grad omega_w) carries Cmu^2. Without it f_beta* ~ 1.7 over the outer
!    half of a channel and k is destroyed too fast (Cf was 0.0038 instead of ~0.0068).
!
!  Reference:
!    [1] C.L. Rumsey, T.B. Gatski, J. Aircraft 38(5), 2001.  [2] T.B. Gatski, C.G. Speziale, J. Fluid Mech. 254, 1993.
!    [3] Turbulence Modeling Resource, EASMko2003: https://turbmodels.larc.nasa.gov/easmko.html
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
  use easm_closure, only: easm_cell

  implicit none

  real(dp), parameter :: CMU = 0.0895_dp
  real(dp), parameter :: SIGMA_K = 1.0_dp
  real(dp), parameter :: GAMMA_W = 0.53_dp
  real(dp), parameter :: BETA_W = 0.83_dp
  real(dp), parameter :: KAPPA_EASM = 0.41_dp
  real(dp), parameter :: SIGMA_W = sqrt(CMU)*(BETA_W-GAMMA_W)/KAPPA_EASM**2

  integer, parameter :: iK = 1, iOm = 2

  real(dp), dimension(:), allocatable :: cmustar, pk, pok, fbs, dif, suv, spv ! pok = P/k with k of the time of evaluation
  logical, dimension(:), allocatable :: fixed
  real(dp), dimension(:), allocatable :: fixval

  private

  public :: modify_viscosity_easm_rumsey_gatski_komega
  public :: modify_viscosity_inlet_easm_rumsey_gatski_komega

contains


subroutine modify_viscosity_easm_rumsey_gatski_komega
  implicit none

  if (.not.allocated(cmustar)) then
    allocate( cmustar(numCells), pk(numCells), pok(numCells), fbs(numCells), suv(numCells), spv(numCells), dif(numTotal) )
    allocate( fixed(numCells), fixval(numCells) )
    cmustar = CMU
  endif

  call algebraic_stress   ! closure: C*mu, aij, production
  call calcsc_tke
  call calcsc_omega
  call modify_mu_eff

end subroutine


subroutine algebraic_stress
!
! C*mu and excess anisotropy from the closure, production P = T_ij dU_i/dx_j, damping function f_beta*.
!
  implicit none
  integer :: inp
  real(dp), dimension(3,3) :: g
  real(dp) :: tau, aG, eta2, chik, mut

  do inp=1,numCells

    g(1,:) = dUdxi(:,inp)
    g(2,:) = dVdxi(:,inp)
    g(3,:) = dWdxi(:,inp)

    tau = 1.0_dp/max(ed(inp),small)

    call easm_cell( g, tau, cmustar(inp), aij(:,inp) )

    ! P = 2 mut S_ij S_ij - rho k a_ij dU_i/dx_j   with   mut = rho C*mu k tau
    eta2 = 0.25_dp*sum( (g + transpose(g))**2 )
    aG = aij(1,inp)*g(1,1) + aij(4,inp)*g(2,2) + aij(6,inp)*g(3,3) &
       + aij(2,inp)*(g(1,2)+g(2,1)) + aij(3,inp)*(g(1,3)+g(3,1)) + aij(5,inp)*(g(2,3)+g(3,2))
    mut = den(inp)*cmustar(inp)*te(inp)*tau
    pk(inp) = 2.0_dp*mut*eta2 - den(inp)*te(inp)*aG

    chik = CMU**2*sum( dTEdxi(:,inp)*dEDdxi(:,inp) )/max(ed(inp),small)**3 ! Wilcox chi_k, omega_w = omega/Cmu
    if ( chik > 0.0_dp ) then
      fbs(inp) = ( 1.0_dp + 680.0_dp*chik**2 )/( 1.0_dp + 400.0_dp*chik**2 )
    else
      fbs(inp) = 1.0_dp
    endif

    pk(inp) = min( pk(inp), 20.0_dp*fbs(inp)*den(inp)*ed(inp)*te(inp) )
    gen(inp) = pk(inp)
    pok(inp) = pk(inp)/max(te(inp),small)   ! P/k stays consistent when k changes in the k solve

  enddo

end subroutine


subroutine calcsc_tke
  implicit none
  integer :: inp

  do inp=1,numCells
    suv(inp) = max(pk(inp),zero)*vol(inp)
    spv(inp) = ( fbs(inp)*den(inp)*ed(inp) - min(pk(inp),zero)/(te(inp)+small) )*vol(inp)
  enddo

  dif = viscos + (vis-viscos)*SIGMA_K

  call transport_scalar( te, dTEdxi, teo, teoo, iK, 'k', suv, spv, dif, wallDirichlet=.true. )

end subroutine


subroutine calcsc_omega
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall

  do inp=1,numCells
    ! gamma (omega/k) P: P/k from the start of the iteration (k has been changed by the k solve, P has not)
    suv(inp) = GAMMA_W*ed(inp)*max(pok(inp),zero)*vol(inp)
    spv(inp) = ( BETA_W*den(inp)*ed(inp) - GAMMA_W*min(pok(inp),zero) )*vol(inp)
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
      fixval(ijp) = 10.0_dp*6.0_dp*viscos/den(ijp)/(BETA_W*dnw(iWall)**2)
    enddo
  enddo

  dif = viscos + (vis-viscos)*SIGMA_W

  call transport_scalar( ed, dEDdxi, edo, edoo, iOm, 'omega', suv, spv, dif, fixed=fixed, fixval=fixval )

end subroutine


subroutine modify_mu_eff
!
! mut = rho C*mu k/omega; wall viscosity = molecular (low-Re), wall shear stress and y+ for output.
!
  implicit none
  integer :: inp, ib, i, iface, ijp, ijb, iWall
  real(dp) :: urf, visold, nxf, nyf, nzf, vtp, xtp, ytp, ztp, ut2, utau

  urf = TurbModel%urfVis

  do inp=1,numCells
    visold = vis(inp)
    vis(inp) = viscos + den(inp)*cmustar(inp)*te(inp)/max(ed(inp),small)
    vis(inp) = urf*vis(inp) + (1.0_dp-urf)*visold
  enddo

  call updateBoundary( vis )

  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      ijb = iBndValueStart(ib) + i
      iWall = iWall + 1
      call wall_normal_tangent( iface, ijp, nxf, nyf, nzf, vtp, xtp, ytp, ztp )
      ut2 = abs( (u(ijb)-u(ijp))*xtp + (v(ijb)-v(ijp))*ytp + (w(ijb)-w(ijp))*ztp )
      tau(iWall) = viscos*ut2/dnw(iWall)
      utau = sqrt( tau(iWall)/den(ijb) )
      ypl(iWall) = den(ijb)*utau*dnw(iWall)/viscos
      visw(iWall) = viscos
      vis(ijb) = viscos
    enddo
  enddo

end subroutine


subroutine modify_viscosity_inlet_easm_rumsey_gatski_komega
  implicit none
  integer :: i, ib, ijb

  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        vis(ijb) = viscos + den(ijb)*CMU*te(ijb)/(ed(ijb)+small)
      end do
    endif
  enddo

end subroutine

end module
