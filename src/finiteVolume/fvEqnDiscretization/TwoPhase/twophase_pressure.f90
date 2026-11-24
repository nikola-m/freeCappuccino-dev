module twophase_pressure
!
!********************************************************************************************************************************132
!
!  Purpose:
!
!    Pressure-velocity coupling for incompressible two-phase flows with variable density (SIMPLE type pressure correction).
!
!  Discussion:
!
!    The stock routines calcp_simple/calcp_piso enforce div(rho u) = 0, which is the continuity equation of a flow with constant
!    density. With a density jump across the interface the incompressibility constraint is div(u) = 0, i.e. it is the volume
!    flux Q_f = u_f.S_f which has to be conservative; the mass flux of the momentum equations is then flmass = rho_f Q_f.
!    This routine therefore solves the equation for the pressure correction p' with the volume flux:
!
!      sum_f Q*_f  +  sum_f ( -K_f Df ) ( p'_N - p'_P )  =  0,        K_f = interpolated  V_P/A_P  (A_P is the momentum diagonal),
!
!    where the face flux of the predicted velocity is obtained by the Rhie-Chow interpolation with the balanced-force correction
!    of Francois et al. (2006): the surface tension and gravity forces enter the face flux in the same compact form in which they
!    appear in the pressure difference,
!
!      Q*_f = U_f.S_f - K_f Df ( p_N - p_P - dp/dx_f.(x_N-x_P) - fbal_f ),
!
!    with fbal_f = (force in compact form) - (interpolated cell force) projected on the cell-center distance vector (clsvof module).
!    In a static drop (or in hydrostatic equilibrium) this makes the discrete pressure gradient and the force balance exactly,
!    and the parasitic currents are only due to the errors of the curvature.
!
!  Limitations:
!
!    Closed domains (wall, symmetry, empty), orthogonal or mildly non-orthogonal meshes (no non-orthogonal correctors).
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use linear_solvers
  use nablap
  use velocity, only: updateVelocityAtBoundary
  use pressure, only: lSolverP, maxiterP, tolAbsP, tolRelP, urfP
  use clsvof, only: flvol, fbal, cl_contErr

  implicit none

  private

  public :: calcp_twophase

contains


!***********************************************************************
!
subroutine calcp_twophase
!
!***********************************************************************
!
  implicit none

  integer :: i, k, ijp, ijn
  real(dp) :: cap, kj, fxp, fxn, ui, vi, wi, dpxi, dpyi, dpzi, xpn, ypn, zpn, ppref
  real(dp), dimension(:), allocatable :: dvol

  a = 0.0_dp
  su = 0.0_dp

  ! > Assemble the pressure correction equation and find the volume flux by Rhie-Chow interpolation

  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    fxn = facint(i)
    fxp = 1.0_dp - fxn

    xpn = xc(ijn) - xc(ijp)
    ypn = yc(ijn) - yc(ijp)
    zpn = zc(ijn) - zc(ijp)

    ! Volume weighted inverse momentum diagonal: velocity correction is u' = - V/A_P grad(p')
    kj = vol(ijp)*apu(ijp)*fxp + vol(ijn)*apu(ijn)*fxn
    cap = -kj*Df(i)

    ui = u(ijp)*fxp + u(ijn)*fxn
    vi = v(ijp)*fxp + v(ijn)*fxn
    wi = w(ijp)*fxp + w(ijn)*fxn

    dpxi = ( dPdxi(1,ijn)*fxn + dPdxi(1,ijp)*fxp )*xpn
    dpyi = ( dPdxi(2,ijn)*fxn + dPdxi(2,ijp)*fxp )*ypn
    dpzi = ( dPdxi(3,ijn)*fxn + dPdxi(3,ijp)*fxp )*zpn

    flvol(i) = ui*arx(i) + vi*ary(i) + wi*arz(i) + cap*( p(ijn) - p(ijp) - dpxi - dpyi - dpzi - fbal(i) )

    ! Off-diagonal elements
    k = icell_jcell_csr_index(i)
    a(k) = cap

    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! Elements on main diagonal
    k = diag(ijp)
    a(k) = a(k) - cap

    k = diag(ijn)
    a(k) = a(k) - cap

    ! Sources
    su(ijp) = su(ijp) - flvol(i)
    su(ijn) = su(ijn) + flvol(i)

  end do

  ! Closed domain: no flux through the boundary faces (wall, symmetry, empty), zero Neumann condition for p'.
  flvol(numInnerFaces+1:numFaces) = 0.0_dp

  if(ltest) write(6,'(19x,a,1pe10.3)') ' Initial sum  =',sum(su)

  pp = 0.0_dp

  call csrsolve( lSolverP, pp, su, resor(4), maxiterP, tolAbsP, tolRelP, 'p' )

  ! Correct the volume fluxes
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    k = icell_jcell_csr_index(i)
    flvol(i) = flvol(i) + a(k)*( pp(ijn) - pp(ijp) )
  enddo

  ! Reference pressure correction
  ppref = pp(pRefCell)

  ! Source(u_i) = sum_f {-p'}_f * S_fi, correction is u = u* + sum_f {-p'}_f * S_fi / ap.
  call gradp_and_sources(pp)

  u(1:numCells) = u(1:numCells) + su*apu
  v(1:numCells) = v(1:numCells) + sv*apv
  w(1:numCells) = w(1:numCells) + sw*apw
  p(1:numCells) = p(1:numCells) + urfp*( pp(1:numCells) - ppref )

  call updateVelocityAtBoundary

  ! The mass flux of the momentum equations (flmass) stays the one of the interface advection (consistent with the density).

  ! Volume continuity error: max_P |sum_f Q_f| dt / V_P
  allocate( dvol(numCells) )
  dvol = 0.0_dp
  do i=1,numInnerFaces
    dvol(owner(i)) = dvol(owner(i)) + flvol(i)
    dvol(neighbour(i)) = dvol(neighbour(i)) - flvol(i)
  enddo
  cl_contErr = maxval( abs(dvol)/vol(1:numCells) )*timestep
  deallocate( dvol )

  write(6,'(19x,a,1pe10.3)') ' Max. volume continuity error [dt div(u)] = ', cl_contErr

end subroutine


end module
