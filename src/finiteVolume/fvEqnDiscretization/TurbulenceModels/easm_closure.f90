module easm_closure
!********************************************************************************************************************************132 
!
!  Purpose:
!
!    Explicit algebraic Reynolds stress closure of Gatski & Speziale (1993) in the form used by Rumsey & Gatski (2003):
!    the weak-equilibrium solution of the algebraic stress equation with the SSG pressure-strain model, truncated
!    to the three-term tensor basis  { S, S W - W S, S S - 1/3 {S S} I }.
!
!  Discussion:
!
!    With  tau = 1/omega  (or k/epsilon),  S_ij = 1/2 (dU_i/dx_j + dU_j/dx_i),  W_ij = 1/2 (dU_i/dx_j - dU_j/dx_i):
!
!      T_ij = -rho <u_i u_j> = 2 mut { S_ij + a2 a4 (S_ik W_kj - W_ik S_kj) - 2 a3 a4 (S_ik S_kj - 1/3 S_kl S_lk delta_ij) }
!                              - 2/3 rho k delta_ij,
!      mut  = rho C*mu k tau,
!      a1 = (4/3 - C2)/2,  a2 = (2 - C4)/2,  a3 = (2 - C3)/2,
!      a4 = tau / ( g1 + 2 g0 C*mu eta^2 tau^2 ),   g0 = C1_1/2,  g1 = C1_0/2 + (Ceps2 - Ceps1)/(Ceps1 - 1),
!
!    and C*mu is a root of the cubic  (-C*mu)^3 + p (-C*mu)^2 + q (-C*mu) + r = 0,
!
!      p = -g1/(eta^2 tau^2 g0),   q = ( g1^2 - 2 eta^2 tau^2 g0 a1 - 2/3 eta^2 tau^2 a3^2 + 2 R^2 eta^2 tau^2 a2^2 ) / (2 eta^2 tau^2 g0)^2,
!      r = g1 a1 / (2 eta^2 tau^2 g0)^2,   eta^2 = S_ij S_ij,   R^2 = W_ij W_ij / (S_ij S_ij).
!
!    As prescribed by Rumsey & Gatski the root with the LOWEST REAL PART is taken (the Cardano / trigonometric solution
!    below follows the Turbulence Modeling Resource, EASMko2003 page). C*mu is floored at 5e-4.
!
!    The routine returns the anisotropy in excess of the Boussinesq one, in the convention of the momentum source
!    ('aij', see asm_stress_terms):  T_ij = 2 mut S_ij - 2/3 rho k delta_ij - rho k a_ij,
!
!      a_ij = -2 C*mu tau { a2 a4 (S W - W S) - 2 a3 a4 (S S - 1/3 {S S} I) }.
!
!    Realizability: the excess part is scaled by the largest factor in [0,1] which keeps <u_i u_i>/(2k) >= 0.02.
!
!  Reference:
!    [1] T.B. Gatski, C.G. Speziale, On explicit algebraic stress models for complex turbulent flows, J. Fluid Mech. 254, 1993.
!    [2] C.L. Rumsey, T.B. Gatski, Recent Turbulence Model Advances Applied to Multielement Airfoil Computations,
!        J. Aircraft 38(5), 2001.  Model description and constants: Turbulence Modeling Resource, EASMko2003,
!        https://turbmodels.larc.nasa.gov/easmko.html
!
  use types

  implicit none

  ! Closure constants (SSG pressure strain)
  real(dp), parameter :: C10 = 3.4_dp, C11 = 1.8_dp, C2 = 0.36_dp, C3 = 1.25_dp, C4 = 0.4_dp
  real(dp), parameter :: CEPS1 = 1.44_dp, CEPS2 = 1.83_dp
  real(dp), parameter :: CMUSTAR_MIN = 5.0e-4_dp

  real(dp), parameter :: A1 = 0.5_dp*(4.0_dp/3.0_dp - C2)
  real(dp), parameter :: A2 = 0.5_dp*(2.0_dp - C4)
  real(dp), parameter :: A3 = 0.5_dp*(2.0_dp - C3)
  real(dp), parameter :: G0 = 0.5_dp*C11
  real(dp), parameter :: G1 = 0.5_dp*C10 + (CEPS2-CEPS1)/(CEPS1-1.0_dp)

  private
  public :: easm_cell, easm_cmustar, CEPS1, CEPS2

contains

pure function easm_cmustar( eta2, w2, tau ) result(cmus)
!
! Cubic equation for C*mu. eta2 = S_ij S_ij,  w2 = W_ij W_ij,  tau = time scale.
!
  implicit none
  real(dp), intent(in) :: eta2, w2, tau
  real(dp) :: cmus
  real(dp) :: e2t2, r2, p, q, r, aa, bb, d, t1, t2, th, tk(3), sq, pi

  pi = acos(-1.0_dp)
  e2t2 = eta2*tau**2

  if ( e2t2 < 1.0e-6_dp ) then
    cmus = G1*A1/( G1**2 + 2.0_dp*w2*tau**2*A2**2 )
  else
    r2 = w2/eta2
    p = -G1/(e2t2*G0)
    q = ( G1**2 - 2.0_dp*e2t2*G0*A1 - 2.0_dp/3.0_dp*e2t2*A3**2 + 2.0_dp*r2*e2t2*A2**2 )/(2.0_dp*e2t2*G0)**2
    r = G1*A1/(2.0_dp*e2t2*G0)**2

    aa = q - p**2/3.0_dp
    bb = ( 2.0_dp*p**3 - 9.0_dp*p*q + 27.0_dp*r )/27.0_dp
    d = bb**2/4.0_dp + aa**3/27.0_dp

    if ( d > 0.0_dp ) then
      sq = sqrt(d)
      t1 = sign( abs(-bb/2.0_dp + sq)**(1.0_dp/3.0_dp), -bb/2.0_dp + sq )
      t2 = sign( abs(-bb/2.0_dp - sq)**(1.0_dp/3.0_dp), -bb/2.0_dp - sq )
      cmus = -min( -p/3.0_dp + t1 + t2, -p/3.0_dp - t1/2.0_dp - t2/2.0_dp )
    elseif ( aa > -1.0e-30_dp ) then
      cmus = p/3.0_dp                       ! triple root
    else
      th = acos( max( min( -(bb/2.0_dp)/sqrt(-aa**3/27.0_dp), 1.0_dp ), -1.0_dp ) )
      tk(1) = -p/3.0_dp + 2.0_dp*sqrt(-aa/3.0_dp)*cos( th/3.0_dp )
      tk(2) = -p/3.0_dp + 2.0_dp*sqrt(-aa/3.0_dp)*cos( (th + 2.0_dp*pi)/3.0_dp )
      tk(3) = -p/3.0_dp + 2.0_dp*sqrt(-aa/3.0_dp)*cos( (th + 4.0_dp*pi)/3.0_dp )
      cmus = -minval(tk)
    endif
  endif

  cmus = max( cmus, CMUSTAR_MIN )

end function


subroutine easm_cell( g, tau, cmus, aex )
!
! In:  g(i,j) = dU_i/dx_j,  tau = time scale.
! Out: cmus = C*mu,  aex(1:6) = excess anisotropy (xx,xy,xz,yy,yz,zz), momentum-source convention.
!
  implicit none
  real(dp), dimension(3,3), intent(in) :: g
  real(dp), intent(in) :: tau
  real(dp), intent(out) :: cmus
  real(dp), dimension(6), intent(out) :: aex

  real(dp), dimension(3,3) :: s, w, sw, ss, bt
  real(dp) :: eta2, w2, trss, a4, fac, alpha, ev, b
  integer :: i, j, k

  do j=1,3
    do i=1,3
      s(i,j) = 0.5_dp*( g(i,j) + g(j,i) )
      w(i,j) = 0.5_dp*( g(i,j) - g(j,i) )
    enddo
  enddo

  eta2 = sum(s*s)
  w2 = sum(w*w)

  cmus = easm_cmustar( eta2, w2, tau )
  a4 = tau/( G1 + 2.0_dp*G0*cmus*eta2*tau**2 )

  ! SW - WS and SS
  do j=1,3
    do i=1,3
      sw(i,j) = 0.0_dp
      ss(i,j) = 0.0_dp
      do k=1,3
        sw(i,j) = sw(i,j) + s(i,k)*w(k,j) - w(i,k)*s(k,j)
        ss(i,j) = ss(i,j) + s(i,k)*s(k,j)
      enddo
    enddo
  enddo
  trss = ss(1,1) + ss(2,2) + ss(3,3)

  ! Bracket {...} of the constitutive relation (excess part, divided by 2 mut)
  do j=1,3
    do i=1,3
      bt(i,j) = A2*a4*sw(i,j) - 2.0_dp*A3*a4*ss(i,j)
    enddo
    bt(j,j) = bt(j,j) + 2.0_dp*A3*a4*trss/3.0_dp
  enddo

  ! a_ij = -2 C*mu tau bt_ij  (T = 2 mut (S + bt) - 2/3 rho k delta  and  T = 2 mut S - 2/3 rho k delta - rho k a)
  fac = -2.0_dp*cmus*tau

  ! Realizability: <u_i u_i>/k = 2/3 - 2 C*mu tau S_ii + a_ii  >=  0.04   (i.e. <u_i u_i>/(2k) >= 0.02)
  alpha = 1.0_dp
  do i=1,3
    b = 2.0_dp/3.0_dp - 2.0_dp*cmus*tau*s(i,i)       ! Boussinesq part of <u_i u_i>/k
    ev = fac*bt(i,i)                                  ! excess part a_ii
    if ( b + ev < 0.04_dp ) then
      if ( ev < -1.0e-14_dp .and. b > 0.04_dp ) then
        alpha = min( alpha, (b - 0.04_dp)/(-ev) )
      else
        alpha = 0.0_dp
      endif
    endif
  enddo
  alpha = max( alpha, 0.0_dp )

  aex(1) = alpha*fac*bt(1,1)
  aex(2) = alpha*fac*bt(1,2)
  aex(3) = alpha*fac*bt(1,3)
  aex(4) = alpha*fac*bt(2,2)
  aex(5) = alpha*fac*bt(2,3)
  aex(6) = alpha*fac*bt(3,3)

end subroutine

end module
