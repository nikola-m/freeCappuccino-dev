module asimple_wallfun
!
!  Purpose:
!
!    Generalized wall functions for the all-speeds SIMPLE algorithm (module simple_allspeeds). For every 'wall' face the
!    module finds the friction velocity u_tau from the velocity of the wall-adjacent cell and sets the effective wall
!    viscosity visw(iWall) = tau_w d_n / U_t, which is exactly the quantity used by the momentum equations
!    (velocity.f90: vsol = max(viscos,visw)*area/d_n) and by the turbulence models. Hence no other code is changed.
!
!  Discussion:
!
!    For each model the law of the wall is written as an explicit function F(u_tau) = U_t, which is monotone in
!    u_tau, so that u_tau is found by bisection (robust, no derivative, ~50 evaluations per wall face):
!
!    'standard'    Launder-Spalding: u* = Cmu^(1/4) sqrt(k_P), y* = rho u* d_n / mu, tau_w = rho u* kappa U_t / ln(E y*) for
!                  y* > ctrans, viscous sublayer otherwise. Needs k (turbulence model); same as the existing code.
!    'reichardt'   Reichardt's law, valid in the viscous sublayer, buffer layer and log layer, explicit in y+:
!                    u+ = ln(1+kappa y+)/kappa + 7.8 [ 1 - exp(-y+/11) - (y+/11) exp(-y+/3) ]
!                  with y+ = rho_w u_tau d_n / mu_w. Independent of k, valid also for laminar near-wall cells.
!    'spl'         Generalized wall function with pressure gradient in the spirit of Shih, Povinelli & Liu (NASA/TM-1999-209626,
!                  2003): the velocity is the sum U = U1 + U2 of the zero-pressure-gradient law U1 = u_tau u+(y+) (Reichardt)
!                  and the pressure-gradient part U2 = u_p [ alpha ln(u_p d_n / nu) + beta ], with the pressure velocity
!                  u_p = ( (nu/rho) |dp/ds| )^(1/3), alpha = 5.0, beta = 8.0, active in ADVERSE pressure gradient (dp/ds>0,
!                  s along the wall velocity); U2 is clipped at zero. The weights u_tau/(u_tau+u_p) and u_p/(u_tau+u_p) of the
!                  two parts are implied by the sum. See NOTE below.
!    'compressible' Compressible, heat-transfer wall function: van Driest transformed Reichardt law with the Crocco-Busemann
!                  temperature-velocity relation and the wall heat flux from Kader's temperature law.
!                    T/T_w = rho_w/rho = 1 - B u+ - A^2 u+^2 ,   A^2 = r u_tau^2/(2 cp T_w),  B = q_w/(rho_w cp u_tau T_w), r = Pr^(1/3)
!                  (q_w > 0 heat flows from the wall into the fluid). The van Driest effective velocity is
!                    u_vD+ = int_0^u+ sqrt(rho/rho_w) du+ ,  rho_w/rho = 1 + Bc u+ - A^2 u+^2 ,  Bc = -B  (in the code: B := Bc)
!                          = (1/A) [ asin((2A^2 u+ - Bc)/S) + asin(Bc/S) ] ,  S = sqrt(Bc^2+4A^2),
!                  inverted analytically for u+; u_vD+ is given by Reichardt's law in the wall variables (rho_w, mu_w). The wall
!                  heat flux comes from Kader's law T+ = Pr y+ exp(-G) + [2.12 ln(1+y+) + beta(Pr)] exp(-1/G),
!                  G = 0.01 (Pr y+)^4 / (1 + 5 Pr^3 y+), q_w = rho_w cp u_tau (T_w - T_P)/T+ . It is used in the velocity
!                  transformation and stored in qwall(:) for output; it is NOT fed to the energy equation boundary condition.
!
!    NOTE: (1) The pressure-gradient constants of 'spl' are taken from the NASA report as available to the author in text
!    extraction (equations garbled) and not verified against the original paper; the form of U2 for the favourable gradient
!    is not implemented (U2 = 0 there). (2) The model of Ong & Chan (2018) is NOT implemented, its formulation was not
!    accessible. 'compressible' is a classical van Driest + Crocco-Busemann + Kader construction, not their model.
!    (3) The momentum equations use max(viscos, visw), so a wall shear stress below the laminar one (separation) is not possible.
!
  use types
  use parameters
  use geometry
  use variables
  use energy_equation, only: cp_const, sutherland, mu_Tref, mu_S

  implicit none

  private

  character(len=12) :: asimple_wallFunc = 'none'   ! none, standard, reichardt, spl, compressible
  real(dp), allocatable, save :: qwall(:)         ! wall heat flux of the 'compressible' wall function [W/m^2]

  real(dp), parameter :: kappa = 0.41_dp
  real(dp), parameter :: cmu25 = 0.5477225575051661_dp   ! 0.09^(1/4)
  real(dp), parameter :: spl_alpha = 5.0_dp, spl_beta = 8.0_dp

  public :: asimple_wallFunc, qwall, update_wall_viscosity

contains


subroutine update_wall_viscosity
!
!  Loop over wall faces, evaluate the selected wall function, fill visw, ypl and tau.
!
  implicit none

  integer :: ib, i, iface, ijp, ijb, iWall
  real(dp) :: are, nxf, nyf, nzf, dn, upb, vpb, wpb, vn, utx, uty, utz, Ut
  real(dp) :: rho, mu, kk, ustar, ystar, utau, tauw, dpds, up_p, Tw, TP, Pr
  real(dp) :: rhow, muw, qw
  real(dp) :: ypmin, ypmax, ypsum, vrmax, vrsum
  integer :: nw

  if ( trim(asimple_wallFunc) == 'none' ) return
  if ( .not. allocated(qwall) ) then
    allocate( qwall(nwal) );  qwall = 0.0_dp
  endif

  iWall = 0
  ypmin = huge(1.0_dp);  ypmax = 0.0_dp;  ypsum = 0.0_dp;  vrmax = 0.0_dp;  vrsum = 0.0_dp;  nw = 0
  do ib = 1, numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i = 1, nfaces(ib)

      iface = startFace(ib) + i
      ijp = owner(iface)
      ijb = iBndValueStart(ib) + i
      iWall = iWall + 1

      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
      nxf = arx(iface)/are;  nyf = ary(iface)/are;  nzf = arz(iface)/are
      dn  = dnw(iWall)

      ! Tangential velocity relative to the wall
      upb = u(ijp)-u(ijb);  vpb = v(ijp)-v(ijb);  wpb = w(ijp)-w(ijb)
      vn  = upb*nxf + vpb*nyf + wpb*nzf
      utx = upb - vn*nxf;  uty = vpb - vn*nyf;  utz = wpb - vn*nzf
      Ut  = sqrt( utx**2 + uty**2 + utz**2 )

      rho = den(ijp)
      mu  = viscosity_of_T( T(ijp) )
      if ( .not. sutherland ) mu = viscos

      if ( Ut < 1.0e-12_dp ) then
        visw(iWall) = max( viscos, mu )
        tau(iWall)  = 0.0_dp;  ypl(iWall) = 0.0_dp
        cycle
      endif

      select case ( trim(asimple_wallFunc) )

      case ( 'standard' )
        kk = max( te(ijp), 0.0_dp )
        ustar = cmu25*sqrt(kk)
        ystar = rho*ustar*dn/mu
        if ( ystar > ctrans ) then
          tauw = rho*ustar*kappa*Ut/log( ELOG*ystar )
        else
          tauw = mu*Ut/dn
        endif
        utau = sqrt( tauw/rho )

      case ( 'reichardt' )
        utau = solve_utau( 0, Ut, dn, rho, mu, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp )
        tauw = rho*utau**2

      case ( 'spl' )
        ! Streamwise pressure gradient (along the tangential velocity direction)
        dpds = ( dPdxi(1,ijp)*utx + dPdxi(2,ijp)*uty + dPdxi(3,ijp)*utz )/Ut
        up_p = 0.0_dp
        if ( dpds > 0.0_dp ) up_p = ( mu/rho**2*dpds )**( 1.0_dp/3.0_dp )
        utau = solve_utau( 1, Ut, dn, rho, mu, up_p, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp )
        tauw = rho*utau**2

      case ( 'compressible' )
        Tw  = T(ijb);  TP = T(ijp)
        muw = viscosity_of_T( Tw )
        if ( .not. sutherland ) muw = viscos
        rhow = max( p(ijb), small )/( Rair*Tw )
        Pr = pranl
        utau = solve_utau( 2, Ut, dn, rhow, muw, 0.0_dp, Tw, TP, Pr, cp_const )
        tauw = rhow*utau**2
        qw = wall_heat_flux( utau, dn, rhow, muw, Tw, TP, Pr, cp_const )
        qwall(iWall) = qw
        rho = rhow;  mu = muw

      case default
        write(*,'(a)') 'Fatal error: asimple_wallFunc must be none, standard, reichardt, spl or compressible.'
        stop
      end select

      ! Effective wall viscosity: tau_w = visw Ut / dn
      visw(iWall) = max( viscos, tauw*dn/Ut )
      tau(iWall)  = tauw
      ypl(iWall)  = rho*sqrt( tauw/rho )*dn/mu

      nw = nw + 1
      ypmin = min( ypmin, ypl(iWall) );  ypmax = max( ypmax, ypl(iWall) );  ypsum = ypsum + ypl(iWall)
      vrmax = max( vrmax, visw(iWall)/viscos );  vrsum = vrsum + visw(iWall)/viscos

    enddo
  enddo

  if ( ltest .and. nw > 0 ) write(*,'(2x,a,3(1pe10.3),a,2(1pe10.3))') 'wall function: y+ min/mean/max =', &
                                  ypmin, ypsum/nw, ypmax, ',  visw/viscos mean/max =', vrsum/nw, vrmax

end subroutine


function viscosity_of_T( Tl ) result(mu)
  implicit none
  real(dp), intent(in) :: Tl
  real(dp) :: mu
  if ( sutherland ) then
    mu = viscos*( Tl/mu_Tref )**1.5_dp*( mu_Tref + mu_S )/( Tl + mu_S )
  else
    mu = viscos
  endif
end function


function reichardt( yp ) result(up)
!
!  Reichardt's law of the wall, u+(y+).
!
  implicit none
  real(dp), intent(in) :: yp
  real(dp) :: up
  up = log( 1.0_dp + kappa*yp )/kappa + 7.8_dp*( 1.0_dp - exp(-yp/11.0_dp) - (yp/11.0_dp)*exp(-yp/3.0_dp) )
end function


function kader_Tplus( yp, Pr ) result(Tp)
!
!  Kader's temperature law of the wall, T+(y+,Pr).
!
  implicit none
  real(dp), intent(in) :: yp, Pr
  real(dp) :: Tp, G, beta
  G = 0.01_dp*( Pr*yp )**4/( 1.0_dp + 5.0_dp*Pr**3*yp )
  beta = ( 3.85_dp*Pr**( 1.0_dp/3.0_dp ) - 1.3_dp )**2 + 2.12_dp*log(Pr)
  Tp = Pr*yp*exp(-G) + ( 2.12_dp*log( 1.0_dp + yp ) + beta )*exp( -1.0_dp/max(G,1.0e-30_dp) )
end function


function wall_heat_flux( utau, dn, rhow, muw, Tw, TP, Pr, cp ) result(qw)
!
!  Wall heat flux [W/m^2], positive from wall into the fluid:  q_w = rho_w cp u_tau (T_w - T_P)/T+.
!
  implicit none
  real(dp), intent(in) :: utau, dn, rhow, muw, Tw, TP, Pr, cp
  real(dp) :: qw, yp
  yp = rhow*utau*dn/muw
  qw = rhow*cp*utau*( Tw - TP )/max( kader_Tplus( yp, Pr ), 1.0e-12_dp )
end function


function velocity_law( imodel, utau, dn, rho, mu, up_p, Tw, TP, Pr, cp ) result(Ut)
!
!  Tangential velocity predicted by the law of the wall for a given friction velocity.
!
  implicit none
  integer,  intent(in) :: imodel
  real(dp), intent(in) :: utau, dn, rho, mu, up_p, Tw, TP, Pr, cp
  real(dp) :: Ut, yp, upl, U2, nu, A2, A, B, S, uvd, phi0, qw

  yp  = rho*utau*dn/mu
  upl = reichardt( yp )

  select case ( imodel )

  case ( 0 )   ! Reichardt
    Ut = utau*upl

  case ( 1 )   ! Reichardt + pressure gradient part
    nu = mu/rho
    U2 = 0.0_dp
    if ( up_p > 0.0_dp ) U2 = max( 0.0_dp, up_p*( spl_alpha*log( up_p*dn/nu ) + spl_beta ) )
    Ut = utau*upl + U2

  case ( 2 )   ! van Driest transformed, Crocco-Busemann
    qw = wall_heat_flux( utau, dn, rho, mu, Tw, TP, Pr, cp )
    B  = - qw/( rho*cp*utau*Tw )                  ! rho_w/rho = 1 + B u+ - A2 u+^2   (B = -q_w/(rho_w cp u_tau T_w))
    A2 = Pr**( 1.0_dp/3.0_dp )*utau**2/( 2.0_dp*cp*Tw )
    uvd = upl                                      ! van Driest effective velocity u_vD+
    if ( A2 < 1.0e-10_dp ) then
      upl = uvd + B*uvd**2/4.0_dp                  ! small A expansion of the inverse
    else
      A = sqrt(A2);  S = sqrt( B*B + 4.0_dp*A2 );  phi0 = asin( B/S )
      ! u_vD = (1/A)[ asin((2 A2 u+ - B)/S) + phi0 ]  =>  u+ = ( B + S sin(A u_vD - phi0) )/(2 A2),  valid while A u_vD - phi0 <= pi/2
      upl = ( B + S*sin( min( A*uvd - phi0, 0.5_dp*acos(-1.0_dp) ) ) )/( 2.0_dp*A2 )
    endif
    Ut = utau*upl

  case default
    Ut = utau*upl
  end select

end function


function solve_utau( imodel, Ut, dn, rho, mu, up_p, Tw, TP, Pr, cp ) result(utau)
!
!  Bisection (geometric) for u_tau such that velocity_law(u_tau) = Ut. The law is monotone increasing in u_tau.
!
  implicit none
  integer,  intent(in) :: imodel
  real(dp), intent(in) :: Ut, dn, rho, mu, up_p, Tw, TP, Pr, cp
  real(dp) :: utau, lo, hi, mid
  integer :: it

  lo = 1.0e-9_dp*Ut
  hi = 10.0_dp*max( Ut, sqrt( mu*Ut/( rho*dn ) ) )

  if ( velocity_law( imodel, lo, dn, rho, mu, up_p, Tw, TP, Pr, cp ) > Ut ) then
    utau = lo     ! pressure-gradient part alone exceeds Ut: shear at its lower bound (see NOTE 3)
    return
  endif

  do it = 1, 80
    mid = sqrt( lo*hi )
    if ( velocity_law( imodel, mid, dn, rho, mu, up_p, Tw, TP, Pr, cp ) > Ut ) then
      hi = mid
    else
      lo = mid
    endif
    if ( hi/lo < 1.0_dp + 1.0e-10_dp ) exit
  enddo
  utau = sqrt( lo*hi )

end function

end module asimple_wallfun
