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
!    'spl'         Generalized wall function of Shih, Povinelli & Liu (J. Turbulence 4, 2003, N15), eq. (9)-(11), (16):
!                    U_t/u_c = (tau_w/(rho u_tau^2)) (u_tau/u_c) f1(Y_tau+) + (nu/rho) (dP_w/dx_i)/u_c^3 f2(Y_c+)
!                  with u_c = u_tau + u_p, u_p = ((nu/rho)|dP_w/dx_i|)^(1/3), Y_tau+ = u_tau y/nu, Y_c+ = u_c y/nu, and the piecewise
!                  polynomial fits f1, f2 of Tables 1 and 2 (viscous, buffer and inertial sublayer; alpha = 5, beta = 8, kappa = 0.41,
!                  C = 5). The vector form (16) gives the wall shear for every sign of the tangential pressure gradient (also near
!                  separation, where u_tau -> 0): for a given u_tau the shear direction is that of R = U_t - (nu/rho) G f2/u_c^2 and
!                  the consistency |R| = u_tau f1(Y_tau+) determines u_tau (bisection). G is the tangential part of the cell
!                  gradient of pressure (eq. (20) of the paper). Near zero shear the equation can have several roots
!                  (u_tau = 0 is one when U_t = u_p f2(Y_p+)); the largest one (attached branch) is taken.
!    'unified'     Ong & Chan, J. Turbulence 19 (2018) 414-430, 'A unified wall function for compressible turbulence modelling'
!                  (the Shih et al. law above with heat transfer and compressibility). The tangential velocity is split
!                  u_t = u1 + u2: u2 = (1/2) u2+ u_p (grad p/|grad p|) with u2+ from  (y_p+)^2 = u2+ + exp(-2 beta/alpha)
!                  [exp(u2+/alpha) - 1 - u2+/alpha]  (eq. 6, Newton), and u_tau from the Spalding-type law (eq. 7 with the
!                  compressible log-law of eq. 10, N = 5):
!                    y_tau+ = u1+ + exp(-kappa C) [ E(u1+) - sum_{n=0..N} (kappa u1+)^n/n! ],   E = exp(kappa u1+) (incompressible) or
!                    E = exp( (kappa/sqrt(g)) [ asin((2 g u1+ - bq)/Q) - phi ] )   (compressible, isothermal wall; eq. 17), with
!                  g = r u_tau^2/(2 cp T_w), bq = u_tau ((T/T_w - 1)/|u1| + r |u1|/(2 cp T_w)), Q = sqrt(bq^2 + 4 g), phi = asin(-bq/Q)
!                  (isothermal wall, T_w = T of the wall face; eq. 17 and 35) or, for an adiabatic wall, T_w = T + r u1^2/(2 cp) and
!                  exp((kappa/sqrt(g)) asin(sqrt(g) u1+)) (eq. 23). The wall is isothermal when bcDFraction(7,ib) = 1 (Dirichlet
!                  T of 0/T), adiabatic otherwise. With compressible = F the law reduces to the incompressible form (g = bq = 0).
!                  The shear stress is tau_w = rho_w u_tau^2 u1/|u1|, the wall heat flux (isothermal) q_w = bq T_w cp rho_w u_tau/Pr.
!                  Differences to the paper's algorithm: u_tau is found by bisection (the function is monotone) instead of
!                  Steffensen's method; and the derivative of eq. (29) is the analytical one, dg/du2+ = 1 + exp(-2 beta/alpha)
!                  (exp(u2+/alpha) - 1)/alpha (the printed equation contains an additional term 2 (y_p+)^2/u2+).
!    'compressible' van Driest transformed Reichardt law with the Crocco-Busemann temperature-velocity relation and the wall
!                  heat flux from Kader's temperature law (classical construction, not the model of Ong & Chan):
!                    T/T_w = rho_w/rho = 1 - B u+ - A^2 u+^2 ,   A^2 = r u_tau^2/(2 cp T_w),  B = q_w/(rho_w cp u_tau T_w), r = Pr^(1/3)
!                  (q_w > 0 heat flows from the wall into the fluid). The van Driest effective velocity is
!                    u_vD+ = int_0^u+ sqrt(rho/rho_w) du+ ,  rho_w/rho = 1 + Bc u+ - A^2 u+^2 ,  Bc = -B  (in the code: B := Bc)
!                          = (1/A) [ asin((2A^2 u+ - Bc)/S) + asin(Bc/S) ] ,  S = sqrt(Bc^2+4A^2),
!                  inverted analytically for u+; u_vD+ is given by Reichardt's law in the wall variables (rho_w, mu_w). The wall
!                  heat flux comes from Kader's law T+ = Pr y+ exp(-G) + [2.12 ln(1+y+) + beta(Pr)] exp(-1/G),
!                  G = 0.01 (Pr y+)^4 / (1 + 5 Pr^3 y+), q_w = rho_w cp u_tau (T_w - T_P)/T+ .
!                  The heat flux of 'compressible' and 'unified' is stored in qwall(:) for output; it is NOT fed to the energy
!                  equation (which keeps the laminar conduction flux k Df (T_w - T_P) at isothermal walls).
!    NOTE: (1) The momentum equations use max(viscos, visw), so a wall shear stress below the laminar one, or of the opposite
!    sign to the velocity (reversed shear), cannot be represented; the shear is projected on the direction of the velocity.
!    (2) The pressure gradient of the wall cell is the cell-centre gradient (as in the paper), no wall-normal pressure
!      variation.
!    (3) The laws assume attached or mildly separated boundary layers: shocks impinging on the wall are outside their range.
!
  use types
  use parameters
  use geometry
  use variables
  use energy_equation, only: cp_const, sutherland, mu_Tref, mu_S

  implicit none

  private

  character(len=12) :: asimple_wallFunc = 'none'   ! none, standard, reichardt, spl (Shih et al.), unified (Ong & Chan), compressible
  real(dp), allocatable, save :: qwall(:)         ! wall heat flux of the 'compressible' wall function [W/m^2]

  real(dp), parameter :: kappa = 0.41_dp
  real(dp), parameter :: cmu25 = 0.5477225575051661_dp   ! 0.09^(1/4)
  real(dp), parameter :: Ccon  = 5.0_dp                 ! log-law constant C of Shih et al. and Ong & Chan
  real(dp), parameter :: wf_alpha = 5.0_dp, wf_beta = 8.0_dp   ! Stratford constants of the pressure gradient part
  integer,  parameter :: Nser = 5                       ! number of terms of the series in the Spalding-type law (Ong & Chan)

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
  real(dp) :: rhow, muw, qw, Gt(3), Rv(3), uv(3), pw, Twall, Tc, rr, dxx, dyy, dzz, projf
  logical :: isothermal
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
        ! Shih, Povinelli & Liu: tangential part of the pressure gradient of the wall cell and the velocity vector
        Gt = dPdxi(1:3,ijp) - ( dPdxi(1,ijp)*nxf + dPdxi(2,ijp)*nyf + dPdxi(3,ijp)*nzf )*(/ nxf, nyf, nzf /)
        uv = (/ utx, uty, utz /)
        call shih_solve( uv, Gt, dn, rho, mu, utau, Rv )
        projf = ( Rv(1)*utx + Rv(2)*uty + Rv(3)*utz )/( max( sqrt(sum(Rv**2)), 1.0e-300_dp )*Ut )
        tauw = rho*utau**2*projf          ! shear stress projected on the direction of the velocity (may be < 0: reversed shear)

      case ( 'unified' )
        ! Ong & Chan: wall pressure by linear extrapolation from the cell centre, isothermal or adiabatic wall
        Gt = dPdxi(1:3,ijp) - ( dPdxi(1,ijp)*nxf + dPdxi(2,ijp)*nyf + dPdxi(3,ijp)*nzf )*(/ nxf, nyf, nzf /)
        uv = (/ utx, uty, utz /)
        dxx = xf(iface)-xc(ijp);  dyy = yf(iface)-yc(ijp);  dzz = zf(iface)-zc(ijp)
        pw = p(ijp) + dPdxi(1,ijp)*dxx + dPdxi(2,ijp)*dyy + dPdxi(3,ijp)*dzz
        isothermal = ( bcDFraction(7,ib) == 1.0_dp )
        Tc = 300.0_dp;  Twall = 300.0_dp
        if ( compressible ) then
          Tc = T(ijp);  Twall = T(ijb)
        else
          pw = 101325.0_dp                   ! not used (incompressible): density of the cell, viscosity viscos
          isothermal = .false.
        endif
        call unified_solve( isothermal, compressible, uv, Gt, dn, pw, Tc, Twall, rho, utau, Rv, rhow, muw, qw )
        projf = ( Rv(1)*utx + Rv(2)*uty + Rv(3)*utz )/Ut
        tauw = rhow*utau**2*projf
        qwall(iWall) = qw
        rho = rhow;  mu = muw

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
        write(*,'(a)') 'Fatal error: asimple_wallFunc must be none, standard, reichardt, spl, unified or compressible.'
        stop
      end select

      ! Effective wall viscosity: tau_w = visw Ut / dn
      visw(iWall) = max( viscos, tauw*dn/Ut )
      tau(iWall)  = tauw
      ypl(iWall)  = rho*utau*dn/mu

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
  real(dp) :: Ut, yp, upl, A2, A, B, S, uvd, phi0, qw

  yp  = rho*utau*dn/mu
  upl = reichardt( yp )

  select case ( imodel )

  case ( 0 )   ! Reichardt
    Ut = utau*upl

  case ( 2 )   ! van Driest transformed, Crocco-Busemann
    qw = wall_heat_flux( utau, dn, rho, mu, Tw, TP, Pr, cp )
    B  = - qw/( rho*cp*utau*Tw )                  ! rho_w/rho = 1 + B u+ - A2 u+^2   (B = -q_w/(rho_w cp u_tau T_w))
    A2 = Pr**( 1.0_dp/3.0_dp )*utau**2/( 2.0_dp*cp*Tw )
    uvd = upl                                      ! van Driest effective velocity u_vD+
    if ( A2 < 1.0e-10_dp ) then
      upl = uvd + B*uvd**2/4.0_dp                  ! small A expansion of the inverse
    else
      A = sqrt(A2);  S = sqrt( B*B + 4.0_dp*A2 );  phi0 = asin( B/S )
      ! u_vD = (1/A)[ asin((2 A2 u+ - B)/S) + phi0 ]  =>  u+ = ( B + S sin(A u_vD - phi0) )/(2 A2),  valid while A u_vD - phi0
      !   <= pi/2
      upl = ( B + S*sin( min( A*uvd - phi0, 0.5_dp*acos(-1.0_dp) ) ) )/( 2.0_dp*A2 )
    endif
    Ut = utau*upl

  case default
    Ut = utau*upl
  end select

end function



function f1_shih( Y ) result(f)
!
!  Piecewise polynomial fit of the friction part of the generalized wall function, Shih et al. (2003), eq. (10), Table 1.
!
  implicit none
  real(dp), intent(in) :: Y
  real(dp) :: f
  if ( Y <= 5.0_dp ) then
    f = Y*( 1.0_dp + Y*( 1.0e-2_dp - 2.9e-3_dp*Y ) )
  elseif ( Y <= 30.0_dp ) then
    f = -0.872_dp + Y*( 1.465_dp + Y*( -7.02e-2_dp + Y*( 1.66e-3_dp - 1.495e-5_dp*Y ) ) )
  elseif ( Y <= 140.0_dp ) then
    f = 8.6_dp + Y*( 0.1864_dp + Y*( -2.006e-3_dp + Y*( 1.144e-5_dp - 2.551e-8_dp*Y ) ) )
  else
    f = log(Y)/kappa + Ccon
  endif
end function


function f2_shih( Y ) result(f)
!
!  Piecewise polynomial fit of the pressure gradient part, Shih et al. (2003), eq. (11), Table 2.
!
  implicit none
  real(dp), intent(in) :: Y
  real(dp) :: f
  if ( Y <= 4.0_dp ) then
    f = Y*Y*( 0.5_dp - 7.31e-3_dp*Y )
  elseif ( Y <= 15.0_dp ) then
    f = -15.138_dp + Y*( 8.4688_dp + Y*( -0.81976_dp + Y*( 3.7292e-2_dp - 6.3866e-4_dp*Y ) ) )
  elseif ( Y <= 30.0_dp ) then
    f = 11.925_dp + Y*( 0.93400_dp + Y*( -2.7805e-2_dp + Y*( 4.6262e-4_dp - 3.1442e-6_dp*Y ) ) )
  else
    f = wf_alpha*log(Y) + wf_beta
  endif
end function


subroutine shih_residual( utau, Utv, Gt, dn, rho, mu, up_p, F, Rv )
!
!  Residual of the generalized wall function for a given friction velocity:  F = u_tau f1(Y_tau+) - |R|,
!  R = U_t - (nu/rho) G f2(Y_c+)/u_c^2   (eq. 9 and 16 of Shih et al. 2003); the wall shear is tau_w = rho u_tau^2 R/|R|.
!
  implicit none
  real(dp), intent(in)  :: utau, Utv(3), Gt(3), dn, rho, mu, up_p
  real(dp), intent(out) :: F, Rv(3)
  real(dp) :: nu, uc, Yt, Yc, coef

  nu  = mu/rho
  uc  = utau + up_p
  Yt  = utau*dn/nu
  Yc  = uc*dn/nu
  coef = 0.0_dp
  if ( up_p > 0.0_dp ) coef = ( mu/rho**2 )*f2_shih( Yc )/uc**2
  Rv = Utv - coef*Gt
  F  = utau*f1_shih( Yt ) - sqrt( sum(Rv**2) )

end subroutine


subroutine shih_solve( Utv, Gt, dn, rho, mu, utau, Rv )
!
!  Friction velocity of the generalized wall function (Shih et al. 2003) by geometric bisection; Rv is the vector whose
!  direction is the direction of the wall shear stress.
!
  implicit none
  real(dp), intent(in)  :: Utv(3), Gt(3), dn, rho, mu
  real(dp), intent(out) :: utau, Rv(3)
  real(dp) :: Ut, up_p, lo, hi, mid, F, Flo, R1(3), xa, xb, rscan
  integer :: it, k
  integer, parameter :: nscan = 80

  Ut   = sqrt( sum(Utv**2) )
  up_p = ( mu/rho**2*sqrt( sum(Gt**2) ) )**( 1.0_dp/3.0_dp )

  lo = 1.0e-10_dp*max( Ut, up_p )
  hi = 10.0_dp*max( Ut, up_p, sqrt( mu*Ut/( rho*dn ) ) )
  call shih_residual( hi, Utv, Gt, dn, rho, mu, up_p, F, Rv )
  do it = 1, 30
    if ( F > 0.0_dp ) exit
    hi = 2.0_dp*hi
    call shih_residual( hi, Utv, Gt, dn, rho, mu, up_p, F, Rv )
  enddo
  call shih_residual( lo, Utv, Gt, dn, rho, mu, up_p, Flo, R1 )
  if ( Flo > 0.0_dp ) then          ! |R| = 0 at u_tau -> 0: zero shear stress (separation point)
    utau = lo;  Rv = Utv;  return
  endif

  ! F(u_tau) can have several roots near zero shear (u_tau = 0 is always one at U_t = u_p f2); the LARGEST one is taken (attached
  ! branch, continuous with the zero pressure gradient limit): scan downwards from hi for the first sign change.
  xb = hi
  rscan = ( lo/hi )**( 1.0_dp/real(nscan,dp) )
  do k = 1, nscan
    xa = hi*rscan**k
    call shih_residual( xa, Utv, Gt, dn, rho, mu, up_p, F, Rv )
    if ( F < 0.0_dp ) then
      lo = xa;  hi = xb
      exit
    endif
    xb = xa
  enddo

  do it = 1, 100
    mid = sqrt( lo*hi )
    call shih_residual( mid, Utv, Gt, dn, rho, mu, up_p, F, Rv )
    if ( F > 0.0_dp ) then
      hi = mid
    else
      lo = mid
    endif
    if ( hi/lo < 1.0_dp + 1.0e-12_dp ) exit
  enddo
  utau = sqrt( lo*hi )
  call shih_residual( utau, Utv, Gt, dn, rho, mu, up_p, F, Rv )
  if ( sqrt( sum(Rv**2) ) < 1.0e-300_dp ) Rv = Utv

end subroutine


function unified_f( utau, ityp, u1m, dn, nu, Tw, TP, cpp, rr ) result(f)
!
!  Residual of the unified law of the wall of Ong & Chan, eq. (31) (isothermal: ityp = 1) / eq. (42) (adiabatic: ityp = 2) /
!  incompressible, eq. (7) (ityp = 0):  f = -Y_tau+ + u1+ + exp(-kappa C) [ E(u1+) - sum_{n=0}^{N} (kappa u1+)^n / n! ].
!  Monotone decreasing in u_tau.
!
  implicit none
  real(dp), intent(in) :: utau, u1m, dn, nu, Tw, TP, cpp, rr
  integer,  intent(in) :: ityp
  real(dp) :: f, u1p, gam, bq, Q, phi, arg, ex, ser, term, ek
  integer :: n

  u1p = u1m/utau
  ek  = exp( -kappa*Ccon )

  select case ( ityp )
  case ( 1 )
    gam = rr*utau**2/( 2.0_dp*cpp*Tw )
    bq  = utau*( ( TP/Tw - 1.0_dp )/u1m + rr*u1m/( 2.0_dp*cpp*Tw ) )
    if ( gam > 1.0e-10_dp ) then
      Q   = sqrt( bq*bq + 4.0_dp*gam )
      phi = asin( max( -1.0_dp, min( 1.0_dp, -bq/Q ) ) )
      arg = max( -1.0_dp, min( 1.0_dp, ( 2.0_dp*gam*u1p - bq )/Q ) )
      ex  = kappa/sqrt(gam)*( asin(arg) - phi )
    elseif ( abs(bq) > 1.0e-10_dp .and. 1.0_dp + bq*u1p > 0.0_dp ) then
      ex  = 2.0_dp*kappa/bq*( sqrt( 1.0_dp + bq*u1p ) - 1.0_dp )       ! heat transfer only, eq. (19)
    else
      ex  = kappa*u1p
    endif
  case ( 2 )
    gam = rr*utau**2/( 2.0_dp*cpp*Tw )
    if ( gam > 1.0e-14_dp ) then
      ex = kappa/sqrt(gam)*asin( max( -1.0_dp, min( 1.0_dp, sqrt(gam)*u1p ) ) )
    else
      ex = kappa*u1p
    endif
  case default
    ex = kappa*u1p
  end select

  ser = 0.0_dp;  term = 1.0_dp
  do n = 0, Nser
    if ( n > 0 ) term = term*kappa*u1p/real(n,dp)
    ser = ser + term
  enddo

  f = -utau*dn/nu + u1p + ek*( exp( min( ex, 700.0_dp ) ) - ser )

end function


subroutine unified_solve( isoth, comp, Utv, Gt, dn, pw, TP, Twb, rhoP, utau, e1, rhow, muw, qw )
!
!  Unified wall function of Ong & Chan (2018), algorithms of Section 3.1 (isothermal) and 3.2 (adiabatic), see the header.
!  Input:  tangential velocity Utv and tangential pressure gradient Gt of the wall cell, wall distance dn, wall pressure pw,
!          cell temperature TP, wall face temperature Twb (isothermal wall), cell density rhoP (incompressible: constant).
!  Output: friction velocity utau, unit vector e1 of the wall shear direction, wall density and viscosity, wall heat flux
!          qw (W/m^2, positive from the wall into the fluid; zero for adiabatic walls and incompressible flow).
!
  implicit none
  logical,  intent(in)  :: isoth, comp
  real(dp), intent(in)  :: Utv(3), Gt(3), dn, pw, TP, Twb, rhoP
  real(dp), intent(out) :: utau, e1(3), rhow, muw, qw
  real(dp) :: Ut, Gm, ghat(3), up_p, yp, s, g, gp, ec, rr, cpp, Pr, Tw, nu, u2(3), u1(3), u1m, lo, hi, mid, f
  real(dp) :: gam, bq
  integer :: it, ipass, npass, ityp

  Pr  = pranl
  rr  = Pr**( 1.0_dp/3.0_dp )
  cpp = cp_const
  Ut  = sqrt( sum(Utv**2) )
  Gm  = sqrt( sum(Gt**2) )
  ghat = 0.0_dp
  if ( Gm > 0.0_dp ) ghat = Gt/Gm
  ec = exp( -2.0_dp*wf_beta/wf_alpha )
  qw = 0.0_dp

  ! Wall state
  Tw = TP
  ityp = 0
  npass = 1
  if ( comp ) then
    if ( isoth ) then
      Tw = Twb;  ityp = 1
    else
      Tw = TP;   ityp = 2;  npass = 2        ! T_w = T + r u1^2/(2 cp) depends on u1, which depends on the properties at T_w
    endif
  endif

  do ipass = 1, npass

    if ( comp ) then
      rhow = pw/( Rair*Tw )
      muw  = viscosity_of_T( Tw )
    else
      rhow = rhoP
      muw  = viscos
    endif
    nu = muw/rhow

    ! Pressure gradient part u2: Newton-Raphson for u2+ of eq. (6)
    up_p = ( muw/rhow**2*Gm )**( 1.0_dp/3.0_dp )
    u2 = 0.0_dp
    if ( up_p > 1.0e-12_dp*max( Ut, 1.0e-12_dp ) ) then
      yp = up_p*dn/nu
      s  = yp*yp
      if ( yp > 1.0_dp ) s = min( s, 2.0_dp*wf_alpha*log(yp) + 2.0_dp*wf_beta + 1.0_dp )   ! start to the right of the root
      do it = 1, 100
        g  = s - yp*yp + ec*( exp( s/wf_alpha ) - 1.0_dp - s/wf_alpha )
        gp = 1.0_dp + ec*( exp( s/wf_alpha ) - 1.0_dp )/wf_alpha
        s  = s - g/gp
        if ( abs(g) < 1.0e-12_dp*max( 1.0_dp, yp*yp ) ) exit
      enddo
      u2 = 0.5_dp*s*up_p*ghat
    endif
    u1  = Utv - u2
    u1m = sqrt( sum(u1**2) )

    if ( ityp == 2 ) Tw = TP + rr*u1m**2/( 2.0_dp*cpp )

  enddo

  if ( ityp == 2 ) then
    rhow = pw/( Rair*Tw );  muw = viscosity_of_T( Tw );  nu = muw/rhow
  endif

  ! Friction velocity: bisection, f is monotone decreasing in u_tau
  if ( u1m < 1.0e-12_dp*max( Ut, 1.0e-12_dp ) .or. u1m < 1.0e-300_dp ) then
    utau = 0.0_dp
    e1 = Utv/max( Ut, 1.0e-300_dp )
    return
  endif
  e1 = u1/u1m

  lo = 1.0e-8_dp*u1m
  hi = 10.0_dp*max( u1m, sqrt( nu*u1m/dn ) )
  do it = 1, 100
    mid = sqrt( lo*hi )
    f = unified_f( mid, ityp, u1m, dn, nu, Tw, TP, cpp, rr )
    if ( f > 0.0_dp ) then
      lo = mid
    else
      hi = mid
    endif
    if ( hi/lo < 1.0_dp + 1.0e-12_dp ) exit
  enddo
  utau = sqrt( lo*hi )

  ! Heat flux of the isothermal wall, eq. (12) with k_w u_tau/nu_w = cp rho_w u_tau/Pr:  q (into the wall) = bq T_w cp
  !   rho_w u_tau/Pr
  if ( ityp == 1 ) then
    gam = rr*utau**2/( 2.0_dp*cpp*Tw )
    bq  = utau*( ( TP/Tw - 1.0_dp )/u1m + rr*u1m/( 2.0_dp*cpp*Tw ) )
    qw  = -bq*Tw*cpp*rhow*utau/Pr
  endif

end subroutine

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
