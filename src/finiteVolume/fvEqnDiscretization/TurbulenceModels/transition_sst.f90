module transition_sst
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Langtry-Menter four-equation transitional SST model (gamma - Re_theta_t), 2009, on top of 'turb_transport'.
!
!  Discussion:
!
!    Four scalars, TurbModel%Scalar(1:4) = k, omega, gamma (intermittency), Re_theta_t~ (transition onset momentum
!    thickness Reynolds number, transported):
!
!      k and omega:  SST-2003 (Menter, Kuntz, Langtry) with   Pk~ = gamma_eff Pk,   Dk~ = min( max(gamma_eff,0.1), 1 ) Dk,
!                    F1 = max( F1, F3 ),  F3 = exp( -(Ry/120)^8 ),  Ry = rho y sqrt(k)/mu
!
!      Dgamma/Dt = Fl ca1 rho S sqrt(gamma Fonset) (1 - ce1 gamma) - ca2 rho W gamma Fturb (ce2 gamma - 1)
!                  + div[ (mu + mut/sigma_f) grad gamma ]
!
!      DReth/Dt  = c_theta rho/t ( Re_theta_t - Reth~ ) ( 1 - F_theta_t ) + div[ sigma_theta (mu + mut) grad Reth~ ],  t = 500 mu/(rho U^2)
!
!    Re_theta_t is the free-stream correlation of Langtry & Menter (function of turbulence intensity and pressure gradient
!    parameter, solved by a short local fixed-point iteration). The separation-induced transition is included
!    (gamma_sep). Wall: zero-flux for gamma and Reth~. Mesh must resolve the viscous sublayer (y+ ~ 1).
!
!    The correlations are the published ones (Langtry & Menter, AIAA J. 47(12), 2009, Eqs. 1-47 and Appendix), entered
!    from the paper text. Checked only on a T3A-like flat plate (examples/turbulentAirfoil/mkplate.py): laminar Cf ~ Blasius,
!    transition onset near Re_x = 1.4e5, turbulent Cf of the flat plate correlation by Re_x ~ 4e5; NOT compared with the ERCOFTAC data.
!
!    Mean flow side: nothing changes, the effective viscosity is the SST one (modify_mu_eff of k_omega_SST,
!    automatic wall treatment).
!
!  Reference:
!    [1] R.B. Langtry, F.R. Menter, Correlation-Based Transition Modeling for Unstructured Parallelized CFD Codes,
!        AIAA Journal, Vol. 47, No. 12, 2009, pp. 2894-2906.
!    [2] F.R. Menter, R. Langtry, S. Volker, Transition modelling for general purpose CFD codes, Flow Turbulence Combust. 77, 2006.
!    [3] F.R. Menter, M. Kuntz, R. Langtry, Ten Years of Industrial Experience with the SST Turbulence Model, 2003.
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
  use k_omega_SST, only: modify_mu_eff

  implicit none

  ! SST constants (same as in k_omega_SST)
  real(dp), parameter :: BETTAST = 0.09_dp
  real(dp), parameter :: SIGMK1 = 0.85_dp,  SIGMK2 = 1.0_dp
  real(dp), parameter :: SIGMOM1 = 0.5_dp,  SIGMOM2 = 0.856_dp
  real(dp), parameter :: BETAI1 = 0.075_dp, BETAI2 = 0.0828_dp
  real(dp), parameter :: ALPHA1 = 5.0_dp/9.0_dp, ALPHA2 = 0.44_dp
  real(dp), parameter :: CMU25 = sqrt(sqrt(BETTAST))

  ! Transition model constants
  real(dp), parameter :: CE1 = 1.0_dp,  CA1 = 2.0_dp
  real(dp), parameter :: CE2 = 50.0_dp, CA2 = 0.06_dp
  real(dp), parameter :: SIGMA_F = 1.0_dp
  real(dp), parameter :: C_THETA = 0.03_dp
  real(dp), parameter :: SIGMA_THETA = 2.0_dp
  real(dp), parameter :: S1 = 2.0_dp

  integer, parameter :: iK = 1, iOm = 2, iGam = 3, iRet = 4

  real(dp), dimension(:), allocatable :: gam, gamo, gamoo       ! intermittency
  real(dp), dimension(:), allocatable :: ret, reto, retoo       ! transported Re_theta_t
  real(dp), dimension(:,:), allocatable :: dGamdxi, dRetdxi
  real(dp), dimension(:), allocatable :: f1, gameff, fthetat, retcorr, pk_sst, dif, suv, spv
  logical, dimension(:), allocatable :: fixed
  real(dp), dimension(:), allocatable :: fixval
  integer :: stampg = -1, stampr = -1

  private

  public :: modify_viscosity_transition_sst
  public :: modify_viscosity_inlet_transition_sst

contains


subroutine modify_viscosity_transition_sst
  implicit none

  if (.not.allocated(gam)) call allocate_and_initialise

  call shift_time_levels( gam, gamo, gamoo, stampg )
  call shift_time_levels( ret, reto, retoo, stampr )
  call transition_quantities
  call calcsc_gamma
  call calcsc_ret
  call calcsc_tke
  call calcsc_omega
  call modify_mu_eff

end subroutine


subroutine allocate_and_initialise
  implicit none
  integer :: inp, ib, i, ijb, iface

  allocate( gam(numTotal), gamo(numTotal), gamoo(numTotal), ret(numTotal), reto(numTotal), retoo(numTotal) )
  allocate( dGamdxi(3,numTotal), dRetdxi(3,numTotal) )
  allocate( f1(numCells), gameff(numCells), fthetat(numCells), retcorr(numCells), pk_sst(numCells) )
  allocate( dif(numTotal), suv(numCells), spv(numCells), fixed(numCells), fixval(numCells) )

  gam = 1.0_dp
  do inp=1,numCells
    ret(inp) = retheta_freestream( te(inp), sqrt(u(inp)**2+v(inp)**2+w(inp)**2), 0.0_dp )
  enddo
  do ib=1,numBoundaries
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijb = iBndValueStart(ib) + i
      ret(ijb) = ret(owner(iface))
    enddo
  enddo
  gamo = gam; gamoo = gam
  reto = ret; retoo = ret

  write(*,'(a)') ' '
  write(*,'(a)') '  **Allocated and initialised intermittency and Re_theta_t fields for the transition SST model.'
  write(*,'(a)') ' '

end subroutine


pure function retheta_freestream( k, umag, lam ) result(ret_eq)
!
! Free-stream Re_theta_t correlation of Langtry-Menter (Tu in percent, lam = pressure gradient parameter).
!
  implicit none
  real(dp), intent(in) :: k, umag, lam
  real(dp) :: ret_eq, tu, f

  tu = max( 100.0_dp*sqrt(2.0_dp*max(k,0.0_dp)/3.0_dp)/max(umag,1e-10_dp), 0.027_dp )

  if ( lam <= 0.0_dp ) then
    f = 1.0_dp - ( -12.986_dp*lam - 123.66_dp*lam**2 - 405.689_dp*lam**3 )*exp( -(tu/1.5_dp)**1.5_dp )
  else
    f = 1.0_dp + 0.275_dp*( 1.0_dp - exp(-35.0_dp*lam) )*exp( -tu/0.5_dp )
  endif

  if ( tu <= 1.3_dp ) then
    ret_eq = ( 1173.51_dp - 589.428_dp*tu + 0.2196_dp/tu**2 )*f
  else
    ret_eq = 331.50_dp*( tu - 0.5658_dp )**(-0.671_dp)*f
  endif

  ret_eq = max( ret_eq, 20.0_dp )

end function


pure function retheta_crit( ret ) result(rec)
!
! Critical Reynolds number as function of Re_theta_t.
!
  implicit none
  real(dp), intent(in) :: ret
  real(dp) :: rec

  if ( ret <= 1870.0_dp ) then
    rec = ret - ( 396.035e-2_dp - 120.656e-4_dp*ret + 868.230e-6_dp*ret**2 - 696.506e-9_dp*ret**3 + 174.105e-12_dp*ret**4 )
  else
    rec = ret - ( 593.11_dp + 0.482_dp*(ret - 1870.0_dp) )
  endif

end function


pure function flength( ret ) result(fl)
!
! Transition length function (before the sublayer correction).
!
  implicit none
  real(dp), intent(in) :: ret
  real(dp) :: fl

  if ( ret < 400.0_dp ) then
    fl = 398.189e-1_dp - 119.270e-4_dp*ret - 132.567e-6_dp*ret**2
  elseif ( ret < 596.0_dp ) then
    fl = 263.404_dp - 123.939e-2_dp*ret + 194.548e-5_dp*ret**2 - 101.695e-8_dp*ret**3
  elseif ( ret < 1200.0_dp ) then
    fl = 0.5_dp - ( ret - 596.0_dp )*3.0e-4_dp
  else
    fl = 0.3188_dp
  endif

end function


subroutine transition_quantities
!
! Cell-wise: SST blending F1 (with F3), free-stream Re_theta_t (with pressure gradient parameter), F_theta_t,
! separation-induced intermittency, effective intermittency. Also unmodified SST production.
!
  implicit none
  integer :: inp, iter
  real(dp) :: nu, d, umag, dUds, lam, theta, retf, rt, rv, rew, rey
  real(dp) :: fwake, delta, deltaBL, thetaBL, ft1, rec, fsep, fre, gsep
  real(dp) :: cdkw, arg1, ksi, domegapl

  do inp=1,numCells

    nu = viscos/den(inp)
    d = walldistance(inp)
    umag = sqrt(u(inp)**2+v(inp)**2+w(inp)**2)

    ! SST F1, plus F3 of the transition model
    domegapl = max( 2*SIGMOM2*den(inp)/(ed(inp)+small)*sum(dTEdxi(:,inp)*dEDdxi(:,inp)), 1e-20_dp )
    ksi = min( max( sqrt(te(inp))/(BETTAST*d*ed(inp)+small), 500.0_dp*nu/(d**2*ed(inp)+small) ), &
               4.0_dp*den(inp)*te(inp)*SIGMOM2/(domegapl*d**2) )
    rey = den(inp)*d*sqrt(te(inp))/viscos
    f1(inp) = max( tanh(ksi**4), exp( -(rey/120.0_dp)**8 ) )

    ! Free-stream Re_theta_t with local pressure gradient parameter (fixed-point iteration)
    if ( umag > 1e-8_dp ) then
      dUds = ( u(inp)*( u(inp)*dUdxi(1,inp)+v(inp)*dVdxi(1,inp)+w(inp)*dWdxi(1,inp) ) &
             + v(inp)*( u(inp)*dUdxi(2,inp)+v(inp)*dVdxi(2,inp)+w(inp)*dWdxi(2,inp) ) &
             + w(inp)*( u(inp)*dUdxi(3,inp)+v(inp)*dVdxi(3,inp)+w(inp)*dWdxi(3,inp) ) )/umag**2
      retf = retheta_freestream( te(inp), umag, 0.0_dp )
      do iter=1,10
        theta = retf*nu/umag
        lam = max( min( theta**2/nu*dUds, 0.1_dp ), -0.1_dp )
        retf = retheta_freestream( te(inp), umag, lam )
      enddo
    else
      retf = retheta_freestream( te(inp), umag, 0.0_dp )
    endif
    retcorr(inp) = retf

    ! F_theta_t - switches off the source term in the boundary layer
    rew = den(inp)*ed(inp)*d**2/viscos
    fwake = exp( -(rew/1.0e5_dp)**2 )
    thetaBL = ret(inp)*viscos/(den(inp)*max(umag,small))
    deltaBL = 7.5_dp*thetaBL
    delta = 50.0_dp*vorticity(inp)*d/max(umag,small)*deltaBL
    ft1 = fwake*exp( -(d/max(delta,small))**4 )
    fthetat(inp) = min( max( ft1, 1.0_dp - ((gam(inp)-1.0_dp/CE2)/(1.0_dp-1.0_dp/CE2))**2 ), 1.0_dp )

    ! Separation-induced transition
    rv = den(inp)*d**2*magStrain(inp)/viscos
    rt = den(inp)*te(inp)/(viscos*ed(inp)+small)
    rec = retheta_crit( ret(inp) )
    fre = exp( -(rt/20.0_dp)**4 )
    gsep = min( S1*max( 0.0_dp, rv/(3.235_dp*max(rec,1.0_dp)) - 1.0_dp )*fre, 2.0_dp )*fthetat(inp)
    gameff(inp) = max( gam(inp), gsep )

  enddo

end subroutine


subroutine calcsc_gamma
!
! Intermittency equation, with transition onset function.
!
  implicit none
  integer :: inp
  real(dp) :: d, rv, rt, rec, fon1, fon2, fon3, fonset, fturb, fsub, rw, flen, pg1, nu

  do inp=1,numCells

    nu = viscos/den(inp)
    d = walldistance(inp)

    rv = den(inp)*d**2*magStrain(inp)/viscos
    rt = den(inp)*te(inp)/(viscos*ed(inp)+small)
    rec = retheta_crit( ret(inp) )

    fon1 = rv/(2.193_dp*max(rec,1.0_dp))
    fon2 = min( max( fon1, fon1**4 ), 2.0_dp )
    fon3 = max( 1.0_dp - (rt/2.5_dp)**3, 0.0_dp )
    fonset = max( fon2 - fon3, 0.0_dp )

    fturb = exp( -(rt/4.0_dp)**4 )

    rw = den(inp)*d**2*ed(inp)/(500.0_dp*viscos)
    fsub = exp( -(rw/0.4_dp)**2 )
    flen = flength(ret(inp))*(1.0_dp-fsub) + 40.0_dp*fsub

    ! Pgamma1 (1 - ce1 gamma) : explicit part with lagged gamma under the root, implicit sink.
    pg1 = flen*CA1*den(inp)*magStrain(inp)*sqrt( max(gam(inp)*fonset, 0.0_dp) )
    suv(inp) = pg1*vol(inp)
    spv(inp) = CE1*pg1*vol(inp)

    ! -E_gamma = ca2 rho W Fturb gamma - ce2 ca2 rho W Fturb gamma^2
    suv(inp) = suv(inp) + CA2*den(inp)*vorticity(inp)*fturb*gam(inp)*vol(inp)
    spv(inp) = spv(inp) + CE2*CA2*den(inp)*vorticity(inp)*fturb*gam(inp)*vol(inp)

  enddo

  dif = viscos + (vis-viscos)/SIGMA_F

  call transport_scalar( gam, dGamdxi, gamo, gamoo, iGam, 'gamma', suv, spv, dif )

  gam(1:numCells) = min( gam(1:numCells), 1.0_dp )
  call updateBoundary( gam )

end subroutine


subroutine calcsc_ret
!
! Transported transition onset Reynolds number; source drives it to the free-stream correlation outside the boundary layer.
!
  implicit none
  integer :: inp
  real(dp) :: umag, rot, tscale

  do inp=1,numCells
    umag = sqrt(u(inp)**2+v(inp)**2+w(inp)**2)
    rot = den(inp)**2*umag**2/(500.0_dp*viscos)  ! rho/t, t = 500 mu/(rho U^2)
    tscale = C_THETA*rot*(1.0_dp-fthetat(inp))*vol(inp)
    suv(inp) = tscale*retcorr(inp)
    spv(inp) = tscale
  enddo

  dif = SIGMA_THETA*( viscos + (vis-viscos) )

  call transport_scalar( ret, dRetdxi, reto, retoo, iRet, 'Re_theta_t', suv, spv, dif )

end subroutine


subroutine calcsc_tke
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall
  real(dp) :: genp, genn, gmod, dmod, sigk

  ! Standard SST production with limiter
  do inp=1,numCells
    gen(inp) = abs(vis(inp)-viscos)*magStrain(inp)**2
    gen(inp) = min( gen(inp), 0.9_dp*den(inp)*te(inp)*ed(inp) )
  enddo

  ! Wall-adjacent cells: production from the wall shear stress (as in k_omega_SST, automatic wall treatment)
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      gen(ijp) = abs(tau(iWall))*cmu25*sqrt(te(ijp))/(dnw(iWall)*cappa)
    enddo
  enddo
  pk_sst = gen(1:numCells)

  do inp=1,numCells
    genp = max(gen(inp),zero)
    genn = min(gen(inp),zero)
    gmod = gameff(inp)                                   ! Pk~ = gamma_eff Pk
    dmod = min( max(gameff(inp),0.1_dp), 1.0_dp )        ! Dk~ = min(max(gamma_eff,0.1),1) Dk
    suv(inp) = gmod*genp*vol(inp)
    spv(inp) = ( dmod*BETTAST*ed(inp)*den(inp) )*vol(inp) - genn*vol(inp)/(te(inp)+small)
  enddo

  do inp=1,numCells
    sigk = f1(inp)*SIGMK1 + (1.0_dp-f1(inp))*SIGMK2
    dif(inp) = viscos + (vis(inp)-viscos)*sigk
  enddo
  call boundary_diffusivity( iK )

  call transport_scalar( te, dTEdxi, teo, teoo, iK, 'k', suv, spv, dif )

end subroutine


subroutine calcsc_omega
  implicit none
  integer :: inp, ib, i, iface, ijp, iWall
  real(dp) :: genp, genn, vist, alphasst, bettasst, domega, sigom, wlog, wvis

  do inp=1,numCells
    genp = max(pk_sst(inp),zero)
    genn = min(pk_sst(inp),zero)
    vist = (vis(inp)-viscos)/densit

    alphasst = f1(inp)*ALPHA1 + (1.0_dp-f1(inp))*ALPHA2
    bettasst = f1(inp)*BETAI1 + (1.0_dp-f1(inp))*BETAI2

    suv(inp) = alphasst*genp*vol(inp)/(vist+small)

    domega = max( 2*(1.0_dp-f1(inp))*den(inp)*SIGMOM2/(ed(inp)+small)*sum(dTEdxi(:,inp)*dEDdxi(:,inp)), 0.0_dp )
    suv(inp) = suv(inp) + domega*vol(inp)

    spv(inp) = bettasst*den(inp)*ed(inp)*vol(inp) - alphasst*genn*vol(inp)/(vist*ed(inp)+small)
  enddo

  fixed = .false.
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      wlog = sqrt(te(ijp))/(cmu25*cappa*dnw(iWall))
      wvis = 6.0_dp*(viscos/den(ijp))/(betai1*dnw(iWall)**2)
      fixed(ijp) = .true.
      fixval(ijp) = sqrt(wvis**2 + wlog**2)    ! automatic wall treatment
    enddo
  enddo

  do inp=1,numCells
    sigom = f1(inp)*SIGMOM1 + (1.0_dp-f1(inp))*SIGMOM2
    dif(inp) = viscos + (vis(inp)-viscos)*sigom
  enddo
  call boundary_diffusivity( iOm )

  call transport_scalar( ed, dEDdxi, edo, edoo, iOm, 'omega', suv, spv, dif, fixed=fixed, fixval=fixval )

end subroutine


subroutine boundary_diffusivity( ifi )
!
! Boundary face diffusivity of k or omega: viscos + (vis_b - viscos)*sigma, with sigma of the owner cell.
!
  implicit none
  integer, intent(in) :: ifi
  integer :: ib, i, iface, ijb, ijp
  real(dp) :: sig

  do ib=1,numBoundaries
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      ijb = iBndValueStart(ib) + i
      if ( ifi == iK ) then
        sig = f1(ijp)*SIGMK1 + (1.0_dp-f1(ijp))*SIGMK2
      else
        sig = f1(ijp)*SIGMOM1 + (1.0_dp-f1(ijp))*SIGMOM2
      endif
      dif(ijb) = viscos + (vis(ijb)-viscos)*sig
    enddo
  enddo

end subroutine


subroutine modify_viscosity_inlet_transition_sst
!
! Update effective viscosity at inlet; set inlet gamma = 1 and Re_theta_t from the free-stream correlation.
!
  implicit none
  integer :: i, ib, ijb

  if (.not.allocated(gam)) call allocate_and_initialise

  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        vis(ijb) = viscos + den(ijb)*te(ijb)/(ed(ijb)+small)
        gam(ijb) = 1.0_dp
        ret(ijb) = retheta_freestream( te(ijb), sqrt(u(ijb)**2+v(ijb)**2+w(ijb)**2), 0.0_dp )
      end do
    endif
  enddo

end subroutine

end module
