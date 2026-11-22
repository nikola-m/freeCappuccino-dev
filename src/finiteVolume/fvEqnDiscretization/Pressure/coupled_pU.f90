module coupled_pU
!
!  Purpose:
!
!    Block-coupled pressure-velocity solver. The unknowns (u,v,w,p) of each cell form one
!    4x4 block, momentum and continuity equations are assembled in a single block-CSR matrix
!    and solved simultaneously with FGMRES preconditioned by block aggregation AMG (block_amg).
!
!  Discretization (cell-centred, collocated, Rhie-Chow):
!
!    Momentum   : convection/diffusion/unsteady/sources/BCs from calcuvw_assemble (the same code as in
!                 the segregated solver). Pressure gradient is implicit, Gauss theorem with linear
!                 face interpolation:  sum_f p_f S_f ,  p_f = (1-lambda) p_P + lambda p_N.
!                 Under-relaxation is implicit (diagonal divided by urfU).
!
!    Continuity : sum_f mdot_f = 0, with Rhie-Chow face flux
!
!        mdot_f = rho_f (Uf . S) - rho_f (Vol/ap)_f Df [ (p_N - p_P) - (grad p)_f . d_PN ]
!
!                 where Uf is the linearly interpolated cell velocity (implicit), the compact pressure
!                 difference is implicit and the wide-stencil gradient term is lagged (deferred).
!                 For a linear pressure field the flux reduces to rho (Uf.S).  (Vol/ap) uses the
!                 un-relaxed ap, hence the converged solution does not depend on urfU.
!
!    Boundaries (implicit where possible):
!      inlet     : prescribed mass flux flmass(iface), u from BC;  p_b = p_P (zero gradient)
!      wall      : zero mass flux, p_b = p_P + grad p . d (linear extrapolation, explicit part lagged)
!      symmetry  : zero mass flux, p_b = p_P
!      outlet    : zero gradient u, mass flux scaled to satisfy global continuity (adjustMassFlow), p_b = p_P
!                  Pressure level is fixed by replacing the continuity row of cell pRefCell with p = p_ref.
!      pressure  : fixed p_b; mass flux with Rhie-Chow term between cell centre and boundary face:
!                  mdot_b = rho U_P . S - rho (Vol/ap)_P Df_b [ p_b - p_P - grad p_P . d_Pb ]
!
!  ALL SPEEDS (compressible = T, ideal gas):  rho = p/(R T), temperature lagged in the p-U step, total enthalpy
!  equation solved afterwards (module allspeeds).  The continuity equation is
!
!        (tc0 rho - tc1 rho_o + tc2 rho_oo) Vol/dt + sum_f mdot_f = 0,         rho = C p,  C = 1/(R T)
!
!        mdot_f = rho_up Vn_f ,    Vn_f = Uf.S - (Vol/ap)_f Df [ (p_N - p_P) - (grad p)_f . d ]   (volume flux, Rhie-Chow)
!
!  with the density of the upwind cell. mdot is linearised (Newton) around the current iterate (*):
!
!        mdot_f = rho_up* Vn_f + Vn_f* C_up ( p_up - p_up* ) + beta_f ( rho_c - rho_up )* Vn_f*
!
!  The second term is the pressure-density coupling (a convection-like term of p in the continuity equation, it is what
!  makes the system hyperbolic for M > 1); the last term is a deferred central-density correction
!  with beta_f = cpl_rhoBlend max(0, 1-M_f^2): second order for subsonic flow, pure upwind density when supersonic.
!  Mass fluxes after the solve satisfy exactly the linear system that was solved.
!
!    Boundaries (compressible):
!      inlet    : mdot_b = C_b p_P (U_b.S), U_b and T_b prescribed (implicit in p_P).
!      outlet   : all variables extrapolated (zero gradient), mdot_b = rho_P U_P.S linearised; no flux scaling, no
!                 pressure pinning (supersonic or non-reflecting-like exit).  Use 'pressure' for a subsonic exit.
!      pressure : fixed p_b, rho_b = p_b/(R T_b); Rhie-Chow flux between cell centre and face, as in the incompressible case.
!      wall, symmetry : zero mass flux.
!    Further boundary condition kinds are attached to these patch types by the namelist arrays cplbc_name / cplbc_kind
!    (module implicit_bc_compressible): stagnation, totalpressure, massflow, compinlet, compoutlet, open, riemann.  Their
!    boundary state is a Picard update from the lagged cell state (cplbc_urf) at the start of every outer iteration.
!    The pressure level is fixed by the pressure boundary; only closed steady domains are pinned (cell pRefCell).
!
!  The momentum equations of the compressible form contain the explicit bulk-viscosity term -(2/3) mu div U (velocity.f90).
!  Not supported (stops with a message): periodic boundaries, CN time stepping.  Turbulence is not coupled to the
!  compressible form of the equations (constant viscosity).
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use velocity, only: calcuvw_assemble, updateVelocityAtBoundary, urfU
  use nablap, only: gradp_and_sources, pscheme
  use allspeeds
  use implicit_bc_compressible, only: bckind, bc_setup, bc_prepare_face, bc_apply_face, bc_post_solve_face, &
                                      bc_pressure_implicit, bc_is_riemann, bc_riemann_momentum, BCK_DEFAULT
  use block_amg, only: amg_hierarchy, amg_setup, amg_fgmres, amg_verbose, amg_smoother, amg_cycle, &
                       amg_theta, amg_coarse_size, amg_nu1, amg_nu2, amg_ccf, amg_max_levels

  implicit none

  private

  ! Switches and parameters set through namelist - see read_input.f90
  logical  :: coupled      = .false.   ! Use block-coupled p-U solver instead of SIMPLE/PISO
  integer  :: cpl_maxit    = 40        ! Max. no. of FGMRES iterations per outer iteration
  integer  :: cpl_restart  = 20        ! FGMRES restart length
  real(dp) :: cpl_tolRel   = 0.05_dp   ! Relative residual drop of linear solver per outer iteration
  real(dp) :: cpl_tolAbs   = 1.0e-14_dp
  real(dp) :: cpl_urfP     = 1.0_dp    ! Under-relaxation of pressure update (1 = none)
  integer  :: cpl_rebuild  = 10        ! Rebuild AMG aggregates every cpl_rebuild outer iterations (values are refreshed always)
  real(dp) :: cpl_contScale = 1.0_dp   ! Scaling of continuity rows (affects only the linear solver norm)
  real(dp) :: cpl_rhoBlend = 1.0_dp    ! All speeds: weight of the deferred central density correction (0 = pure upwind density)
  real(dp) :: cpl_pseudoCFL = 0.0_dp   ! All speeds: pseudo-time term (C Vol/dtau)(p-p*) in continuity, dtau = CFL Vol/sum_f(|Vn|+c|S|)/2; 0 = off
  logical  :: cpl_sutherland = .false.  ! All speeds: Sutherland viscosity law, viscos = mu at T = cpl_Tref
  real(dp) :: cpl_Tref = 288.15_dp     ! All speeds: Sutherland reference temperature
  real(dp) :: cpl_Ssuth = 110.4_dp     ! All speeds: Sutherland temperature S [K]
  integer  :: cpl_cflRamp = 0          ! All speeds: the pseudo CFL grows geometrically from 1 to cpl_pseudoCFL during the first cpl_cflRamp iterations

  public :: coupled, calc_coupled_pU
  public :: cpl_maxit, cpl_restart, cpl_tolRel, cpl_tolAbs, cpl_urfP, cpl_rebuild, cpl_contScale
  public :: cpl_rhoBlend, cpl_pseudoCFL, cpl_cflRamp, cpl_sutherland, cpl_Tref, cpl_Ssuth
  ! AMG parameters, re-exported so that read_input only needs this module
  public :: amg_verbose, amg_smoother, amg_cycle, amg_theta, amg_coarse_size, amg_nu1, amg_nu2, amg_ccf, amg_max_levels

  type(amg_hierarchy), save :: amg
  real(dp), allocatable, save :: bval(:,:,:)       ! (4,4,nnz) block matrix values, CSR order of sparse_matrix
  real(dp), allocatable, save :: rhsb(:,:)         ! (4,numCells)
  real(dp), allocatable, save :: solb(:,:)         ! (4,numCells)
  real(dp), allocatable, save :: Ef(:)             ! deferred Rhie-Chow part of inner face fluxes
  real(dp), allocatable, save :: Eb(:)             ! deferred part for 'pressure' boundary faces (indexed by boundary face)
  real(dp), allocatable, save :: Vs(:)             ! all speeds: lagged volume flux through inner faces
  real(dp), allocatable, save :: Edc(:)             ! all speeds: deferred central density correction of the mass flux
  real(dp), allocatable, save :: Vsb(:)            ! all speeds: lagged U_P.S on boundary faces
  real(dp), allocatable, save :: wsum(:)           ! all speeds: sum of face wave speeds (for the pseudo time step)
  logical, save :: first_call = .true.
  logical, save :: has_pressure_bnd = .false.
  logical, save :: has_outlet = .false.
  integer, save :: nouter = 0
  real(dp), save :: pref_value = 0.0_dp


contains


subroutine setup_coupled
  implicit none
  integer :: ib

  if ( numPeriodic > 0 ) then
    write(*,'(a)') 'Fatal error: block-coupled p-U solver does not support periodic boundaries yet.'
    stop
  endif
  if ( cn ) then
    write(*,'(a)') 'Fatal error: block-coupled p-U solver supports only BDF1/2 or steady state (no Crank-Nicolson).'
    stop
  endif
  if ( compressible .and. bdf3 ) then
    write(*,'(a)') 'Fatal error: compressible block-coupled p-U solver supports only BDF1/2 or steady state.'
    stop
  endif
  if ( pscheme /= 'linear' .and. pscheme /= 'central' ) then
    write(*,'(a)') 'Fatal error: block-coupled p-U solver supports only linear pressure interpolation (pscheme).'
    stop
  endif

  do ib=1,numBoundaries
    if ( bctype(ib) == 'pressure' ) has_pressure_bnd = .true.
    if ( bctype(ib) == 'outlet' )   has_outlet = .true.
  enddo

  allocate( bval(4,4,nnz), rhsb(4,numCells), solb(4,numCells) )
  allocate( Ef(numInnerFaces), Eb(numBoundaryFaces) )
  Ef = 0.0_dp
  Eb = 0.0_dp
  pref_value = p(pRefCell)

  if ( compressible ) then
    allocate( Vs(numInnerFaces), Edc(numInnerFaces), Vsb(numBoundaryFaces), wsum(numCells) )
    Vs = 0.0_dp; Edc = 0.0_dp; Vsb = 0.0_dp; wsum = 0.0_dp
    if ( cpl_sutherland ) then; sutherland = .true.; mu_Tref = cpl_Tref; mu_S = cpl_Ssuth; endif
    call as_init
    call bc_setup                      ! compressible boundary conditions of the namelist (implicit_bc_compressible)
    call initial_mass_fluxes
  endif

  first_call = .false.

end subroutine


subroutine scale_outlet_flux
!
!  Same as adjustMassFlow of the segregated solvers - outlet mass fluxes are the extrapolated velocity flux
!  scaled to match the total inlet mass flow flomas - but guarded for the start of a simulation, when the
!  velocity at the outlet may be zero (then the mass flux distribution set in init is kept).
!
  implicit none
  integer :: i, ib, ijp, ijb, iface
  real(dp) :: flowo, fac

  flowo = 0.0_dp
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'outlet' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      flowo = flowo + den(ijp)*( u(ijp)*arx(iface)+v(ijp)*ary(iface)+w(ijp)*arz(iface) )
    enddo
  enddo

  if ( flowo > 1.0e-6_dp*abs(flomas) ) then

    fac = flomas/flowo

    do ib=1,numBoundaries
      if ( bctype(ib) /= 'outlet' ) cycle
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        u(ijb) = u(ijp)*fac
        v(ijb) = v(ijp)*fac
        w(ijb) = w(ijp)*fac
        flmass(iface) = den(ijp)*( u(ijb)*arx(iface)+v(ijb)*ary(iface)+w(ijb)*arz(iface) )
      enddo
    enddo

  endif

end subroutine


subroutine calc_coupled_pU
!
!  One outer (Picard) iteration of the block-coupled solver.
!
  implicit none

  integer :: i, k, ijp, ijn, ib, iface, ijb, inp, kd, k1, k2, its
  real(dp) :: dgu, dgv, dgw, sumoff
  real(dp) :: lam, dene, Dpu, cap, capb, are2, dpbx, sdb
  real(dp) :: sx, sy, sz, dx, dy, dz, gdotd, extrap
  real(dp) :: res0, resf, r4(4), rsum(4)
  real(dp) :: ui, vi, wi
  real(dp) :: capv, Evv, rhou, Cu, rhoc, Mf2, beta, Cc, rb, Cb, vnb, dp_new, cf2, dtau, cflnow
  logical :: rebuild, need_pin

  if ( first_call ) call setup_coupled

  nouter = nouter + 1

  if ( compressible ) then
    call as_time_levels
  elseif ( .not. has_pressure_bnd .and. has_outlet .and. .not. const_mflux ) then
    ! Outlet mass fluxes scaled for global mass conservation (only if there is no fixed-pressure boundary).
    call scale_outlet_flux
  endif

  !
  ! 1. Assemble momentum equations without the pressure gradient (also updates boundary values and gradients)
  !
  call calcuvw_assemble( .true. )

  bval = 0.0_dp
  rhsb = 0.0_dp

  do inp=1,numCells

    kd = diag(inp)

    ! Sum of off-diagonal coefficients (the diagonal is excluded explicitly, no cancellation involved).
    sumoff = sum( a( ia(inp) : kd-1 ) ) + sum( a( kd+1 : ia(inp+1)-1 ) )

    ! Un-relaxed diagonal, its reciprocal is used by Rhie-Chow interpolation as in the segregated solver.
    dgu = spu(inp) - sumoff
    dgv = spv(inp) - sumoff
    dgw = sp(inp)  - sumoff
    apu(inp) = 1.0_dp/(dgu+small)
    apv(inp) = 1.0_dp/(dgv+small)
    apw(inp) = 1.0_dp/(dgw+small)

    ! Implicit under-relaxation
    dgu = dgu/urfU(1)
    dgv = dgv/urfU(2)
    dgw = dgw/urfU(3)

    do k = ia(inp), ia(inp+1)-1
      if ( k == kd ) then
        bval(1,1,k) = dgu
        bval(2,2,k) = dgv
        bval(3,3,k) = dgw
      else
        bval(1,1,k) = a(k)
        bval(2,2,k) = a(k)
        bval(3,3,k) = a(k)
      endif
    enddo

    rhsb(1,inp) = su(inp) + (1.0_dp-urfU(1))*dgu*u(inp)
    rhsb(2,inp) = sv(inp) + (1.0_dp-urfU(2))*dgv*v(inp)
    rhsb(3,inp) = sw(inp) + (1.0_dp-urfU(3))*dgw*w(inp)

  enddo

  !
  ! 2. Inner faces: implicit pressure gradient in momentum, Rhie-Chow mass flux in continuity
  !
  if ( compressible ) wsum = 0.0_dp
  do i=1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)

    sx = arx(i); sy = ary(i); sz = arz(i)

    k1 = icell_jcell_csr_index(i)    ! (P,N)
    k2 = jcell_icell_csr_index(i)    ! (N,P)

    ! --- Momentum: + sum_f p_f S_f  on owner,  - sum_f p_f S_f on neighbour
    bval(1,4,diag(ijp)) = bval(1,4,diag(ijp)) + (1.0_dp-lam)*sx
    bval(2,4,diag(ijp)) = bval(2,4,diag(ijp)) + (1.0_dp-lam)*sy
    bval(3,4,diag(ijp)) = bval(3,4,diag(ijp)) + (1.0_dp-lam)*sz

    bval(1,4,k1) = bval(1,4,k1) + lam*sx
    bval(2,4,k1) = bval(2,4,k1) + lam*sy
    bval(3,4,k1) = bval(3,4,k1) + lam*sz

    bval(1,4,diag(ijn)) = bval(1,4,diag(ijn)) - lam*sx
    bval(2,4,diag(ijn)) = bval(2,4,diag(ijn)) - lam*sy
    bval(3,4,diag(ijn)) = bval(3,4,diag(ijn)) - lam*sz

    bval(1,4,k2) = bval(1,4,k2) - (1.0_dp-lam)*sx
    bval(2,4,k2) = bval(2,4,k2) - (1.0_dp-lam)*sy
    bval(3,4,k2) = bval(3,4,k2) - (1.0_dp-lam)*sz

    ! --- Continuity
    Dpu  = lam*vol(ijn)*apu(ijn) + (1.0_dp-lam)*vol(ijp)*apu(ijp)

    dx = xc(ijn)-xc(ijp)
    dy = yc(ijn)-yc(ijp)
    dz = zc(ijn)-zc(ijp)

    ! Interpolated pressure gradient projected on the cell-centre distance vector
    gdotd = (lam*dPdxi(1,ijn)+(1.0_dp-lam)*dPdxi(1,ijp))*dx + &
            (lam*dPdxi(2,ijn)+(1.0_dp-lam)*dPdxi(2,ijp))*dy + &
            (lam*dPdxi(3,ijn)+(1.0_dp-lam)*dPdxi(3,ijp))*dz

    if ( .not. compressible ) then
      dene = den(ijp)*(1.0_dp-lam) + den(ijn)*lam
      cap  = -dene*Dpu*Df(i)                 ! < 0
    else
      !
      ! All speeds: upwind density, Newton linearisation of rho_up * Vn, deferred central density correction
      !
      capv = -Dpu*Df(i)                      ! volume flux coefficient of (p_N - p_P), < 0
      Evv  = -capv*gdotd                     ! deferred part of the volume flux

      ui = u(ijp) + lam*(u(ijn)-u(ijp))
      vi = v(ijp) + lam*(v(ijn)-v(ijp))
      wi = w(ijp) + lam*(w(ijn)-w(ijp))

      Vs(i) = ui*sx + vi*sy + wi*sz + capv*( p(ijn)-p(ijp) ) + Evv

      if ( Vs(i) >= 0.0_dp ) then
        rhou = den(ijp);  Cu = 1.0_dp/(Rair*T(ijp)+small)
      else
        rhou = den(ijn);  Cu = 1.0_dp/(Rair*T(ijn)+small)
      endif

      rhoc = den(ijp)*(1.0_dp-lam) + den(ijn)*lam
      Mf2  = Vs(i)**2/( sx*sx+sy*sy+sz*sz )/speed_of_sound_sq( T(ijp)*(1.0_dp-lam)+T(ijn)*lam )
      beta = cpl_rhoBlend*max( 0.0_dp, 1.0_dp-Mf2 )
      Edc(i) = beta*(rhoc-rhou)*Vs(i)

      ! face wave speed  (|Vn| + c |S|)/2 for the pseudo time step
      cf2 = 0.5_dp*( abs(Vs(i)) + sqrt( speed_of_sound_sq( T(ijp)*(1.0_dp-lam)+T(ijn)*lam )*( sx*sx+sy*sy+sz*sz ) ) )
      wsum(ijp) = wsum(ijp) + cf2
      wsum(ijn) = wsum(ijn) + cf2

      dene = rhou
      cap  = rhou*capv                       ! < 0
    endif

    Ef(i) = -cap*gdotd                       ! explicit (deferred) part of the face flux

    ! owner row  ( + mdot_f )
    bval(4,1,diag(ijp)) = bval(4,1,diag(ijp)) + dene*(1.0_dp-lam)*sx
    bval(4,2,diag(ijp)) = bval(4,2,diag(ijp)) + dene*(1.0_dp-lam)*sy
    bval(4,3,diag(ijp)) = bval(4,3,diag(ijp)) + dene*(1.0_dp-lam)*sz
    bval(4,1,k1) = bval(4,1,k1) + dene*lam*sx
    bval(4,2,k1) = bval(4,2,k1) + dene*lam*sy
    bval(4,3,k1) = bval(4,3,k1) + dene*lam*sz
    bval(4,4,diag(ijp)) = bval(4,4,diag(ijp)) - cap
    bval(4,4,k1)        = bval(4,4,k1)        + cap
    rhsb(4,ijp) = rhsb(4,ijp) - Ef(i)

    ! neighbour row ( - mdot_f )
    bval(4,1,diag(ijn)) = bval(4,1,diag(ijn)) - dene*lam*sx
    bval(4,2,diag(ijn)) = bval(4,2,diag(ijn)) - dene*lam*sy
    bval(4,3,diag(ijn)) = bval(4,3,diag(ijn)) - dene*lam*sz
    bval(4,1,k2) = bval(4,1,k2) - dene*(1.0_dp-lam)*sx
    bval(4,2,k2) = bval(4,2,k2) - dene*(1.0_dp-lam)*sy
    bval(4,3,k2) = bval(4,3,k2) - dene*(1.0_dp-lam)*sz
    bval(4,4,diag(ijn)) = bval(4,4,diag(ijn)) - cap
    bval(4,4,k2)        = bval(4,4,k2)        + cap
    rhsb(4,ijn) = rhsb(4,ijn) + Ef(i)

    if ( compressible ) then
      ! pressure-density coupling  Vs C_up (p_up - p_up*)  and the deferred density correction
      if ( Vs(i) >= 0.0_dp ) then
        bval(4,4,diag(ijp)) = bval(4,4,diag(ijp)) + Vs(i)*Cu
        bval(4,4,k2)        = bval(4,4,k2)        - Vs(i)*Cu
      else
        bval(4,4,k1)        = bval(4,4,k1)        + Vs(i)*Cu
        bval(4,4,diag(ijn)) = bval(4,4,diag(ijn)) - Vs(i)*Cu
      endif
      rhsb(4,ijp) = rhsb(4,ijp) + Vs(i)*rhou - Edc(i)
      rhsb(4,ijn) = rhsb(4,ijn) - Vs(i)*rhou + Edc(i)
    endif

  enddo

  !
  ! 3. Boundary faces
  !
  Eb = 0.0_dp

  do ib=1,numBoundaries

    do i=1,nfaces(ib)

      iface = startFace(ib) + i
      ijp   = owner(iface)
      ijb   = iBndValueStart(ib) + i
      kd    = diag(ijp)
      sx = arx(iface); sy = ary(iface); sz = arz(iface)

      dx = xf(iface)-xc(ijp)
      dy = yf(iface)-yc(ijp)
      dz = zf(iface)-zc(ijp)

      ! Compressible boundary conditions of the namelist: boundary state from the current cell state
      if ( compressible .and. bckind(ib) /= BCK_DEFAULT ) call bc_prepare_face( ib, i, iface, ijp, ijb )

      ! --- Momentum: pressure force on boundary face
      if ( compressible .and. bckind(ib) /= BCK_DEFAULT .and. bc_is_riemann( iface ) ) then

        ! Subsonic characteristic boundary: pressure force implicit in U_P and p_P
        call bc_riemann_momentum( iface, ijp, ijb, kd, bval, rhsb )

      elseif ( bctype(ib) == 'pressure' .and. .not. ( compressible .and. bc_pressure_implicit( iface ) ) ) then

        ! Fixed boundary pressure - known, goes to the right hand side
        rhsb(1,ijp) = rhsb(1,ijp) - p(ijb)*sx
        rhsb(2,ijp) = rhsb(2,ijp) - p(ijb)*sy
        rhsb(3,ijp) = rhsb(3,ijp) - p(ijb)*sz

      else

        ! Zero gradient, p_b = p_P, (walls: + linear extrapolation with lagged gradient)
        bval(1,4,kd) = bval(1,4,kd) + sx
        bval(2,4,kd) = bval(2,4,kd) + sy
        bval(3,4,kd) = bval(3,4,kd) + sz

        if ( bctype(ib) == 'wall' ) then
          extrap = dPdxi(1,ijp)*dx + dPdxi(2,ijp)*dy + dPdxi(3,ijp)*dz
          rhsb(1,ijp) = rhsb(1,ijp) - extrap*sx
          rhsb(2,ijp) = rhsb(2,ijp) - extrap*sy
          rhsb(3,ijp) = rhsb(3,ijp) - extrap*sz
        endif

      endif

      if ( compressible .and. bctype(ib) /= 'empty' ) then
        wsum(ijp) = wsum(ijp) + 0.5_dp*( abs( u(ijp)*sx + v(ijp)*sy + w(ijp)*sz ) &
                  + sqrt( speed_of_sound_sq( T(ijp) )*( sx*sx+sy*sy+sz*sz ) ) )
      endif

      ! --- Continuity: mass flux through boundary face
      if ( compressible .and. bckind(ib) /= BCK_DEFAULT ) then

        ! Boundary condition types of implicit_bc_compressible
        call bc_apply_face( ib, i, iface, ijp, ijb, kd, bval, rhsb, Eb, Vsb )

      elseif ( compressible .and. bctype(ib) == 'inlet' ) then

        ! mdot_b = rho_b U_b.S, rho_b = p_P/(R T_b): linear in p_P
        Cb = 1.0_dp/(Rair*T(ijb)+small)
        bval(4,4,kd) = bval(4,4,kd) + Cb*( u(ijb)*sx + v(ijb)*sy + w(ijb)*sz )

      elseif ( compressible .and. bctype(ib) == 'outlet' ) then

        ! Zero gradient: mdot_b = rho_P U_P.S, linearised around the current iterate
        Cc  = 1.0_dp/(Rair*T(ijp)+small)
        vnb = u(ijp)*sx + v(ijp)*sy + w(ijp)*sz
        Vsb(iface-numInnerFaces) = vnb
        bval(4,1,kd) = bval(4,1,kd) + den(ijp)*sx
        bval(4,2,kd) = bval(4,2,kd) + den(ijp)*sy
        bval(4,3,kd) = bval(4,3,kd) + den(ijp)*sz
        bval(4,4,kd) = bval(4,4,kd) + Cc*vnb
        rhsb(4,ijp)  = rhsb(4,ijp)  + vnb*den(ijp)

      elseif ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' ) then

        ! Prescribed (inlet) or scaled (outlet) mass flux. Known.
        rhsb(4,ijp) = rhsb(4,ijp) - flmass(iface)

      elseif ( bctype(ib) == 'pressure' ) then

        if ( compressible ) then
          rb = den(ijb)                    ! p_b/(R T_b)
        else
          rb = den(ijp)
        endif

        are2 = sx*sx+sy*sy+sz*sz
        sdb  = sx*dx+sy*dy+sz*dz
        capb = -rb*vol(ijp)*apu(ijp)*are2/sdb     ! < 0

        extrap = dPdxi(1,ijp)*dx + dPdxi(2,ijp)*dy + dPdxi(3,ijp)*dz
        Eb(iface-numInnerFaces) = capb*( p(ijb) - extrap )          ! known part of the flux, no p_P, no u_P

        bval(4,1,kd) = bval(4,1,kd) + rb*sx
        bval(4,2,kd) = bval(4,2,kd) + rb*sy
        bval(4,3,kd) = bval(4,3,kd) + rb*sz
        bval(4,4,kd) = bval(4,4,kd) - capb
        rhsb(4,ijp)  = rhsb(4,ijp) - Eb(iface-numInnerFaces)

      endif
      ! wall, symmetry, empty, etc.: zero mass flux.

    enddo

  enddo

  !
  ! 3b. All speeds, unsteady term  (tc0 rho - tc1 rho_o + tc2 rho_oo) Vol/dt,  rho = C p
  !
  if ( compressible .and. ltransient ) then
    do inp=1,numCells
      Cc = 1.0_dp/(Rair*T(inp)+small)
      bval(4,4,diag(inp)) = bval(4,4,diag(inp)) + tc0*Cc*vol(inp)/timestep
      rhsb(4,inp) = rhsb(4,inp) + ( tc1*den_o(inp) - tc2*den_oo(inp) )*vol(inp)/timestep
    enddo
  endif

  !
  ! 3c. All speeds, pseudo-time stabilisation  (C Vol/dtau)(p - p*)  - vanishes at convergence
  !
  if ( compressible .and. cpl_pseudoCFL > 0.0_dp ) then
    cflnow = cpl_pseudoCFL
    if ( cpl_cflRamp > 0 .and. nouter < cpl_cflRamp .and. cpl_pseudoCFL > 1.0_dp ) &
      cflnow = cpl_pseudoCFL**( real(nouter-1,dp)/real(cpl_cflRamp,dp) )
    do inp=1,numCells
      dtau = cflnow*vol(inp)/max( wsum(inp), small )
      Cc = 1.0_dp/(Rair*T(inp)+small)
      bval(4,4,diag(inp)) = bval(4,4,diag(inp)) + Cc*vol(inp)/dtau
      rhsb(4,inp) = rhsb(4,inp) + Cc*vol(inp)/dtau*p(inp)
    enddo
  endif

  !
  ! 4. Pressure level for problems without a fixed-pressure boundary
  !
  need_pin = .not. has_pressure_bnd
  if ( compressible ) need_pin = need_pin .and. .not. has_outlet .and. .not. ltransient
  if ( need_pin ) then
    do k = ia(pRefCell), ia(pRefCell+1)-1
      bval(4,:,k) = 0.0_dp
    enddo
    bval(4,4,diag(pRefCell)) = 1.0_dp
    rhsb(4,pRefCell) = pref_value
  endif

  ! Optional scaling of the continuity rows
  if ( cpl_contScale /= 1.0_dp ) then
    do inp=1,numCells
      if ( need_pin .and. inp == pRefCell ) cycle
      do k = ia(inp), ia(inp+1)-1
        bval(4,:,k) = bval(4,:,k)*cpl_contScale
      enddo
      rhsb(4,inp) = rhsb(4,inp)*cpl_contScale
    enddo
  endif

  !
  ! 5. Residuals of the assembled (nonlinear) system for the current iterate: L1 norms, like in segregated solvers
  !
  rsum = 0.0_dp
  do inp=1,numCells
    r4 = rhsb(:,inp)
    do k = ia(inp), ia(inp+1)-1
      ijn = ja(k)
      r4(1) = r4(1) - ( bval(1,1,k)*u(ijn) + bval(1,4,k)*p(ijn) )
      r4(2) = r4(2) - ( bval(2,2,k)*v(ijn) + bval(2,4,k)*p(ijn) )
      r4(3) = r4(3) - ( bval(3,3,k)*w(ijn) + bval(3,4,k)*p(ijn) )
      r4(4) = r4(4) - ( bval(4,1,k)*u(ijn) + bval(4,2,k)*v(ijn) + bval(4,3,k)*w(ijn) + bval(4,4,k)*p(ijn) )
    enddo
    rsum = rsum + abs(r4)
  enddo
  resor(1) = rsum(1)
  resor(2) = rsum(2)
  resor(3) = rsum(3)
  resor(4) = rsum(4)/max(cpl_contScale,small)

  !
  ! 6. Solve the coupled system with FGMRES + block AMG
  !
  do inp=1,numCells
    solb(1,inp) = u(inp)
    solb(2,inp) = v(inp)
    solb(3,inp) = w(inp)
    solb(4,inp) = p(inp)
  enddo

  rebuild = ( nouter == 1 .or. mod(nouter-1, max(cpl_rebuild,1)) == 0 )
  call amg_setup( amg, numCells, ia, ja, bval, rebuild )
  call amg_fgmres( amg, solb, rhsb, cpl_tolRel, cpl_tolAbs, cpl_maxit, cpl_restart, res0, resf, its )

  write(*,'(2x,a,i0,a,es9.2,a,es9.2)') 'Coupled p-U: FGMRES(AMG) iterations = ', its, &
        ', initial residual = ', res0, ', final residual = ', resf
  write(6,'(2x,a,i0,a,es9.2,a,es9.2)') 'Coupled p-U: FGMRES(AMG) iterations = ', its, &
        ', initial residual = ', res0, ', final residual = ', resf
  write(6,'(2x,a,4(es10.3,1x))') 'Coupled p-U: L1 residuals u,v,w,continuity: ', resor(1:4)

  !
  ! 7. Update fields
  !
  do inp=1,numCells
    u(inp) = solb(1,inp)
    v(inp) = solb(2,inp)
    w(inp) = solb(3,inp)
    p(inp) = p(inp) + cpl_urfP*( solb(4,inp) - p(inp) )
  enddo

  call updateVelocityAtBoundary

  !
  ! 8. Mass fluxes consistent with the solved continuity equation
  !
  do i=1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)

    Dpu  = lam*vol(ijn)*apu(ijn) + (1.0_dp-lam)*vol(ijp)*apu(ijp)

    ui = u(ijp) + lam*(u(ijn)-u(ijp))
    vi = v(ijp) + lam*(v(ijn)-v(ijp))
    wi = w(ijp) + lam*(w(ijn)-w(ijp))

    if ( .not. compressible ) then

      dene = den(ijp)*(1.0_dp-lam) + den(ijn)*lam
      cap  = -dene*Dpu*Df(i)

      flmass(i) = dene*( ui*arx(i) + vi*ary(i) + wi*arz(i) ) + cap*( p(ijn)-p(ijp) ) + Ef(i)

    else

      ! The mass flux of the linearised continuity equation that was solved (den is still rho* here)
      if ( Vs(i) >= 0.0_dp ) then
        rhou = den(ijp);  Cu = 1.0_dp/(Rair*T(ijp)+small);  dp_new = p(ijp)
      else
        rhou = den(ijn);  Cu = 1.0_dp/(Rair*T(ijn)+small);  dp_new = p(ijn)
      endif
      cap  = -rhou*Dpu*Df(i)

      flmass(i) = rhou*( ui*arx(i) + vi*ary(i) + wi*arz(i) ) + cap*( p(ijn)-p(ijp) ) + Ef(i) &
                + Vs(i)*( Cu*dp_new - rhou ) + Edc(i)

    endif

  enddo

  do ib=1,numBoundaries

    if ( compressible .and. bckind(ib) /= BCK_DEFAULT ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp   = owner(iface)
        ijb   = iBndValueStart(ib) + i
        call bc_post_solve_face( ib, i, iface, ijp, ijb, Eb, Vsb )
      enddo

    elseif ( bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp   = owner(iface)
        ijb   = iBndValueStart(ib) + i
        sx = arx(iface); sy = ary(iface); sz = arz(iface)
        dx = xf(iface)-xc(ijp); dy = yf(iface)-yc(ijp); dz = zf(iface)-zc(ijp)
        if ( compressible ) then
          rb = den(ijb)
        else
          rb = den(ijp)
        endif
        capb = -rb*vol(ijp)*apu(ijp)*(sx*sx+sy*sy+sz*sz)/(sx*dx+sy*dy+sz*dz)
        flmass(iface) = rb*( u(ijp)*sx + v(ijp)*sy + w(ijp)*sz ) - capb*p(ijp) + Eb(iface-numInnerFaces)

        ! Velocity at the boundary face consistent with the flux (needed by momentum convection/diffusion terms):
        ! U_b = U_P - S (Vol/ap)/(S.d) [ p_b - p_P - grad p . d ]
        dpbx = ( p(ijb) - p(ijp) - ( dPdxi(1,ijp)*dx + dPdxi(2,ijp)*dy + dPdxi(3,ijp)*dz ) ) &
               * vol(ijp)*apu(ijp)/(sx*dx+sy*dy+sz*dz)
        u(ijb) = u(ijp) - sx*dpbx
        v(ijb) = v(ijp) - sy*dpbx
        w(ijb) = w(ijp) - sz*dpbx
      enddo

    elseif ( compressible .and. bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp   = owner(iface)
        ijb   = iBndValueStart(ib) + i
        flmass(iface) = p(ijp)/(Rair*T(ijb)+small)*( u(ijb)*arx(iface) + v(ijb)*ary(iface) + w(ijb)*arz(iface) )
      enddo

    elseif ( compressible .and. bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp   = owner(iface)
        ijb   = iBndValueStart(ib) + i
        u(ijb) = u(ijp); v(ijb) = v(ijp); w(ijb) = w(ijp)
        Cc = 1.0_dp/(Rair*T(ijp)+small)
        flmass(iface) = den(ijp)*( u(ijp)*arx(iface) + v(ijp)*ary(iface) + w(ijp)*arz(iface) ) &
                      + Vsb(iface-numInnerFaces)*( Cc*p(ijp) - den(ijp) )
      enddo

    endif

  enddo

  ! Gradient of the new pressure field, for use by the other equations and the next iteration
  call gradp_and_sources( p )

  ! All speeds: density from the new pressure (temperature is lagged until the energy equation is solved)
  if ( compressible ) call update_density

end subroutine


subroutine initial_mass_fluxes
!
!  All speeds: make the initial mass fluxes and boundary values consistent with rho = p/(R T).
!
  implicit none
  integer :: i, ib, ijp, ijn, ijb, iface
  real(dp) :: lam, ui, vi, wi

  do i=1,numInnerFaces
    ijp = owner(i); ijn = neighbour(i); lam = facint(i)
    ui = u(ijp) + lam*(u(ijn)-u(ijp))
    vi = v(ijp) + lam*(v(ijn)-v(ijp))
    wi = w(ijp) + lam*(w(ijn)-w(ijp))
    flmass(i) = ( den(ijp)*(1.0_dp-lam)+den(ijn)*lam )*( ui*arx(i) + vi*ary(i) + wi*arz(i) )
  enddo

  do ib=1,numBoundaries
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      ijb = iBndValueStart(ib) + i
      if ( bctype(ib) == 'outlet' ) then
        u(ijb) = u(ijp); v(ijb) = v(ijp); w(ijb) = w(ijp)
      endif
      if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' .or. bctype(ib) == 'pressure' ) then
        flmass(iface) = den(ijb)*( u(ijb)*arx(iface) + v(ijb)*ary(iface) + w(ijb)*arz(iface) )
      else
        flmass(iface) = 0.0_dp
      endif
    enddo
  enddo

end subroutine

end module coupled_pU
