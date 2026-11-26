module simple_allspeeds
!
!  Purpose:
!
!    Alternative implementation of the SIMPLE pressure-velocity coupling for ALL SPEEDS (incompressible to supersonic flow),
!    pressure-based, collocated, in the spirit of Przulj & Basara. It sits next to the existing SIMPLE/PISO/coupled
!    solvers (switch asimple = T in input.nml) and does not change them: it is meant for testing and comparison.
!
!  Discussion:
!
!    One iteration of the algorithm (calc_simple_allspeeds):
!
!      0. time levels (BDF1/BDF2 coefficients tc0,tc1,tc2 of the energy module), generalized wall functions -> visw
!      1. momentum predictor (the existing assembly calcuvw: convection, stress-tensor diffusion, BDF, pressure gradient,
!         bulk viscosity term, under-relaxation urfU)
!      2. pressure-correction equation, assembled here, solved by the AMG of module asimple_amg (or any csrsolve solver)
!      3. correction of mass fluxes, velocity, pressure and density
!      4. the energy equation is solved by main (calc_energy of the energy module), which closes the loop through rho(p,T)
!
!    Pressure-correction equation. The mass flux of a face is  m_f = rho_f V_f , where V_f is the Rhie-Chow volume flux
!
!         V_f = Uf.S - (Vol/ap)_f Df [ (p_N - p_P) - (grad p)_f . d_PN ]
!
!    and the face density is a Mach number blend of the central and the upwind value (beta = asimple_beta max(0,1-M_f^2)):
!
!         rho_f = beta ( (1-lam) rho_P + lam rho_N ) + (1-beta) rho_up
!
!    Perturbations: m' = rho_f V' + rho_f' V*,  with  V' = -(Vol/ap)_f Df (p'_N - p'_P),  rho' = p'/(R T) (T lagged), hence
!
!         m'_f = cap (p'_N - p'_P) + cvP p'_P + cvN p'_N ,   cap = -rho_f (Vol/ap)_f Df ,
!         cvP = V* [ beta (1-lam) + (1-beta) H(V*) ] /(R T_P) ,   cvN = V* [ beta lam + (1-beta)(1-H(V*)) ] /(R T_N)
!
!    (H is the Heaviside function). Continuity of the cell, with the BDF density term:
!
!         (tc0 rho - tc1 rho_o + tc2 rho_oo) Vol/dt + sum_f m_f = 0
!         ->  tc0 Vol/(R T dt) p'_P + sum_f m'_f = -( tc0 rho* - tc1 rho_o + tc2 rho_oo ) Vol/dt - sum_f m*_f
!
!    The first term is the 'compressible' (transient) part of Przulj & Basara, the cv terms are the convective density-pressure
!    coupling: elliptic for M<1 (cap dominates), hyperbolic (upwind) for M>1 (cap is small, cv upwinded). With compressible = F
!    all cv terms are zero and the algorithm reduces to the classical incompressible SIMPLE of freeCappuccino (a useful check).
!
!    Boundaries:  inlet - fixed mass flux rho_b U_b.S (explicit);  outlet - zero gradient, m_b = rho_P U_P.S with p' coupling
!    cvb = (U_P.S)+/(R T_P) on the diagonal (compressible) or global mass flow scaling (incompressible);  pressure - the same
!    Rhie-Chow boundary flux and correction as in calcp_simple;  wall, symmetry - no flux.  Periodic boundaries are not supported.
!
!    Correction:  u = u* + sum_f(-p' S)/ap ,  p = p* + urfP p' ,  rho = p/(R T)  (T is updated by the energy equation afterwards).
!
!  Reference:
!
!    Przulj V., Basara B., Bounded convection schemes for unstructured grids, AIAA 2001-2593; Basara B., Przulj V. et al.
!    Ferziger J.H., Peric M., Computational Methods for Fluid Dynamics.
!    Demirdzic I., Lilek Z., Peric M., A collocated finite volume method for predicting flows at all speeds, IJNMF 16 (1993).
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
  use gradients
  use linear_solvers
  use faceflux_mass
  use nablap, only: gradp_and_sources
  use velocity, only: calcuvw, calcU, updateVelocityAtBoundary
  use allspeeds, only: as_init, as_time_levels, as_ready, update_density, den_o, den_oo, tc0, tc1, tc2, speed_of_sound_sq
  use asimple_amg, only: asimple_amg_solve, pamg_theta, pamg_coarse_size, pamg_nu1, pamg_nu2, pamg_max_levels, &
                         pamg_restart, pamg_verbose
  use asimple_wallfun, only: asimple_wallFunc, update_wall_viscosity
  use implicit_bc_compressible, only: bckind, bcmode, bc_setup, bc_prepare_face, mdot_face, edir, U_set, &
                                      BCK_DEFAULT, BCK_COMP_INLET, BCK_MASSFLOW, BCK_TOTAL_PRESSURE, BCK_STAGNATION, &
                                      BCK_RIEMANN, BCK_OPEN, MODE_EXTRAP

  implicit none

  private

  ! Switches and parameters set through namelist - see read_input.f90
  logical  :: asimple = .false.                 ! Use all-speeds SIMPLE instead of SIMPLE/PISO/coupled, True/False.
  real(dp) :: asimple_urfP = 0.3_dp             ! Under-relaxation of the pressure update.
  real(dp) :: asimple_beta = 1.0_dp             ! Weight of the central part of the face density (0 = pure upwind density).
  character(len=30) :: asimple_lSolverP = 'amg' ! Linear solver for p': 'amg' (own AMG + FGMRES) or any csrsolve solver name.
  integer  :: asimple_maxitP = 40               ! Max. no. of iterations of the p' solver.
  real(dp) :: asimple_tolRelP = 0.05_dp         ! Relative residual drop in the p' solver.
  real(dp) :: asimple_tolAbsP = 1.0e-14_dp      ! Absolute residual level of the p' solver.
  real(dp) :: asimple_pseudoCFL = 0.0_dp        ! Steady compressible runs: pseudo-time CFL of the term (Vol/(R T dtau)) p' (0 = off).
  integer  :: asimple_cflRamp = 0               ! The pseudo CFL grows geometrically from 1 to asimple_pseudoCFL in this many iterations.

  public :: asimple, asimple_urfP, asimple_beta, asimple_lSolverP, asimple_maxitP, asimple_tolRelP, asimple_tolAbsP
  public :: asimple_pseudoCFL, asimple_cflRamp
  public :: calc_simple_allspeeds
  ! Re-exported so that read_input only needs this module
  public :: pamg_theta, pamg_coarse_size, pamg_nu1, pamg_nu2, pamg_max_levels, pamg_restart, pamg_verbose, asimple_wallFunc

  real(dp), allocatable, save :: capf(:)        ! cap  of inner faces
  real(dp), allocatable, save :: cvP(:), cvN(:) ! convective density-pressure coefficients of inner faces
  real(dp), allocatable, save :: cvB(:)         ! same for boundary faces (outlet), indexed by iface-numInnerFaces
  real(dp), allocatable, save :: capB(:)        ! cap of pressure-type boundary faces with a compressible BC kind (< 0)
  logical, save :: bc_kinds = .false.           ! compressible BC kinds (cplbc_* namelist arrays) are present
  real(dp), allocatable, save :: wsum(:)        ! sum of face wave speeds (|V_f| + c |S_f|) of a cell, for the pseudo time step
  integer, save :: nouter = 0                   ! Outer iteration counter (pseudo CFL ramp)
  logical, save :: first_call = .true.
  integer, save :: nits_amg = 0

contains


subroutine calc_simple_allspeeds
!
!  One outer iteration of the all-speeds SIMPLE algorithm. Called from main instead of calcuvw + calcp_simple.
!
  implicit none

  if ( first_call ) call setup_simple_allspeeds

  ! Time levels, BDF coefficients, wall functions
  if ( compressible ) call as_time_levels
  call update_wall_viscosity

  ! Boundary state (p_b, T_b, rho_b, face mode) of the compressible BC kinds from the current cell state
  if ( bc_kinds ) call prepare_bc_kinds

  ! Momentum predictor
  if ( calcU ) call calcuvw

  ! Pressure correction and corrections of the flow field
  call calcp_allspeeds

end subroutine


subroutine setup_simple_allspeeds
  implicit none

  if ( numPeriodic > 0 ) then
    write(*,'(a)') 'Fatal error: all-speeds SIMPLE does not support periodic boundaries.'
    stop
  endif
  if ( cn ) then
    write(*,'(a)') 'Fatal error: all-speeds SIMPLE supports BDF1/BDF2 or steady state, not Crank-Nicolson.'
    stop
  endif
  if ( compressible .and. bdf3 ) then
    write(*,'(a)') 'Fatal error: compressible all-speeds SIMPLE supports BDF1/BDF2 or steady state.'
    stop
  endif

  allocate( capf(numInnerFaces), cvP(numInnerFaces), cvN(numInnerFaces) )
  allocate( cvB(max(numBoundaryFaces,1)), capB(max(numBoundaryFaces,1)), wsum(numCells) )
  capf = 0.0_dp;  cvP = 0.0_dp;  cvN = 0.0_dp;  cvB = 0.0_dp;  capB = 0.0_dp

  ! Compressible boundary conditions of the namelist (cplbc_name/cplbc_kind...), the same kinds as in the coupled solver
  if ( compressible ) then
    call bc_setup
    bc_kinds = any( bckind(1:numBoundaries) /= BCK_DEFAULT )
  endif

  if ( compressible .and. .not. as_ready ) call as_init

  ! BDF2 in the first time step of a new run: time_shift has copied the (never set) 'uo' into 'uoo', so the first step
  ! would use zero velocities at level n-1. Make the first step BDF1-like, as the energy module does.
  if ( ltransient .and. bdf2 .and. .not. lread ) then
    uoo(1:numCells) = uo(1:numCells)
    voo(1:numCells) = vo(1:numCells)
    woo(1:numCells) = wo(1:numCells)
  endif

  first_call = .false.

  write(*,'(a)') ' '
  write(*,'(a)') ' All-speeds SIMPLE (Przulj-Basara type) pressure-velocity coupling is active.'
  write(*,'(2a)') '   p-correction solver : ', trim(asimple_lSolverP)
  write(*,'(2a)') '   wall function       : ', trim(asimple_wallFunc)
  write(*,'(a,l1)') '   compressible        : ', compressible
  write(*,'(a)') ' '

end subroutine


subroutine calcp_allspeeds
!
!  Assemble and solve the pressure-correction equation, correct fluxes, velocity, pressure and density.
!
  implicit none

  integer :: i, k, ib, iface, inp, ijp, ijn, ijb
  real(dp) :: ppref, cap, can, fmcor
  real(dp) :: rhoc, rhoup, rhof, Vf, are, Tf, c2, Mf, beta, CP, CN, hP, tr, Vb, res0, res
  integer :: its, npc
  real(dp) :: cfl

  ! Clear matrix coefficient and default rhs vector arrays
  a  = 0.0_dp
  su = 0.0_dp
  wsum = 0.0_dp

  !
  ! > Inner faces: Rhie-Chow mass flux, matrix coefficients
  !
  do i = 1, numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    ! Central-density Rhie-Chow flux and coefficient cap = -rho_c (Vol/ap)_f Df
    call facefluxmass2( i, ijp, ijn, arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i) )

    if ( compressible ) then

      ! Volume flux, face density: Mach number blend of central and upwind values
      rhoc = den(ijp)*(1.0_dp-facint(i)) + den(ijn)*facint(i)
      Vf   = flmass(i)/rhoc
      rhoup = den(ijp);  hP = 1.0_dp
      if ( Vf < 0.0_dp ) then
        rhoup = den(ijn);  hP = 0.0_dp
      endif

      are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 )
      Tf  = T(ijp)*(1.0_dp-facint(i)) + T(ijn)*facint(i)
      c2  = speed_of_sound_sq( Tf )
      Mf  = abs(Vf)/( are*sqrt(c2) + small )
      beta = asimple_beta*max( 0.0_dp, 1.0_dp - Mf**2 )

      rhof = beta*rhoc + (1.0_dp-beta)*rhoup

      wsum(ijp) = wsum(ijp) + abs(Vf) + are*sqrt(c2)
      wsum(ijn) = wsum(ijn) + abs(Vf) + are*sqrt(c2)

      ! Flux and diffusion coefficient with the blended density
      flmass(i) = rhof*Vf
      cap = cap*rhof/rhoc
      can = cap

      ! Density-pressure coupling: rho' = p'/(R T)
      CP = 1.0_dp/( Rair*T(ijp) + small )
      CN = 1.0_dp/( Rair*T(ijn) + small )
      cvP(i) = Vf*( beta*(1.0_dp-facint(i)) + (1.0_dp-beta)*hP )*CP
      cvN(i) = Vf*( beta*facint(i)          + (1.0_dp-beta)*(1.0_dp-hP) )*CN

    else
      cvP(i) = 0.0_dp
      cvN(i) = 0.0_dp
    endif

    capf(i) = cap

    ! m' = cap (p'_N - p'_P) + cvP p'_P + cvN p'_N
    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = cap + cvN(i)

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap - cvP(i)

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - cap + cvP(i)

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap - cvN(i)

    ! > Sources:
    su(ijp) = su(ijp) - flmass(i)
    su(ijn) = su(ijn) + flmass(i)

  end do


  !
  ! > Boundaries
  !

  ! Incompressible: scaled mass flow at outlet (global continuity)
  if ( .not.compressible .and. .not.const_mflux ) call adjustMassFlow

  cvB = 0.0_dp

  do ib = 1, numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i = 1, nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        if ( compressible ) flmass(iface) = den(ijb)*( u(ijb)*arx(iface) + v(ijb)*ary(iface) + w(ijb)*arz(iface) )
        if ( bc_kinds ) then
          if ( bckind(ib) == BCK_MASSFLOW ) flmass(iface) = -mdot_face(iface-numInnerFaces)   ! prescribed mass flow
        endif
        su(ijp) = su(ijp) - flmass(iface)
      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i = 1, nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        if ( compressible ) then
          Vb = u(ijp)*arx(iface) + v(ijp)*ary(iface) + w(ijp)*arz(iface)
          flmass(iface) = den(ijp)*Vb
          cvB(iface-numInnerFaces) = max( Vb, 0.0_dp )/( Rair*T(ijp) + small )
          k = diag(ijp)
          a(k) = a(k) + cvB(iface-numInnerFaces)
        endif
        su(ijp) = su(ijp) - flmass(iface)
      end do

    elseif ( bctype(ib) == 'pressure' ) then

      do i = 1, nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        if ( bc_kinds ) then
          if ( bcmode(iface-numInnerFaces) == MODE_EXTRAP ) then
            ! Supersonic exit of a pressure patch: everything extrapolated, mdot_b = rho_P U_P.S, as the outlet
            Vb = u(ijp)*arx(iface) + v(ijp)*ary(iface) + w(ijp)*arz(iface)
            flmass(iface) = den(ijp)*Vb
            cvB(iface-numInnerFaces) = max( Vb, 0.0_dp )/( Rair*T(ijp) + small )
            k = diag(ijp)
            a(k) = a(k) + cvB(iface-numInnerFaces)
            su(ijp) = su(ijp) - flmass(iface)
            pp(ijb) = 0.0_dp
            cycle
          elseif ( bckind(ib) /= BCK_DEFAULT ) then
            call pressure_bnd_flux_kind( ijp, ijb, iface, cap, flmass(iface) )
            capB(iface-numInnerFaces) = cap
            k = diag(ijp)
            a(k) = a(k) - cap
            su(ijp) = su(ijp) - flmass(iface)
            pp(ijb) = 0.0_dp
            cycle
          endif
        endif

        call facefluxmassPressBnd( ijp, ijb, xf(iface), yf(iface), zf(iface), &
                                   arx(iface), ary(iface), arz(iface), cap, flmass(iface) )

        k = diag(ijp)
        a(k) = a(k) - cap
        su(ijp) = su(ijp) - flmass(iface)

        ! Boundary value of the pressure correction
        pp(ijb) = 0.0_dp
      end do

    endif

  enddo


  !
  ! > Transient (density) term of the compressible continuity equation, BDF1/BDF2.
  !
  if ( compressible .and. ltransient ) then
    do inp = 1, numCells
      tr = vol(inp)/timestep
      k = diag(inp)
      a(k) = a(k) + tc0*tr/( Rair*T(inp) + small )
      su(inp) = su(inp) - ( tc0*den(inp) - tc1*den_o(inp) + tc2*den_oo(inp) )*tr
    end do
  endif

  !
  ! > Pseudo-time term for steady compressible flow: (Vol/(R T dtau)) p', dtau = CFL Vol / ( 0.5 sum_f (|V_f| + c |S_f|) )
  !
  if ( compressible .and. asimple_pseudoCFL > 0.0_dp ) then
    nouter = nouter + 1
    cfl = asimple_pseudoCFL
    if ( asimple_cflRamp > 0 .and. nouter < asimple_cflRamp ) &
      cfl = asimple_pseudoCFL**( real(nouter-1,dp)/real(asimple_cflRamp,dp) )
    do inp = 1, numCells
      k = diag(inp)
      a(k) = a(k) + 0.5_dp*wsum(inp)/( cfl*( Rair*T(inp) + small ) )
    end do
  endif

  if ( ltest ) write(6,'(19x,a,1pe10.3)') ' Initial sum  =', sum(su)


  ! Multiple pressure corrections loop (non-orthogonality) ***********************************************
  npc = max( npcor, 1 )
  do ipcorr = 1, npc

    pp(1:numCells) = 0.0_dp

    ! Solve for pressure correction
    if ( trim(asimple_lSolverP) == 'amg' ) then
      call asimple_amg_solve( numCells, ia, ja, a, su, pp, asimple_tolRelP, asimple_tolAbsP, asimple_maxitP, res0, res, its )
      resor(4) = res0
      nits_amg = its
      write(*,'(2x,a,1pe10.3,a,1pe10.3,a,i0)') 'AMG-FGMRES:  Solving for p, Initial residual = ', res0, &
                                                ', Final residual = ', res, ', No Iterations ', its
    else
      call csrsolve( asimple_lSolverP, pp, su, resor(4), asimple_maxitP, asimple_tolAbsP, asimple_tolRelP, 'p' )
    endif

    !
    ! Correct mass fluxes
    !
    do iface = 1, numInnerFaces
      ijp = owner(iface)
      ijn = neighbour(iface)
      flmass(iface) = flmass(iface) + capf(iface)*( pp(ijn)-pp(ijp) ) + cvP(iface)*pp(ijp) + cvN(iface)*pp(ijn)
    enddo

    do ib = 1, numBoundaries

      if ( bctype(ib) == 'pressure' ) then

        do i = 1, nfaces(ib)
          iface = startFace(ib) + i
          ijp = owner(iface)
          ijb = iBndValueStart(ib) + i
          if ( bc_kinds ) then
            if ( bcmode(iface-numInnerFaces) == MODE_EXTRAP ) then
              flmass(iface) = flmass(iface) + cvB(iface-numInnerFaces)*pp(ijp)
              cycle
            elseif ( bckind(ib) /= BCK_DEFAULT ) then
              flmass(iface) = flmass(iface) - capB(iface-numInnerFaces)*pp(ijp)
              cycle
            endif
          endif
          call facefluxmassCorrPressBnd( ijp, ijb, xf(iface), yf(iface), zf(iface), &
                                         arx(iface), ary(iface), arz(iface), flmass(iface) )
        end do

      elseif ( bctype(ib) == 'outlet' .and. compressible ) then

        do i = 1, nfaces(ib)
          iface = startFace(ib) + i
          ijp = owner(iface)
          flmass(iface) = flmass(iface) + cvB(iface-numInnerFaces)*pp(ijp)
        end do

      endif

    enddo

    !
    ! Correct velocity components, pressure and density
    !
    ppref = pp(pRefCell)
    do ib = 1, numBoundaries
      if ( bctype(ib) == 'pressure' ) then
        ppref = 0.0_dp
        exit
      endif
    enddo
    ! Compressible flow: the level of pressure is set by the equation of state and density, no reference
    if ( compressible ) ppref = 0.0_dp

    ! su,sv,sw = sum_f {-p'}_f S_fi
    call gradp_and_sources( pp )

    u(1:numCells) = u(1:numCells) + su*apu
    v(1:numCells) = v(1:numCells) + sv*apv
    w(1:numCells) = w(1:numCells) + sw*apw
    p(1:numCells) = p(1:numCells) + asimple_urfP*( pp(1:numCells) - ppref )

    if ( compressible ) call update_density

    ! Explicit correction of velocity at symmetry boundaries
    call updateVelocityAtBoundary

    ! Boundary velocity/density of the compressible BC kinds consistent with the corrected fluxes
    if ( bc_kinds ) call post_bc_kinds

    if ( ipcorr /= npc ) then
    !
    ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    !
      su = 0.0_dp
      do i = 1, numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)
        call fluxmc2( ijp, ijn, arx(i), ary(i), arz(i), fmcor )
        flmass(i) = flmass(i) + fmcor
        su(ijp) = su(ijp) - fmcor
        su(ijn) = su(ijn) + fmcor
      enddo
    endif

  enddo
  ! END: multiple pressure corrections loop

  ! Write continuity error report:
  call continuityErrors

  ! Mass flow through the open boundaries (inflow positive): inlet, outlet and pressure patches
  if ( ltest ) call report_boundary_mass_flow

  ! Correct driving force for a constant mass flow rate simulation:
  if ( const_mflux ) call constant_mass_flow_forcing

end subroutine


subroutine prepare_bc_kinds
!
!  Boundary state (p_b, T_b, rho_b) and face mode of every face with a compressible BC kind (implicit_bc_compressible),
!  Picard-lagged from the cell state with the under-relaxation cplbc_urf.
!
  implicit none
  integer :: ib, i, iface, ijp, ijb

  do ib = 1, numBoundaries
    if ( bckind(ib) == BCK_DEFAULT ) cycle
    do i = 1, nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      ijb = iBndValueStart(ib) + i
      call bc_prepare_face( ib, i, iface, ijp, ijb )
    enddo
  enddo

end subroutine


subroutine pressure_bnd_flux_kind( ijp, ijb, iface, cap, fluxmass )
!
!  Mass flux of a pressure-type face with known p_b: Rhie-Chow between cell centre and face with the boundary density
!  rho_b = p_b/(R T_b) (not the cell density as in the stock routine); cap < 0 is the coefficient of the p' equation.
!
  implicit none
  integer, intent(in) :: ijp, ijb, iface
  real(dp), intent(out) :: cap, fluxmass
  real(dp) :: xpn, ypn, zpn, capp, dpcor, sx, sy, sz

  sx = arx(iface); sy = ary(iface); sz = arz(iface)
  xpn = xf(iface)-xc(ijp); ypn = yf(iface)-yc(ijp); zpn = zf(iface)-zc(ijp)
  capp = vol(ijp)*apu(ijp)/( sx*xpn + sy*ypn + sz*zpn )
  dpcor = p(ijb) - p(ijp) - ( dPdxi(1,ijp)*xpn + dPdxi(2,ijp)*ypn + dPdxi(3,ijp)*zpn )

  u(ijb) = u(ijp) - sx*capp*dpcor
  v(ijb) = v(ijp) - sy*capp*dpcor
  w(ijb) = w(ijp) - sz*capp*dpcor

  fluxmass = den(ijb)*( u(ijb)*sx + v(ijb)*sy + w(ijb)*sz )
  cap = -den(ijb)*( sx*sx + sy*sy + sz*sz )*capp

end subroutine


subroutine post_bc_kinds
!
!  After the corrections: boundary pressure/density/velocity of the compressible BC kinds used in the next momentum
!  assembly. Mass flow inlet and velocity inlet: p_b = p_P, U_b from the mass flow. Pressure-type kinds with inflow:
!  velocity along the inflow direction with the magnitude from the mass flux (free stream tangential part for 'riemann').
!
  implicit none
  integer :: ib, i, iface, ijp, ijb, ib2
  real(dp) :: sx, sy, sz, are, nx, ny, nz, ain, q, rb, uninf, utx, uty, utz

  do ib = 1, numBoundaries

    if ( bckind(ib) == BCK_DEFAULT ) cycle

    do i = 1, nfaces(ib)

      iface = startFace(ib) + i
      ijp = owner(iface)
      ijb = iBndValueStart(ib) + i
      ib2 = iface - numInnerFaces
      sx = arx(iface); sy = ary(iface); sz = arz(iface)

      select case ( bckind(ib) )

      case ( BCK_COMP_INLET )
        p(ijb) = p(ijp)
        den(ijb) = p(ijb)/( Rair*T(ijb) + small )

      case ( BCK_MASSFLOW )
        p(ijb) = p(ijp)
        den(ijb) = p(ijb)/( Rair*T(ijb) + small )
        ain = -( edir(1,ib2)*sx + edir(2,ib2)*sy + edir(3,ib2)*sz )
        q = mdot_face(ib2)/( den(ijb)*ain )
        u(ijb) = edir(1,ib2)*q; v(ijb) = edir(2,ib2)*q; w(ijb) = edir(3,ib2)*q
        flmass(iface) = -mdot_face(ib2)

      case ( BCK_TOTAL_PRESSURE, BCK_STAGNATION, BCK_OPEN, BCK_RIEMANN )
        if ( bcmode(ib2) == MODE_EXTRAP ) then
          u(ijb) = u(ijp); v(ijb) = v(ijp); w(ijb) = w(ijp)
        elseif ( flmass(iface) < 0.0_dp ) then
          rb = den(ijb)
          are = sqrt( sx*sx + sy*sy + sz*sz )
          nx = sx/are; ny = sy/are; nz = sz/are
          if ( bckind(ib) == BCK_RIEMANN ) then
            uninf = U_set(1,ib2)*nx + U_set(2,ib2)*ny + U_set(3,ib2)*nz
            utx = U_set(1,ib2) - uninf*nx; uty = U_set(2,ib2) - uninf*ny; utz = U_set(3,ib2) - uninf*nz
            q = -flmass(iface)/( rb*are )
            u(ijb) = -q*nx + utx; v(ijb) = -q*ny + uty; w(ijb) = -q*nz + utz
          else
            ain = -( edir(1,ib2)*sx + edir(2,ib2)*sy + edir(3,ib2)*sz )
            q = -flmass(iface)/( rb*ain )
            u(ijb) = edir(1,ib2)*q; v(ijb) = edir(2,ib2)*q; w(ijb) = edir(3,ib2)*q
          endif
        endif

      end select

    enddo
  enddo

end subroutine


subroutine report_boundary_mass_flow
!
!  Mass flow through inlet, outlet and pressure boundaries [kg/s], inflow positive; for diagnostics (ltest = T).
!
  implicit none
  integer :: ib, i, iface
  real(dp) :: fin, fout, fpr, dmdt

  fin = 0.0_dp;  fout = 0.0_dp;  fpr = 0.0_dp
  do ib = 1, numBoundaries
    do i = 1, nfaces(ib)
      iface = startFace(ib) + i
      select case ( trim(bctype(ib)) )
      case ( 'inlet' );    fin  = fin  - flmass(iface)
      case ( 'outlet' );   fout = fout - flmass(iface)
      case ( 'pressure' ); fpr  = fpr  - flmass(iface)
      end select
    enddo
  enddo
  dmdt = 0.0_dp
  if ( compressible .and. ltransient ) &
    dmdt = sum( ( tc0*den(1:numCells) - tc1*den_o(1:numCells) + tc2*den_oo(1:numCells) )*vol )/timestep
  write(*,'(2x,a,3(1pe11.3),a,2(1pe11.3))') 'mass flow in (inlet/outlet/pressure) =', fin, fout, fpr, &
                                             ',  net inflow, d(mass)/dt =', fin+fout+fpr, dmdt

end subroutine

end module simple_allspeeds
