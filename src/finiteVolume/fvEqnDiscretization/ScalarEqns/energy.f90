module energy_equation
!
! Purpose:
!   Energy equation assembly and solution for compressible (and incompressible) flows.
!
! Description:
!   - Assembles the conservative energy equation for one of four energy variables E:
!
!       TOTAL_ENTHALPY  H = h + K  :  d(rho H)/dt + div(rho U H) - dp/dt              = div(k grad T) + div(tau.U)
!       TOTAL_ENERGY    E = e + K  :  d(rho E)/dt + div(rho U (E + p/rho))            = div(k grad T) + div(tau.U)
!       ENTHALPY        h          :  d(rho h)/dt + div(rho U h) - dp/dt - U.grad p   = div(k grad T) + tau:grad U
!       INTERNAL_ENERGY e          :  d(rho e)/dt + div(rho U e) + p div U            = div(k grad T) + tau:grad U
!
!     with K = |U|^2/2, h = h(T) = integral of cp(T) dT, e = h - R T, ideal gas.
!   - Solves for the energy variable, inverts energy -> temperature with Newton-Raphson (handles cp = cp(T)).
!   - Convection: implicit first order upwind, higher order schemes (CDS, LUDS, SOU, QUICK, MUSCL with limiters)
!     by deferred correction. The implicit upwind part acts on the energy variable, the explicit correction is
!     the difference between the higher order flux and the upwind flux evaluated with lagged values. For the
!     TOTAL_ENERGY formulation the convected quantity is the total enthalpy E + p/rho (pressure work is in the flux,
!     it is not added again as a source), see the remark below - this is the subtle point that is easy to miss.
!   - Conduction: the energy variable is not a temperature, therefore the implicit diffusion coefficient is
!     k/c_x (c_x = cp for the enthalpy forms, cv for the energy forms) applied to E, and the difference between the
!     true heat flux k grad T and its implicit counterpart evaluated with lagged values is added explicitly
!     (this also removes the gradient of K from the total forms, and handles variable cp). At convergence the
!     flux is exactly k grad T.
!   - Viscous terms: div(tau.U) face-wise (total forms), tau:grad U cell-wise (static forms). Optional, addViscDiss.
!   - Unsteady term: BDF1 / BDF2 (coefficients tc0,tc1,tc2) or steady. Conservative (reference) or non-conservative
!     form of the convective / unsteady operator (conservativeEn). In the conservative form the mass fluxes must
!     satisfy the continuity equation (they do after the coupled p-U solve).
!   - Thermodynamic properties: cp(T) models (constant, polynomial, NASA 7-coefficient, Shomate), cv = cp - R,
!     Sutherland viscosity, k = cp mu / Pr.
!
! USAGE:
!
!   use energy_equation
!
!   ! Set convection scheme (also through the namelist: cSchemeEn, limiterEn)
!   call set_convection_scheme(MUSCL, LIMITER_VANLEER)
!
!   ! Or for maximum accuracy
!   call set_convection_scheme(QUICK)
!
!   ! For shock-capturing
!   call set_convection_scheme(MUSCL, LIMITER_SUPERBEE)
!
!   ! One call does it all (assemble, solve, temperature, boundary values, properties, density)
!   call calc_energy
!
!   ! or step by step
!   call assemble_energy_equation( TOTAL_ENTHALPY )
!   call solve_energy_equation( lSolverEn, maxiterEn, tolAbsEn, tolRelEn )
!   call compute_temperature_from_energy( TOTAL_ENTHALPY )
!   call update_thermodynamic_properties
!
!
! Deferred correction approach:
!
!   ! IMPLICIT: first order upwind (stable), written for the non-conservative / conservative form
!   A(P,N) += min(mdot,0)        A(N,P) += -max(mdot,0)           (matrix convention of the code: off-diagonals,
!   diagonal = sp - sum of off-diagonals; the conservative form adds  sp(P) += mdot, sp(N) -= mdot )
!
!   ! EXPLICIT: higher order correction
!   deferred = gdsEn * mdot * (phi_f_HO - phi_f_UDS)
!   su(P) -= deferred ;  su(N) += deferred
!
!   Remark: the implicit part contains both the upwind face value and the mass flux, i.e. mdot and its split
!   into outflow/inflow parts max(mdot,0), min(mdot,0) have to be consistent with the sign convention of the
!   matrix (the reference file uses the opposite sign for the off-diagonal).
!
! Differences with respect to the reference file energy_equation_HO_schemes.f90 (bugs fixed):
!   - the implicit upwind coefficients of the reference file cancel (off-diagonal = -(mdot - cap + can) = 0);
!   - the diffusion coefficient is k/c_x for the energy variable, not k;
!   - -dp/dt had the wrong sign, and p div U was added to the total energy formulation although the pressure work is
!     already in the convective flux (double counting); viscous dissipation tau:grad U was added to the total forms where
!     div(tau.U) belongs;
!   - SOU used 0.5 grad.(x_f - x_C), which is a quarter of the Warming-Beam correction; here SOU = phi_C + (phi_C - phi_U)/2
!     and LUDS = phi_C + grad(phi_C).(x_f - x_C);
!   - the gradient ratio r of the limited schemes is the one of Jasak (2 grad(phi_C).d / (phi_D - phi_C) - 1);
!   - Fortran: declarations inside CASE blocks, undefined parameters.
!
! No OpenMP.
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use sparse_matrix
  use linear_solvers

  implicit none

  private

  !
  ! Energy variables
  !
  integer, parameter, public :: TOTAL_ENERGY    = 1    ! Solve for E = e + 0.5*u^2
  integer, parameter, public :: INTERNAL_ENERGY = 2    ! Solve for e
  integer, parameter, public :: ENTHALPY        = 3    ! Solve for h
  integer, parameter, public :: TOTAL_ENTHALPY  = 4    ! Solve for H = h + 0.5*u^2

  ! Convection scheme types
  integer, parameter, public :: UDS   = 1    ! First-order upwind
  integer, parameter, public :: CDS   = 2    ! Central differencing
  integer, parameter, public :: SOU   = 3    ! Second-order upwind
  integer, parameter, public :: QUICK = 4    ! Quadratic upwind
  integer, parameter, public :: MUSCL = 5    ! Monotonic upwind
  integer, parameter, public :: LUDS  = 6    ! Linear upwind

  ! Limiter types for MUSCL
  integer, parameter, public :: LIMITER_NONE     = 0
  integer, parameter, public :: LIMITER_VANLEER  = 1
  integer, parameter, public :: LIMITER_MINMOD   = 2
  integer, parameter, public :: LIMITER_SUPERBEE = 3
  integer, parameter, public :: LIMITER_VANALBADA = 4

  ! Specific heat models
  integer, parameter :: CP_CONSTANT   = 1
  integer, parameter :: CP_POLYNOMIAL = 2
  integer, parameter :: CP_NASA       = 3
  integer, parameter :: CP_SHOMATE    = 4

  !
  ! Discretization and solution parameters - modified through the input.nml file
  !
  logical  :: calcEn = .False.                      ! To activate the solution of this field in the main function.
  real(dp) :: urfEn  = 0.9_dp                       ! Under-relaxation factor.
  real(dp) :: gdsEn = 1.0_dp                        ! Deferred correction factor.
  character ( len=30 ) :: cSchemeEn = 'linearUpwind'! Convection scheme: uds, cds, luds (linearUpwind), sou, quick, muscl.
  character ( len=30 ) :: limiterEn = 'vanLeer'     ! Limiter of the MUSCL scheme: none, vanLeer, minmod, superbee, vanAlbada.
  character ( len=12 ) :: dSchemeEn = 'skewness'    ! Diffusion scheme (kept for compatibility, not used).
  integer :: nrelaxEn = 0                           ! Non-orthogonal correction type (kept for compatibility, not used).
  character ( len=12 ) :: lSolverEn = 'bicgstab'    ! Linear algebraic solver.
  integer  :: maxiterEn = 10                        ! Max number of iterations in linear solver.
  real(dp) :: tolAbsEn = 1e-13                      ! Absolute residual level.
  real(dp) :: tolRelEn = 0.025_dp                   ! Relative drop in residual to exit linear solver.
  real(dp) :: sigtEn = 0.85_dp                      ! sigma_t (kept for compatibility, not used)
  logical :: solveTotalEnthalpy = .True.            ! What we solve here - default is total enthalpy, change it
  logical :: solveTotalEnergy = .False.             ! in input.nml file. Priority: enthalpy, internal energy,
  logical :: solveInternalEnergy = .False.          ! total energy, total enthalpy.
  logical :: solveEnthalpy = .False.                !
  logical :: addViscDiss = .False.                  ! Add viscous terms (dissipation / work of viscous stresses)? T/F
  logical :: conservativeEn = .False.               ! Conservative form of the convective term (adds m to the diagonal: diverges in cold-start pseudo-time
                                                   ! iterations while the mass flux is unbalanced; same converged solution), or non-conservative (default).

  ! Specific heat
  character ( len=12 ) :: cpModel = 'constant'      ! constant, polynomial, nasa, shomate
  real(dp) :: cp_const = 1006.0_dp                  ! Heat capacity at constant pressure [J/kgK] (air)
  real(dp) :: cp_Tref = 0.0_dp                      ! Reference temperature of the enthalpy, h(T_ref) = 0 [K] (constant, polynomial)
  real(dp) :: cp_poly(0:3) = (/ 1006.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)  ! cp = c0 + c1 T + c2 T^2 + c3 T^3
  ! NASA 7-coefficient polynomials (default: N2, GRI-Mech 3.0), cp/R_nasa = a1 + a2 T + ... , R_nasa = nasa_R
  real(dp) :: nasa_lo(7) = (/ 3.298677_dp, 1.408240e-3_dp, -3.963222e-6_dp, 5.641515e-9_dp, -2.444854e-12_dp, &
                              -1.020900e3_dp, 3.950372_dp /)
  real(dp) :: nasa_hi(7) = (/ 2.926640_dp, 1.487977e-3_dp, -5.684761e-7_dp, 1.009704e-10_dp, -6.753351e-15_dp, &
                              -9.227977e2_dp, 5.980528_dp /)
  real(dp) :: nasa_R = 296.80_dp                    ! specific gas constant of the NASA polynomials [J/kgK] (N2)
  ! Shomate equation (default: N2, NIST, 298-6000 K), cp [J/molK] = A + B t + C t^2 + D t^3 + E/t^2, t = T/1000
  real(dp) :: shomate(6) = (/ 19.50583_dp, 19.88705_dp, -8.598535_dp, 1.369784_dp, 0.527601_dp, -4.935202_dp /)
  real(dp) :: shomate_M = 0.0280134_dp              ! molar mass [kg/mol]

  ! Sutherland viscosity law, mu = viscos (T/T_ref)^1.5 (T_ref + S)/(T + S)
  logical  :: sutherland = .False.
  real(dp) :: mu_Tref = 288.15_dp
  real(dp) :: mu_S = 110.4_dp

  real(dp), parameter :: T_floor = 50.0_dp, T_ceil = 6000.0_dp   ! physical bounds of the temperature

  ! Current settings
  integer :: energy_formulation = TOTAL_ENTHALPY
  integer :: convection_scheme = LUDS
  integer :: limiter_type = LIMITER_VANLEER

  ! Fields owned by this module
  real(dp), allocatable, save :: E(:), E_o(:), E_oo(:)       ! energy variable and its time levels
  real(dp), allocatable, save :: dEdxi(:,:), dCvdxi(:,:)
  real(dp), allocatable, save :: Cv(:)                       ! convected quantity (E, or E + p/rho for TOTAL_ENERGY)
  real(dp), allocatable, save :: Kf(:)                       ! kinetic energy
  real(dp), allocatable, save :: cpf(:), cxf(:), kth(:)      ! cp, c_x, thermal conductivity
  real(dp), allocatable, save :: den_o(:), den_oo(:), p_o(:), p_oo(:)
  real(dp), save :: tc0 = 0.0_dp, tc1 = 0.0_dp, tc2 = 0.0_dp ! BDF coefficients: (tc0 x - tc1 xo + tc2 xoo)/dt
  integer, save :: last_itime = -1, nsteps_done = 0
  logical, save :: energy_ready = .false.

  public :: calc_energy, assemble_energy_equation, solve_energy_equation, compute_temperature_from_energy, &
            update_thermodynamic_properties, set_convection_scheme, energy_init, energy_time_levels, &
            compute_cp, compute_cv, compute_enthalpy, compute_internal_energy, gamma_of_T, speed_of_sound_sq, &
            energy_formulation_selected

  public :: calcEn, urfEn, gdsEn, cSchemeEn, limiterEn, dSchemeEn, nrelaxEn, lSolverEn, maxiterEn, tolAbsEn, tolRelEn, &
            sigtEn, solveEnthalpy, solveInternalEnergy, solveTotalEnergy, solveTotalEnthalpy, addViscDiss, conservativeEn, &
            cpModel, cp_const, cp_Tref, cp_poly, nasa_lo, nasa_hi, nasa_R, shomate, shomate_M, &
            sutherland, mu_Tref, mu_S
  public :: den_o, den_oo, p_o, p_oo, tc0, tc1, tc2, energy_ready

contains


!***********************************************************************
!
subroutine calc_energy
!
! Purpose:
!   Assemble and solve the energy equation, recover temperature, update boundary values,
!   thermodynamic properties, density (compressible flow) and viscosity (Sutherland).
!
  implicit none

  integer :: formulation

  if ( .not. energy_ready ) call energy_init

  formulation = energy_formulation_selected()

  ! Time levels are shifted once per time step (the coupled p-U solver does it before its own step)
  call energy_time_levels

  call assemble_energy_equation( formulation )

  call solve_energy_equation( lSolverEn, maxiterEn, tolAbsEn, tolRelEn )

  call compute_temperature_from_energy( formulation )

  call update_thermodynamic_properties

  ! Copy for output (the 'Energy' field of the VTK files)
  if ( allocated(En) ) En = E

end subroutine calc_energy


!***********************************************************************
!
function energy_formulation_selected() result(formulation)
!
! Description:
!   Energy variable chosen through the namelist.
!
  implicit none
  integer :: formulation

  if ( solveEnthalpy ) then
    formulation = ENTHALPY
  elseif ( solveInternalEnergy ) then
    formulation = INTERNAL_ENERGY
  elseif ( solveTotalEnergy ) then
    formulation = TOTAL_ENERGY
  else
    formulation = TOTAL_ENTHALPY
  endif

end function energy_formulation_selected


!***********************************************************************
!
subroutine energy_init
!
! Description:
!   Allocation, scheme selection from the input strings, consistent initial state.
!
  implicit none

  allocate( E(numTotal), E_o(numTotal), E_oo(numTotal), Cv(numTotal), Kf(numTotal) )
  allocate( dEdxi(3,numTotal), dCvdxi(3,numTotal) )
  allocate( cpf(numTotal), cxf(numTotal), kth(numTotal) )
  allocate( den_o(numTotal), den_oo(numTotal), p_o(numTotal), p_oo(numTotal) )

  dEdxi = 0.0_dp
  dCvdxi = 0.0_dp

  call select_scheme_from_input

  energy_formulation = energy_formulation_selected()

  call update_thermodynamic_properties
  call kinetic_energy
  call energy_from_temperature( energy_formulation )

  E_o = E;  E_oo = E
  den_o = den;  den_oo = den
  p_o = p;  p_oo = p

  energy_ready = .true.

end subroutine energy_init


!***********************************************************************
!
subroutine select_scheme_from_input
!
! Description:
!   Convection scheme and limiter names of the namelist -> integer codes.
!
  implicit none
  character(len=30) :: s

  s = adjustl( cSchemeEn )
  call lowercase( s )
  select case ( trim(s) )
  case ('uds','upwind')
    convection_scheme = UDS
  case ('cds','central','linear')
    convection_scheme = CDS
  case ('sou')
    convection_scheme = SOU
  case ('quick')
    convection_scheme = QUICK
  case ('muscl')
    convection_scheme = MUSCL
  case default                      ! linearUpwind, luds and anything else of the segregated solvers
    convection_scheme = LUDS
  end select

  s = adjustl( limiterEn )
  call lowercase( s )
  select case ( trim(s) )
  case ('none')
    limiter_type = LIMITER_NONE
  case ('minmod')
    limiter_type = LIMITER_MINMOD
  case ('superbee')
    limiter_type = LIMITER_SUPERBEE
  case ('vanalbada')
    limiter_type = LIMITER_VANALBADA
  case default
    limiter_type = LIMITER_VANLEER
  end select

end subroutine select_scheme_from_input


subroutine lowercase( s )
  implicit none
  character(len=*), intent(inout) :: s
  integer :: i, ic
  do i=1,len_trim(s)
    ic = iachar( s(i:i) )
    if ( ic >= 65 .and. ic <= 90 ) s(i:i) = achar( ic+32 )
  enddo
end subroutine lowercase


!***********************************************************************
!
subroutine energy_time_levels
!
! Description:
!   Once per time step (detected by the change of itime): shift the time levels of density, pressure and the energy
!   variable and set the BDF coefficients. The first step of a run is BDF1. Steady runs: coefficients are zero.
!   Time levels of E are taken from the current temperature and velocity (end of the previous time step).
!
  implicit none

  if ( .not. energy_ready ) call energy_init

  if ( .not. ltransient ) then
    tc0 = 0.0_dp; tc1 = 0.0_dp; tc2 = 0.0_dp
    return
  endif

  if ( itime == last_itime ) return
  last_itime = itime

  call kinetic_energy
  call energy_from_temperature( energy_formulation )

  den_oo = den_o;  p_oo = p_o;  E_oo = E_o
  den_o  = den;    p_o  = p;    E_o  = E

  nsteps_done = nsteps_done + 1

  if ( bdf2 .and. nsteps_done > 1 ) then
    tc0 = 1.5_dp; tc1 = 2.0_dp; tc2 = 0.5_dp
  else
    tc0 = 1.0_dp; tc1 = 1.0_dp; tc2 = 0.0_dp
  endif

end subroutine energy_time_levels


!***********************************************************************
!
subroutine kinetic_energy
  implicit none
  Kf = 0.5_dp*( u**2 + v**2 + w**2 )
end subroutine kinetic_energy


!***********************************************************************
!
subroutine energy_from_temperature( formulation )
!
! Description:
!   Energy variable, cells and boundary faces, from temperature and (boundary) velocity.
!
  implicit none
  integer, intent(in) :: formulation
  integer :: i

  do i=1,numTotal
    select case ( formulation )
    case ( TOTAL_ENERGY )
      E(i) = compute_internal_energy( T(i) ) + Kf(i)
    case ( INTERNAL_ENERGY )
      E(i) = compute_internal_energy( T(i) )
    case ( ENTHALPY )
      E(i) = compute_enthalpy( T(i) )
    case default
      E(i) = compute_enthalpy( T(i) ) + Kf(i)
    end select
  enddo

end subroutine energy_from_temperature


!***********************************************************************
!
function energy_of_T( Tl, Kl, formulation ) result(Er)
!
! Description:
!   Energy variable for a given temperature and kinetic energy.
!
  implicit none
  real(dp), intent(in) :: Tl, Kl
  integer, intent(in) :: formulation
  real(dp) :: Er

  select case ( formulation )
  case ( TOTAL_ENERGY )
    Er = compute_internal_energy( Tl ) + Kl
  case ( INTERNAL_ENERGY )
    Er = compute_internal_energy( Tl )
  case ( ENTHALPY )
    Er = compute_enthalpy( Tl )
  case default
    Er = compute_enthalpy( Tl ) + Kl
  end select

end function energy_of_T


!***********************************************************************
!
subroutine assemble_energy_equation( formulation )
!
! Description:
!   Assembles the energy equation for the chosen energy variable (see the header of the module).
!   Matrix (a, su, sp) of the sparse_matrix module, CSR, unknown E.
!
  implicit none

  integer, intent(in) :: formulation

  integer :: i, k, ib, inp, ijp, ijn, ijb, iface
  real(dp) :: lam, sx, sy, sz, dx, dy, dz, Dfb
  real(dp) :: m, apot, dpdt
  real(dp) :: gam_f, k_f, d_diff
  real(dp) :: phi_ho, phi_uds, deferred
  real(dp) :: true_flux, impl_flux
  real(dp) :: urfr, urfm, off_diagonal_terms
  real(dp) :: mu_f, dux(3), dvx(3), dwx(3), divu, tx, ty, tz, uf, vf, wf
  real(dp) :: Eb, are2, sdb, de, Ew
  real(dp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, Phi, Psi, divUc
  logical :: total_form, enthalpy_form, isothermal

  energy_formulation = formulation
  if ( .not. energy_ready ) call energy_init

  total_form    = ( formulation == TOTAL_ENTHALPY .or. formulation == TOTAL_ENERGY )
  enthalpy_form = ( formulation == TOTAL_ENTHALPY .or. formulation == ENTHALPY )

  ! State: kinetic energy, energy variable consistent with T and the current velocity, properties
  call kinetic_energy
  call energy_from_temperature( formulation )
  call update_thermodynamic_properties

  ! Convected quantity
  if ( formulation == TOTAL_ENERGY ) then
    Cv = E + p/( den + small )
  else
    Cv = E
  endif

  ! Gradients
  call grad( E, dEdxi )
  if ( formulation == TOTAL_ENERGY ) then
    call grad( Cv, dCvdxi )
  else
    dCvdxi = dEdxi
  endif
  call grad( T, dTdxi )
  if ( addViscDiss .or. formulation == INTERNAL_ENERGY .or. formulation == ENTHALPY ) then
    call grad( u, dUdxi )
    call grad( v, dVdxi )
    call grad( w, dWdxi )
  endif

  a  = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

  !---------------------------------------------------------------------------
  ! TRANSIENT TERM: d(rho E)/dt, BDF1 / BDF2 (conservative) or rho dE/dt (non-conservative)
  !---------------------------------------------------------------------------
  if ( ltransient ) then
    do inp=1,numCells
      apot = den(inp)*vol(inp)/timestep
      sp(inp) = sp(inp) + tc0*apot
      if ( conservativeEn ) then
        su(inp) = su(inp) + vol(inp)/timestep*( tc1*den_o(inp)*E_o(inp) - tc2*den_oo(inp)*E_oo(inp) )
      else
        su(inp) = su(inp) + apot*( tc1*E_o(inp) - tc2*E_oo(inp) )
      endif
    enddo
  endif

  !---------------------------------------------------------------------------
  ! SOURCE TERMS: pressure work and viscous dissipation (cell-wise forms)
  !---------------------------------------------------------------------------
  do inp=1,numCells

    ! dp/dt
    if ( ltransient ) then
      dpdt = ( tc0*p(inp) - tc1*p_o(inp) + tc2*p_oo(inp) )/timestep
    else
      dpdt = 0.0_dp
    endif

    select case ( formulation )

    case ( TOTAL_ENTHALPY )
      su(inp) = su(inp) + vol(inp)*dpdt

    case ( ENTHALPY )
      ! Dp/Dt = dp/dt + U.grad p
      su(inp) = su(inp) + vol(inp)*( dpdt + u(inp)*dPdxi(1,inp) + v(inp)*dPdxi(2,inp) + w(inp)*dPdxi(3,inp) )

    case ( INTERNAL_ENERGY )
      ! - p div U
      divUc = dUdxi(1,inp) + dVdxi(2,inp) + dWdxi(3,inp)
      su(inp) = su(inp) - vol(inp)*p(inp)*divUc

    end select

    ! Viscous dissipation tau:grad U (static forms)
    if ( addViscDiss .and. .not. total_form ) then
      dudx = dUdxi(1,inp); dudy = dUdxi(2,inp); dudz = dUdxi(3,inp)
      dvdx = dVdxi(1,inp); dvdy = dVdxi(2,inp); dvdz = dVdxi(3,inp)
      dwdx = dWdxi(1,inp); dwdy = dWdxi(2,inp); dwdz = dWdxi(3,inp)
      Psi = ( dudx + dvdy + dwdz )**2
      Phi = 2.0_dp*( dudx**2 + dvdy**2 + dwdz**2 ) + ( dudy + dvdx )**2 + ( dudz + dwdx )**2 + ( dwdy + dvdz )**2
      su(inp) = su(inp) + vol(inp)*vis(inp)*( Phi - twothirds*Psi )
    endif

  enddo

  !---------------------------------------------------------------------------
  ! INNER FACES: convection + diffusion
  !---------------------------------------------------------------------------
  do i=1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)

    sx = arx(i); sy = ary(i); sz = arz(i)
    dx = xc(ijn)-xc(ijp); dy = yc(ijn)-yc(ijp); dz = zc(ijn)-zc(ijp)

    !-------------------------------------------------------------------
    ! CONVECTION TERM: div(rho*U*E) or div(rho*U*(E + p/rho)), mass flux from the p-U solver
    !-------------------------------------------------------------------
    m = flmass(i)

    ! IMPLICIT PART: first-order upwind (UDS)
    k = icell_jcell_csr_index(i)
    a(k) = a(k) + min( m, 0.0_dp )
    k = jcell_icell_csr_index(i)
    a(k) = a(k) - max( m, 0.0_dp )

    if ( conservativeEn ) then
      ! diagonal = sp - sum of off-diagonals, adds sum of outflow fluxes only if sum(mdot) is added to sp
      sp(ijp) = sp(ijp) + m
      sp(ijn) = sp(ijn) - m
    endif

    ! EXPLICIT PART: higher-order correction (deferred), phi_f_uds of the convected variable
    if ( m >= 0.0_dp ) then
      phi_uds = Cv(ijp)
      deferred = m*( Cv(ijp) - E(ijp) )          ! = m p/rho for TOTAL_ENERGY, zero otherwise
    else
      phi_uds = Cv(ijn)
      deferred = m*( Cv(ijn) - E(ijn) )
    endif

    if ( convection_scheme /= UDS ) then
      phi_ho = compute_face_value_ho( ijp, ijn, i, m, convection_scheme, limiter_type )
      deferred = deferred + gdsEn*m*( phi_ho - phi_uds )
    endif

    su(ijp) = su(ijp) - deferred
    su(ijn) = su(ijn) + deferred

    ! Non-conservative form of the p/rho part of the TOTAL_ENERGY convected quantity:  - (p/rho)_C sum(mdot)
    if ( .not. conservativeEn .and. formulation == TOTAL_ENERGY ) then
      su(ijp) = su(ijp) + m*p(ijp)/( den(ijp)+small )
      su(ijn) = su(ijn) - m*p(ijn)/( den(ijn)+small )
    endif

    !-------------------------------------------------------------------
    ! DIFFUSION TERM: div(k grad T), implicit k/c_x on E, explicit remainder
    !-------------------------------------------------------------------
    k_f   = kth(ijp) + ( kth(ijn)-kth(ijp) )*lam
    gam_f = ( kth(ijp)/cxf(ijp) )*( 1.0_dp-lam ) + ( kth(ijn)/cxf(ijn) )*lam
    d_diff = gam_f*Df(i)

    k = icell_jcell_csr_index(i)
    a(k) = a(k) - d_diff
    k = jcell_icell_csr_index(i)
    a(k) = a(k) - d_diff

    ! True heat flux into P through the face:  k [ Df (T_N - T_P) + grad T_f . (S - Df d) ]
    true_flux = k_f*( Df(i)*( T(ijn)-T(ijp) ) &
              + ( ( dTdxi(1,ijp)*(1.0_dp-lam)+dTdxi(1,ijn)*lam )*( sx - Df(i)*dx ) &
                + ( dTdxi(2,ijp)*(1.0_dp-lam)+dTdxi(2,ijn)*lam )*( sy - Df(i)*dy ) &
                + ( dTdxi(3,ijp)*(1.0_dp-lam)+dTdxi(3,ijn)*lam )*( sz - Df(i)*dz ) ) )

    ! Implicit counterpart evaluated with lagged values
    impl_flux = d_diff*( E(ijn)-E(ijp) )

    su(ijp) = su(ijp) + ( true_flux - impl_flux )
    su(ijn) = su(ijn) - ( true_flux - impl_flux )

    !-------------------------------------------------------------------
    ! Work of the viscous stresses  div(tau.U)  (total forms), explicit
    !-------------------------------------------------------------------
    if ( addViscDiss .and. total_form ) then
      mu_f = vis(ijp) + ( vis(ijn)-vis(ijp) )*lam
      dux = dUdxi(:,ijp)*(1.0_dp-lam) + dUdxi(:,ijn)*lam
      dvx = dVdxi(:,ijp)*(1.0_dp-lam) + dVdxi(:,ijn)*lam
      dwx = dWdxi(:,ijp)*(1.0_dp-lam) + dWdxi(:,ijn)*lam
      divu = dux(1)+dvx(2)+dwx(3)
      tx = mu_f*( 2*dux(1)*sx + (dux(2)+dvx(1))*sy + (dux(3)+dwx(1))*sz ) - twothirds*mu_f*divu*sx
      ty = mu_f*( (dvx(1)+dux(2))*sx + 2*dvx(2)*sy + (dvx(3)+dwx(2))*sz ) - twothirds*mu_f*divu*sy
      tz = mu_f*( (dwx(1)+dux(3))*sx + (dwx(2)+dvx(3))*sy + 2*dwx(3)*sz ) - twothirds*mu_f*divu*sz
      uf = u(ijp)*(1.0_dp-lam) + u(ijn)*lam
      vf = v(ijp)*(1.0_dp-lam) + v(ijn)*lam
      wf = w(ijp)*(1.0_dp-lam) + w(ijn)*lam
      su(ijp) = su(ijp) + tx*uf + ty*vf + tz*wf
      su(ijn) = su(ijn) - ( tx*uf + ty*vf + tz*wf )
    endif

  enddo

  !---------------------------------------------------------------------------
  ! BOUNDARY FACES
  !---------------------------------------------------------------------------
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' .or. bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        sx = arx(iface); sy = ary(iface); sz = arz(iface)
        dx = xf(iface)-xc(ijp); dy = yf(iface)-yc(ijp); dz = zf(iface)-zc(ijp)
        Dfb = ( sx*sx+sy*sy+sz*sz )/( sx*dx+sy*dy+sz*dz )
        m = flmass(iface)

        ! Boundary value of the energy variable and of the convected quantity (from boundary T and U)
        Eb = energy_of_T( T(ijb), Kf(ijb), formulation )
        if ( formulation == TOTAL_ENERGY ) then
          Ew = Eb + p(ijb)/( den(ijb)+small )        ! convected: total enthalpy at the boundary
        else
          Ew = Eb
        endif

        !-----------------------------------------------------------------
        ! Convection
        !   inflow  (m < 0): known boundary value (T_b, U_b of the boundary condition)
        !   outflow (m >= 0): zero gradient, boundary value is the cell value
        !-----------------------------------------------------------------
        if ( m < 0.0_dp ) then
          su(ijp) = su(ijp) - m*Ew
          if ( .not. conservativeEn ) sp(ijp) = sp(ijp) - m
        else
          if ( conservativeEn ) sp(ijp) = sp(ijp) + m
          if ( formulation == TOTAL_ENERGY ) su(ijp) = su(ijp) - m*p(ijp)/( den(ijp)+small )
        endif
        if ( .not. conservativeEn .and. formulation == TOTAL_ENERGY ) su(ijp) = su(ijp) + m*p(ijp)/( den(ijp)+small )

        !-----------------------------------------------------------------
        ! Diffusion at inlets: Dirichlet temperature
        !-----------------------------------------------------------------
        if ( bctype(ib) == 'inlet' ) then
          gam_f  = kth(ijb)/cxf(ijb)
          d_diff = gam_f*Dfb
          sp(ijp) = sp(ijp) + d_diff
          su(ijp) = su(ijp) + d_diff*Eb
          ! explicit remainder: k Dfb (T_b - T_P) - gam Dfb (E_b - E_P), both with lagged values
          su(ijp) = su(ijp) + kth(ijb)*Dfb*( T(ijb)-T(ijp) ) - d_diff*( Eb-E(ijp) )
        endif

        !-----------------------------------------------------------------
        ! Viscous work with boundary velocity (inlet/outlet/pressure faces)
        !-----------------------------------------------------------------
        if ( addViscDiss .and. total_form ) then
          mu_f = vis(ijb)
          dux = dUdxi(:,ijp); dvx = dVdxi(:,ijp); dwx = dWdxi(:,ijp)
          divu = dux(1)+dvx(2)+dwx(3)
          tx = mu_f*( 2*dux(1)*sx + (dux(2)+dvx(1))*sy + (dux(3)+dwx(1))*sz ) - twothirds*mu_f*divu*sx
          ty = mu_f*( (dvx(1)+dux(2))*sx + 2*dvx(2)*sy + (dvx(3)+dwx(2))*sz ) - twothirds*mu_f*divu*sy
          tz = mu_f*( (dwx(1)+dux(3))*sx + (dwx(2)+dvx(3))*sy + 2*dwx(3)*sz ) - twothirds*mu_f*divu*sz
          su(ijp) = su(ijp) + tx*u(ijb) + ty*v(ijb) + tz*w(ijb)
        endif

      enddo

    elseif ( bctype(ib) == 'wall' ) then

      isothermal = ( bcDFraction(7,ib) == 1.0_dp )

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        sx = arx(iface); sy = ary(iface); sz = arz(iface)
        dx = xf(iface)-xc(ijp); dy = yf(iface)-yc(ijp); dz = zf(iface)-zc(ijp)
        Dfb = ( sx*sx+sy*sy+sz*sz )/( sx*dx+sy*dy+sz*dz )

        ! Work of the wall shear on the fluid: (traction).U_wall, non-zero only for a moving wall.
        ! traction = mu Df_b (U_b - U_P)_tangential
        if ( addViscDiss .and. total_form .and. ( abs(u(ijb))+abs(v(ijb))+abs(w(ijb)) ) > 0.0_dp ) then
          are2 = sx*sx+sy*sy+sz*sz
          tx = u(ijb)-u(ijp); ty = v(ijb)-v(ijp); tz = w(ijb)-w(ijp)
          sdb = ( tx*sx+ty*sy+tz*sz )/are2
          tx = tx - sdb*sx; ty = ty - sdb*sy; tz = tz - sdb*sz
          su(ijp) = su(ijp) + vis(ijb)*Dfb*( tx*u(ijb) + ty*v(ijb) + tz*w(ijb) )
        endif

        if ( isothermal ) then
          ! Wall energy with the kinetic energy of the cell, so that only the temperature difference drives the flux:
          ! flux = (k/c_x) Df_b [ E(T_w, K_P) - E_P ],  = k Df_b (T_w - T_P) at convergence.
          Ew = energy_of_T( T(ijb), Kf(ijp), formulation )
          gam_f  = kth(ijb)/cxf(ijb)
          de = gam_f*Dfb
          sp(ijp) = sp(ijp) + de
          su(ijp) = su(ijp) + de*Ew
          su(ijp) = su(ijp) + kth(ijb)*Dfb*( T(ijb)-T(ijp) ) - de*( Ew - E(ijp) )
        endif
        ! adiabatic: nothing (zero gradient of T)

      enddo

    endif

  enddo

  !---------------------------------------------------------------------------
  ! Diagonal, under-relaxation
  !---------------------------------------------------------------------------
  urfr = 1.0_dp/urfEn
  urfm = 1.0_dp-urfEn

  do inp=1,numCells
    off_diagonal_terms = sum( a( ia(inp) : ia(inp+1)-1 ) ) - a(diag(inp))
    a(diag(inp)) = sp(inp) - off_diagonal_terms
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = su(inp) + urfm*a(diag(inp))*E(inp)
  enddo

end subroutine assemble_energy_equation


!***********************************************************************
!
subroutine solve_energy_equation( solver_type, max_iter, tolAbs, tolRel )
!
! Description:
!   Solves the energy equation using the specified linear solver (see linear_solvers).
!
  implicit none

  character(len=*), intent(in) :: solver_type
  integer, intent(in) :: max_iter
  real(dp), intent(in) :: tolAbs, tolRel

  call csrsolve( solver_type, E, su, resor(7), max_iter, tolAbs, tolRel, 'Energy' )

end subroutine solve_energy_equation


!***********************************************************************
!
subroutine compute_temperature_from_energy( formulation )
!
! Description:
!   Inverts energy -> temperature, cells and boundary values.
!   Handles temperature-dependent specific heat: cp = cp(T), Newton-Raphson iteration for each cell.
!
!   Relations:
!     h(T) = integral(cp(T')dT', T_ref, T)
!     e(T) = h(T) - R T = integral(cv(T')dT', T_ref, T)
!
  implicit none

  integer, intent(in) :: formulation

  integer :: i, ib, iter, ijp, ijb
  integer, parameter :: max_iter_nr = 20
  real(dp) :: T_guess, T_new, ke
  real(dp) :: target_val, f_val, df_dT, tmin, tmax
  real(dp), parameter :: tol_nr = 1.0e-9_dp
  logical :: static_h

  do i=1,numCells

    T_guess = T(i)
    T_new = T_guess
    ke = Kf(i)

    select case ( formulation )
    case ( TOTAL_ENERGY )
      target_val = E(i) - ke
      static_h = .false.
    case ( INTERNAL_ENERGY )
      target_val = E(i)
      static_h = .false.
    case ( ENTHALPY )
      target_val = E(i)
      static_h = .true.
    case default
      target_val = E(i) - ke
      static_h = .true.
    end select

    !-----------------------------------------------------------------------
    ! NEWTON-RAPHSON ITERATION: find T such that h(T) = h_target (or e(T) = e_target)
    !-----------------------------------------------------------------------
    do iter = 1, max_iter_nr

      if ( static_h ) then
        f_val = compute_enthalpy( T_guess ) - target_val
        df_dT = compute_cp( T_guess )
      else
        f_val = compute_internal_energy( T_guess ) - target_val
        df_dT = compute_cv( T_guess )
      endif

      T_new = T_guess - f_val/df_dT
      T_new = max( T_new, T_floor )
      T_new = min( T_new, T_ceil )

      if ( abs(T_new - T_guess) < tol_nr*T_guess ) exit

      T_guess = T_new

    enddo

    T(i) = T_new

  enddo

  ! Boundary values: zero gradient on outlet, pressure, symmetry, empty (inlets keep their values);
  ! adiabatic walls: zero gradient.
  call updateBoundary( T )

  do ib=1,numBoundaries
    if ( bctype(ib) == 'wall' .and. bcDFraction(7,ib) /= 1.0_dp ) then
      do i=1,nfaces(ib)
        ijp = owner(startFace(ib)+i)
        ijb = iBndValueStart(ib)+i
        T(ijb) = T(ijp)
      enddo
    endif
  enddo

  tmin = minval( T(1:numCells) )
  tmax = maxval( T(1:numCells) )
  write(6,'(2x,es11.4,a,es11.4)') tmin,' <= Temperature <= ',tmax

end subroutine compute_temperature_from_energy


!***********************************************************************
!
function compute_cp( Tl ) result(cp)
!
! Description:
!   Computes specific heat at constant pressure as a function of temperature [J/kgK].
!
!   Options:
!     - Constant cp (simple)
!     - Polynomial
!     - NASA 7-coefficient polynomial (accurate)
!     - Shomate equation
!
  implicit none

  real(dp), intent(in) :: Tl
  real(dp) :: cp
  real(dp) :: ts

  select case ( cp_model_id() )

  case ( CP_POLYNOMIAL )
    cp = cp_poly(0) + cp_poly(1)*Tl + cp_poly(2)*Tl**2 + cp_poly(3)*Tl**3

  case ( CP_NASA )
    ! cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    if ( Tl < 1000.0_dp ) then
      cp = nasa_R*( nasa_lo(1) + nasa_lo(2)*Tl + nasa_lo(3)*Tl**2 + nasa_lo(4)*Tl**3 + nasa_lo(5)*Tl**4 )
    else
      cp = nasa_R*( nasa_hi(1) + nasa_hi(2)*Tl + nasa_hi(3)*Tl**2 + nasa_hi(4)*Tl**3 + nasa_hi(5)*Tl**4 )
    endif

  case ( CP_SHOMATE )
    ! cp = A + B*t + C*t^2 + D*t^3 + E/t^2  where t = T/1000 [J/molK] -> [J/kgK]
    ts = Tl/1000.0_dp
    cp = ( shomate(1) + shomate(2)*ts + shomate(3)*ts**2 + shomate(4)*ts**3 + shomate(5)/ts**2 )/shomate_M

  case default
    cp = cp_const

  end select

end function compute_cp


function cp_model_id() result(id)
  implicit none
  integer :: id
  select case ( trim(cpModel) )
  case ('polynomial')
    id = CP_POLYNOMIAL
  case ('nasa')
    id = CP_NASA
  case ('shomate')
    id = CP_SHOMATE
  case default
    id = CP_CONSTANT
  end select
end function cp_model_id


!***********************************************************************
!
function compute_cv( Tl ) result(cv)
!
! Description:
!   Computes specific heat at constant volume, cv = cp - R
!
  implicit none
  real(dp), intent(in) :: Tl
  real(dp) :: cv

  cv = compute_cp( Tl ) - Rair

end function compute_cv


!***********************************************************************
!
function gamma_of_T( Tl ) result(g)
!
! Description:
!   Ratio of specific heats, cp/cv.
!
  implicit none
  real(dp), intent(in) :: Tl
  real(dp) :: g, cpl

  cpl = compute_cp( Tl )
  g = cpl/( cpl - Rair )

end function gamma_of_T


!***********************************************************************
!
function speed_of_sound_sq( Tl ) result(c2)
!
! Description:
!   c^2 = gamma R T  (ideal gas)
!
  implicit none
  real(dp), intent(in) :: Tl
  real(dp) :: c2

  c2 = gamma_of_T( Tl )*Rair*Tl

end function speed_of_sound_sq


!***********************************************************************
!
function compute_enthalpy( Tl ) result(h)
!
! Description:
!   Computes enthalpy from temperature, h = integral(cp(T')dT', T_ref, T)  [J/kg]
!
  implicit none
  real(dp), intent(in) :: Tl
  real(dp) :: h
  real(dp) :: ts, t0

  select case ( cp_model_id() )

  case ( CP_POLYNOMIAL )
    t0 = cp_Tref
    h = cp_poly(0)*( Tl - t0 ) + cp_poly(1)*( Tl**2 - t0**2 )/2.0_dp + &
        cp_poly(2)*( Tl**3 - t0**3 )/3.0_dp + cp_poly(3)*( Tl**4 - t0**4 )/4.0_dp

  case ( CP_NASA )
    ! Integrated NASA polynomial (a6 is the enthalpy constant)
    if ( Tl < 1000.0_dp ) then
      h = nasa_R*( nasa_lo(1)*Tl + nasa_lo(2)*Tl**2/2.0_dp + nasa_lo(3)*Tl**3/3.0_dp + &
                   nasa_lo(4)*Tl**4/4.0_dp + nasa_lo(5)*Tl**5/5.0_dp + nasa_lo(6) )
    else
      h = nasa_R*( nasa_hi(1)*Tl + nasa_hi(2)*Tl**2/2.0_dp + nasa_hi(3)*Tl**3/3.0_dp + &
                   nasa_hi(4)*Tl**4/4.0_dp + nasa_hi(5)*Tl**5/5.0_dp + nasa_hi(6) )
    endif

  case ( CP_SHOMATE )
    ! h - h(298.15) = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F   [kJ/mol] (cp in J/molK, t = T/1000)
    ts = Tl/1000.0_dp
    h = ( shomate(1)*ts + shomate(2)*ts**2/2.0_dp + shomate(3)*ts**3/3.0_dp + shomate(4)*ts**4/4.0_dp &
        - shomate(5)/ts + shomate(6) )*1000.0_dp/shomate_M

  case default
    h = cp_const*( Tl - cp_Tref )

  end select

end function compute_enthalpy


!***********************************************************************
!
function compute_internal_energy( Tl ) result(e)
!
! Description:
!   Computes internal energy from temperature, e = h - R T = integral(cv(T')dT', T_ref, T)
!
  implicit none
  real(dp), intent(in) :: Tl
  real(dp) :: e

  e = compute_enthalpy( Tl ) - Rair*Tl

end function compute_internal_energy


!***********************************************************************
!
subroutine update_thermodynamic_properties
!
! Description:
!   Updates all thermodynamic and transport properties based on temperature, cells and boundary faces:
!
!   - cp(T), cv(T): specific heats (c_x = cp for the enthalpy forms, cv for the energy forms)
!   - mu(T): dynamic viscosity (Sutherland's law, optional; laminar: vis = mu)
!   - k(T) = cp mu / Pr: thermal conductivity
!   - den: from the equation of state (rho = p/(R*T)) for compressible flow
!
  implicit none

  integer :: i
  real(dp) :: cp_val, mu_val
  logical :: enth

  enth = ( energy_formulation == TOTAL_ENTHALPY .or. energy_formulation == ENTHALPY )

  if ( sutherland ) then
    vis = viscos*( T/mu_Tref )**1.5_dp*( mu_Tref + mu_S )/( T + mu_S )
  endif

  do i=1,numTotal
    cp_val = compute_cp( T(i) )
    cpf(i) = cp_val
    if ( enth ) then
      cxf(i) = cp_val
    else
      cxf(i) = cp_val - Rair
    endif
    ! Thermal conductivity: Prandtl number based
    mu_val = vis(i)
    kth(i) = cp_val*mu_val/pranl
  enddo

  if ( compressible ) den = p/( Rair*T + small )

end subroutine update_thermodynamic_properties


!***********************************************************************
!
subroutine set_convection_scheme( scheme, limiter )
!
! Description:
!   Sets the convection scheme and limiter for the energy equation (overrides the namelist strings).
!
  implicit none

  integer, intent(in) :: scheme
  integer, intent(in), optional :: limiter

  ! the input strings (read at initialisation) must not overwrite an explicit call
  if ( .not. energy_ready ) call energy_init
  convection_scheme = scheme
  if ( present(limiter) ) limiter_type = limiter

  write(*,'(a)') '========================================='
  write(*,'(a)') 'Energy Equation Convection Settings:'
  write(*,'(a)', advance='no') '  Scheme: '

  select case ( scheme )
  case ( UDS )
    write(*,*) 'First-Order Upwind (UDS)'
  case ( CDS )
    write(*,*) 'Central Differencing (CDS)'
  case ( SOU )
    write(*,*) 'Second-Order Upwind (SOU)'
  case ( QUICK )
    write(*,*) 'QUICK'
  case ( MUSCL )
    write(*,*) 'MUSCL'
    write(*,'(a)', advance='no') '  Limiter: '
    select case ( limiter_type )
    case ( LIMITER_NONE )
      write(*,*) 'None'
    case ( LIMITER_VANLEER )
      write(*,*) 'Van Leer'
    case ( LIMITER_MINMOD )
      write(*,*) 'Minmod'
    case ( LIMITER_SUPERBEE )
      write(*,*) 'Superbee'
    case ( LIMITER_VANALBADA )
      write(*,*) 'Van Albada'
    end select
  case ( LUDS )
    write(*,*) 'Linear Upwind (LUDS)'
  end select
  write(*,'(a)') '========================================='

end subroutine set_convection_scheme


!***********************************************************************
!
function compute_face_value_ho( ijp, ijn, iface, mass_flux, scheme, limiter ) result(phi_f)
!
! Description:
!   Computes higher-order face value of the convected quantity Cv (gradient dCvdxi).
!   Implements: CDS, LUDS, SOU, QUICK, MUSCL with limiters.
!   C = upwind cell, D = downwind cell, U = upwind-upwind cell (an existing cell if there is one aligned
!   with the C-D line, otherwise the virtual node phi_U = phi_D - 2 grad(phi_C).d_CD).
!
  implicit none

  integer, intent(in) :: ijp, ijn, iface, scheme, limiter
  real(dp), intent(in) :: mass_flux
  real(dp) :: phi_f

  integer :: ic, id, iju
  real(dp) :: phi_c, phi_d, phi_u
  real(dp) :: dcd(3), gdc, r_factor, limiter_val

  if ( mass_flux >= 0.0_dp ) then
    ic = ijp; id = ijn
  else
    ic = ijn; id = ijp
  endif

  phi_c = Cv(ic)
  phi_d = Cv(id)

  select case ( scheme )

  case ( UDS )
    !-----------------------------------------------------------------------
    ! First-order upwind
    !-----------------------------------------------------------------------
    phi_f = phi_c

  case ( CDS )
    !-----------------------------------------------------------------------
    ! Central differencing (second-order)
    !-----------------------------------------------------------------------
    phi_f = Cv(ijp) + ( Cv(ijn) - Cv(ijp) )*facint(iface)

  case ( LUDS )
    !-----------------------------------------------------------------------
    ! Linear Upwind Differencing Scheme
    !-----------------------------------------------------------------------
    phi_f = phi_c + dCvdxi(1,ic)*( xf(iface)-xc(ic) ) + dCvdxi(2,ic)*( yf(iface)-yc(ic) ) + dCvdxi(3,ic)*( zf(iface)-zc(ic) )

  case default
    !-----------------------------------------------------------------------
    ! SOU, QUICK, MUSCL: need the upwind-upwind value
    !-----------------------------------------------------------------------
    dcd(1) = xc(id)-xc(ic); dcd(2) = yc(id)-yc(ic); dcd(3) = zc(id)-zc(ic)
    gdc = dCvdxi(1,ic)*dcd(1) + dCvdxi(2,ic)*dcd(2) + dCvdxi(3,ic)*dcd(3)

    iju = find_upwind_cell( ic, id )
    if ( iju > 0 ) then
      phi_u = Cv(iju)
    else
      phi_u = phi_d - 2.0_dp*gdc              ! virtual upwind node
    endif

    select case ( scheme )

    case ( SOU )
      ! Second-order upwind (Warming-Beam): phi_f = phi_C + (phi_C - phi_U)/2
      phi_f = phi_c + 0.5_dp*( phi_c - phi_u )

    case ( QUICK )
      ! QUICK: phi_f = (6 phi_C + 3 phi_D - phi_U)/8
      phi_f = 0.75_dp*phi_c + 0.375_dp*phi_d - 0.125_dp*phi_u

    case ( MUSCL )
      ! MUSCL with limiter, kappa = 1/3 (third-order)
      r_factor = compute_gradient_ratio_unstructured( phi_c, phi_d, gdc )
      limiter_val = apply_limiter( r_factor, limiter )
      phi_f = phi_c + 0.25_dp*limiter_val*( ( 1.0_dp - 1.0_dp/3.0_dp )*( phi_c - phi_u ) + &
                                            ( 1.0_dp + 1.0_dp/3.0_dp )*( phi_d - phi_c ) )

    case default
      phi_f = phi_c

    end select

  end select

end function compute_face_value_ho


!***********************************************************************
!
function compute_gradient_ratio_unstructured( phi_c, phi_d, gdc ) result(r)
!
! Description:
!   Gradient ratio for the limited schemes on unstructured meshes, Jasak (1996) / Darwish & Moukalled (2003):
!
!      r = (phi_C - phi_U)/(phi_D - phi_C) = 2 grad(phi)_C . d_CD / (phi_D - phi_C) - 1
!
!   where d_CD is the vector from the upwind cell centre C to the downwind cell centre D.
!
  implicit none

  real(dp), intent(in) :: phi_c, phi_d, gdc
  real(dp) :: r, den_r

  den_r = phi_d - phi_c

  if ( abs(den_r) > 1.0e-12_dp*( abs(phi_c) + abs(phi_d) ) + small ) then
    r = 2.0_dp*gdc/den_r - 1.0_dp
  else
    r = 1.0_dp          ! uniform field
  endif

  ! Clip extreme values that can occur on highly skewed meshes
  r = max( -10.0_dp, min( 10.0_dp, r ) )

end function compute_gradient_ratio_unstructured


!***********************************************************************
!
function apply_limiter( r, limiter ) result(phi_r)
!
! Description:
!   Applies flux limiter function
!
  implicit none

  real(dp), intent(in) :: r
  integer, intent(in) :: limiter
  real(dp) :: phi_r

  select case ( limiter )

  case ( LIMITER_NONE )
    phi_r = 1.0_dp

  case ( LIMITER_VANLEER )
    ! Van Leer: phi(r) = (r + |r|) / (1 + |r|)
    phi_r = ( r + abs(r) )/( 1.0_dp + abs(r) )

  case ( LIMITER_MINMOD )
    ! Minmod: phi(r) = max(0, min(1, r))
    phi_r = max( 0.0_dp, min( 1.0_dp, r ) )

  case ( LIMITER_SUPERBEE )
    ! Superbee: phi(r) = max(0, min(2r, 1), min(r, 2))
    phi_r = max( 0.0_dp, min( 2.0_dp*r, 1.0_dp ), min( r, 2.0_dp ) )

  case ( LIMITER_VANALBADA )
    ! Van Albada: phi(r) = (r + r^2) / (1 + r^2), zero for r < 0
    phi_r = max( 0.0_dp, ( r + r*r )/( 1.0_dp + r*r ) )

  case default
    phi_r = 1.0_dp

  end select

  ! Ensure TVD condition: 0 <= phi(r) <= 2
  phi_r = max( 0.0_dp, min( 2.0_dp, phi_r ) )

end function apply_limiter


!***********************************************************************
!
function find_upwind_cell( ic, id ) result(iju)
!
! Description:
!   Finds the upwind-upwind cell for higher-order schemes: among the neighbours of the upwind cell C the one that is
!   best aligned with the line D -> C (cosine larger than 0.7), 0 if there is none (boundary cell, unstructured mesh).
!
  implicit none

  integer, intent(in) :: ic, id
  integer :: iju

  integer :: k, j
  real(dp) :: dcd(3), dcu(3), ncd, ncu, cosa, cosmax

  iju = 0
  cosmax = 0.7_dp

  dcd(1) = xc(ic)-xc(id); dcd(2) = yc(ic)-yc(id); dcd(3) = zc(ic)-zc(id)    ! direction D -> C
  ncd = sqrt( dcd(1)**2 + dcd(2)**2 + dcd(3)**2 )

  do k = ia(ic), ia(ic+1)-1
    j = ja(k)
    if ( j == ic .or. j == id ) cycle
    dcu(1) = xc(j)-xc(ic); dcu(2) = yc(j)-yc(ic); dcu(3) = zc(j)-zc(ic)    ! direction C -> U, wanted: -(D -> C)
    ncu = sqrt( dcu(1)**2 + dcu(2)**2 + dcu(3)**2 )
    cosa = -( dcd(1)*dcu(1) + dcd(2)*dcu(2) + dcd(3)*dcu(3) )/( ncd*ncu + small )
    if ( cosa > cosmax ) then
      cosmax = cosa
      iju = j
    endif
  enddo

end function find_upwind_cell


end module energy_equation
