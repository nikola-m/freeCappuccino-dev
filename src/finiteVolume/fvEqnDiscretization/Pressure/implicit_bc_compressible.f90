module implicit_bc_compressible
!
! Purpose:
!   Implicit boundary conditions of the all-speeds (compressible) block-coupled pressure-velocity solver.
!
! Description:
!   Implements implicit treatment of the boundary faces of the compressible coupled solver, following the organisation
!   of the reference file implicit_boundary_conditions_compressible.f90 (Darwish, Manganai & Moukalled (2018)
!   "Implicit boundary conditions for coupled solvers", Computers and Fluids 168, pp. 54-66):
!   one routine per boundary condition type, each linearises the boundary mass flux in the cell pressure p_P,
!   the cell velocity U_P (when the boundary value is extrapolated) and uses the same 4x4 block of the owner cell.
!
!   BOUNDARY CONDITION TYPES (selected per patch by the namelist, see read_input.f90):
!
!     'compinlet'     compressible_inlet_bc        patch type 'inlet'     U_b, T_b prescribed, p_b = p_P.
!     'massflow'      mass_flow_inlet_bc           patch type 'inlet'     mdot prescribed, U_b follows, T_b prescribed.
!     'totalpressure' total_pressure_inlet_bc      patch type 'pressure'  p_0 prescribed, T_b prescribed.
!     'stagnation'    stagnation_bc                patch type 'pressure'  p_0 and T_0 prescribed (isentropic relations).
!     'riemann'       riemann_invariant_bc         patch type 'pressure'  characteristic far field (p, T, U of 0/ files).
!     'compoutlet'    compressible_outlet_bc       patch type 'pressure'  subsonic: static p prescribed; supersonic:
!                                                                        everything extrapolated.
!     'open'          open_boundary_bc             patch type 'pressure'  switches between outlet (static p) and
!                                                                        total pressure inlet (entrainment) per face.
!
!   Patches without an entry keep the default behaviour of coupled_pU ('inlet', 'outlet', 'pressure').
!
!   The patch types ('inlet', 'pressure') are the ones the rest of the code (momentum, energy, boundary updates)
!   understands, the type of the namelist selects only the way the boundary state (p_b, T_b, rho_b, U_b) and the
!   mass flux are computed. Pressure-like conditions (total pressure, stagnation, Riemann, outlet, open) share one
!   assembly: the boundary pressure p_b is evaluated from the current cell state before each coupled iteration (the
!   boundary condition relation), after which the boundary face is treated as a fixed pressure face with the
!   Rhie-Chow flux between the cell centre and the face, exactly as the 'pressure' patch of coupled_pU:
!
!        mdot_b = rho_b U_P.S - rho_b (Vol/ap)_P Df_b [ p_b - p_P - grad p_P . d_Pb ]        (linear in U_P and p_P)
!
!   so that the velocity-pressure coupling at the boundary is implicit, and the boundary relation (p_0 - dynamic
!   pressure, Riemann invariants, ...) is a lagged Picard update (under-relaxation: cplbc_urf).
!   Supersonic outflow faces ('compoutlet', 'open', 'riemann') switch to the extrapolated form
!
!        mdot_b = rho_P U_P.S   linearised in (U_P, p_P),   p_b = p_P
!
!   Block structure at the boundary (row of the owner cell, diagonal block):
!   +-----------------+------+
!   | A_uu  A_uv  A_uw| G_px |  <- momentum: boundary pressure force p_b S (implicit p_P or known p_b)
!   +-----------------+------+
!   | D_u   D_v   D_w | D_pp |  <- continuity: boundary mass flux
!   +-----------------+------+
!
  use types
  use parameters
  use geometry
  use sparse_matrix, only: apu
  use variables
  use energy_equation, only: gamma_of_T, compute_cp, speed_of_sound_sq

  implicit none

  private

  ! Kinds of boundary conditions
  integer, parameter, public :: BCK_DEFAULT       = 0
  integer, parameter, public :: BCK_COMP_INLET    = 1
  integer, parameter, public :: BCK_MASSFLOW      = 2
  integer, parameter, public :: BCK_TOTAL_PRESSURE = 3
  integer, parameter, public :: BCK_STAGNATION    = 4
  integer, parameter, public :: BCK_RIEMANN       = 5
  integer, parameter, public :: BCK_COMP_OUTLET   = 6
  integer, parameter, public :: BCK_OPEN          = 7

  ! Face modes (how the mass flux of a boundary face is assembled)
  integer, parameter, public :: MODE_NONE     = 0   ! default treatment of coupled_pU
  integer, parameter, public :: MODE_PRESSURE = 1   ! p_b known, Rhie-Chow flux
  integer, parameter, public :: MODE_EXTRAP   = 2   ! zero gradient, mdot = rho_P U_P.S linearised
  integer, parameter, public :: MODE_FIXFLUX  = 3   ! prescribed mass flux
  integer, parameter, public :: MODE_VELINLET = 4   ! prescribed U_b, T_b, p_b = p_P
  integer, parameter, public :: MODE_RIEMANN  = 5   ! subsonic characteristic boundary: p_b, u_nb linear in (U_P, p_P)

  integer, parameter :: maxbc = 20

  ! Namelist input (see read_input.f90): per patch name, type and values
  character(len=15), save :: cplbc_name(maxbc) = ''
  character(len=24), save :: cplbc_kind(maxbc) = ''
  real(dp), save :: cplbc_p0(maxbc)   = 0.0_dp      ! total (stagnation) pressure, or static pressure for the outlet kinds [Pa]
  real(dp), save :: cplbc_T0(maxbc)   = 0.0_dp      ! total (stagnation) temperature [K]
  real(dp), save :: cplbc_mdot(maxbc) = 0.0_dp      ! mass flow rate through the patch [kg/s] (massflow), positive = inflow
  real(dp), save :: cplbc_dirx(maxbc) = 0.0_dp      ! inflow direction (default: opposite to the patch normal)
  real(dp), save :: cplbc_diry(maxbc) = 0.0_dp
  real(dp), save :: cplbc_dirz(maxbc) = 0.0_dp
  real(dp), save :: cplbc_urf = 0.3_dp              ! under-relaxation of the boundary pressure update

  integer, allocatable, save :: bckind(:)           ! kind of boundary condition per patch
  integer, allocatable, save :: bcmode(:)           ! mode per boundary face (index iface - numInnerFaces)
  real(dp), allocatable, save :: p_set(:), T_set(:), U_set(:,:)   ! values of the 0/ files per boundary face
  real(dp), allocatable, save :: edir(:,:)          ! inflow direction (unit vector) per boundary face
  real(dp), allocatable, save :: mdot_face(:)       ! prescribed mass flow through a face (massflow kind), >= 0 = inflow
  real(dp), allocatable, save :: pbc0(:), Tbc0(:)   ! p_0, T_0 per patch
  logical, save :: bc_active = .false.

  public :: bckind, bcmode, bc_active
  public :: cplbc_name, cplbc_kind, cplbc_p0, cplbc_T0, cplbc_mdot, cplbc_dirx, cplbc_diry, cplbc_dirz, cplbc_urf
  public :: bc_setup, bc_prepare_face, bc_apply_face, bc_post_solve_face, bc_pressure_implicit
  public :: bc_is_riemann, bc_riemann_momentum
  public :: apply_compressible_inlet_bc, apply_compressible_outlet_bc, apply_total_pressure_inlet_bc, &
            apply_mass_flow_inlet_bc, apply_stagnation_bc, apply_riemann_invariant_bc, apply_open_boundary_bc

contains


!***********************************************************************
!
subroutine bc_setup
!
! Description:
!   Matches the patch names of the namelist with the patches of the mesh, checks the patch type, stores the values of
!   the initial fields at the boundary (free stream state of the characteristic conditions), the inflow directions
!   and the prescribed mass flow per face.
!
  implicit none

  integer :: ib, ik, i, iface, ijb, ib2, nbf, nentries
  real(dp) :: are, area, dnorm, nx, ny, nz, mdot_pt
  character(len=24) :: kn

  nbf = numBoundaryFaces
  allocate( bckind(numBoundaries), bcmode(nbf), pbc0(numBoundaries), Tbc0(numBoundaries) )
  allocate( p_set(nbf), T_set(nbf), U_set(3,nbf), edir(3,nbf), mdot_face(nbf) )
  bckind = BCK_DEFAULT
  bcmode = MODE_NONE
  pbc0 = 0.0_dp; Tbc0 = 0.0_dp
  p_set = 0.0_dp; T_set = 0.0_dp; U_set = 0.0_dp; edir = 0.0_dp; mdot_face = 0.0_dp

  nentries = 0
  do ik=1,maxbc
    if ( len_trim( cplbc_name(ik) ) > 0 ) nentries = nentries + 1
  enddo
  if ( nentries > 0 .and. .not. compressible ) then
    write(*,'(a)') 'Fatal error: cplbc_name/cplbc_kind boundary conditions need compressible = T.'
    stop
  endif

  do ik=1,maxbc

    if ( len_trim( cplbc_name(ik) ) == 0 ) cycle

    ib = 0
    do i=1,numBoundaries
      if ( trim(adjustl(bcname(i))) == trim(adjustl(cplbc_name(ik))) ) ib = i
    enddo

    if ( ib == 0 ) then
      write(*,'(a,a)') 'Fatal error: compressible boundary condition for a non-existing patch: ', trim(cplbc_name(ik))
      stop
    endif

    kn = normalise_kind( cplbc_kind(ik) )

    select case ( trim(kn) )
    case ('compinlet','compressibleinlet')
      bckind(ib) = BCK_COMP_INLET
    case ('massflow','massflowinlet')
      bckind(ib) = BCK_MASSFLOW
    case ('totalpressure','totalpressureinlet')
      bckind(ib) = BCK_TOTAL_PRESSURE
    case ('stagnation')
      bckind(ib) = BCK_STAGNATION
    case ('riemann','riemanninvariant')
      bckind(ib) = BCK_RIEMANN
    case ('compoutlet','compressibleoutlet')
      bckind(ib) = BCK_COMP_OUTLET
    case ('open','openboundary')
      bckind(ib) = BCK_OPEN
    case default
      write(*,'(a,a,a,a)') 'Fatal error: unknown compressible boundary condition kind ', trim(cplbc_kind(ik)), &
                           ' for patch ', trim(cplbc_name(ik))
      stop
    end select

    ! Patch type check
    select case ( bckind(ib) )
    case ( BCK_COMP_INLET, BCK_MASSFLOW )
      if ( bctype(ib) /= 'inlet' ) then
        write(*,'(a,a,a)') 'Fatal error: patch ', trim(cplbc_name(ik)), ' must be of type inlet for this boundary condition.'
        stop
      endif
    case default
      if ( bctype(ib) /= 'pressure' ) then
        write(*,'(a,a,a)') 'Fatal error: patch ', trim(cplbc_name(ik)), ' must be of type pressure for this boundary condition.'
        stop
      endif
    end select

    pbc0(ib) = cplbc_p0(ik)
    Tbc0(ib) = cplbc_T0(ik)

    if ( bckind(ib) == BCK_TOTAL_PRESSURE .or. bckind(ib) == BCK_STAGNATION ) then
      if ( pbc0(ib) <= 0.0_dp ) then
        write(*,'(a,a)') 'Fatal error: cplbc_p0 must be positive for patch ', trim(cplbc_name(ik))
        stop
      endif
    endif
    if ( bckind(ib) == BCK_STAGNATION .and. Tbc0(ib) <= 0.0_dp ) then
      write(*,'(a,a)') 'Fatal error: cplbc_T0 must be positive for patch ', trim(cplbc_name(ik))
      stop
    endif

    ! Area of the patch (mass flow distribution)
    area = 0.0_dp
    do i=1,nfaces(ib)
      iface = startFace(ib)+i
      area = area + sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
    enddo
    mdot_pt = cplbc_mdot(ik)

    dnorm = sqrt( cplbc_dirx(ik)**2 + cplbc_diry(ik)**2 + cplbc_dirz(ik)**2 )

    do i=1,nfaces(ib)

      iface = startFace(ib)+i
      ijb   = iBndValueStart(ib)+i
      ib2   = iface - numInnerFaces

      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
      nx = arx(iface)/are; ny = ary(iface)/are; nz = arz(iface)/are

      ! Values of the initial fields at the boundary = free stream / set values
      p_set(ib2) = p(ijb)
      T_set(ib2) = T(ijb)
      U_set(1,ib2) = u(ijb); U_set(2,ib2) = v(ijb); U_set(3,ib2) = w(ijb)

      ! The static pressure kinds take it from cplbc_p0 if given
      if ( ( bckind(ib) == BCK_COMP_OUTLET .or. bckind(ib) == BCK_OPEN ) .and. pbc0(ib) > 0.0_dp ) p_set(ib2) = pbc0(ib)

      ! Inflow direction: user vector or opposite to the outward normal
      if ( dnorm > 0.0_dp ) then
        edir(1,ib2) = cplbc_dirx(ik)/dnorm; edir(2,ib2) = cplbc_diry(ik)/dnorm; edir(3,ib2) = cplbc_dirz(ik)/dnorm
        if ( edir(1,ib2)*nx + edir(2,ib2)*ny + edir(3,ib2)*nz >= -0.05_dp ) then
          write(*,'(a,a)') 'Fatal error: inflow direction is not directed into the domain on patch ', trim(cplbc_name(ik))
          stop
        endif
      else
        edir(1,ib2) = -nx; edir(2,ib2) = -ny; edir(3,ib2) = -nz
      endif

      ! Mass flow per face proportional to its area
      if ( bckind(ib) == BCK_MASSFLOW ) mdot_face(ib2) = mdot_pt*are/area

    enddo

    bc_active = .true.

    write(*,'(a,a,a,a,a,es10.3,a,es10.3,a,es10.3)') '  Compressible BC: patch ', trim(cplbc_name(ik)), '  kind ', trim(kn), &
          '   p0 = ', cplbc_p0(ik), '  T0 = ', cplbc_T0(ik), '  mdot = ', cplbc_mdot(ik)

  enddo

  ! Initial boundary state consistent with the kinds
  call bc_initial_state

end subroutine bc_setup


function normalise_kind( s ) result(r)
!
!  Lower case and without underscores / blanks:  Total_Pressure_Inlet -> totalpressureinlet
!
  implicit none
  character(len=*), intent(in) :: s
  character(len=24) :: r
  integer :: i, ic, n
  r = ''
  n = 0
  do i=1,len_trim(s)
    ic = iachar( s(i:i) )
    if ( s(i:i) == '_' .or. s(i:i) == ' ' .or. s(i:i) == '-' ) cycle
    if ( ic >= 65 .and. ic <= 90 ) ic = ic + 32
    n = n + 1
    if ( n > len(r) ) exit
    r(n:n) = achar(ic)
  enddo
end function normalise_kind


!***********************************************************************
!
subroutine bc_initial_state
!
! Description:
!   Initial boundary velocity and mass flux of the mass flow inlets (velocity from the prescribed mass flow).
!
  implicit none
  integer :: ib, i, iface, ijb, ib2
  real(dp) :: ain, rb

  do ib=1,numBoundaries
    if ( bckind(ib) /= BCK_MASSFLOW ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib)+i
      ijb   = iBndValueStart(ib)+i
      ib2   = iface - numInnerFaces
      ain = -( edir(1,ib2)*arx(iface) + edir(2,ib2)*ary(iface) + edir(3,ib2)*arz(iface) )
      rb  = p(ijb)/( Rair*T(ijb)+small )
      u(ijb) = edir(1,ib2)*mdot_face(ib2)/( rb*ain )
      v(ijb) = edir(2,ib2)*mdot_face(ib2)/( rb*ain )
      w(ijb) = edir(3,ib2)*mdot_face(ib2)/( rb*ain )
      flmass(iface) = -mdot_face(ib2)
    enddo
  enddo

end subroutine bc_initial_state


!***********************************************************************
!
function bc_pressure_implicit( iface ) result(impl)
!
! Description:
!   Is the boundary pressure in the momentum equation implicit (p_b = p_P, extrapolating faces) or known?
!
  implicit none
  integer, intent(in) :: iface
  logical :: impl
  impl = ( bcmode(iface-numInnerFaces) == MODE_EXTRAP )
end function bc_pressure_implicit


!***********************************************************************
!
subroutine bc_prepare_face( ib, i, iface, ijp, ijb )
!
! Description:
!   Boundary state (p_b, T_b, rho_b) and face mode from the boundary condition relation and the current cell state
!   (Picard lagging). Called before the face is assembled.
!
  implicit none

  integer, intent(in) :: ib, i, iface, ijp, ijb

  integer :: ib2
  real(dp) :: are, nx, ny, nz, un, q, pb, Tb, cP, gam, cpl, M, rhoP
  real(dp) :: uP, vP, wP
  real(dp) :: p0, T0
  real(dp) :: R_out, R_in, unb, cb, Mb, gm1, Sent, rho_b, pinf, Tinf, rinf, cinf, uninf, unp
  logical :: inflow

  ib2 = iface - numInnerFaces

  if ( bckind(ib) == BCK_DEFAULT ) return
  if ( bckind(ib) == BCK_COMP_INLET ) then
    ! U_b, T_b prescribed, p_b = p_P; the density of the face follows the cell pressure with the Picard relaxation
    bcmode(ib2) = MODE_VELINLET
    p(ijb) = p(ijb) + cplbc_urf*( p(ijp) - p(ijb) )
    den(ijb) = p(ijb)/( Rair*T(ijb) + small )
    return
  endif
  if ( bckind(ib) == BCK_MASSFLOW ) then
    bcmode(ib2) = MODE_FIXFLUX
    return
  endif

  are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
  nx = arx(iface)/are; ny = ary(iface)/are; nz = arz(iface)/are

  uP = u(ijp); vP = v(ijp); wP = w(ijp)
  rhoP = den(ijp)
  gam = gamma_of_T( T(ijp) )
  cP  = sqrt( speed_of_sound_sq( T(ijp) ) )
  un  = uP*nx + vP*ny + wP*nz               ! outward normal velocity
  inflow = ( un < 0.0_dp )

  ! speed of the inflow along the inflow direction
  ! total-pressure type kinds use the speed of the cell (as OpenFOAM's totalPressure): with the normal component alone
  ! the flow tangential to a far-field patch would be assigned p_b = p_0
  q = sqrt( uP*uP + vP*vP + wP*wP )

  p0 = pbc0(ib); T0 = Tbc0(ib)

  pb = p(ijp); Tb = T(ijp)
  bcmode(ib2) = MODE_PRESSURE

  select case ( bckind(ib) )

  case ( BCK_TOTAL_PRESSURE )
    !-----------------------------------------------------------------
    ! total pressure inlet: p_b = p_0 (1+(g-1)/2 M^2)^(-g/(g-1)) (inflow; M from the inflow speed and T_b), p_b = p_0 (outflow); T_b prescribed
    !-----------------------------------------------------------------
    if ( inflow ) then
      Tb = T_set(ib2)
      M  = q/sqrt( speed_of_sound_sq( Tb ) )
      pb = p0/( 1.0_dp + 0.5_dp*(gam-1.0_dp)*M*M )**( gam/(gam-1.0_dp) )      ! isentropic: p_0 = p (1+(g-1)/2 M^2)^(g/(g-1))
    else
      pb = p0
      Tb = T(ijp)
    endif
    pb = max( pb, 0.05_dp*p0 )

  case ( BCK_STAGNATION )
    !-----------------------------------------------------------------
    ! stagnation inlet:  T_b = T_0 - U^2/(2 cp),  p_b = p_0 (T_b/T_0)^(gamma/(gamma-1))
    !-----------------------------------------------------------------
    if ( inflow ) then
      cpl = compute_cp( T(ijp) )
      Tb = max( T0 - 0.5_dp*q*q/cpl, 0.2_dp*T0 )
      pb = p0*( Tb/T0 )**( gam/(gam-1.0_dp) )
    else
      pb = p0                                          ! backflow: p_0 acts as static pressure
      Tb = T(ijp)
    endif

  case ( BCK_COMP_OUTLET )
    !-----------------------------------------------------------------
    ! compressible outlet: subsonic - static pressure prescribed; supersonic - extrapolated
    !-----------------------------------------------------------------
    M = abs(un)/cP
    if ( .not. inflow .and. M >= 1.0_dp ) then
      bcmode(ib2) = MODE_EXTRAP
      pb = p(ijp); Tb = T(ijp)
    else
      pb = p_set(ib2)
      if ( inflow ) then
        Tb = T_set(ib2)
      else
        Tb = T(ijp)
      endif
    endif

  case ( BCK_OPEN )
    !-----------------------------------------------------------------
    ! open boundary: outflow - outlet (static p), inflow - total pressure inlet with entrainment velocity
    !-----------------------------------------------------------------
    M = abs(un)/cP
    if ( .not. inflow .and. M >= 1.0_dp ) then
      bcmode(ib2) = MODE_EXTRAP
      pb = p(ijp); Tb = T(ijp)
    elseif ( inflow ) then
      pb = max( p_set(ib2) - 0.5_dp*rhoP*un*un, 0.05_dp*p_set(ib2) )
      Tb = T_set(ib2)
    else
      pb = p_set(ib2)
      Tb = T(ijp)
    endif

  case ( BCK_RIEMANN )
    !-----------------------------------------------------------------
    ! characteristic far field: R+ = u_n + 2c/(gamma-1) from the interior (outgoing),
    ! R- = u_n - 2c/(gamma-1) from the free stream (incoming), entropy and tangential velocity from the upwind side.
    !-----------------------------------------------------------------
    gm1 = gam - 1.0_dp
    pinf = p_set(ib2); Tinf = T_set(ib2)
    rinf = pinf/( Rair*Tinf )
    cinf = sqrt( gam*Rair*Tinf )
    uninf = U_set(1,ib2)*nx + U_set(2,ib2)*ny + U_set(3,ib2)*nz
    unp = un

    R_out = unp + 2.0_dp*cP/gm1
    R_in  = uninf - 2.0_dp*cinf/gm1
    unb = 0.5_dp*( R_out + R_in )
    cb  = max( 0.25_dp*gm1*( R_out - R_in ), 0.1_dp*cinf )
    Mb  = unb/cb

    if ( Mb <= -1.0_dp ) then
      pb = pinf; Tb = Tinf                           ! supersonic inflow: all from the free stream
    elseif ( Mb >= 1.0_dp ) then
      bcmode(ib2) = MODE_EXTRAP                      ! supersonic outflow: all from the interior
      pb = p(ijp); Tb = T(ijp)
    else
      if ( unb < 0.0_dp ) then
        Sent = pinf/rinf**gam                        ! entropy of the free stream
      else
        Sent = p(ijp)/rhoP**gam                      ! entropy of the interior
      endif
      rho_b = ( cb*cb/( gam*Sent ) )**( 1.0_dp/gm1 )
      pb = rho_b*cb*cb/gam
      Tb = cb*cb/( gam*Rair )
      bcmode(ib2) = MODE_RIEMANN                     ! subsonic: implicit characteristic relations (see riemann_*)
    endif

  end select

  ! Under-relaxation of the boundary pressure, boundary temperature and density
  pb = p(ijb) + cplbc_urf*( pb - p(ijb) )
  p(ijb) = pb
  T(ijb) = Tb
  den(ijb) = pb/( Rair*Tb + small )

end subroutine bc_prepare_face


!***********************************************************************
!
subroutine bc_apply_face( ib, i, iface, ijp, ijb, kd, bval, rhsb, Eb, Vsb )
!
! Description:
!   Dispatches the continuity row of the boundary face to the routine of its boundary condition type.
!
  implicit none

  integer, intent(in) :: ib, i, iface, ijp, ijb, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:), Vsb(:)

  integer :: ib2

  ib2 = iface - numInnerFaces

  select case ( bckind(ib) )

  case ( BCK_COMP_INLET )
    call apply_compressible_inlet_bc( ijb, iface, ijp, kd, bval, rhsb )

  case ( BCK_MASSFLOW )
    call apply_mass_flow_inlet_bc( ijb, iface, ijp, mdot_face(ib2), kd, bval, rhsb )

  case ( BCK_TOTAL_PRESSURE )
    call apply_total_pressure_inlet_bc( ijb, iface, ijp, pbc0(ib), kd, bval, rhsb, Eb )

  case ( BCK_STAGNATION )
    call apply_stagnation_bc( ijb, iface, ijp, pbc0(ib), Tbc0(ib), kd, bval, rhsb, Eb )

  case ( BCK_RIEMANN )
    call apply_riemann_invariant_bc( ijb, iface, ijp, kd, bval, rhsb, Eb, Vsb )

  case ( BCK_COMP_OUTLET )
    call apply_compressible_outlet_bc( ijb, iface, ijp, kd, bval, rhsb, Eb, Vsb )

  case ( BCK_OPEN )
    call apply_open_boundary_bc( ijb, iface, ijp, kd, bval, rhsb, Eb, Vsb )

  end select

end subroutine bc_apply_face


!***********************************************************************
!
subroutine pressure_face_assembly( iface, ijp, ijb, kd, bval, rhsb, Eb )
!
! Description:
!   Fixed (known) boundary pressure, Rhie-Chow mass flux between the cell centre and the face:
!
!      mdot_b = rho_b U_P.S - rho_b (Vol/ap)_P Df_b [ p_b - p_P - grad p_P . d_Pb ]
!
!   linear in U_P and p_P. Pressure force p_b S of the momentum equation is added by the caller (known p_b).
!
  implicit none

  integer, intent(in) :: iface, ijp, ijb, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:)

  real(dp) :: sx, sy, sz, dx, dy, dz, rb, are2, sdb, capb, extrap

  sx = arx(iface); sy = ary(iface); sz = arz(iface)
  dx = xf(iface)-xc(ijp); dy = yf(iface)-yc(ijp); dz = zf(iface)-zc(ijp)

  rb   = den(ijb)                              ! p_b/(R T_b)
  are2 = sx*sx+sy*sy+sz*sz
  sdb  = sx*dx+sy*dy+sz*dz
  capb = -rb*vol(ijp)*apu(ijp)*are2/sdb        ! < 0

  extrap = dPdxi(1,ijp)*dx + dPdxi(2,ijp)*dy + dPdxi(3,ijp)*dz
  Eb(iface-numInnerFaces) = capb*( p(ijb) - extrap )          ! known part of the flux, no p_P, no u_P

  bval(4,1,kd) = bval(4,1,kd) + rb*sx
  bval(4,2,kd) = bval(4,2,kd) + rb*sy
  bval(4,3,kd) = bval(4,3,kd) + rb*sz
  bval(4,4,kd) = bval(4,4,kd) - capb
  rhsb(4,ijp)  = rhsb(4,ijp) - Eb(iface-numInnerFaces)

end subroutine pressure_face_assembly


!***********************************************************************
!
subroutine extrapolated_face_assembly( iface, ijp, kd, bval, rhsb, Vsb )
!
! Description:
!   Zero gradient (supersonic outflow): mdot_b = rho_P U_P.S, rho_P = C p_P, linearised around the current iterate
!
!      mdot_b = rho* U.S + (U*.S) C ( p_P - p_P* ) + rho* (U - U*).S
!
!   Momentum: p_b = p_P (implicit, added by the caller).
!
  implicit none

  integer, intent(in) :: iface, ijp, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Vsb(:)

  real(dp) :: sx, sy, sz, Cc, vnb

  sx = arx(iface); sy = ary(iface); sz = arz(iface)

  Cc  = 1.0_dp/( Rair*T(ijp) + small )
  vnb = u(ijp)*sx + v(ijp)*sy + w(ijp)*sz
  Vsb(iface-numInnerFaces) = vnb

  bval(4,1,kd) = bval(4,1,kd) + den(ijp)*sx
  bval(4,2,kd) = bval(4,2,kd) + den(ijp)*sy
  bval(4,3,kd) = bval(4,3,kd) + den(ijp)*sz
  bval(4,4,kd) = bval(4,4,kd) + Cc*vnb
  rhsb(4,ijp)  = rhsb(4,ijp)  + vnb*den(ijp)

end subroutine extrapolated_face_assembly


!***********************************************************************
!
subroutine apply_compressible_inlet_bc( ijb, iface, ijp, kd, bval, rhsb )
!
! Description:
!   Implicit compressible inlet boundary condition
!
!   For compressible flows, density varies with pressure/temperature: rho = p/(R*T) (ideal gas), which creates an
!   additional coupling: mdot = rho(p) u.S. With the velocity U_b and temperature T_b prescribed (Dirichlet) and
!   the boundary pressure extrapolated (zero gradient, p_b = p_P):
!
!      mdot_b = C_b p_P (U_b.S),     C_b = 1/(R T_b)
!
!   is linear in the cell pressure - exact, no linearisation needed (dmdot/dp = C_b U_b.S, dmdot/dT neglected
!   because T_b is prescribed). Since the velocity is prescribed there is no Rhie-Chow term at the face.
!   Momentum: p_b = p_P (implicit).
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:)

  real(dp) :: Cb

  ! Density of the face from the lagged cell pressure: the implicit form (mdot = C p_P U.S) makes the inflow grow with
  ! the cell pressure, a destabilising feedback in the pseudo-time iteration (diverges in the nozzle test)
  Cb = den(ijb)
  rhsb(4,ijp) = rhsb(4,ijp) - Cb*( u(ijb)*arx(iface) + v(ijb)*ary(iface) + w(ijb)*arz(iface) )

end subroutine apply_compressible_inlet_bc


!***********************************************************************
!
subroutine apply_mass_flow_inlet_bc( ijb, iface, ijp, mdot_prescribed, kd, bval, rhsb )
!
! Description:
!   Implicit mass flow inlet boundary condition
!
!   Prescribed: mass flow rate mdot (inflow, per face: patch mass flow times area fraction)
!   Unknown: velocity magnitude and pressure
!
!   Approach:
!   1. Velocity direction known (normal to the inlet or the user direction)
!   2. Magnitude adjusted to match the mass flow:  |U_b| = mdot / (rho_b A_in),  rho_b = p_P/(R T_b),
!      updated after every coupled solve (bc_post_solve_face)
!   3. Pressure extrapolated: p_b = p_P (implicit in the momentum equation)
!
!   The mass flux through the face is known, it goes to the right hand side of the continuity equation:
!
!      mdot_b = - mdot_prescribed   (outflow positive)
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(in) :: mdot_prescribed
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:)

  rhsb(4,ijp) = rhsb(4,ijp) + mdot_prescribed

end subroutine apply_mass_flow_inlet_bc


!***********************************************************************
!
subroutine apply_total_pressure_inlet_bc( ijb, iface, ijp, p_total, kd, bval, rhsb, Eb )
!
! Description:
!   Implicit total pressure inlet (for compressible flows)
!
!   Prescribed: total pressure p_0 = p + 0.5 rho |U|^2, temperature T_b
!   Computed: static pressure and velocity
!
!   Boundary relation (bc_prepare_face):  p_b = p_0 - 0.5 rho_P q^2,  q = U_P . e_in  (inflow),  p_b = p_0 (outflow).
!   Implicit coupling of pressure and velocity at the boundary is the Rhie-Chow flux between the cell centre
!   and the face, linear in (U_P, p_P) (pressure_face_assembly): the flow adjusts itself to the prescribed
!   total pressure through the momentum balance of the first cell. The inflow velocity direction is e_in,
!   its magnitude follows from the mass flux (bc_post_solve_face).
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(in) :: p_total
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:)

  if ( bcmode(iface-numInnerFaces) == MODE_EXTRAP ) then
    return
  endif

  call pressure_face_assembly( iface, ijp, ijb, kd, bval, rhsb, Eb )

end subroutine apply_total_pressure_inlet_bc


!***********************************************************************
!
subroutine apply_stagnation_bc( ijb, iface, ijp, p0_bc, T0_bc, kd, bval, rhsb, Eb )
!
! Description:
!   Stagnation (total) pressure and temperature inlet
!
!   For compressible flows:
!   p_0 = p (1 + (gamma-1)/2 M^2)^(gamma/(gamma-1))
!   T_0 = T (1 + (gamma-1)/2 M^2)
!
!   Evaluated from the current cell state (bc_prepare_face):  T_b = T_0 - q^2/(2 cp),
!   p_b = p_0 (T_b/T_0)^(gamma/(gamma-1)),  rho_b = p_b/(R T_b). The dependence of p_b on the velocity is the
!   Newton-type Picard update of the boundary (lagged), the implicit part is the Rhie-Chow flux of the pressure face
!   (see apply_total_pressure_inlet_bc). T_b enters the energy equation as the inflow temperature: the
!   total enthalpy of the inflow is c_p T_0 for the converged solution.
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(in) :: p0_bc, T0_bc
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:)

  call pressure_face_assembly( iface, ijp, ijb, kd, bval, rhsb, Eb )

end subroutine apply_stagnation_bc


!***********************************************************************
!
subroutine apply_riemann_invariant_bc( ijb, iface, ijp, kd, bval, rhsb, Eb, Vsb )
!
! Description:
!   Characteristic-based boundary condition using Riemann invariants
!
!   For compressible flows, characteristics determine info propagation:
!   - Incoming characteristics: prescribed from the free stream (0/ files)
!   - Outgoing characteristics: extrapolated from the interior
!
!   Riemann invariants (1D normal to the boundary, n outward):
!   R+ = u_n + 2c/(gamma-1)  (outgoing)
!   R- = u_n - 2c/(gamma-1)  (incoming)
!   S  = p/rho^gamma         (entropy wave: free stream for inflow, interior for outflow)
!
!   u_n = (R+ + R-)/2, c = (gamma-1)(R+ - R-)/4, T_b = c^2/(gamma R), rho_b = (c^2/(gamma S))^(1/(gamma-1)),
!   p_b = rho_b c^2/gamma. Supersonic faces: all from the free stream (inflow) or from the interior (outflow).
!   Tangential velocity: free stream (inflow) or interior (outflow).
!
!   The boundary pressure is evaluated in bc_prepare_face; the continuity row is the Rhie-Chow pressure face (subsonic) or
!   the extrapolated flux (supersonic outflow).
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:), Vsb(:)

  if ( bcmode(iface-numInnerFaces) == MODE_EXTRAP ) then
    call extrapolated_face_assembly( iface, ijp, kd, bval, rhsb, Vsb )
  elseif ( bcmode(iface-numInnerFaces) == MODE_RIEMANN ) then
    call riemann_continuity( iface, ijp, ijb, kd, bval, rhsb )
  else
    call pressure_face_assembly( iface, ijp, ijb, kd, bval, rhsb, Eb )
  endif

end subroutine apply_riemann_invariant_bc


!***********************************************************************
!
subroutine riemann_characteristic_state( iface, ijp, ijb, rhoc, pinf, uninf, A, nx, ny, nz )
!
! Description:
!   Coefficients of the linearised characteristic relations at a subsonic boundary face (n outward, rho c from the cell):
!
!      u_nb = 1/2 [ (u_nP + u_n,inf) + (p_P - p_inf)/(rho c) ]
!      p_b  = 1/2 [ (p_P + p_inf) + rho c (u_nP - u_n,inf) ]
!
!   (exact acoustic Riemann problem between the cell state and the free stream, linear in U_P and p_P).
!
  implicit none
  integer, intent(in) :: iface, ijp, ijb
  real(dp), intent(out) :: rhoc, pinf, uninf, A, nx, ny, nz
  integer :: ib2
  ib2 = iface - numInnerFaces
  A = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
  nx = arx(iface)/A; ny = ary(iface)/A; nz = arz(iface)/A
  rhoc  = den(ijp)*sqrt( speed_of_sound_sq( T(ijp) ) )
  pinf  = p_set(ib2)
  uninf = U_set(1,ib2)*nx + U_set(2,ib2)*ny + U_set(3,ib2)*nz
end subroutine riemann_characteristic_state


!***********************************************************************
!
subroutine riemann_continuity( iface, ijp, ijb, kd, bval, rhsb )
!
! Description:
!   Continuity row of a subsonic characteristic boundary: mdot_b = rho_b u_nb A, linear in U_P and p_P.
!
  implicit none
  integer, intent(in) :: iface, ijp, ijb, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:)
  real(dp) :: rhoc, pinf, uninf, A, nx, ny, nz, rb

  call riemann_characteristic_state( iface, ijp, ijb, rhoc, pinf, uninf, A, nx, ny, nz )
  rb = den(ijb)
  bval(4,1,kd) = bval(4,1,kd) + 0.5_dp*rb*arx(iface)
  bval(4,2,kd) = bval(4,2,kd) + 0.5_dp*rb*ary(iface)
  bval(4,3,kd) = bval(4,3,kd) + 0.5_dp*rb*arz(iface)
  bval(4,4,kd) = bval(4,4,kd) + 0.5_dp*rb*A/rhoc
  rhsb(4,ijp)  = rhsb(4,ijp)  - 0.5_dp*rb*A*( uninf - pinf/rhoc )

end subroutine riemann_continuity


!***********************************************************************
!
function bc_is_riemann( iface ) result(r)
  implicit none
  integer, intent(in) :: iface
  logical :: r
  r = ( bcmode(iface-numInnerFaces) == MODE_RIEMANN )
end function bc_is_riemann


!***********************************************************************
!
subroutine bc_riemann_momentum( iface, ijp, ijb, kd, bval, rhsb )
!
! Description:
!   Pressure force of a subsonic characteristic boundary face in the momentum equations, implicit in U_P and p_P:
!   p_b S = S/2 [ p_P + p_inf + rho c ( n.U_P - u_n,inf ) ]  (all terms on the left hand side of the equations).
!
  implicit none
  integer, intent(in) :: iface, ijp, ijb, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:)
  real(dp) :: rhoc, pinf, uninf, A, nx, ny, nz, sv(3), nv(3), uold(3)
  integer :: ii, jj

  call riemann_characteristic_state( iface, ijp, ijb, rhoc, pinf, uninf, A, nx, ny, nz )
  uold = (/ u(ijp), v(ijp), w(ijp) /)
  sv = (/ arx(iface), ary(iface), arz(iface) /)
  nv = (/ nx, ny, nz /)
  do ii=1,3
    bval(ii,4,kd) = bval(ii,4,kd) + 0.5_dp*sv(ii)
    ! diagonal part implicit, the cross coupling between the velocity components lagged
    bval(ii,ii,kd) = bval(ii,ii,kd) + 0.5_dp*rhoc*sv(ii)*nv(ii)
    do jj=1,3
      if ( jj /= ii ) rhsb(ii,ijp) = rhsb(ii,ijp) - 0.5_dp*rhoc*sv(ii)*nv(jj)*uold(jj)
    enddo
    rhsb(ii,ijp) = rhsb(ii,ijp) - 0.5_dp*sv(ii)*( pinf - rhoc*uninf )
  enddo

end subroutine bc_riemann_momentum


!***********************************************************************
!
subroutine apply_compressible_outlet_bc( ijb, iface, ijp, kd, bval, rhsb, Eb, Vsb )
!
! Description:
!   Implicit compressible outlet boundary condition
!
!   Two cases based on the Mach number of the cell next to the boundary:
!   1. Subsonic (M < 1): prescribe the static pressure (Rhie-Chow flux between cell centre and face)
!   2. Supersonic (M >= 1): all variables extrapolated, mdot = rho_P U_P.S linearised in (U_P, p_P)
!
!   Implicit treatment accounts for:
!   - Pressure-density coupling through EOS  (drho/dp = 1/(R T))
!   - Velocity-pressure coupling through the momentum equation (Rhie-Chow)
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:), Vsb(:)

  if ( bcmode(iface-numInnerFaces) == MODE_EXTRAP ) then
    call extrapolated_face_assembly( iface, ijp, kd, bval, rhsb, Vsb )
  else
    call pressure_face_assembly( iface, ijp, ijb, kd, bval, rhsb, Eb )
  endif

end subroutine apply_compressible_outlet_bc


!***********************************************************************
!
subroutine apply_open_boundary_bc( ijb, iface, ijp, kd, bval, rhsb, Eb, Vsb )
!
! Description:
!   Implicit open boundary (inflow/outflow)
!
!   Automatically switches between inlet and outlet based on flow direction (sign of U_P.n)
!
!   Outflow: static pressure prescribed (subsonic) or zero gradient (supersonic); zero gradient for U and T
!   Inflow:  total pressure p_b = p_set - 0.5 rho U_n^2 with the temperature of the 0/T file (entrainment), velocity
!            normal to the boundary
!
  implicit none

  integer, intent(in) :: ijb, iface, ijp, kd
  real(dp), intent(inout) :: bval(:,:,:), rhsb(:,:), Eb(:), Vsb(:)

  if ( bcmode(iface-numInnerFaces) == MODE_EXTRAP ) then
    call extrapolated_face_assembly( iface, ijp, kd, bval, rhsb, Vsb )
  else
    call pressure_face_assembly( iface, ijp, ijb, kd, bval, rhsb, Eb )
  endif

end subroutine apply_open_boundary_bc


!***********************************************************************
!
subroutine bc_post_solve_face( ib, i, iface, ijp, ijb, Eb, Vsb )
!
! Description:
!   After the coupled solve: mass flux of the boundary face consistent with the solved continuity equation, and the
!   boundary velocity (and, for the mass flow inlet, the boundary density) used by the next momentum assembly.
!
  implicit none

  integer, intent(in) :: ib, i, iface, ijp, ijb
  real(dp), intent(in) :: Eb(:), Vsb(:)

  integer :: ib2
  real(dp) :: sx, sy, sz, dx, dy, dz, rb, are2, sdb, capb, dpbx, Cc, ain, q, are, nx, ny, nz, utx, uty, utz, un_inf, unP

  ib2 = iface - numInnerFaces
  sx = arx(iface); sy = ary(iface); sz = arz(iface)
  dx = xf(iface)-xc(ijp); dy = yf(iface)-yc(ijp); dz = zf(iface)-zc(ijp)

  if ( bckind(ib) == BCK_RIEMANN .and. bcmode(ib2) == MODE_RIEMANN ) then
    ! subsonic characteristic boundary: flux and velocity from the solved cell state (linearised Riemann problem)
    call riemann_characteristic_state( iface, ijp, ijb, rb, dpbx, un_inf, are, nx, ny, nz )
    ! rb = rho c, dpbx = p_inf, un_inf = u_n,inf, are = |S| (names reused)
    unP = u(ijp)*nx + v(ijp)*ny + w(ijp)*nz
    q = 0.5_dp*( unP + un_inf + ( p(ijp) - dpbx )/rb )
    flmass(iface) = den(ijb)*q*are
    if ( q < 0.0_dp ) then
      utx = U_set(1,ib2) - un_inf*nx; uty = U_set(2,ib2) - un_inf*ny; utz = U_set(3,ib2) - un_inf*nz
    else
      utx = u(ijp) - unP*nx; uty = v(ijp) - unP*ny; utz = w(ijp) - unP*nz
    endif
    u(ijb) = utx + q*nx; v(ijb) = uty + q*ny; w(ijb) = utz + q*nz
    return
  endif

  select case ( bckind(ib) )

  case ( BCK_COMP_INLET )

    flmass(iface) = den(ijb)*( u(ijb)*sx + v(ijb)*sy + w(ijb)*sz )
    p(ijb) = p(ijp)
    den(ijb) = p(ijb)/( Rair*T(ijb)+small )

  case ( BCK_MASSFLOW )

    ! boundary density from the cell pressure, velocity from the mass flow
    p(ijb) = p(ijp)
    den(ijb) = p(ijb)/( Rair*T(ijb)+small )
    ain = -( edir(1,ib2)*sx + edir(2,ib2)*sy + edir(3,ib2)*sz )
    q = mdot_face(ib2)/( den(ijb)*ain )
    u(ijb) = edir(1,ib2)*q; v(ijb) = edir(2,ib2)*q; w(ijb) = edir(3,ib2)*q
    flmass(iface) = -mdot_face(ib2)

  case default

    if ( bcmode(ib2) == MODE_EXTRAP ) then

      u(ijb) = u(ijp); v(ijb) = v(ijp); w(ijb) = w(ijp)
      Cc = 1.0_dp/( Rair*T(ijp)+small )
      flmass(iface) = den(ijp)*( u(ijp)*sx + v(ijp)*sy + w(ijp)*sz ) + Vsb(ib2)*( Cc*p(ijp) - den(ijp) )

    else

      rb   = den(ijb)
      are2 = sx*sx+sy*sy+sz*sz
      sdb  = sx*dx+sy*dy+sz*dz
      capb = -rb*vol(ijp)*apu(ijp)*are2/sdb
      flmass(iface) = rb*( u(ijp)*sx + v(ijp)*sy + w(ijp)*sz ) - capb*p(ijp) + Eb(ib2)

      if ( flmass(iface) < 0.0_dp .and. ( bckind(ib) == BCK_TOTAL_PRESSURE .or. bckind(ib) == BCK_STAGNATION &
           .or. bckind(ib) == BCK_OPEN .or. bckind(ib) == BCK_RIEMANN ) ) then

        ! Inflow: velocity along the inflow direction, magnitude from the mass flux
        are = sqrt( are2 )
        nx = sx/are; ny = sy/are; nz = sz/are
        utx = 0.0_dp; uty = 0.0_dp; utz = 0.0_dp
        if ( bckind(ib) == BCK_RIEMANN ) then
          ! tangential velocity of the free stream, normal velocity along -n
          un_inf = U_set(1,ib2)*nx + U_set(2,ib2)*ny + U_set(3,ib2)*nz
          utx = U_set(1,ib2) - un_inf*nx; uty = U_set(2,ib2) - un_inf*ny; utz = U_set(3,ib2) - un_inf*nz
          q = -flmass(iface)/( rb*are )
          u(ijb) = -q*nx + utx; v(ijb) = -q*ny + uty; w(ijb) = -q*nz + utz
        else
          ain = -( edir(1,ib2)*sx + edir(2,ib2)*sy + edir(3,ib2)*sz )
          q = -flmass(iface)/( rb*ain )
          u(ijb) = edir(1,ib2)*q; v(ijb) = edir(2,ib2)*q; w(ijb) = edir(3,ib2)*q
        endif

      else

        ! Velocity at the boundary face consistent with the flux (needed by momentum convection/diffusion terms):
        ! U_b = U_P - S (Vol/ap)/(S.d) [ p_b - p_P - grad p . d ]
        dpbx = ( p(ijb) - p(ijp) - ( dPdxi(1,ijp)*dx + dPdxi(2,ijp)*dy + dPdxi(3,ijp)*dz ) ) &
               * vol(ijp)*apu(ijp)/(sx*dx+sy*dy+sz*dz)
        u(ijb) = u(ijp) - sx*dpbx
        v(ijb) = v(ijp) - sy*dpbx
        w(ijb) = w(ijp) - sz*dpbx

      endif

    endif

  end select

end subroutine bc_post_solve_face

end module implicit_bc_compressible
