module turb_transport
!
! Purpose:
!   Shared finite volume engine for turbulence-model scalar transport equations
!
!     d(rho phi)/dt + div(rho U phi) - div( Gamma grad phi ) = Su - Sp*phi
!
!   Every transport equation of a turbulence model (k, epsilon, omega, zeta, gamma, Re_theta, Reynolds
!   stress components, ...) has this form - only the diffusivity Gamma and the source terms differ.
!   A model module therefore only has to:
!     1) evaluate the volume-integrated sources 'suv' (explicit) and 'spv' (implicit, >=0, multiplies phi),
!     2) evaluate the diffusivity field 'dif' (cell centres and boundary faces),
!     3) call 'transport_scalar'.
!   The engine takes care of: gradient, convection (any scheme of 'facefluxsc'), diffusion with
!   non-orthogonal correction, BDF/BDF2/Crank-Nicolson time terms, inlet/outlet/wall/periodic boundaries,
!   under-relaxation, solution of the linear system, update of boundary values, clipping.
!
!   Solution parameters (urf, convection scheme, linear solver, tolerances) come from TurbModel%Scalar(ifi).
!
!   Optional features, needed by some models:
!     * lconv=.false., ltrans=.false. - purely elliptic (Helmholtz-type) equation, e.g. the 'f' equation in zeta-f
!       and the alpha equation in the elliptic blending models.
!     * fixed/fixval - cells in which the value is prescribed (wall-function treatment of epsilon/omega).
!     * wallDirichlet - walls hold a Dirichlet value phi(ijb) with diffusive flux (otherwise zero flux, which is
!       what k-omega SST uses for k).
!     * positive=.false. - the variable may be negative (e.g. Reynolds shear stresses).
!
!   Time levels of fields owned by a model (not the global te/ed) are handled by 'shift_time_levels'.
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use sparse_matrix
  use linear_solvers
  use TurbModelData, only: TurbModel
  use scalar_fluxes, only: facefluxsc

  implicit none

  private

  public :: transport_scalar, shift_time_levels, report_range, wall_normal_tangent, copy_owner_to_boundary
  public :: wall_viscosity_loglaw

contains


subroutine transport_scalar( fi, dfidxi, fio, fioo, ifi, fname, suv, spv, dif, &
                             lconv, ltrans, wallDirichlet, fixed, fixval, positive )
!
! Assemble and solve one scalar transport equation, see module header.
!
  implicit none

  real(dp), dimension(numTotal), intent(inout) :: fi            ! Field, size numCells + numBoundaryFaces
  real(dp), dimension(3,numCells), intent(inout) :: dfidxi      ! Its gradient. CAUTION: 'grad' also writes boundary
                                                                  ! values, so the actual array must be (3,numTotal)!
  real(dp), dimension(numTotal), intent(in) :: fio, fioo        ! Old time levels (n, n-1); ignored in steady state
  integer, intent(in) :: ifi                                    ! Index of scalar in TurbModel%Scalar(:)
  character(len=*), intent(in) :: fname                         ! Name for the solver log
  real(dp), dimension(numCells), intent(in) :: suv              ! Volume integrated explicit source
  real(dp), dimension(numCells), intent(in) :: spv              ! Volume integrated implicit source coefficient (>= 0)
  real(dp), dimension(numTotal), intent(in) :: dif              ! Diffusivity (dynamic) at cells and boundary faces
  logical, intent(in), optional :: lconv, ltrans, wallDirichlet, positive
  logical, dimension(numCells), intent(in), optional :: fixed
  real(dp), dimension(numCells), intent(in), optional :: fixval

!
! Locals
!
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iper, l, iftwin
  real(dp) :: cap, can, suadd, dcoef, apotime, urfrs, urfms, off_diagonal_terms, fm, fimin
  logical :: conv, trans, wallDir, pos
  real(dp) :: gam, urf, tolAbs, tolRel
  integer :: maxiter
  character(len=12) :: lSolver
  character(len=30) :: cScheme

  conv = .true.;  if (present(lconv)) conv = lconv
  trans = ltransient;  if (present(ltrans)) trans = ltrans .and. ltransient
  wallDir = .false.;  if (present(wallDirichlet)) wallDir = wallDirichlet
  pos = .true.;  if (present(positive)) pos = positive

  gam     = TurbModel%Scalar(ifi)%gds
  urf     = TurbModel%Scalar(ifi)%urf
  lSolver = TurbModel%Scalar(ifi)%lSolver
  maxiter = TurbModel%Scalar(ifi)%maxiter
  tolAbs  = TurbModel%Scalar(ifi)%tolAbs
  tolRel  = TurbModel%Scalar(ifi)%tolRel
  cScheme = TurbModel%Scalar(ifi)%cScheme

  call grad(fi,dfidxi)

  a = 0.0_dp
  su(1:numCells) = suv
  sp(1:numCells) = spv

  !
  ! Unsteady term
  !
  if ( trans ) then
    do inp=1,numCells
      apotime = den(inp)*vol(inp)/timestep
      if( bdf .or. cn ) then
        su(inp) = su(inp) + apotime*fio(inp)
        sp(inp) = sp(inp) + apotime
      elseif( bdf2 ) then
        su(inp) = su(inp) + apotime*( 2*fio(inp) - 0.5_dp*fioo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime
      endif
    enddo
  endif

  !
  ! Inner faces
  !
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)

    dcoef = dif(ijp) + (dif(ijn)-dif(ijp))*facint(i)
    fm = 0.0_dp;  if (conv) fm = flmass(i)

    call facefluxsc( i, ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     fm, facint(i), gam, cScheme, fi, dfidxi, dcoef, cap, can, suadd )

    k = icell_jcell_csr_index(i)
    a(k) = can
    k = jcell_icell_csr_index(i)
    a(k) = cap

    k = diag(ijp)
    a(k) = a(k) - can
    k = diag(ijn)
    a(k) = a(k) - cap

    su(ijp) = su(ijp) + suadd
    su(ijn) = su(ijn) - suadd
  enddo

  !
  ! Boundaries
  !
  iPer = 0
  l = numInnerFaces

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' .or. bctype(ib) == 'pressure' .or. &
         ( bctype(ib) == 'wall' .and. wallDir ) ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        dcoef = dif(ijb)
        fm = 0.0_dp;  if (conv) fm = flmass(iface)

        call facefluxsc( ijp, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         fm, dfidxi, dcoef, cap, can, suadd )

        sp(ijp) = sp(ijp) - can
        su(ijp) = su(ijp) - can*fi(ijb) + suadd
      end do

    elseif ( bctype(ib) == 'periodic' ) then

      iPer = iPer + 1

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        iftwin = startFaceTwin(iPer) + i
        ijn = owner(iftwin)

        dcoef = half*( dif(ijp) + dif(ijn) )
        fm = 0.0_dp;  if (conv) fm = flmass(iface)

        call facefluxsc( i, ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         fm, gam, fi, dfidxi, dcoef, cap, can, suadd )

        l = l + 1

        k = icell_jcell_csr_index(l)
        a(k) = can
        k = jcell_icell_csr_index(l)
        a(k) = cap

        k = diag(ijp)
        a(k) = a(k) - can
        k = diag(ijn)
        a(k) = a(k) - cap

        su(ijp) = su(ijp) + suadd
        su(ijn) = su(ijn) - suadd
      end do

    endif

  enddo

  !
  ! Crank-Nicolson: half of the spatial operator at the old time level goes to the right hand side
  !
  if ( cn .and. trans ) then

    a = 0.5_dp*a

    do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      k = icell_jcell_csr_index(i)
      su(ijp) = su(ijp) - a(k)*fio(ijn)
      k = jcell_icell_csr_index(i)
      su(ijn) = su(ijn) - a(k)*fio(ijp)
    enddo

    do ijp=1,numCells
      apotime = den(ijp)*vol(ijp)/timestep
      off_diagonal_terms = sum( a( ia(ijp) : ia(ijp+1)-1 ) ) - a(diag(ijp))
      su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*fio(ijp)
      sp(ijp) = sp(ijp) + apotime
    enddo

  endif

  !
  ! Cells with prescribed value (wall functions): row replaced by identity.
  !
  if ( present(fixed) .and. present(fixval) ) then
    do inp=1,numCells
      if ( fixed(inp) ) then
        a( ia(inp):ia(inp+1)-1 ) = 0.0_dp
        sp(inp) = 1.0_dp
        su(inp) = fixval(inp)
        fi(inp) = fixval(inp)
      endif
    enddo
  endif

  !
  ! Diagonal assembly and implicit under-relaxation
  !
  urfrs = 1.0_dp/urf
  urfms = 1.0_dp - urf

  do inp = 1,numCells

    a(diag(inp)) = sp(inp)
    do k = ia(inp),ia(inp+1)-1
      if (k == diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) - a(k)
    enddo

    if ( present(fixed) ) then
      if ( fixed(inp) ) cycle ! no relaxation of prescribed values
    endif

    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)

  enddo

  call csrsolve(lSolver, fi, su, TurbModel%Scalar(ifi)%resor, maxiter, tolAbs, tolRel, trim(fname) )

  call updateBoundary( fi )

  call report_range( fi, fname )

  if ( pos ) then
    fimin = minval(fi(1:numCells))
    if ( fimin < 0.0_dp ) fi(1:numCells) = max( fi(1:numCells), small )
  endif

end subroutine


subroutine shift_time_levels( phi, phio, phioo, stamp )
!
! Time levels for fields that are owned by a turbulence model (global te, ed are shifted in time_shift.f90).
! To be called at the beginning of the model update: it detects a new time step by comparing 'itime' with 'stamp'.
!
  implicit none
  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(numTotal), intent(inout) :: phio, phioo
  integer, intent(inout) :: stamp

  if ( stamp < 0 ) then      ! first call - initialise both levels
    phio = phi
    phioo = phi
    stamp = itime
  elseif ( itime /= stamp ) then
    phioo = phio
    phio = phi
    stamp = itime
  endif

end subroutine


subroutine copy_owner_to_boundary( phi )
!
! Boundary face values of a cell field taken from the owner cell (zero-gradient). Useful for diffusivities.
!
  implicit none
  real(dp), dimension(numTotal), intent(inout) :: phi
  integer :: ib, i, ijb, iface

  do ib=1,numBoundaries
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijb = iBndValueStart(ib) + i
      phi(ijb) = phi(owner(iface))
    enddo
  enddo

end subroutine


subroutine wall_viscosity_loglaw( cmu25 )
!
! Wall viscosity, y* and wall shear stress for models with first cell either in the viscous sublayer or the log layer
! (hard switch at y* = ctrans):
!    y*  = rho Cmu^1/4 sqrt(k) y / mu
!    y* > ctrans : standard wall function, mu_w = y* mu kappa / ln(E y*),  tau_w = kappa rho Vt Cmu^1/4 sqrt(k) / ln(E y*)
!    otherwise   : viscous sublayer, tau_w = mu |U_t - U_wall| / y, mu_w = mu
! Sets visw, vis at wall faces, ypl, tau.
!
  implicit none
  real(dp), intent(in) :: cmu25
  integer :: ib, i, iface, ijp, ijb, iWall
  real(dp) :: nxf, nyf, nzf, vtp, xtp, ytp, ztp, ut2, utau, viscw

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

      ypl(iWall) = den(ijp)*cmu25*sqrt(te(ijp))*dnw(iWall)/viscos
      viscw = zero

      if ( ypl(iWall) > ctrans ) then
        viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
        tau(iWall) = cappa*den(ijp)*vtp*cmu25*sqrt(te(ijp))/log(Elog*ypl(iWall))
      else
        tau(iWall) = viscos*ut2/dnw(iWall)
        utau = sqrt( tau(iWall)/den(ijb) )
        ypl(iWall) = den(ijb)*utau*dnw(iWall)/viscos
      endif

      visw(iWall) = max(viscos,viscw)
      vis(ijb) = visw(iWall)
    enddo

  enddo

end subroutine


subroutine report_range( fi, fname )
  implicit none
  real(dp), dimension(numTotal), intent(in) :: fi
  character(len=*), intent(in) :: fname
  write(*,'(2x,es11.4,3a,es11.4)') minval(fi(1:numCells)),' <= ',trim(fname),' <= ',maxval(fi(1:numCells))
end subroutine


subroutine wall_normal_tangent( iface, ijp, nx, ny, nz, vt, tx, ty, tz )
!
! Wall unit normal, magnitude of the wall-parallel velocity of the first cell and its unit direction.
! Falls back to the flow direction-free tangent (0) when the cell velocity is parallel to the normal.
!
  implicit none
  integer, intent(in) :: iface, ijp
  real(dp), intent(out) :: nx, ny, nz, vt, tx, ty, tz
  real(dp) :: are, vn

  are = sqrt(arx(iface)**2 + ary(iface)**2 + arz(iface)**2)
  nx = arx(iface)/are;  ny = ary(iface)/are;  nz = arz(iface)/are
  vn = u(ijp)*nx + v(ijp)*ny + w(ijp)*nz
  tx = u(ijp) - vn*nx;  ty = v(ijp) - vn*ny;  tz = w(ijp) - vn*nz
  vt = sqrt(tx*tx + ty*ty + tz*tz)
  if ( vt > small ) then
    tx = tx/vt;  ty = ty/vt;  tz = tz/vt
  else
    tx = 0.0_dp;  ty = 0.0_dp;  tz = 0.0_dp
  endif

end subroutine

end module
