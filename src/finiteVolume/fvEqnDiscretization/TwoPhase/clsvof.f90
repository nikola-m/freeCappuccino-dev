module clsvof
!
!********************************************************************************************************************************132
!
!  Purpose:
!
!    Coupled Level Set and Volume-of-Fluid (CLSVOF) method for incompressible two-phase flows on unstructured meshes.
!
!  Discussion:
!
!    Two fields describe the interface between phase 1 (vof = 1, lsf > 0) and phase 2 (vof = 0, lsf < 0):
!
!      vof - volume fraction of phase 1, advected by a geometric (PLIC) flux calculation: conserves the volume of the phases,
!      lsf - level set (signed distance to the interface), gives the normal, the curvature and smooth rho and mu.
!
!    One time step of the interface (subroutine clsvof_advance) is, following Sussman & Puckett (2000):
!
!      1. Volume flux of the step (from the flow solver, AB2 extrapolated, or analytic for the pure advection tests)
!      2. Advect vof with the unsplit PLIC flux: for each face the swept prism (face extruded by Q*dt/|S| into the donor cell) is
!         clipped by the PLIC plane of the donor cell, the clipped volume is the volume of phase 1 which crosses the face.
!      3. Advect lsf (SSP-RK3, TVD limited, upwind-biased) with the same volume flux
!      4. Coupling: the normal of the interface is taken from lsf, the plane constant is found from the new vof (bisection on
!         the clipped cell volume), and lsf is reset to the exact distance to the PLIC planes in a narrow band of cells
!         around the interface (Sussman & Puckett 2000, Sec. 2.4). Level set is thus volume-consistent and a distance function.
!      5. rho and mu from the smoothed Heaviside function H(lsf), curvature kappa = -div( grad lsf / |grad lsf| ),
!         CSF surface force f = sigma kappa grad H together with the gravity force (rho - densit) g, and the balanced-force
!         face correction fbal used in the Rhie-Chow flux of the pressure equation (see twophase_pressure).
!
!    Sign convention: lsf > 0 in phase 1, normal n = grad(lsf)/|grad(lsf)| points into phase 1.
!    Plane of a cell: phase 1 is the side { n.x - d >= 0 } of the plane (pnx,pny,pnz,pd).
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Reference:
!
!    Sussman M., Puckett E.G., A Coupled Level Set and Volume-of-Fluid Method for Computing 3D and Axisymmetric
!      Incompressible Two-phase Flows, J. Comput. Phys. 162 (2000) 301-337.
!    Sussman M., Smereka P., Osher S., A Level Set Approach for Computing Solutions to Incompressible Two-Phase Flow,
!      J. Comput. Phys. 114 (1994) 146-159.
!    Francois M.M. et al., A balanced-force algorithm for continuous and sharp interfacial surface tension models within a
!      volume tracking framework, J. Comput. Phys. 213 (2006) 141-173.
!
  use types
  use parameters
  use geometry
  use variables, only: den, vis, visw, u, v, w
  use gradients, only: grad

  implicit none

  !
  ! Parameters - modified trough input.nml file
  !
  logical  :: twoPhase = .false.               ! Activates the two-phase (CLSVOF) solver.
  real(dp) :: cl_rho1 = 1000.0_dp              ! Density of phase 1 [kg/m3].
  real(dp) :: cl_rho2 = 1.0_dp                 ! Density of phase 2 [kg/m3].
  real(dp) :: cl_mu1 = 1.0e-3_dp               ! Dynamic viscosity of phase 1 [Pa s].
  real(dp) :: cl_mu2 = 1.8e-5_dp               ! Dynamic viscosity of phase 2 [Pa s].
  real(dp) :: cl_sigma = 0.0_dp                ! Surface tension coefficient [N/m].
  real(dp) :: cl_eps = 1.5_dp                  ! Half-width of the smoothed Heaviside function in cell sizes.
  integer  :: cl_band = 4                      ! Half-width of the narrow band where lsf is a distance function [cells].
  integer  :: cl_nReinit = 0                   ! Sussman-Smereka-Osher PDE reinitialization sweeps outside of the band.
  real(dp) :: cl_cfl = 0.25_dp                 ! Max. Courant number of one advection substep (subcycling above it).
  real(dp) :: cl_fmin = 1.0e-8_dp              ! Cells with cl_fmin < vof < 1-cl_fmin are interface cells.
  logical  :: cl_ab2 = .true.                  ! Second order Adams-Bashforth extrapolation of the volume flux.
  logical  :: cl_harmonic = .false.            ! Harmonic instead of arithmetic mixing of the viscosity.
  character( len=12 ) :: cl_ic = 'sphere'      ! Initial interface: sphere, cylinder (axis z), plane (phase 1 below z=center(3)).
  integer  :: cl_inside = 1                    ! Phase inside the sphere/cylinder: 1 (drop) or 2 (bubble).
  real(dp) :: cl_center(3) = 0.0_dp            ! Center of the sphere/cylinder.
  real(dp) :: cl_radius = 0.1_dp               ! Radius of the sphere/cylinder.
  character( len=12 ) :: cl_flow = 'solve'     ! solve (Navier-Stokes), or prescribed: translate, rotate, vortex.
  real(dp) :: cl_vel(3) = 0.0_dp               ! Velocity of the translation test.
  real(dp) :: cl_omega = 1.0_dp                ! Angular velocity of the rotation test (about z through cl_vcenter).
  real(dp) :: cl_vcenter(2) = 0.5_dp           ! Center of rotation in x-y plane.
  real(dp) :: cl_T = 8.0_dp                    ! Period of the reversing vortex-in-a-box test.
  real(dp) :: cl_inflowF = 0.0_dp              ! Volume fraction of fluid entering through open boundaries (pure advection tests).

  ! Fields
  real(dp), dimension(:), allocatable :: vof          ! Volume fraction of phase 1, [numTotal]
  real(dp), dimension(:), allocatable :: lsf          ! Level set function, [numTotal]
  real(dp), dimension(:), allocatable :: kappa        ! Curvature, [numTotal]
  real(dp), dimension(:), allocatable :: hs           ! Smoothed Heaviside function of lsf, [numTotal]
  real(dp), dimension(:), allocatable :: flvol        ! Volume flux through faces (along face area vector), [numFaces]
  real(dp), dimension(:), allocatable :: fbal         ! Balanced-force correction of pressure difference, [numInnerFaces]
  real(dp), dimension(:,:), allocatable :: fst        ! Body force per unit volume: surface tension + gravity, [3,numCells]

  logical :: cl_solveFlow = .true.                    ! False if the flow field is prescribed.
  real(dp) :: cl_time = 0.0_dp                        ! Time of the interface.
  real(dp) :: cl_volume0 = 0.0_dp                     ! Initial volume of phase 1.
  real(dp) :: cl_clipVol = 0.0_dp                     ! Cumulated volume of phase 1 clipped to keep 0 <= vof <= 1.
  real(dp) :: cl_contErr = 0.0_dp                     ! Max. volume continuity error of the last pressure correction.

  ! Private data
  real(dp), dimension(:), allocatable, private :: flvolo, qeff                  ! Flux of previous step, effective flux
  real(dp), dimension(:), allocatable, private :: q1acc                        ! Volume of phase 1 advected through faces in the step
  real(dp), dimension(:), allocatable, private :: fvx, fvy, fvz                 ! Velocity at face centers, normal component = qeff/|S|
  logical, dimension(:), allocatable, private :: planevalid                      ! Plane of the cell is valid
  real(dp), dimension(:), allocatable, private :: pnx, pny, pnz, pd              ! PLIC planes
  real(dp), dimension(:), allocatable, private :: hloc, rcell, vtet, rfrag       ! Cell size, radius, volume, plane fragment radius
  real(dp), dimension(:,:), allocatable, private :: gphi                         ! Gradient of lsf
  integer, dimension(:), allocatable, private :: cf_ia, cf_ja                    ! Cell to faces connectivity (CSR)
  integer, dimension(:), allocatable, private :: mixed, bfs_queue, bfs_depth, bfs_stamp
  integer, private :: nmixed = 0, stamp = 0
  real(dp), private :: dt_old = 0.0_dp
  logical, private :: first_step = .true.
  integer, private :: cl_unit = 0

  private

  public :: twoPhase, cl_rho1, cl_rho2, cl_mu1, cl_mu2, cl_sigma, cl_eps, cl_band, cl_nReinit, cl_cfl, cl_fmin, cl_ab2
  public :: cl_harmonic, cl_ic, cl_inside, cl_center, cl_radius, cl_flow, cl_vel, cl_omega, cl_vcenter, cl_T, cl_inflowF
  public :: vof, lsf, kappa, hs, flvol, fbal, fst, cl_solveFlow, cl_time, cl_contErr, cl_clipVol
  public :: clsvof_init, clsvof_advance, clsvof_monitor, clsvof_write_restart


contains


!***********************************************************************
!
subroutine clsvof_init
!
!***********************************************************************
!
!  Allocates the fields, builds the connectivity used by the geometric routines, creates the initial interface
!  (or reads the restart file) and sets the initial density, viscosity, curvature and forces.
!
  implicit none

  integer :: inp, i, ib, ijp, ijn, ijb
  real(dp) :: uf, vf, wf, lam

  if( .not. twoPhase ) return

  cl_solveFlow = ( trim(cl_flow) == 'solve' )

  if( .not.ltransient ) then
    write(6,'(a)') ' *** CLSVOF: two-phase flow is solved only as a transient problem (ltransient = T)!'
    stop
  endif

  call check_boundaries

  allocate( vof(numTotal), lsf(numTotal), kappa(numTotal), hs(numTotal) )
  allocate( flvol(numFaces), flvolo(numFaces), qeff(numFaces), fbal(numInnerFaces), fst(3,numCells) )
  allocate( fvx(numFaces), fvy(numFaces), fvz(numFaces), planevalid(numCells), q1acc(numFaces) )
  q1acc = 0.0_dp
  allocate( pnx(numCells), pny(numCells), pnz(numCells), pd(numCells) )
  allocate( hloc(numCells), rcell(numCells), vtet(numCells), rfrag(numCells) )
  allocate( gphi(3,numTotal) )
  allocate( mixed(numCells), bfs_queue(numCells), bfs_depth(numCells), bfs_stamp(numCells) )

  vof = 0.0_dp
  lsf = 0.0_dp
  kappa = 0.0_dp
  hs = 0.0_dp
  flvol = 0.0_dp
  flvolo = 0.0_dp
  qeff = 0.0_dp
  fbal = 0.0_dp
  fst = 0.0_dp
  pnx = 0.0_dp; pny = 0.0_dp; pnz = 1.0_dp; pd = 0.0_dp
  gphi = 0.0_dp
  bfs_stamp = 0
  fvx = 0.0_dp; fvy = 0.0_dp; fvz = 0.0_dp
  planevalid = .false.

  ! Wall viscosity of the momentum equations is max(viscos,visw): use the smaller viscosity as the floor
  viscos = min( cl_mu1, cl_mu2 )

  call build_cell_faces
  call cell_sizes

  if( lread ) then
    call clsvof_read_restart
  else
    call init_interface
  endif

  ! Volume flux of the initial velocity field
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)
    uf = u(ijp) + (u(ijn)-u(ijp))*lam
    vf = v(ijp) + (v(ijn)-v(ijp))*lam
    wf = w(ijp) + (w(ijn)-w(ijp))*lam
    flvol(i) = uf*arx(i) + vf*ary(i) + wf*arz(i)
  enddo
  do ib=1,numBoundaries
    do i=1,nfaces(ib)
      ijb = iBndValueStart(ib) + i
      flvol( startFace(ib)+i ) = u(ijb)*arx(startFace(ib)+i) + v(ijb)*ary(startFace(ib)+i) + w(ijb)*arz(startFace(ib)+i)
    enddo
  enddo
  flvolo = flvol

  call update_properties
  call init_mass_flux

  cl_volume0 = sum( vof(1:numCells)*vol(1:numCells) )

  ! History file
  open( newunit=cl_unit, file='clsvof_history.dat', status='replace' )
  write(cl_unit,'(a)') '# time  V1  (V1-V10)/V10  clip  xc1 yc1 zc1  xc2 yc2 zc2  U1x U1y U1z  U2x U2y U2z  '// &
                       'maxU  area  contErr  nInterface'

  write(6,'(a)') ' '
  write(6,'(a)') '  Coupled Level Set - Volume of Fluid (CLSVOF) two-phase solver'
  write(6,'(a,es12.4)') '    Density ratio rho1/rho2       : ', cl_rho1/cl_rho2
  write(6,'(a,es12.4)') '    Viscosity ratio mu1/mu2       : ', cl_mu1/cl_mu2
  write(6,'(a,es12.4)') '    Surface tension [N/m]         : ', cl_sigma
  write(6,'(a,es12.4)') '    Initial volume of phase 1     : ', cl_volume0
  write(6,'(a,i8)')     '    Interface cells               : ', nmixed
  write(6,'(a)') ' '

end subroutine


!***********************************************************************
!
subroutine check_boundaries
!
!***********************************************************************
!
  implicit none
  integer :: ib

  do ib=1,numBoundaries
    select case( trim(bctype(ib)) )
    case( 'wall', 'symmetry', 'empty' )
      ! supported
    case default
      if( cl_solveFlow ) then
        write(6,'(3a)') ' *** CLSVOF: boundary type ''', trim(bctype(ib)), &
                        ''' is not supported by the two-phase flow solver (wall, symmetry, empty only)!'
        stop
      endif
    end select
  enddo

end subroutine


!***********************************************************************
!
subroutine build_cell_faces
!
!***********************************************************************
!
!  Cell to face connectivity in CSR format (inner and boundary faces).
!
  implicit none
  integer :: iface, inp, k
  integer, dimension(:), allocatable :: cnt

  allocate( cf_ia(numCells+1) )
  allocate( cnt(numCells) )
  cnt = 0

  do iface=1,numInnerFaces
    cnt(owner(iface)) = cnt(owner(iface)) + 1
    cnt(neighbour(iface)) = cnt(neighbour(iface)) + 1
  enddo
  do iface=numInnerFaces+1,numFaces
    cnt(owner(iface)) = cnt(owner(iface)) + 1
  enddo

  cf_ia(1) = 1
  do inp=1,numCells
    cf_ia(inp+1) = cf_ia(inp) + cnt(inp)
  enddo

  allocate( cf_ja(cf_ia(numCells+1)-1) )

  cnt = 0
  do iface=1,numFaces
    inp = owner(iface)
    k = cf_ia(inp) + cnt(inp)
    cf_ja(k) = iface
    cnt(inp) = cnt(inp) + 1
    if( iface <= numInnerFaces ) then
      inp = neighbour(iface)
      k = cf_ia(inp) + cnt(inp)
      cf_ja(k) = iface
      cnt(inp) = cnt(inp) + 1
    endif
  enddo

  deallocate( cnt )

end subroutine


!***********************************************************************
!
subroutine cell_sizes
!
!***********************************************************************
!
!  Local cell size (smallest distance to a neighbour cell center - robust for one-cell-thick 2D meshes),
!  radius of the cell (largest distance of a vertex from cell center), and the cell volume
!  as the sum of the volumes of the tetrahedra used by the clipping routines.
!
  implicit none
  integer :: i, inp, k, iface, m, nn, n1
  real(dp) :: d, pc(3), pf(3), pa(3), pb(3), vsum

  hloc = great
  do i=1,numInnerFaces
    d = sqrt( (xc(neighbour(i))-xc(owner(i)))**2 + (yc(neighbour(i))-yc(owner(i)))**2 + (zc(neighbour(i))-zc(owner(i)))**2 )
    hloc(owner(i)) = min( hloc(owner(i)), d )
    hloc(neighbour(i)) = min( hloc(neighbour(i)), d )
  enddo

  do inp=1,numCells
    pc = (/ xc(inp), yc(inp), zc(inp) /)
    rcell(inp) = 0.0_dp
    vsum = 0.0_dp
    do k=cf_ia(inp), cf_ia(inp+1)-1
      iface = cf_ja(k)
      nn = nnodes(iface)
      pf = (/ xf(iface), yf(iface), zf(iface) /)
      do m=1,nn
        n1 = node(m,iface)
        pa = (/ x(n1), y(n1), z(n1) /)
        n1 = node(mod(m,nn)+1,iface)
        pb = (/ x(n1), y(n1), z(n1) /)
        rcell(inp) = max( rcell(inp), sqrt(sum((pa-pc)**2)) )
        vsum = vsum + tet_volume( pc, pf, pa, pb )
      enddo
    enddo
    vtet(inp) = vsum
    rfrag(inp) = 0.5_dp*hloc(inp)
  enddo

end subroutine


!***********************************************************************
!
! GEOMETRY KERNELS - volumes of polyhedra clipped by a plane
!
!***********************************************************************

pure function tet_volume( a, b, c, d ) result(vv)
  implicit none
  real(dp), dimension(3), intent(in) :: a, b, c, d
  real(dp) :: vv
  real(dp), dimension(3) :: e1, e2, e3

  e1 = b - a
  e2 = c - a
  e3 = d - a
  vv = abs( e1(1)*(e2(2)*e3(3)-e2(3)*e3(2)) - e1(2)*(e2(1)*e3(3)-e2(3)*e3(1)) + e1(3)*(e2(1)*e3(2)-e2(2)*e3(1)) )/6.0_dp

end function


pure function tet_pos_volume( p, s ) result(vp)
!
!  Volume of the part of the tetrahedron with vertices p(:,1:4) where the linear function with the vertex values s(1:4)
!  is positive. Exact; no division by differences of values of the same sign (no cancellation).
!
  implicit none
  real(dp), dimension(3,4), intent(in) :: p
  real(dp), dimension(4), intent(in) :: s
  real(dp) :: vp

  integer :: i, j, k, np, nn, ia, ib, ic, id
  integer, dimension(4) :: ipos, ineg
  real(dp) :: t
  real(dp), dimension(3,3) :: q
  real(dp), dimension(3) :: pac, pad, pbc, pbd

  np = 0
  nn = 0
  do i=1,4
    if( s(i) > 0.0_dp ) then
      np = np + 1
      ipos(np) = i
    else
      nn = nn + 1
      ineg(nn) = i
    endif
  enddo

  select case( np )

  case( 0 )
    vp = 0.0_dp

  case( 4 )
    vp = tet_volume( p(:,1), p(:,2), p(:,3), p(:,4) )

  case( 1 )
    ia = ipos(1)
    do k=1,3
      j = ineg(k)
      t = s(ia)/( s(ia) - s(j) )
      q(:,k) = p(:,ia) + t*( p(:,j) - p(:,ia) )
    enddo
    vp = tet_volume( p(:,ia), q(:,1), q(:,2), q(:,3) )

  case( 3 )
    ia = ineg(1)
    do k=1,3
      j = ipos(k)
      t = s(ia)/( s(ia) - s(j) )
      q(:,k) = p(:,ia) + t*( p(:,j) - p(:,ia) )
    enddo
    vp = tet_volume( p(:,1), p(:,2), p(:,3), p(:,4) ) - tet_volume( p(:,ia), q(:,1), q(:,2), q(:,3) )

  case( 2 )
    ia = ipos(1)
    ib = ipos(2)
    ic = ineg(1)
    id = ineg(2)
    pac = p(:,ia) + s(ia)/( s(ia)-s(ic) )*( p(:,ic) - p(:,ia) )
    pad = p(:,ia) + s(ia)/( s(ia)-s(id) )*( p(:,id) - p(:,ia) )
    pbc = p(:,ib) + s(ib)/( s(ib)-s(ic) )*( p(:,ic) - p(:,ib) )
    pbd = p(:,ib) + s(ib)/( s(ib)-s(id) )*( p(:,id) - p(:,ib) )
    ! Wedge decomposed into three tetrahedra
    vp = tet_volume( p(:,ia), pac, pad, p(:,ib) ) + tet_volume( pac, pad, p(:,ib), pbc ) + tet_volume( pad, p(:,ib), pbc, pbd )

  end select

end function


function cell_pos_volume( inp, nx, ny, nz, d ) result(vp)
!
!  Volume of the part of cell inp on the side { n.x - d > 0 } of the plane.
!
  implicit none
  integer, intent(in) :: inp
  real(dp), intent(in) :: nx, ny, nz, d
  real(dp) :: vp

  integer :: k, iface, nn, m, n1
  real(dp), dimension(3,4) :: p
  real(dp), dimension(4) :: s

  vp = 0.0_dp
  p(:,1) = (/ xc(inp), yc(inp), zc(inp) /)
  s(1) = nx*xc(inp) + ny*yc(inp) + nz*zc(inp) - d

  do k=cf_ia(inp), cf_ia(inp+1)-1
    iface = cf_ja(k)
    nn = nnodes(iface)
    p(:,2) = (/ xf(iface), yf(iface), zf(iface) /)
    s(2) = nx*xf(iface) + ny*yf(iface) + nz*zf(iface) - d
    do m=1,nn
      n1 = node(m,iface)
      p(:,3) = (/ x(n1), y(n1), z(n1) /)
      s(3) = nx*x(n1) + ny*y(n1) + nz*z(n1) - d
      n1 = node(mod(m,nn)+1,iface)
      p(:,4) = (/ x(n1), y(n1), z(n1) /)
      s(4) = nx*x(n1) + ny*y(n1) + nz*z(n1) - d
      vp = vp + tet_pos_volume( p, s )
    enddo
  enddo

end function


subroutine plane_constant( inp, nx, ny, nz, fv, d )
!
!  Plane constant d such that the volume of cell inp on the side { n.x - d > 0 } is fv*V (bisection).
!
  implicit none
  integer, intent(in) :: inp
  real(dp), intent(in) :: nx, ny, nz, fv
  real(dp), intent(out) :: d

  integer :: k, iface, m, it
  real(dp) :: lo, hi, mid, sp, target, vp

  lo = great
  hi = -great
  do k=cf_ia(inp), cf_ia(inp+1)-1
    iface = cf_ja(k)
    do m=1,nnodes(iface)
      sp = nx*x(node(m,iface)) + ny*y(node(m,iface)) + nz*z(node(m,iface))
      lo = min( lo, sp )
      hi = max( hi, sp )
    enddo
  enddo

  target = fv*vtet(inp)

  ! The volume on the positive side decreases with d, from V (d=lo) to 0 (d=hi).
  do it=1,50
    mid = 0.5_dp*( lo + hi )
    vp = cell_pos_volume( inp, nx, ny, nz, mid )
    if( vp > target ) then
      lo = mid
    else
      hi = mid
    endif
  enddo

  d = 0.5_dp*( lo + hi )

end subroutine


function plane_area( inp, nx, ny, nz, d ) result(ar)
!
!  Area of the plane inside the cell = - d(Vpos)/dd, by central difference (Vpos is a piecewise cubic function of d).
!
  implicit none
  integer, intent(in) :: inp
  real(dp), intent(in) :: nx, ny, nz, d
  real(dp) :: ar, dl

  dl = 1.0e-4_dp*hloc(inp)
  ar = ( cell_pos_volume( inp, nx, ny, nz, d-dl ) - cell_pos_volume( inp, nx, ny, nz, d+dl ) )/( 2*dl )

end function


!***********************************************************************
!
subroutine init_interface
!
!***********************************************************************
!
!  Initial vof from the exact normal and distance of the shape: the plane tangent at the closest point of the surface
!  is clipped with the cell. Initial lsf is then rebuilt from the planes, as in every time step.
!
  implicit none

  integer :: inp, ib
  real(dp) :: phi, nx, ny, nz, dd, fv

  nmixed = 0

  do inp=1,numCells

    call shape_distance( xc(inp), yc(inp), zc(inp), phi, nx, ny, nz )

    if( phi > rcell(inp) ) then
      vof(inp) = 1.0_dp
    elseif( phi < -rcell(inp) ) then
      vof(inp) = 0.0_dp
    else
      dd = nx*xc(inp) + ny*yc(inp) + nz*zc(inp) - phi
      fv = cell_pos_volume( inp, nx, ny, nz, dd )/vtet(inp)
      if( fv > 1.0_dp-cl_fmin ) fv = 1.0_dp
      if( fv < cl_fmin ) fv = 0.0_dp
      vof(inp) = fv
      pnx(inp) = nx
      pny(inp) = ny
      pnz(inp) = nz
      pd(inp)  = dd
      planevalid(inp) = .true.
    endif

  enddo

  call set_boundary_values( vof )

  call list_mixed_cells
  call lsf_from_planes

end subroutine


subroutine shape_distance( xx, yy, zz, phi, nx, ny, nz )
!
!  Signed distance (positive in phase 1) and unit normal (towards phase 1) of the initial shape.
!
  implicit none
  real(dp), intent(in) :: xx, yy, zz
  real(dp), intent(out) :: phi, nx, ny, nz

  real(dp) :: rx, ry, rz, r, sg

  sg = 1.0_dp
  if( cl_inside == 2 ) sg = -1.0_dp  ! phase 2 inside: lsf = r - R

  select case( trim(cl_ic) )

  case( 'sphere' )
    rx = xx - cl_center(1)
    ry = yy - cl_center(2)
    rz = zz - cl_center(3)
    r = max( sqrt( rx**2 + ry**2 + rz**2 ), small )
    phi = sg*( cl_radius - r )
    nx = -sg*rx/r
    ny = -sg*ry/r
    nz = -sg*rz/r

  case( 'cylinder' )
    rx = xx - cl_center(1)
    ry = yy - cl_center(2)
    r = max( sqrt( rx**2 + ry**2 ), small )
    phi = sg*( cl_radius - r )
    nx = -sg*rx/r
    ny = -sg*ry/r
    nz = 0.0_dp

  case( 'plane' )
    phi = cl_center(3) - zz
    nx = 0.0_dp
    ny = 0.0_dp
    nz = -1.0_dp

  case default
    write(6,'(3a)') ' *** CLSVOF: unknown initial interface cl_ic = ', trim(cl_ic), ' (sphere, cylinder, plane)'
    stop

  end select

end subroutine


!***********************************************************************
!
subroutine set_boundary_values( fld )
!
!***********************************************************************
!
!  Zero gradient values of a field at boundary faces.
!
  implicit none
  real(dp), dimension(numTotal), intent(inout) :: fld
  integer :: i

  do i=1,numBoundaryFaces
    fld(numCells+i) = fld( owner(numInnerFaces+i) )
  enddo

end subroutine


subroutine list_mixed_cells
  implicit none
  integer :: inp

  nmixed = 0
  do inp=1,numCells
    if( vof(inp) > cl_fmin .and. vof(inp) < 1.0_dp-cl_fmin ) then
      nmixed = nmixed + 1
      mixed(nmixed) = inp
    endif
  enddo

end subroutine


!***********************************************************************
!
subroutine lsf_from_planes
!
!***********************************************************************
!
!  Level set as the distance to the PLIC planes. Interface cells get the exact signed distance of the cell center to
!  the plane of the cell. Cells within cl_band layers of an interface cell get the smallest distance to the plane
!  fragments (disks of radius rfrag in the planes) of the interface cells which have them in their layers;
!  further cells get the capped value +-(cl_band+1.5)*h.
!
  implicit none

  integer :: j, jc, inp, k, iface, ic, head, tail, depth, kk
  real(dp) :: sg, dist, absd, cap, sj, tx, ty, tz, tt, ex, qx, qy, qz
  real(dp), dimension(:), allocatable :: absmin

  allocate( absmin(numCells) )

  do inp=1,numCells
    cap = ( cl_band + 1.5_dp )*hloc(inp)
    if( vof(inp) > 0.5_dp ) then
      lsf(inp) = cap
    else
      lsf(inp) = -cap
    endif
    absmin(inp) = cap
  enddo

  do j=1,nmixed

    jc = mixed(j)
    stamp = stamp + 1
    head = 1
    tail = 1
    bfs_queue(1) = jc
    bfs_depth(1) = 0
    bfs_stamp(jc) = stamp

    ! Foot point of the cell center on the plane and the signed distance
    sj = pnx(jc)*xc(jc) + pny(jc)*yc(jc) + pnz(jc)*zc(jc) - pd(jc)
    qx = xc(jc) - sj*pnx(jc)
    qy = yc(jc) - sj*pny(jc)
    qz = zc(jc) - sj*pnz(jc)

    do while( head <= tail )

      ic = bfs_queue(head)
      depth = bfs_depth(head)
      head = head + 1

      ! Distance of the cell center ic to the fragment of the plane of cell jc
      sg = pnx(jc)*xc(ic) + pny(jc)*yc(ic) + pnz(jc)*zc(ic) - pd(jc)
      tx = xc(ic) - sg*pnx(jc) - qx
      ty = yc(ic) - sg*pny(jc) - qy
      tz = zc(ic) - sg*pnz(jc) - qz
      tt = sqrt( tx**2 + ty**2 + tz**2 )
      ex = max( 0.0_dp, tt - rfrag(jc) )
      absd = sqrt( sg**2 + ex**2 )

      if( absd < absmin(ic) ) then
        absmin(ic) = absd
        lsf(ic) = sign( absd, sg )
      endif

      if( depth < cl_band ) then
        do k=cf_ia(ic), cf_ia(ic+1)-1
          iface = cf_ja(k)
          if( iface > numInnerFaces ) cycle
          kk = owner(iface)
          if( kk == ic ) kk = neighbour(iface)
          if( bfs_stamp(kk) /= stamp ) then
            bfs_stamp(kk) = stamp
            tail = tail + 1
            bfs_queue(tail) = kk
            bfs_depth(tail) = depth + 1
          endif
        enddo
      endif

    enddo

  enddo

  ! Exact values in interface cells
  do j=1,nmixed
    jc = mixed(j)
    lsf(jc) = pnx(jc)*xc(jc) + pny(jc)*yc(jc) + pnz(jc)*zc(jc) - pd(jc)
  enddo

  call set_boundary_values( lsf )

  deallocate( absmin )

  if( cl_nReinit > 0 ) call reinitialize_pde

end subroutine


!***********************************************************************
!
subroutine reinitialize_pde
!
!***********************************************************************
!
!  Sussman, Smereka & Osher (1994) redistancing   lsf_t + S(lsf0) ( |grad lsf| - 1 ) = 0   outside the cells which
!  carry the exact distance (interface cells and their face neighbours; Russo & Smereka fixing, so the interface does not move).
!  |grad lsf| is a Godunov-type Gauss gradient with the face value taken from the upwind cell (the cell closer to the interface).
!
  implicit none

  integer :: it, i, ijp, ijn, inp, ib
  real(dp) :: dtau, fc, gx, gy, gz, gm, s0
  real(dp), dimension(:), allocatable :: phi0, phin
  logical, dimension(:), allocatable :: fixed

  allocate( phi0(numTotal), phin(numTotal), fixed(numCells) )
  phi0 = lsf
  fixed = .false.

  do inp=1,nmixed
    fixed( mixed(inp) ) = .true.
  enddo
  do i=1,numInnerFaces
    if( fixed(owner(i)) .or. fixed(neighbour(i)) ) then
      if( vof(owner(i)) > cl_fmin .and. vof(owner(i)) < 1.0_dp-cl_fmin ) fixed(neighbour(i)) = .true.
      if( vof(neighbour(i)) > cl_fmin .and. vof(neighbour(i)) < 1.0_dp-cl_fmin ) fixed(owner(i)) = .true.
    endif
  enddo

  do it=1,cl_nReinit

    gphi(:,1:numCells) = 0.0_dp

    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      s0 = sign( 1.0_dp, 0.5_dp*( phi0(ijp) + phi0(ijn) ) )
      ! Upwind value: the characteristics leave the interface, i.e. the face value is taken from the cell with smaller |lsf|
      if( abs(lsf(ijp)) <= abs(lsf(ijn)) ) then
        fc = lsf(ijp)
      else
        fc = lsf(ijn)
      endif
      gphi(1,ijp) = gphi(1,ijp) + fc*arx(i);  gphi(1,ijn) = gphi(1,ijn) - fc*arx(i)
      gphi(2,ijp) = gphi(2,ijp) + fc*ary(i);  gphi(2,ijn) = gphi(2,ijn) - fc*ary(i)
      gphi(3,ijp) = gphi(3,ijp) + fc*arz(i);  gphi(3,ijn) = gphi(3,ijn) - fc*arz(i)
    enddo
    do i=1,numBoundaryFaces
      ijp = owner(numInnerFaces+i)
      fc = lsf(ijp)
      gphi(1,ijp) = gphi(1,ijp) + fc*arx(numInnerFaces+i)
      gphi(2,ijp) = gphi(2,ijp) + fc*ary(numInnerFaces+i)
      gphi(3,ijp) = gphi(3,ijp) + fc*arz(numInnerFaces+i)
    enddo

    phin = lsf
    do inp=1,numCells
      if( fixed(inp) ) cycle
      gm = sqrt( sum( (gphi(:,inp)/vol(inp))**2 ) )
      s0 = phi0(inp)/sqrt( phi0(inp)**2 + hloc(inp)**2 )
      dtau = 0.3_dp*hloc(inp)
      phin(inp) = lsf(inp) - dtau*s0*( gm - 1.0_dp )
    enddo
    lsf = phin
    call set_boundary_values( lsf )

  enddo

  deallocate( phi0, phin, fixed )

end subroutine


!***********************************************************************
!
! TIME ADVANCEMENT OF THE INTERFACE
!
!***********************************************************************

subroutine clsvof_advance
!
!  Advances the interface from time level n to n+1 and sets the properties for the momentum and pressure equations
!  of the time step. Called once per time step before the iteration loop.
!
  implicit none

  integer :: isub, nsub
  real(dp) :: dt, r, cflmax, tmid

  if( .not. twoPhase ) return

  dt = timestep

  ! 1. Volume flux for the advection
  if( cl_solveFlow ) then
    if( cl_ab2 .and. .not.first_step .and. dt_old > 0.0_dp ) then
      r = dt/dt_old
      qeff = ( 1.0_dp + 0.5_dp*r )*flvol - 0.5_dp*r*flvolo
    else
      qeff = flvol
    endif
    flvolo = flvol
  else
    tmid = cl_time + 0.5_dp*dt
    call analytic_flux( tmid, qeff )
    flvol = qeff
  endif

  ! Velocity at faces for the swept regions of the geometric flux
  call face_velocities( qeff, cl_time + 0.5_dp*dt )

  ! Number of substeps for the Courant number limit
  call advective_courant( qeff, dt, cflmax )
  nsub = max( 1, ceiling( cflmax/cl_cfl ) )

  ! 2. Advect the volume fraction with the PLIC planes of the beginning of the step
  q1acc = 0.0_dp
  do isub=1,nsub
    if( isub > 1 ) call update_plane_constants
    call advect_vof( qeff, dt/nsub )
    call list_mixed_cells
  enddo

  ! 3. Advect the level set
  do isub=1,nsub
    call advect_lsf( qeff, dt/nsub )
  enddo

  ! 4. Coupling: planes from the advected level set normal and the new volume fraction, lsf from the planes
  call reconstruct_planes
  call lsf_from_planes

  ! 5. Properties and forces for the momentum and pressure equations
  call update_properties
  call advection_mass_flux( qeff, dt )

  cl_time = cl_time + dt
  dt_old = dt
  first_step = .false.

end subroutine


subroutine advective_courant( q, dt, cflmax )
  implicit none
  real(dp), dimension(numFaces), intent(in) :: q
  real(dp), intent(in) :: dt
  real(dp), intent(out) :: cflmax
  integer :: i
  real(dp), dimension(:), allocatable :: outv

  allocate( outv(numCells) )
  outv = 0.0_dp

  do i=1,numInnerFaces
    if( q(i) > 0.0_dp ) then
      outv(owner(i)) = outv(owner(i)) + q(i)
    else
      outv(neighbour(i)) = outv(neighbour(i)) - q(i)
    endif
  enddo
  do i=numInnerFaces+1,numFaces
    if( q(i) > 0.0_dp ) outv(owner(i)) = outv(owner(i)) + q(i)
  enddo

  cflmax = maxval( outv*dt/vol(1:numCells) )

  deallocate( outv )

end subroutine


subroutine update_plane_constants
!
!  Recompute the plane constant of interface cells for the current vof, with the normals unchanged (between substeps).
!
  implicit none
  integer :: j, inp

  do j=1,nmixed
    inp = mixed(j)
    if( .not. planevalid(inp) ) then
      ! Cell which became an interface cell within this time step: normal from the volume fraction
      call youngs_normal( inp, pnx(inp), pny(inp), pnz(inp) )
      planevalid(inp) = .true.
    endif
    call plane_constant( inp, pnx(inp), pny(inp), pnz(inp), vof(inp), pd(inp) )
  enddo

end subroutine


subroutine reconstruct_planes
!
!  PLIC planes: normal from the gradient of the level set, plane constant from the volume fraction.
!  In interface cells with a vanishing level set gradient the normal of the vof gradient (Youngs) is used.
!
  implicit none
  integer :: j, inp
  real(dp) :: gm, nx, ny, nz

  call set_boundary_values( lsf )
  call grad( lsf, gphi )

  call list_mixed_cells

  planevalid = .false.

  do j=1,nmixed
    inp = mixed(j)
    planevalid(inp) = .true.
    gm = sqrt( sum( gphi(:,inp)**2 ) )
    if( gm > 1.0e-6_dp ) then
      nx = gphi(1,inp)/gm
      ny = gphi(2,inp)/gm
      nz = gphi(3,inp)/gm
    else
      call youngs_normal( inp, nx, ny, nz )
    endif
    pnx(inp) = nx
    pny(inp) = ny
    pnz(inp) = nz
    call plane_constant( inp, nx, ny, nz, vof(inp), pd(inp) )
  enddo

end subroutine


subroutine youngs_normal( inp, nx, ny, nz )
  implicit none
  integer, intent(in) :: inp
  real(dp), intent(out) :: nx, ny, nz
  integer :: k, iface, kk
  real(dp) :: fx, fy, fz, gm, fc

  fx = 0.0_dp
  fy = 0.0_dp
  fz = 0.0_dp
  do k=cf_ia(inp), cf_ia(inp+1)-1
    iface = cf_ja(k)
    if( iface <= numInnerFaces ) then
      kk = owner(iface)
      if( kk == inp ) kk = neighbour(iface)
      fc = 0.5_dp*( vof(inp) + vof(kk) )
      if( owner(iface) == inp ) then
        fx = fx + fc*arx(iface); fy = fy + fc*ary(iface); fz = fz + fc*arz(iface)
      else
        fx = fx - fc*arx(iface); fy = fy - fc*ary(iface); fz = fz - fc*arz(iface)
      endif
    endif
  enddo
  gm = sqrt( fx**2 + fy**2 + fz**2 )
  if( gm > small ) then
    nx = fx/gm; ny = fy/gm; nz = fz/gm
  else
    nx = 0.0_dp; ny = 0.0_dp; nz = 1.0_dp
  endif

end subroutine


!***********************************************************************
!
subroutine advect_vof( q, dt )
!
!***********************************************************************
!
!  Unsplit geometric advection of the volume fraction.
!
!    vof^{n+1} = vof^n - 1/V sum_f q1_f + vof^n dt/V sum_f Q_f
!
!  where q1_f is the volume of phase 1 which leaves the cell through the face f during dt, obtained by clipping the swept prism
!  of the face by the plane of the donor cell. The last term (zero for a divergence-free flux) keeps vof = 0 and vof = 1
!  invariant under the discrete divergence of the flux.
!
  implicit none
  real(dp), dimension(numFaces), intent(in) :: q
  real(dp), intent(in) :: dt

  integer :: i, ijp, ijn, donor, iface
  real(dp) :: qq, frac, q1, vtot, vnew
  real(dp), dimension(:), allocatable :: vold, divq, fnew

  allocate( vold(numCells), divq(numCells), fnew(numCells) )
  vold = vof(1:numCells)
  fnew = vold
  divq = 0.0_dp

  do i=1,numInnerFaces
    qq = q(i)
    ijp = owner(i)
    ijn = neighbour(i)
    divq(ijp) = divq(ijp) + qq
    divq(ijn) = divq(ijn) - qq
    if( abs(qq) < small ) cycle
    if( qq > 0.0_dp ) then
      donor = ijp
    else
      donor = ijn
    endif
    frac = donor_fraction( donor, i, dt )
    q1 = sign( frac*abs(qq)*dt, qq )
    q1acc(i) = q1acc(i) + q1
    fnew(ijp) = fnew(ijp) - q1/vol(ijp)
    fnew(ijn) = fnew(ijn) + q1/vol(ijn)
  enddo

  ! Open boundaries (only used by the pure advection tests; the flow solver supports closed domains)
  do iface=numInnerFaces+1,numFaces
    qq = q(iface)
    ijp = owner(iface)
    divq(ijp) = divq(ijp) + qq
    if( abs(qq) < small ) cycle
    if( qq > 0.0_dp ) then
      frac = donor_fraction( ijp, iface, dt )
      fnew(ijp) = fnew(ijp) - frac*qq*dt/vol(ijp)
      q1acc(iface) = q1acc(iface) + frac*qq*dt
    else
      fnew(ijp) = fnew(ijp) - cl_inflowF*qq*dt/vol(ijp)
      q1acc(iface) = q1acc(iface) + cl_inflowF*qq*dt
    endif
  enddo

  do i=1,numCells
    fnew(i) = fnew(i) + vold(i)*dt*divq(i)/vol(i)
  enddo

  ! Unsplit fluxes of the swept regions overlap/leave gaps near corners: values outside [0,1] are clipped and the volume
  ! defect is redistributed to the neighbour cells, so the total volume is conserved.
  call redistribute_defects( fnew )

  vof(1:numCells) = fnew
  call set_boundary_values( vof )

  deallocate( vold, divq, fnew )

end subroutine


function donor_fraction( donor, iface, dt ) result(frac)
!
!  Fraction of phase 1 in the volume swept through the face iface during dt by the flow which leaves the donor cell.
!  The swept region is the face polygon displaced backwards by u_f dt (a slanted prism, the Lagrangian flux region); its volume
!  is |S.u_f| dt = |Q| dt. It is clipped by the PLIC plane of the donor cell.
!
  implicit none
  integer, intent(in) :: donor, iface
  real(dp), intent(in) :: dt
  real(dp) :: frac

  integer :: nn, m, n1
  real(dp) :: fd, vpos, vsum
  real(dp), dimension(3,3) :: a, b
  real(dp), dimension(3,4) :: pt
  real(dp), dimension(4) :: st
  real(dp), dimension(3) :: disp

  fd = vof(donor)

  if( fd <= cl_fmin ) then
    frac = 0.0_dp
    return
  elseif( fd >= 1.0_dp-cl_fmin ) then
    frac = 1.0_dp
    return
  endif

  disp(1) = -fvx(iface)*dt
  disp(2) = -fvy(iface)*dt
  disp(3) = -fvz(iface)*dt

  nn = nnodes(iface)
  vpos = 0.0_dp
  vsum = 0.0_dp

  do m=1,nn

    a(:,1) = (/ xf(iface), yf(iface), zf(iface) /)
    n1 = node(m,iface)
    a(:,2) = (/ x(n1), y(n1), z(n1) /)
    n1 = node(mod(m,nn)+1,iface)
    a(:,3) = (/ x(n1), y(n1), z(n1) /)
    b(:,1) = a(:,1) + disp
    b(:,2) = a(:,2) + disp
    b(:,3) = a(:,3) + disp

    ! Prism a1 a2 a3 - b1 b2 b3 as three tetrahedra: (a1 a2 a3 b1), (a2 a3 b1 b2), (a3 b1 b2 b3)
    pt(:,1) = a(:,1); pt(:,2) = a(:,2); pt(:,3) = a(:,3); pt(:,4) = b(:,1)
    call plane_values( pt, donor, st )
    vpos = vpos + tet_pos_volume( pt, st )
    vsum = vsum + tet_volume( pt(:,1), pt(:,2), pt(:,3), pt(:,4) )

    pt(:,1) = a(:,2); pt(:,2) = a(:,3); pt(:,3) = b(:,1); pt(:,4) = b(:,2)
    call plane_values( pt, donor, st )
    vpos = vpos + tet_pos_volume( pt, st )
    vsum = vsum + tet_volume( pt(:,1), pt(:,2), pt(:,3), pt(:,4) )

    pt(:,1) = a(:,3); pt(:,2) = b(:,1); pt(:,3) = b(:,2); pt(:,4) = b(:,3)
    call plane_values( pt, donor, st )
    vpos = vpos + tet_pos_volume( pt, st )
    vsum = vsum + tet_volume( pt(:,1), pt(:,2), pt(:,3), pt(:,4) )

  enddo

  if( vsum > small ) then
    frac = min( 1.0_dp, max( 0.0_dp, vpos/vsum ) )
  else
    frac = fd
  endif

end function


subroutine face_velocities( q, t )
!
!  Velocity at face centers, with its normal component corrected to q/|S|.
!
  implicit none
  real(dp), dimension(numFaces), intent(in) :: q
  real(dp), intent(in) :: t
  integer :: i
  real(dp) :: lam, are, nx, ny, nz, un
  real(dp), dimension(3) :: uv

  do i=1,numFaces
    if( cl_solveFlow ) then
      if( i <= numInnerFaces ) then
        lam = facint(i)
        uv(1) = u(owner(i)) + lam*( u(neighbour(i)) - u(owner(i)) )
        uv(2) = v(owner(i)) + lam*( v(neighbour(i)) - v(owner(i)) )
        uv(3) = w(owner(i)) + lam*( w(neighbour(i)) - w(owner(i)) )
      else
        uv = 0.0_dp
      endif
    else
      call analytic_velocity( xf(i), yf(i), zf(i), t, uv )
    endif
    are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 )
    nx = arx(i)/are; ny = ary(i)/are; nz = arz(i)/are
    un = uv(1)*nx + uv(2)*ny + uv(3)*nz
    fvx(i) = uv(1) + ( q(i)/are - un )*nx
    fvy(i) = uv(2) + ( q(i)/are - un )*ny
    fvz(i) = uv(3) + ( q(i)/are - un )*nz
  enddo

end subroutine


subroutine analytic_velocity( xx, yy, zz, t, uv )
!
!  Velocity U = rot(A) of the prescribed flows by central differences of the vector potential.
!
  implicit none
  real(dp), intent(in) :: xx, yy, zz, t
  real(dp), dimension(3), intent(out) :: uv
  real(dp), parameter :: dl = 1.0e-5_dp
  real(dp), dimension(3) :: ap, am

  call vector_potential( xx, yy+dl, zz, t, ap );  call vector_potential( xx, yy-dl, zz, t, am )
  uv(1) = ( ap(3) - am(3) )/( 2*dl )
  call vector_potential( xx, yy, zz+dl, t, ap );  call vector_potential( xx, yy, zz-dl, t, am )
  uv(1) = uv(1) - ( ap(2) - am(2) )/( 2*dl )
  call vector_potential( xx, yy, zz+dl, t, ap );  call vector_potential( xx, yy, zz-dl, t, am )
  uv(2) = ( ap(1) - am(1) )/( 2*dl )
  call vector_potential( xx+dl, yy, zz, t, ap );  call vector_potential( xx-dl, yy, zz, t, am )
  uv(2) = uv(2) - ( ap(3) - am(3) )/( 2*dl )
  uv(3) = ( ap(2) - am(2) )/( 2*dl )
  call vector_potential( xx, yy+dl, zz, t, ap );  call vector_potential( xx, yy-dl, zz, t, am )
  uv(3) = uv(3) - ( ap(1) - am(1) )/( 2*dl )

end subroutine


subroutine redistribute_defects( f )
!
!  Clips f to [0,1] and gives the volume defect of the clipped cells to their face neighbours which can take it
!  (cells with room for phase 1 for an overshoot, cells with phase 1 for an undershoot). A few Jacobi passes;
!  what remains is accumulated in cl_clipVol (volume of phase 1 removed, negative: added).
!
  implicit none
  real(dp), dimension(numCells), intent(inout) :: f

  integer :: pass, inp, k, iface, kk
  real(dp) :: d, cap, tot, give
  real(dp), dimension(:), allocatable :: add
  logical :: any_defect

  allocate( add(numCells) )

  do pass=1,6

    add = 0.0_dp
    any_defect = .false.

    do inp=1,numCells
      if( f(inp) >= 0.0_dp .and. f(inp) <= 1.0_dp ) cycle
      any_defect = .true.
      if( f(inp) > 1.0_dp ) then
        d = ( f(inp) - 1.0_dp )*vol(inp)        ! phase 1 volume to place elsewhere
        f(inp) = 1.0_dp
      else
        d = f(inp)*vol(inp)                      ! negative: phase 1 volume to take away elsewhere
        f(inp) = 0.0_dp
      endif
      tot = 0.0_dp
      do k=cf_ia(inp), cf_ia(inp+1)-1
        iface = cf_ja(k)
        if( iface > numInnerFaces ) cycle
        kk = owner(iface)
        if( kk == inp ) kk = neighbour(iface)
        if( d > 0.0_dp ) then
          tot = tot + max( 0.0_dp, 1.0_dp - f(kk) )*vol(kk)
        else
          tot = tot + max( 0.0_dp, f(kk) )*vol(kk)
        endif
      enddo
      if( tot <= small ) then
        cl_clipVol = cl_clipVol + d
        cycle
      endif
      give = min( abs(d), tot )
      cl_clipVol = cl_clipVol + d - sign( give, d )
      do k=cf_ia(inp), cf_ia(inp+1)-1
        iface = cf_ja(k)
        if( iface > numInnerFaces ) cycle
        kk = owner(iface)
        if( kk == inp ) kk = neighbour(iface)
        if( d > 0.0_dp ) then
          cap = max( 0.0_dp, 1.0_dp - f(kk) )*vol(kk)
        else
          cap = max( 0.0_dp, f(kk) )*vol(kk)
        endif
        add(kk) = add(kk) + sign( give*cap/tot, d )/vol(kk)
      enddo
    enddo

    f = f + add
    if( .not. any_defect ) exit

  enddo

  ! Final clip of what is left
  do inp=1,numCells
    if( f(inp) > 1.0_dp ) then
      cl_clipVol = cl_clipVol + ( f(inp) - 1.0_dp )*vol(inp)
      f(inp) = 1.0_dp
    elseif( f(inp) < 0.0_dp ) then
      cl_clipVol = cl_clipVol + f(inp)*vol(inp)
      f(inp) = 0.0_dp
    endif
    if( f(inp) < cl_fmin ) then
      cl_clipVol = cl_clipVol + f(inp)*vol(inp)
      f(inp) = 0.0_dp
    elseif( f(inp) > 1.0_dp - cl_fmin ) then
      cl_clipVol = cl_clipVol - ( 1.0_dp - f(inp) )*vol(inp)
      f(inp) = 1.0_dp
    endif
  enddo

  deallocate( add )

end subroutine


subroutine plane_values( pt, inp, st )
  implicit none
  real(dp), dimension(3,4), intent(in) :: pt
  integer, intent(in) :: inp
  real(dp), dimension(4), intent(out) :: st
  integer :: k

  do k=1,4
    st(k) = pnx(inp)*pt(1,k) + pny(inp)*pt(2,k) + pnz(inp)*pt(3,k) - pd(inp)
  enddo

end subroutine


!***********************************************************************
!
subroutine advect_lsf( q, dt )
!
!***********************************************************************
!
!  Third order TVD Runge-Kutta (Shu & Osher) for   lsf_t + u.grad(lsf) = 0.
!
  implicit none
  real(dp), dimension(numFaces), intent(in) :: q
  real(dp), intent(in) :: dt

  real(dp), dimension(:), allocatable :: p0, p1, p2, rhs

  allocate( p0(numTotal), p1(numTotal), p2(numTotal), rhs(numCells) )

  p0 = lsf

  call lsf_rhs( p0, q, rhs )
  p1 = p0
  p1(1:numCells) = p0(1:numCells) + dt*rhs
  call set_boundary_values( p1 )

  call lsf_rhs( p1, q, rhs )
  p2 = p0
  p2(1:numCells) = 0.75_dp*p0(1:numCells) + 0.25_dp*( p1(1:numCells) + dt*rhs )
  call set_boundary_values( p2 )

  call lsf_rhs( p2, q, rhs )
  lsf(1:numCells) = p0(1:numCells)/3.0_dp + 2.0_dp/3.0_dp*( p2(1:numCells) + dt*rhs )
  call set_boundary_values( lsf )

  deallocate( p0, p1, p2, rhs )

end subroutine


subroutine lsf_rhs( phi, q, rhs )
!
!  rhs = - 1/V sum_f Q_f ( phi_f - phi_P ); phi_f from the upwind cell C with the TVD (van Leer) limiter of Jasak et al.:
!  phi_f = phi_C + psi(r)/2 (phi_D - phi_C),  r = 2 grad(phi)_C.d_CD / (phi_D - phi_C) - 1.
!
  implicit none
  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(numFaces), intent(in) :: q
  real(dp), dimension(numCells), intent(out) :: rhs

  integer :: i, ijp, ijn, ic, id, iface
  real(dp) :: qq, dphi, gd, r, psi, pf, cap

  call grad( phi, gphi )

  rhs = 0.0_dp

  do i=1,numInnerFaces
    qq = q(i)
    ijp = owner(i)
    ijn = neighbour(i)
    if( qq >= 0.0_dp ) then
      ic = ijp
      id = ijn
    else
      ic = ijn
      id = ijp
    endif
    dphi = phi(id) - phi(ic)
    if( abs(dphi) < 1.0e-14_dp ) then
      pf = phi(ic)
    else
      gd = gphi(1,ic)*(xc(id)-xc(ic)) + gphi(2,ic)*(yc(id)-yc(ic)) + gphi(3,ic)*(zc(id)-zc(ic))
      r = 2.0_dp*gd/dphi - 1.0_dp
      psi = ( r + abs(r) )/( 1.0_dp + abs(r) )
      pf = phi(ic) + 0.5_dp*psi*dphi
    endif
    rhs(ijp) = rhs(ijp) - qq*( pf - phi(ijp) )/vol(ijp)
    rhs(ijn) = rhs(ijn) + qq*( pf - phi(ijn) )/vol(ijn)
  enddo

  ! Inflow through open boundaries (pure advection tests): value of phase 2
  do iface=numInnerFaces+1,numFaces
    qq = q(iface)
    ijp = owner(iface)
    if( qq < 0.0_dp ) then
      cap = -( cl_band + 1.5_dp )*hloc(ijp)
      if( cl_inflowF > 0.5_dp ) cap = -cap
      rhs(ijp) = rhs(ijp) - qq*( cap - phi(ijp) )/vol(ijp)
    endif
  enddo

end subroutine


!***********************************************************************
!
subroutine update_properties
!
!***********************************************************************
!
!  Smoothed Heaviside function, density, viscosity, curvature, body force and the balanced-force face correction.
!
  implicit none

  integer :: i, ib, ijp, ijn, inp, ijb, iw, iface
  real(dp) :: eps, lam, gx, gy, gz, gm, flux, hf, kf, denf, fcomp, fint
  real(dp) :: dx, dy, dz, fxm, fym, fzm
  real(dp), dimension(:,:), allocatable :: gh, nrm

  allocate( gh(3,numCells), nrm(3,numTotal) )

  call set_boundary_values( lsf )

  do inp=1,numCells
    eps = cl_eps*hloc(inp)
    hs(inp) = heaviside( lsf(inp), eps )
    ! Density and viscosity from the volume fraction: the density is then consistent with the mass flux of the momentum
    ! equations, which is built from the same volume fluxes of phase 1 (advection_mass_flux).
    den(inp) = cl_rho2 + ( cl_rho1 - cl_rho2 )*vof(inp)
    if( cl_harmonic ) then
      vis(inp) = 1.0_dp/( vof(inp)/cl_mu1 + (1.0_dp-vof(inp))/cl_mu2 )
    else
      vis(inp) = cl_mu2 + ( cl_mu1 - cl_mu2 )*vof(inp)
    endif
  enddo

  call set_boundary_values( hs )
  call set_boundary_values( den )
  call set_boundary_values( vis )

  ! Wall viscosity, in the order of the wall faces used by the momentum equations
  iw = 0
  do ib=1,numBoundaries
    if( bctype(ib) == 'wall' ) then
      do i=1,nfaces(ib)
        iw = iw + 1
        visw(iw) = vis( owner(startFace(ib)+i) )
      enddo
    endif
  enddo

  ! Curvature  kappa = - div( n ),  n = grad(lsf)/|grad(lsf)|
  call grad( lsf, gphi )

  do inp=1,numCells
    gm = sqrt( sum( gphi(:,inp)**2 ) )
    if( gm > 1.0e-8_dp ) then
      nrm(:,inp) = gphi(:,inp)/gm
    else
      nrm(:,inp) = 0.0_dp
    endif
  enddo
  do i=1,numBoundaryFaces
    nrm(:,numCells+i) = nrm(:,owner(numInnerFaces+i))
  enddo

  kappa(1:numCells) = 0.0_dp
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)
    gx = nrm(1,ijp) + lam*( nrm(1,ijn) - nrm(1,ijp) )
    gy = nrm(2,ijp) + lam*( nrm(2,ijn) - nrm(2,ijp) )
    gz = nrm(3,ijp) + lam*( nrm(3,ijn) - nrm(3,ijp) )
    gm = sqrt( gx**2 + gy**2 + gz**2 )
    if( gm > small ) then
      gx = gx/gm; gy = gy/gm; gz = gz/gm
    endif
    flux = gx*arx(i) + gy*ary(i) + gz*arz(i)
    kappa(ijp) = kappa(ijp) - flux
    kappa(ijn) = kappa(ijn) + flux
  enddo
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    kappa(ijp) = kappa(ijp) - ( nrm(1,ijp)*arx(iface) + nrm(2,ijp)*ary(iface) + nrm(3,ijp)*arz(iface) )
  enddo
  do inp=1,numCells
    kappa(inp) = kappa(inp)/vol(inp)
    ! Under-resolved interface (filaments, drops of a size of a few cells): limit the curvature
    kappa(inp) = max( -2.0_dp/hloc(inp), min( 2.0_dp/hloc(inp), kappa(inp) ) )
  enddo
  call set_boundary_values( kappa )

  ! Gradient of the Heaviside function (Gauss)
  gh = 0.0_dp
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    hf = hs(ijp) + facint(i)*( hs(ijn) - hs(ijp) )
    gh(1,ijp) = gh(1,ijp) + hf*arx(i);  gh(1,ijn) = gh(1,ijn) - hf*arx(i)
    gh(2,ijp) = gh(2,ijp) + hf*ary(i);  gh(2,ijn) = gh(2,ijn) - hf*ary(i)
    gh(3,ijp) = gh(3,ijp) + hf*arz(i);  gh(3,ijn) = gh(3,ijn) - hf*arz(i)
  enddo
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    hf = hs(ijp)
    gh(1,ijp) = gh(1,ijp) + hf*arx(iface)
    gh(2,ijp) = gh(2,ijp) + hf*ary(iface)
    gh(3,ijp) = gh(3,ijp) + hf*arz(iface)
  enddo

  ! Force per unit volume: CSF surface tension + gravity with the reference density densit
  do inp=1,numCells
    fst(1,inp) = cl_sigma*kappa(inp)*gh(1,inp)/vol(inp) + ( den(inp) - densit )*gravx
    fst(2,inp) = cl_sigma*kappa(inp)*gh(2,inp)/vol(inp) + ( den(inp) - densit )*gravy
    fst(3,inp) = cl_sigma*kappa(inp)*gh(3,inp)/vol(inp) + ( den(inp) - densit )*gravz
  enddo

  ! Balanced-force correction: difference between the compact (face) form of the force and the interpolated cell forces,
  ! projected on the line between the cell centers; it enters the Rhie-Chow flux as an additional pressure difference.
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)
    dx = xc(ijn) - xc(ijp)
    dy = yc(ijn) - yc(ijp)
    dz = zc(ijn) - zc(ijp)
    kf = kappa(ijp) + lam*( kappa(ijn) - kappa(ijp) )
    denf = den(ijp) + lam*( den(ijn) - den(ijp) )
    fcomp = cl_sigma*kf*( hs(ijn) - hs(ijp) ) + ( denf - densit )*( gravx*dx + gravy*dy + gravz*dz )
    fxm = fst(1,ijp) + lam*( fst(1,ijn) - fst(1,ijp) )
    fym = fst(2,ijp) + lam*( fst(2,ijn) - fst(2,ijp) )
    fzm = fst(3,ijp) + lam*( fst(3,ijn) - fst(3,ijp) )
    fint = fxm*dx + fym*dy + fzm*dz
    fbal(i) = fcomp - fint
  enddo

  deallocate( gh, nrm )

end subroutine


pure function heaviside( phi, eps ) result(hv)
  implicit none
  real(dp), intent(in) :: phi, eps
  real(dp) :: hv

  if( phi < -eps ) then
    hv = 0.0_dp
  elseif( phi > eps ) then
    hv = 1.0_dp
  else
    hv = 0.5_dp*( 1.0_dp + phi/eps + sin(pi*phi/eps)/pi )
  endif

end function


subroutine init_mass_flux
!
!  Initial mass flux of the momentum equations from the volume flux: flmass = rho_f Q_f.
!
  use variables, only: flmass
  implicit none
  integer :: i, iface

  do i=1,numInnerFaces
    flmass(i) = ( den(owner(i)) + facint(i)*( den(neighbour(i)) - den(owner(i)) ) )*flvol(i)
  enddo
  do iface=numInnerFaces+1,numFaces
    flmass(iface) = den(owner(iface))*flvol(iface)
  enddo

end subroutine


subroutine advection_mass_flux( q, dt )
!
!  Mass flux of the momentum equations consistent with the advection of the volume fraction:
!
!    F_f = rho2 Q_f + (rho1 - rho2) q1_f / dt,
!
!  where q1_f is the volume of phase 1 which crossed the face in the step. The density rho = rho2 + (rho1-rho2) vof then
!  satisfies the discrete mass balance  (rho^{n+1} - rho^n) V + dt sum_f F_f = 0  and the momentum equations are
!  consistent with it (this is essential at high density ratios). The flux is kept constant in the outer iterations.
!
  use variables, only: flmass
  implicit none
  real(dp), dimension(numFaces), intent(in) :: q
  real(dp), intent(in) :: dt
  integer :: i

  do i=1,numFaces
    flmass(i) = cl_rho2*q(i) + ( cl_rho1 - cl_rho2 )*q1acc(i)/dt
  enddo

end subroutine


!***********************************************************************
!
! PRESCRIBED, DIVERGENCE FREE FLOW FIELDS FOR THE ADVECTION TESTS
!
!***********************************************************************

subroutine vector_potential( xx, yy, zz, t, aa )
!
!  Velocity is the curl of the vector potential A, U = rot(A), so the face fluxes obtained as the circulation of A
!  along the face edges are exactly divergence free for every cell, whatever the quadrature.
!
  implicit none
  real(dp), intent(in) :: xx, yy, zz, t
  real(dp), dimension(3), intent(out) :: aa
  real(dp) :: psi

  select case( trim(cl_flow) )

  case( 'translate' )
    ! A = 1/2 U x r
    aa(1) = 0.5_dp*( cl_vel(2)*zz - cl_vel(3)*yy )
    aa(2) = 0.5_dp*( cl_vel(3)*xx - cl_vel(1)*zz )
    aa(3) = 0.5_dp*( cl_vel(1)*yy - cl_vel(2)*xx )

  case( 'rotate' )
    ! U = omega (-(y-yc), (x-xc), 0)
    aa(1) = 0.0_dp
    aa(2) = 0.0_dp
    aa(3) = -0.5_dp*cl_omega*( (xx-cl_vcenter(1))**2 + (yy-cl_vcenter(2))**2 )

  case( 'vortex' )
    ! LeVeque (1996) vortex in a box, reversing at t = T/2:  psi = sin^2(pi x) sin^2(pi y) cos(pi t/T)/pi
    psi = sin(pi*xx)**2*sin(pi*yy)**2*cos(pi*t/cl_T)/pi
    aa(1) = 0.0_dp
    aa(2) = 0.0_dp
    aa(3) = psi

  case default
    write(6,'(3a)') ' *** CLSVOF: unknown cl_flow = ', trim(cl_flow), ' (solve, translate, rotate, vortex)'
    stop

  end select

end subroutine


subroutine analytic_flux( t, q )
  implicit none
  real(dp), intent(in) :: t
  real(dp), dimension(numFaces), intent(out) :: q

  integer :: iface, m, nn, n1, n2, ig
  real(dp), parameter, dimension(3) :: gp = (/ 0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp /)
  real(dp), parameter, dimension(3) :: gw = (/ 0.2777777777777778_dp, 0.4444444444444444_dp, 0.2777777777777778_dp /)
  real(dp) :: circ, xx, yy, zz, ex, ey, ez
  real(dp), dimension(3) :: aa, sn

  do iface=1,numFaces

    nn = nnodes(iface)
    circ = 0.0_dp
    sn = 0.0_dp

    do m=1,nn
      n1 = node(m,iface)
      n2 = node(mod(m,nn)+1,iface)
      ex = x(n2) - x(n1)
      ey = y(n2) - y(n1)
      ez = z(n2) - z(n1)
      do ig=1,3
        xx = x(n1) + gp(ig)*ex
        yy = y(n1) + gp(ig)*ey
        zz = z(n1) + gp(ig)*ez
        call vector_potential( xx, yy, zz, t, aa )
        circ = circ + gw(ig)*( aa(1)*ex + aa(2)*ey + aa(3)*ez )
      enddo
      ! Area vector of the polygon in the order of the nodes, 1/2 sum r_k x r_k+1
      sn(1) = sn(1) + 0.5_dp*( y(n1)*z(n2) - z(n1)*y(n2) )
      sn(2) = sn(2) + 0.5_dp*( z(n1)*x(n2) - x(n1)*z(n2) )
      sn(3) = sn(3) + 0.5_dp*( x(n1)*y(n2) - y(n1)*x(n2) )
    enddo

    ! Circulation is the flux through the face for the orientation given by the order of the nodes
    if( sn(1)*arx(iface) + sn(2)*ary(iface) + sn(3)*arz(iface) < 0.0_dp ) circ = -circ
    q(iface) = circ

  enddo

end subroutine


!***********************************************************************
!
! MONITORING AND RESTART
!
!***********************************************************************

subroutine clsvof_monitor
!
!  One line of clsvof_history.dat per time step: volume of phase 1, its change, clipped volume,
!  centroids and mean velocities of both phases, max. velocity, interface area, continuity error.
!
  implicit none

  integer :: inp, j
  real(dp) :: v1, v2, w1, w2, c1(3), c2(3), u1(3), u2(3), umax, area, vf, m1, m2

  v1 = 0.0_dp; v2 = 0.0_dp
  c1 = 0.0_dp; c2 = 0.0_dp; u1 = 0.0_dp; u2 = 0.0_dp
  umax = 0.0_dp

  do inp=1,numCells
    w1 = vof(inp)*vol(inp)
    w2 = ( 1.0_dp - vof(inp) )*vol(inp)
    v1 = v1 + w1
    v2 = v2 + w2
    c1 = c1 + w1*(/ xc(inp), yc(inp), zc(inp) /)
    c2 = c2 + w2*(/ xc(inp), yc(inp), zc(inp) /)
    u1 = u1 + w1*(/ u(inp), v(inp), w(inp) /)
    u2 = u2 + w2*(/ u(inp), v(inp), w(inp) /)
    umax = max( umax, sqrt( u(inp)**2 + v(inp)**2 + w(inp)**2 ) )
  enddo

  if( v1 > small ) then
    c1 = c1/v1
    u1 = u1/v1
  endif
  if( v2 > small ) then
    c2 = c2/v2
    u2 = u2/v2
  endif

  area = 0.0_dp
  do j=1,nmixed
    inp = mixed(j)
    area = area + plane_area( inp, pnx(inp), pny(inp), pnz(inp), pd(inp) )
  enddo

  write(cl_unit,'(es15.7,2es17.9,es12.4,12es15.7,es14.6,es14.6,es11.3,i7)') &
    cl_time, v1, (v1-cl_volume0)/max(cl_volume0,small), cl_clipVol, c1, c2, u1, u2, umax, area, cl_contErr, nmixed
  flush(cl_unit)

end subroutine


subroutine clsvof_write_restart
  implicit none
  integer :: iu

  if( .not. twoPhase ) return

  open( newunit=iu, file='clsvof.restart', form='unformatted', access='stream', status='replace' )
  write(iu) cl_time, dt_old, cl_volume0, cl_clipVol
  write(iu) vof(1:numCells), lsf(1:numCells), flvol, flvolo
  close(iu)

end subroutine


subroutine clsvof_read_restart
  implicit none
  integer :: iu, ios
  logical :: ex

  inquire( file='clsvof.restart', exist=ex )
  if( .not. ex ) then
    write(6,'(a)') ' *** CLSVOF: lread = T but the file clsvof.restart does not exist!'
    stop
  endif

  open( newunit=iu, file='clsvof.restart', form='unformatted', access='stream', status='old' )
  read(iu) cl_time, dt_old, cl_volume0, cl_clipVol
  read(iu) vof(1:numCells), lsf(1:numCells), flvol, flvolo
  close(iu)

  first_step = .false.
  call set_boundary_values( vof )
  call set_boundary_values( lsf )
  call list_mixed_cells
  call reconstruct_planes

end subroutine


end module clsvof
