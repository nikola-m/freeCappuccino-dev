module moving_mesh
!
!********************************************************************************************************************************132
!
!  Purpose:
!
!    Dynamic (moving and deforming) mesh capability: mesh deformation with radial basis functions (RBF), space conservation law
!    (GCL) fluxes of the faces, recalculation of the finite volume geometry, and the corrections that the flow solver needs
!    on a moving mesh (relative mass fluxes, volume change source of the pressure correction equation).
!
!  Discussion:
!
!    The fluid part of the solver needs to know how the mesh moves, whatever the way the nodes are moved. The data of this module
!    are the nodal coordinates and the cell volumes of the previous time level (xo,yo,zo,volo), the volume swept by every face in
!    the time step (flux_grid) and the corresponding mass flux of the moving grid (flmass_grid).
!
!    One time step of the mesh (subroutine mm_advance, called from the time loop after the time level is shifted) is:
!
!      1. Old nodal coordinates and cell volumes are stored by swapping: xo = x, volo = vol.
!      2. New nodal coordinates by mesh deformation. The displacement of the control points (nodes of the boundary patches:
!         prescribed motion on the moving patches, zero on the others) is interpolated to all nodes with RBF. The displacement is
!         always measured from the initial mesh (xini,yini,zini) so the mesh is a function of time only, which makes restart simple.
!      3. Face fluxes of the moving grid, flux_grid = volume swept by the face in the time step, are calculated exactly for the
!         linear motion of the nodes during the step. The new cell volumes follow as the sum of the fluxes of the faces,
!         vol = volo + sum_f flux_grid (outward from the cell). This is the space conservation law (GCL) imposed in the
!         conservative way: the volume of a cell is its old volume plus what its faces have swept, and it agrees with the volume
!         calculated from the new nodal coordinates to the round-off error (see mm_update_geometry, reported in every time step).
!      4. Mass fluxes of the grid motion, flmass_grid = den/dt * flux_grid.
!      5. Geometry data of the finite volume discretization are recalculated: xc,yc,zc, xf,yf,zf, arx,ary,arz, facint, Df, and
!         the wall and symmetry distances dnw, srdw, dns, srds.
!      6. The mass fluxes through the inner faces are made relative to the moving grid, flmass = flmass - flmass_grid. The
!         solver keeps flmass relative to the mesh of the current time level (it is made absolute again before the next move),
!         so it enters the convective terms of all transport equations without modifications.
!      7. The source of the pressure(-correction) equation gets the volume change of the cells, see mm_pressure_source:
!                su(inp) = su(inp) + den(inp)/dt*( volo(inp) - vol(inp) ),
!         which makes the absolute mass flux (relative flux plus the flux of the grid) divergence free.
!
!    The unsteady term of the transport equations on a moving mesh is den*volo/dt*(phi - phio) when the convection is written as
!    sum_f F_f (phi_f - phi_P), as it is in this code. It is the exact first order form of d(rho V phi)/dt + sum_f F_rel phi_f.
!    Use the function vol_unsteady(inp) in place of vol(inp) in the first order unsteady term.
!
!    Supported: laminar incompressible flow, first order implicit time stepping (bdf = T), SIMPLE and PISO, Gauss gradients,
!    momentum, temperature and passive scalar equations. Everything else stops with a message (see init.f90).
!    The moving patches have to be of type 'wall' (the velocity of the fluid at the wall is the velocity of the moving face),
!    all other boundaries should be at rest.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    Rendall T.C.S., Allen C.B., Unified fluid-structure interpolation and mesh motion using radial basis functions,
!      Int. J. Numer. Meth. Engng. 74 (2008) 1519-1559.
!    Demirdzic I., Peric M., Space conservation law in finite volume calculations of fluid flow,
!      Int. J. Numer. Meth. Fluids 8 (1988) 1037-1050.
!    Ferziger J.H., Peric M., Computational Methods for Fluid Dynamics, 3rd ed., Springer, 2002, Sec. 12.2 (moving grids).
!    Wendland H., Piecewise polynomial, positive definite and compactly supported radial functions of minimal degree,
!      Adv. Comput. Math. 4 (1995) 389-396.
!
  use types
  use parameters
  use geometry
  use variables, only: flmass, den, u, v, w

  implicit none

  !
  ! Parameters - modified trough input.nml file
  !
  logical  :: DynamicMesh = .false.                      ! Turn on/off dynamic mesh flow solution.
  character( len=12 ) :: mm_motion = 'oscillate'         ! Motion of the moving patches: oscillate, translate, rotate, user.
  character( len=30 ), dimension(8) :: mm_patch = ''     ! Names of the patches (type wall) that move.
  real(dp) :: mm_amp(3) = 0.0_dp                         ! oscillate: amplitude vector, d = amp*sin(2*pi*freq*t).
  real(dp) :: mm_freq = 1.0_dp                           ! oscillate: frequency [1/s].
  real(dp) :: mm_vel(3) = 0.0_dp                         ! translate: velocity vector, d = vel*t.
  real(dp) :: mm_center(3) = 0.0_dp                      ! rotate: point on the axis of rotation.
  real(dp) :: mm_axis(3) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)  ! rotate: direction of the axis of rotation.
  real(dp) :: mm_omega = 0.0_dp                          ! rotate: angular velocity [rad/s], angle = omega*t.
  character( len=20 ) :: mm_rbf = 'wendland'             ! RBF: wendland, gaussian, multiquadric, inverse_multiquadric, thin_plate.
  real(dp) :: mm_radius = 0.0_dp                         ! wendland: support radius (<= 0: bounding box diagonal).
  real(dp) :: mm_epsilon = 1.0_dp                        ! Shape parameter of gaussian, multiquadric, inverse_multiquadric [1/m].

  ! Fields
  real(dp), dimension(:), allocatable :: xo, yo, zo       ! Nodal coordinates of the previous time level, [numNodes]
  real(dp), dimension(:), allocatable :: volo             ! Cell volumes of the previous time level, [numCells]
  real(dp), dimension(:), allocatable :: flux_grid        ! Volume swept by a face in the time step (along area vector), [numFaces]
  real(dp), dimension(:), allocatable :: flmass_grid      ! Mass flux of the grid motion, den/dt*flux_grid, [numFaces]

  ! Private data
  real(dp), dimension(:), allocatable, private :: xini, yini, zini        ! Nodal coordinates of the initial mesh
  real(dp), dimension(:), allocatable, private :: vgeo, r8tmp             ! Cell volume from the nodal coordinates, work array
  integer, dimension(:), allocatable, private :: ictrl                    ! Index in the list of control points of a node, 0 if free
  integer, dimension(:), allocatable, private :: ctrl                     ! Node numbers of the control points
  logical, dimension(:), allocatable, private :: cmoving                  ! Control point belongs to a moving patch
  real(dp), dimension(:), allocatable, private :: cdx, cdy, cdz           ! Prescribed displacements of the control points
  real(dp), dimension(:), allocatable, private :: alx, aly, alz           ! RBF weights
  real(dp), dimension(:,:), allocatable, private :: lu                    ! LU factors of the RBF matrix
  integer, dimension(:), allocatable, private :: piv                      ! Pivots of the LU factorization
  integer, dimension(:), allocatable, private :: mv_face, mv_ijb          ! Faces of the moving patches, their boundary value index
  real(dp), dimension(:), allocatable, private :: xfo, yfo, zfo           ! Old centers of the faces of the moving patches
  logical, dimension(:), allocatable, private :: patch_moves              ! Boundary patch is one of the moving patches
  integer, private :: nctrl = 0, nmvface = 0, rbf_kind = 1
  real(dp), private :: rbf_R = 1.0_dp
  real(dp), dimension(3), private :: axis_n = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
  logical, private :: warned_fixed = .false.

  private

  public :: DynamicMesh, mm_motion, mm_patch, mm_amp, mm_freq, mm_vel, mm_center, mm_axis, mm_omega, mm_rbf, mm_radius, mm_epsilon
  public :: xo, yo, zo, volo, flux_grid, flmass_grid
  public :: mm_init, mm_advance, mm_pressure_source, mm_unsupported, mm_write_points, vol_unsteady


contains


!***********************************************************************
!
subroutine mm_unsupported( what )
!
!***********************************************************************
!
!  Stops the program with a message: the feature 'what' cannot be combined with the moving mesh.
!
  implicit none

  character( len=* ), intent(in) :: what

  write(*,'(/,a)')   ' *** Moving mesh (DynamicMesh = T) cannot be used with: '//trim(what)//'.'
  write(*,'(a,/)')   '     Supported: laminar incompressible flow, bdf = T, SIMPLE or PISO, gauss = T, no periodic patches.'
  write(6,'(/,a)')   ' *** Moving mesh (DynamicMesh = T) cannot be used with: '//trim(what)//'.'
  write(6,'(a,/)')   '     Supported: laminar incompressible flow, bdf = T, SIMPLE or PISO, gauss = T, no periodic patches.'
  stop

end subroutine mm_unsupported


!***********************************************************************
!
subroutine mm_init
!
!***********************************************************************
!
!  Checks the setup, allocates the data, finds the control points and the faces of the moving patches, factorizes the RBF matrix,
!  and, for a restarted simulation, moves the mesh to the time level of the restart file.
!
  implicit none

  integer :: ib, i, k, iface, ivertex, inode, ip, nmoving_nodes, nshared
  integer, dimension(:), allocatable :: nodetype  ! 0 free, 1 fixed control, 2 moving control
  logical :: found
  real(dp) :: bbmin(3), bbmax(3), diag, s
  real(dp) :: mem_mb

  if( .not. DynamicMesh ) return

  !
  ! Setup that cannot be combined with the moving mesh
  !
  if( .not. ltransient )                 call mm_unsupported('steady simulations (ltransient = T is required)')
  if( .not. bdf )                        call mm_unsupported('time stepping other than bdf = T (first order backward Euler)')
  if( bdf2 .or. bdf3 .or. cn )           call mm_unsupported('bdf2, bdf3 or Crank-Nicolson (use bdf = T only)')
  if( lturb )                            call mm_unsupported('turbulence models (TurbModel = none is required)')
  if( compressible )                     call mm_unsupported('compressible flow')
  if( nPer > 0 )                         call mm_unsupported('periodic boundaries')
  if( .not. (SIMPLE .or. PISO) )         call mm_unsupported('pressure-velocity coupling other than SIMPLE or PISO')

  select case( trim(mm_rbf) )
    case('wendland');             rbf_kind = 1
    case('gaussian');             rbf_kind = 2
    case('multiquadric');         rbf_kind = 3
    case('inverse_multiquadric'); rbf_kind = 4
    case('thin_plate');           rbf_kind = 5
    case default
      write(*,'(a)') ' *** Moving mesh: unknown RBF type '//trim(mm_rbf)// &
                     ' (wendland, gaussian, multiquadric, inverse_multiquadric, thin_plate).'
      stop
  end select

  select case( trim(mm_motion) )
    case('oscillate','translate','user')
    case('rotate')
      s = sqrt( sum( mm_axis**2 ) )
      if( s < small ) then
        write(*,'(a)') ' *** Moving mesh: the axis of rotation mm_axis has zero length.'
        stop
      endif
      axis_n = mm_axis / s
    case default
      write(*,'(a)') ' *** Moving mesh: unknown motion '//trim(mm_motion)//' (oscillate, translate, rotate, user).'
      stop
  end select

  write(*,'(a)') ' '
  write(*,'(a)') '  Initializing moving mesh...'

  !
  ! Allocate and initialize the data
  !
  allocate( xo(numNodes), yo(numNodes), zo(numNodes) )
  allocate( xini(numNodes), yini(numNodes), zini(numNodes) )
  allocate( volo(numCells), vgeo(numCells), r8tmp(numCells) )
  allocate( flux_grid(numFaces), flmass_grid(numFaces) )
  allocate( ictrl(numNodes), nodetype(numNodes) )

  xini = x;  yini = y;  zini = z
  xo = x;    yo = y;    zo = z
  volo = vol
  flux_grid = 0.0_dp
  flmass_grid = 0.0_dp
  ictrl = 0
  nodetype = 0

  !
  ! Nodes of the moving patches: moving control points. Nodes of all other patches (except 'empty' and 'symmetry', where the
  ! nodes slide along the plane) are fixed control points. Moving patches have priority over the fixed ones.
  !
  nmvface = 0
  allocate( patch_moves(numBoundaries) )
  patch_moves = .false.
  do ib=1,numBoundaries
    found = .false.
    do ip=1,size(mm_patch)
      if( len_trim(mm_patch(ip)) > 0 .and. trim(bcname(ib)) == trim(mm_patch(ip)) ) found = .true.
    enddo
    patch_moves(ib) = found
    if( found ) then
      if( bctype(ib) /= 'wall' ) then
        write(*,'(a)') ' *** Moving mesh: moving patch '//trim(bcname(ib))//' is not of type wall.'
        stop
      endif
      nmvface = nmvface + nfaces(ib)
    endif
  enddo

  do ip=1,size(mm_patch)
    if( len_trim(mm_patch(ip)) == 0 ) cycle
    found = .false.
    do ib=1,numBoundaries
      if( trim(bcname(ib)) == trim(mm_patch(ip)) ) found = .true.
    enddo
    if( .not. found ) then
      write(*,'(a)') ' *** Moving mesh: patch '//trim(mm_patch(ip))//' (mm_patch) does not exist in the mesh.'
      stop
    endif
  enddo

  if( nmvface == 0 ) then
    write(*,'(a)') ' *** Moving mesh: no moving patch, set mm_patch = ''name'', ...'
    stop
  endif

  allocate( mv_face(nmvface), mv_ijb(nmvface), xfo(nmvface), yfo(nmvface), zfo(nmvface) )

  ! First fixed control points
  do ib=1,numBoundaries
    if( bctype(ib) == 'empty' .or. bctype(ib) == 'symmetry' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      do ivertex=1,nnodes(iface)
        nodetype( node(ivertex,iface) ) = 1
      enddo
    enddo
  enddo

  ! Then moving control points, and the list of moving faces
  k = 0
  nshared = 0
  do ib=1,numBoundaries
    found = .false.
    do ip=1,size(mm_patch)
      if( len_trim(mm_patch(ip)) > 0 .and. trim(bcname(ib)) == trim(mm_patch(ip)) ) found = .true.
    enddo
    if( .not. found ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      k = k + 1
      mv_face(k) = iface
      mv_ijb(k) = iBndValueStart(ib) + i
      do ivertex=1,nnodes(iface)
        nodetype( node(ivertex,iface) ) = 2
      enddo
    enddo
  enddo

  ! Nodes shared between a moving patch and a patch of another kind (inlet, outlet, a wall at rest): they move with the patch.
  do ib=1,numBoundaries
    found = .false.
    do ip=1,size(mm_patch)
      if( len_trim(mm_patch(ip)) > 0 .and. trim(bcname(ib)) == trim(mm_patch(ip)) ) found = .true.
    enddo
    if( found .or. bctype(ib) == 'empty' .or. bctype(ib) == 'symmetry' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      do ivertex=1,nnodes(iface)
        if( nodetype( node(ivertex,iface) ) == 2 ) nshared = nshared + 1
      enddo
    enddo
  enddo

  nctrl = count( nodetype > 0 )
  nmoving_nodes = count( nodetype == 2 )

  allocate( ctrl(nctrl), cmoving(nctrl), cdx(nctrl), cdy(nctrl), cdz(nctrl), alx(nctrl), aly(nctrl), alz(nctrl) )

  k = 0
  do inode=1,numNodes
    if( nodetype(inode) > 0 ) then
      k = k + 1
      ctrl(k) = inode
      cmoving(k) = ( nodetype(inode) == 2 )
      ictrl(inode) = k
    endif
  enddo

  deallocate( nodetype )

  !
  ! RBF support radius: diagonal of the bounding box of the initial mesh, if not specified.
  !
  bbmin = (/ minval(x), minval(y), minval(z) /)
  bbmax = (/ maxval(x), maxval(y), maxval(z) /)
  diag = sqrt( sum( (bbmax-bbmin)**2 ) )
  rbf_R = mm_radius
  if( rbf_R <= 0.0_dp ) rbf_R = diag

  !
  ! RBF matrix of the control points (initial positions, so it is constant) and its LU factorization.
  !
  mem_mb = 8.0_dp*real(nctrl,dp)**2/1.0e6_dp
  write(*,'(a,i0,a,i0,a)')   '  Control points: ', nctrl, ' (', nmoving_nodes, ' on moving patches)'
  write(*,'(a,a,a,es10.3)')  '  RBF: ', trim(mm_rbf), ', support radius / shape parameter: ', merge(rbf_R,mm_epsilon,rbf_kind==1)
  write(*,'(a,f9.1,a)')      '  RBF matrix (dense LU): ', mem_mb, ' MB'
  if( nshared > 0 ) write(*,'(a,i0,a)') '  Note: ', nshared, &
       ' nodes of the moving patches are shared with patches at rest; they move with the moving patch.'

  allocate( lu(nctrl,nctrl), piv(nctrl) )
  call mm_rbf_matrix
  call lu_factor( nctrl, lu, piv, k )
  if( k /= 0 ) then
    write(*,'(a,i0)') ' *** Moving mesh: the RBF matrix is singular, zero pivot in column ', k
    write(*,'(a)')    '     Coincident control points, or a RBF that needs a polynomial: try mm_rbf = ''wendland''.'
    stop
  endif

  call mm_check_interpolation

  !
  ! Restarted simulation: the mesh at the time of the restart file (time) and the grid fluxes of the last time step.
  !
  if( lread .and. time > 0.0_dp ) then
    write(*,'(a,es10.3)') '  Restart: moving the mesh to the time level ', time
    ! Mesh of the previous time level, and its geometry. The volumes are taken from the nodal coordinates here.
    call mm_move_nodes( time - timestep )
    call mm_update_geometry
    vol = vgeo
    ! Standard time step to time, but the mass fluxes of the restart file are already relative to this mesh.
    call mm_step( time, timestep, .false. )
  endif

  write(*,'(a)') '  Moving mesh is initialized.'
  write(*,'(a)') ' '

end subroutine mm_init


!***********************************************************************
!
subroutine mm_advance
!
!***********************************************************************
!
!  One time step of the moving mesh. Called from the time loop after the shift of the time levels and after the time
!  has been advanced (time is the time of the new level).
!
  implicit none

  if( .not. DynamicMesh ) return

  call mm_step( time, timestep, .true. )

end subroutine mm_advance


!***********************************************************************
!
subroutine mm_step( tnew, dt, shift_flux )
!
!***********************************************************************
!
!  Steps 1-6 of the description in the header of the module.
!  If shift_flux is true, mass fluxes flmass are made absolute with the grid flux of the previous step before the mesh moves,
!  and relative to the new one afterwards.
!
  implicit none

  real(dp), intent(in) :: tnew, dt
  logical, intent(in) :: shift_flux

  integer :: i, k, iface, ijp, ijn, ijb
  real(dp) :: fx, fy, fz, are2, qn, vn, dmax, volmin, voltot, gclerr, fixedmax
  real(dp) :: denf

  ! mass fluxes become absolute again (the flux of the grid of the previous step is added back)
  if( shift_flux ) then
    do i=1,numInnerFaces
      flmass(i) = flmass(i) + flmass_grid(i)
    enddo
  endif

  !
  ! 1. Swap: old coordinates and volumes. Old face centers of the moving patches (wall velocity).
  !
  xo = x
  yo = y
  zo = z
  volo = vol

  do k=1,nmvface
    xfo(k) = xf(mv_face(k))
    yfo(k) = yf(mv_face(k))
    zfo(k) = zf(mv_face(k))
  enddo

  !
  ! 2. New nodal coordinates: RBF mesh deformation.
  !
  call mm_move_nodes( tnew )

  !
  ! 3. Volume swept by the faces, new cell volumes by the space conservation law.
  !
  call mm_swept_volumes

  vol = volo
  do iface=1,numInnerFaces
    ijp = owner(iface)
    ijn = neighbour(iface)
    vol(ijp) = vol(ijp) + flux_grid(iface)
    vol(ijn) = vol(ijn) - flux_grid(iface)
  enddo
  do iface=numInnerFaces+1,numFaces
    ijp = owner(iface)
    vol(ijp) = vol(ijp) + flux_grid(iface)
  enddo

  volmin = minval( vol )
  if( volmin <= 0.0_dp ) then
    write(*,'(a,es10.3,a)') ' *** Moving mesh: negative cell volume (', volmin, ')!'
    write(*,'(a)')          '     Reduce the motion per time step or the amplitude, or increase the RBF support radius.'
    write(6,'(a,es10.3,a)') ' *** Moving mesh: negative cell volume (', volmin, ')!'
    write(6,'(a)')          '     Reduce the motion per time step or the amplitude, or increase the RBF support radius.'
    stop
  endif

  !
  ! 4. Mass fluxes of the grid motion: den/dt * swept volume.
  !
  do iface=1,numInnerFaces
    ijp = owner(iface)
    ijn = neighbour(iface)
    denf = den(ijp) + facint(iface)*( den(ijn) - den(ijp) )
    flmass_grid(iface) = denf/dt * flux_grid(iface)
  enddo
  do iface=numInnerFaces+1,numFaces
    flmass_grid(iface) = den(owner(iface))/dt * flux_grid(iface)
  enddo

  !
  ! 5. Geometry for the finite volume discretization from the new nodal coordinates.
  !
  call mm_update_geometry

  !
  ! Velocity of the fluid at the moving walls = velocity of the face. Its normal component is the swept volume over dt,
  ! so that there is exactly no flux through the moving wall relative to the grid.
  !
  do k=1,nmvface
    iface = mv_face(k)
    ijb = mv_ijb(k)
    fx = ( xf(iface) - xfo(k) )/dt
    fy = ( yf(iface) - yfo(k) )/dt
    fz = ( zf(iface) - zfo(k) )/dt
    are2 = arx(iface)**2 + ary(iface)**2 + arz(iface)**2
    qn = flux_grid(iface)/dt
    vn = fx*arx(iface) + fy*ary(iface) + fz*arz(iface)
    u(ijb) = fx + arx(iface)*( qn - vn )/are2
    v(ijb) = fy + ary(iface)*( qn - vn )/are2
    w(ijb) = fz + arz(iface)*( qn - vn )/are2
  enddo

  !
  ! 6. Mass fluxes relative to the moving grid.
  !
  if( shift_flux ) then
    do i=1,numInnerFaces
      flmass(i) = flmass(i) - flmass_grid(i)
    enddo
  endif

  !
  ! Report: size of the motion, volumes, space conservation error.
  !
  dmax = maxval( sqrt( (x-xo)**2 + (y-yo)**2 + (z-zo)**2 ) )
  voltot = sum( vol )
  gclerr = maxval( abs( vgeo - vol )/vol )

  ! Boundary faces other than moving walls, symmetry and empty planes should be at rest.
  fixedmax = 0.0_dp
  if( .not. warned_fixed ) then
    do i=1,numBoundaries
      if( patch_moves(i) .or. bctype(i) == 'empty' .or. bctype(i) == 'symmetry' ) cycle
      do k=1,nfaces(i)
        fixedmax = max( fixedmax, abs( flux_grid( startFace(i) + k ) ) )
      enddo
    enddo
  endif
  if( fixedmax > 1.0e-10_dp*maxval(vol) .and. .not. warned_fixed ) then
    write(6,'(a,es10.3)') '  Warning: boundary faces that are not on moving patches sweep volume (max): ', fixedmax
    warned_fixed = .true.
  endif

  write(6,'(a,es10.3,2(a,es10.3),a,es10.3)') '  Moving mesh: max displacement step = ', dmax, '  min vol = ', volmin, &
       '  total vol = ', voltot, '  GCL error = ', gclerr

end subroutine mm_step


!***********************************************************************
!
subroutine mm_displacement( t, px, py, pz, dx, dy, dz )
!
!***********************************************************************
!
!  Displacement of a point of a moving patch, initially at p = (px,py,pz), at time t (relative to the initial mesh).
!
  implicit none

  real(dp), intent(in) :: t, px, py, pz
  real(dp), intent(out) :: dx, dy, dz

  real(dp) :: s, theta, ct, st, rx, ry, rz, kdr, kx, ky, kz

  select case( trim(mm_motion) )

    case('oscillate')

      s = sin( 2.0_dp*pi*mm_freq*t )
      dx = mm_amp(1)*s
      dy = mm_amp(2)*s
      dz = mm_amp(3)*s

    case('translate')

      dx = mm_vel(1)*t
      dy = mm_vel(2)*t
      dz = mm_vel(3)*t

    case('rotate')

      ! Rodrigues' rotation formula about the axis through mm_center with direction axis_n
      theta = mm_omega*t
      ct = cos( theta )
      st = sin( theta )
      rx = px - mm_center(1)
      ry = py - mm_center(2)
      rz = pz - mm_center(3)
      kx = axis_n(1);  ky = axis_n(2);  kz = axis_n(3)
      kdr = kx*rx + ky*ry + kz*rz
      dx = rx*ct + ( ky*rz - kz*ry )*st + kx*kdr*( 1.0_dp - ct ) - rx
      dy = ry*ct + ( kz*rx - kx*rz )*st + ky*kdr*( 1.0_dp - ct ) - ry
      dz = rz*ct + ( kx*ry - ky*rx )*st + kz*kdr*( 1.0_dp - ct ) - rz

    case('user')

      call mm_user_displacement( t, px, py, pz, dx, dy, dz )

    case default

      dx = 0.0_dp;  dy = 0.0_dp;  dz = 0.0_dp

  end select

end subroutine mm_displacement


!***********************************************************************
!
subroutine mm_user_displacement( t, px, py, pz, dx, dy, dz )
!
!***********************************************************************
!
!  Displacement of the nodes of the moving patches for mm_motion = 'user'. Edit this routine to prescribe any motion
!  (translation, rotation, pitching, deformation of the patch). The displacement must be zero at t = 0 and
!  it is measured from the initial position p = (px,py,pz) of the node.
!
!  Example, pitching about the z-axis through the point (0.25, 0, 0) with amplitude 5 degrees and frequency 1 Hz:
!
!    theta = 5.0_dp*pi/180.0_dp*sin( 2.0_dp*pi*t )
!    dx = 0.25_dp + (px-0.25_dp)*cos(theta) - py*sin(theta) - px
!    dy = (px-0.25_dp)*sin(theta) + py*cos(theta) - py
!    dz = 0.0_dp
!
  implicit none

  real(dp), intent(in) :: t, px, py, pz
  real(dp), intent(out) :: dx, dy, dz

  dx = 0.0_dp
  dy = 0.0_dp
  dz = 0.0_dp

end subroutine mm_user_displacement


!***********************************************************************
!
function rbf_phi( r ) result( phi )
!
!***********************************************************************
!
!  Radial basis function of the distance r.
!
  implicit none

  real(dp), intent(in) :: r
  real(dp) :: phi
  real(dp) :: rho

  select case( rbf_kind )

    case(1) ! Wendland C2, compact support rbf_R, positive definite
      rho = r/rbf_R
      if( rho >= 1.0_dp ) then
        phi = 0.0_dp
      else
        phi = (1.0_dp-rho)**4 * ( 4.0_dp*rho + 1.0_dp )
      endif

    case(2) ! Gaussian
      phi = exp( -( mm_epsilon*r )**2 )

    case(3) ! Multiquadric
      phi = sqrt( 1.0_dp + ( mm_epsilon*r )**2 )

    case(4) ! Inverse multiquadric
      phi = 1.0_dp/sqrt( 1.0_dp + ( mm_epsilon*r )**2 )

    case default ! Thin plate spline
      if( r < 1.0e-10_dp ) then
        phi = 0.0_dp
      else
        phi = r**2*log( r )
      endif

  end select

end function rbf_phi


!***********************************************************************
!
subroutine mm_rbf_matrix
!
!***********************************************************************
!
!  RBF matrix of the control points at their initial positions: A(i,j) = phi( |p_i - p_j| ), stored in lu.
!
  implicit none

  integer :: i, j
  real(dp) :: r

  do j=1,nctrl
    do i=1,nctrl
      r = sqrt( ( xini(ctrl(i)) - xini(ctrl(j)) )**2 + ( yini(ctrl(i)) - yini(ctrl(j)) )**2 + &
                ( zini(ctrl(i)) - zini(ctrl(j)) )**2 )
      lu(i,j) = rbf_phi( r )
    enddo
  enddo

end subroutine mm_rbf_matrix


!***********************************************************************
!
subroutine mm_check_interpolation
!
!***********************************************************************
!
!  Quality of the RBF system: interpolates a test displacement (1 on the moving control points, 0 on the fixed ones) back to
!  the control points. A large error means an ill-conditioned matrix (support radius too large for the number of control points).
!
  implicit none

  integer :: i, j
  real(dp) :: r, s, err

  do i=1,nctrl
    if( cmoving(i) ) then
      alx(i) = 1.0_dp
    else
      alx(i) = 0.0_dp
    endif
  enddo

  call lu_solve( nctrl, lu, piv, alx )

  err = 0.0_dp
  do i=1,nctrl
    s = 0.0_dp
    do j=1,nctrl
      r = sqrt( ( xini(ctrl(i)) - xini(ctrl(j)) )**2 + ( yini(ctrl(i)) - yini(ctrl(j)) )**2 + &
                ( zini(ctrl(i)) - zini(ctrl(j)) )**2 )
      s = s + alx(j)*rbf_phi( r )
    enddo
    if( cmoving(i) ) then
      err = max( err, abs( s - 1.0_dp ) )
    else
      err = max( err, abs( s ) )
    endif
  enddo

  write(*,'(a,es10.3)') '  RBF interpolation error at the control points (test displacement 1): ', err
  if( err > 1.0e-6_dp ) then
    write(*,'(a)') '  Warning: the RBF matrix is ill-conditioned. Reduce mm_radius or use fewer control points.'
  endif

end subroutine mm_check_interpolation


!***********************************************************************
!
subroutine mm_move_nodes( t )
!
!***********************************************************************
!
!  Nodal coordinates at time t: initial coordinates plus displacement. The control points get the prescribed displacement,
!  all other nodes the RBF interpolation of the displacement of the control points.
!
  implicit none

  real(dp), intent(in) :: t

  integer :: i, j, inode
  real(dp) :: r, phi, sx, sy, sz

  ! Prescribed displacement of the control points
  do i=1,nctrl
    inode = ctrl(i)
    if( cmoving(i) ) then
      call mm_displacement( t, xini(inode), yini(inode), zini(inode), cdx(i), cdy(i), cdz(i) )
    else
      cdx(i) = 0.0_dp
      cdy(i) = 0.0_dp
      cdz(i) = 0.0_dp
    endif
  enddo

  ! RBF weights
  alx = cdx
  aly = cdy
  alz = cdz
  call lu_solve( nctrl, lu, piv, alx )
  call lu_solve( nctrl, lu, piv, aly )
  call lu_solve( nctrl, lu, piv, alz )

  ! Interpolation to all nodes
  do inode=1,numNodes

    i = ictrl(inode)

    if( i > 0 ) then

      x(inode) = xini(inode) + cdx(i)
      y(inode) = yini(inode) + cdy(i)
      z(inode) = zini(inode) + cdz(i)

    else

      sx = 0.0_dp
      sy = 0.0_dp
      sz = 0.0_dp

      do j=1,nctrl
        r = sqrt( ( xini(inode) - xini(ctrl(j)) )**2 + ( yini(inode) - yini(ctrl(j)) )**2 + ( zini(inode) - zini(ctrl(j)) )**2 )
        phi = rbf_phi( r )
        sx = sx + alx(j)*phi
        sy = sy + aly(j)*phi
        sz = sz + alz(j)*phi
      enddo

      x(inode) = xini(inode) + sx
      y(inode) = yini(inode) + sy
      z(inode) = zini(inode) + sz

    endif

  enddo

end subroutine mm_move_nodes


!***********************************************************************
!
subroutine mm_swept_volumes
!
!***********************************************************************
!
!  Volume swept by every face while the nodes move linearly from (xo,yo,zo) to (x,y,z) in the time step.
!
!  Faces are split in triangles (1,i+1,i+2), i=1..nnodes-2, as in the calculation of the geometry in the geometry module.
!  The volume swept by a triangle is the integral over the time step of the flux of the velocity of the moving surface through it.
!  The velocity is linear over the triangle and constant in time, so the flux is S(t).dc, with dc the displacement of the
!  centroid of the triangle and S(t) the area vector, a quadratic function of time. Simpson's rule is thus exact:
!
!     dV = ( S(0) + 4 S(1/2) + S(1) ).dc / 6.
!
!  Since the volume of a polyhedron bounded by such triangles is a cubic function of time, the sum of the swept volumes of its
!  faces is exactly its volume change: space conservation law (GCL) is satisfied to the round-off.
!
  implicit none

  integer :: iface, i, n1, n2, n3
  real(dp) :: s0x, s0y, s0z, smx, smy, smz, s1x, s1y, s1z, dcx, dcy, dcz, dv
  real(dp) :: ax, ay, az, bx, by, bz, cx, cy, cz

  do iface=1,numFaces

    flux_grid(iface) = 0.0_dp
    n1 = node(1,iface)

    do i=1,nnodes(iface)-2

      n2 = node(i+1,iface)
      n3 = node(i+2,iface)

      ! Old triangle
      call tri_area_vector( xo(n1),yo(n1),zo(n1), xo(n2),yo(n2),zo(n2), xo(n3),yo(n3),zo(n3), s0x,s0y,s0z )

      ! New triangle
      call tri_area_vector( x(n1),y(n1),z(n1), x(n2),y(n2),z(n2), x(n3),y(n3),z(n3), s1x,s1y,s1z )

      ! Triangle at the middle of the time step
      ax = 0.5_dp*( xo(n1) + x(n1) );  ay = 0.5_dp*( yo(n1) + y(n1) );  az = 0.5_dp*( zo(n1) + z(n1) )
      bx = 0.5_dp*( xo(n2) + x(n2) );  by = 0.5_dp*( yo(n2) + y(n2) );  bz = 0.5_dp*( zo(n2) + z(n2) )
      cx = 0.5_dp*( xo(n3) + x(n3) );  cy = 0.5_dp*( yo(n3) + y(n3) );  cz = 0.5_dp*( zo(n3) + z(n3) )
      call tri_area_vector( ax,ay,az, bx,by,bz, cx,cy,cz, smx,smy,smz )

      ! Displacement of the centroid of the triangle
      dcx = ( x(n1)-xo(n1) + x(n2)-xo(n2) + x(n3)-xo(n3) )/3.0_dp
      dcy = ( y(n1)-yo(n1) + y(n2)-yo(n2) + y(n3)-yo(n3) )/3.0_dp
      dcz = ( z(n1)-zo(n1) + z(n2)-zo(n2) + z(n3)-zo(n3) )/3.0_dp

      dv = ( ( s0x + 4.0_dp*smx + s1x )*dcx + ( s0y + 4.0_dp*smy + s1y )*dcy + ( s0z + 4.0_dp*smz + s1z )*dcz )/6.0_dp

      flux_grid(iface) = flux_grid(iface) + dv

    enddo

  enddo

end subroutine mm_swept_volumes


!***********************************************************************
!
pure subroutine tri_area_vector( x1, y1, z1, x2, y2, z2, x3, y3, z3, sx, sy, sz )
!
!***********************************************************************
!
!  Area vector of the triangle (1,2,3): 1/2 (p2-p1) x (p3-p1).
!
  implicit none

  real(dp), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
  real(dp), intent(out) :: sx, sy, sz

  real(dp) :: px, py, pz, qx, qy, qz

  px = x2 - x1;  py = y2 - y1;  pz = z2 - z1
  qx = x3 - x1;  qy = y3 - y1;  qz = z3 - z1

  sx = 0.5_dp*( py*qz - pz*qy )
  sy = 0.5_dp*( pz*qx - px*qz )
  sz = 0.5_dp*( px*qy - py*qx )

end subroutine tri_area_vector


!***********************************************************************
!
subroutine mm_update_geometry
!
!***********************************************************************
!
!  Recalculates the geometry of the finite volume discretization from the nodal coordinates. It follows the calculation of the
!  geometry in the geometry module (read_mesh_native): face area vectors and centers (area weighted centroid of the triangles of
!  the face), cell centroids (Wang 1999), interpolation factor, Df, and distances of the wall and symmetry faces.
!
!  Cell volume from the nodal coordinates is kept in vgeo only for the check of the space conservation law: the volume used by
!  the solver is the one from the sum of the swept volumes (vol).
!
  implicit none

  integer :: iface, i, ib, inp, inn, iWall, iSym
  real(dp), parameter :: third = 1.0_dp/3.0_dp
  real(dp) :: px, py, pz, qx, qy, qz, nx, ny, nz, cx, cy, cz, are, areasum, riSi
  real(dp) :: xpn, ypn, zpn, dpn, djp, djn, nxf, nyf, nzf, xjp, yjp, zjp

  arx = 0.0_dp;  ary = 0.0_dp;  arz = 0.0_dp
  xf = 0.0_dp;   yf = 0.0_dp;   zf = 0.0_dp
  xc = 0.0_dp;   yc = 0.0_dp;   zc = 0.0_dp
  vgeo = 0.0_dp
  r8tmp = 0.0_dp

  do iface=1,numFaces

    inp = owner(iface)
    areasum = 0.0_dp

    do i=1,nnodes(iface)-2

      px = x( node(i+1,iface) ) - x( node(1,iface) )
      py = y( node(i+1,iface) ) - y( node(1,iface) )
      pz = z( node(i+1,iface) ) - z( node(1,iface) )
      qx = x( node(i+2,iface) ) - x( node(1,iface) )
      qy = y( node(i+2,iface) ) - y( node(1,iface) )
      qz = z( node(i+2,iface) ) - z( node(1,iface) )

      call tri_area_vector( 0.0_dp,0.0_dp,0.0_dp, px,py,pz, qx,qy,qz, nx,ny,nz )

      arx(iface) = arx(iface) + nx
      ary(iface) = ary(iface) + ny
      arz(iface) = arz(iface) + nz

      cx = third*( x( node(i+2,iface) ) + x( node(i+1,iface) ) + x( node(1,iface) ) )
      cy = third*( y( node(i+2,iface) ) + y( node(i+1,iface) ) + y( node(1,iface) ) )
      cz = third*( z( node(i+2,iface) ) + z( node(i+1,iface) ) + z( node(1,iface) ) )

      are = sqrt( nx**2 + ny**2 + nz**2 )

      xf(iface) = xf(iface) + are*cx
      yf(iface) = yf(iface) + are*cy
      zf(iface) = zf(iface) + are*cz

      areasum = areasum + are

      riSi = cx*nx + cy*ny + cz*nz

      vgeo(inp) = vgeo(inp) + third*riSi
      xc(inp) = xc(inp) + 0.75_dp*riSi*cx
      yc(inp) = yc(inp) + 0.75_dp*riSi*cy
      zc(inp) = zc(inp) + 0.75_dp*riSi*cz
      r8tmp(inp) = r8tmp(inp) + riSi

      if( iface <= numInnerFaces ) then

        inn = neighbour(iface)
        riSi = -riSi

        vgeo(inn) = vgeo(inn) + third*riSi
        xc(inn) = xc(inn) + 0.75_dp*riSi*cx
        yc(inn) = yc(inn) + 0.75_dp*riSi*cy
        zc(inn) = zc(inn) + 0.75_dp*riSi*cz
        r8tmp(inn) = r8tmp(inn) + riSi

      endif

    enddo

    xf(iface) = xf(iface)/areasum
    yf(iface) = yf(iface)/areasum
    zf(iface) = zf(iface)/areasum

  enddo

  do inp=1,numCells
    xc(inp) = xc(inp)/r8tmp(inp)
    yc(inp) = yc(inp)/r8tmp(inp)
    zc(inp) = zc(inp)/r8tmp(inp)
  enddo

  !
  ! Interpolation factor, the same variant as in the geometry module.
  !
  do iface=1,numInnerFaces

    inp = owner(iface)
    inn = neighbour(iface)

    if( interpolation_coeff_variant == 1 ) then

      xpn = xc(inn) - xc(inp)
      ypn = yc(inn) - yc(inp)
      zpn = zc(inn) - zc(inp)
      dpn = sqrt( xpn**2 + ypn**2 + zpn**2 )

      call find_intersection_point( x(node(1,iface)), y(node(1,iface)), z(node(1,iface)), &
                                    x(node(2,iface)), y(node(2,iface)), z(node(2,iface)), &
                                    x(node(3,iface)), y(node(3,iface)), z(node(3,iface)), &
                                    xc(inp), yc(inp), zc(inp), xc(inn), yc(inn), zc(inn), xjp, yjp, zjp )

      djn = sqrt( (xjp-xc(inp))**2 + (yjp-yc(inp))**2 + (zjp-zc(inp))**2 )
      facint(iface) = djn/dpn

    elseif( interpolation_coeff_variant == 2 ) then

      djp = sqrt( (xf(iface)-xc(inp))**2 + (yf(iface)-yc(inp))**2 + (zf(iface)-zc(inp))**2 )
      djn = sqrt( (xf(iface)-xc(inn))**2 + (yf(iface)-yc(inn))**2 + (zf(iface)-zc(inn))**2 )
      facint(iface) = djp/( djp + djn )

    else

      djp = abs( arx(iface)*(xf(iface)-xc(inp)) + ary(iface)*(yf(iface)-yc(inp)) + arz(iface)*(zf(iface)-zc(inp)) )
      dpn = abs( arx(iface)*(xc(inn)-xc(inp)) + ary(iface)*(yc(inn)-yc(inp)) + arz(iface)*(zc(inn)-zc(inp)) )
      facint(iface) = djp/dpn

    endif

    ! Df = Sf.Sf / (Sf.dpn)
    xpn = xc(inn) - xc(inp)
    ypn = yc(inn) - yc(inp)
    zpn = zc(inn) - zc(inp)
    are = arx(iface)**2 + ary(iface)**2 + arz(iface)**2
    Df(iface) = are/( arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn )

  enddo

  !
  ! Distances of the wall and symmetry faces to the cell centers.
  !
  iWall = 0
  iSym = 0

  do ib=1,numBoundaries

    if( bctype(ib) == 'symmetry' ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        inp = owner(iface)
        iSym = iSym + 1
        are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
        nxf = arx(iface)/are;  nyf = ary(iface)/are;  nzf = arz(iface)/are
        dns(iSym) = (xf(iface)-xc(inp))*nxf + (yf(iface)-yc(inp))*nyf + (zf(iface)-zc(inp))*nzf
        srds(iSym) = are/dns(iSym)
      enddo

    elseif( bctype(ib) == 'wall' ) then

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        inp = owner(iface)
        iWall = iWall + 1
        are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
        nxf = arx(iface)/are;  nyf = ary(iface)/are;  nzf = arz(iface)/are
        dnw(iWall) = (xf(iface)-xc(inp))*nxf + (yf(iface)-yc(inp))*nyf + (zf(iface)-zc(inp))*nzf
        srdw(iWall) = are/dnw(iWall)
      enddo

    endif

  enddo

end subroutine mm_update_geometry


!***********************************************************************
!
subroutine mm_pressure_source( su )
!
!***********************************************************************
!
!  Volume change of the cells in the source of the pressure(-correction) equation:
!
!       su(inp) = su(inp) + den(inp)/dt*( volo(inp) - vol(inp) ).
!
!  Call it after the continuity residual sum_f flmass_rel of the cells (flmass relative to the moving grid) has been
!  assembled in su (su = - sum_f flmass_rel). The absolute mass flux through the faces of a cell, flmass_rel + flmass_grid,
!  then adds up to zero, since sum_f flmass_grid = den/dt*( vol - volo ) (space conservation law).
!  The same term has to be in the continuity error report.
!
  implicit none

  real(dp), dimension(numCells), intent(inout) :: su

  integer :: inp

  if( .not. DynamicMesh ) return

  do inp=1,numCells
    su(inp) = su(inp) + den(inp)/timestep*( volo(inp) - vol(inp) )
  enddo

end subroutine mm_pressure_source


!***********************************************************************
!
pure function vol_unsteady( inp ) result( v )
!
!***********************************************************************
!
!  Cell volume for the first order unsteady term den*V/dt*(phi - phio): volo on the moving mesh, vol on a fixed one.
!
  implicit none

  integer, intent(in) :: inp
  real(dp) :: v

  if( DynamicMesh ) then
    v = volo(inp)
  else
    v = vol(inp)
  endif

end function vol_unsteady


!***********************************************************************
!
subroutine mm_write_points( unit, ib )
!
!***********************************************************************
!
!  Writes the current nodal coordinates in the ASCII Points array of the VTK files in place of those of the template:
!  ib = 0: all nodes (interior.vtu); ib > 0: nodes of the faces of the boundary region ib in the order of ascending node
!  number (the points of the .vtp templates written by examples/coupledPU/vtkmesh.py).
!
  implicit none

  integer, intent(in) :: unit, ib

  integer :: i, iface, ivertex
  logical, dimension(:), allocatable :: used

  if( ib == 0 ) then

    do i=1,numNodes
      write(unit,'(10x,3(es16.8,1x))') x(i), y(i), z(i)
    enddo

  else

    allocate( used(numNodes) )
    used = .false.
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      do ivertex=1,nnodes(iface)
        used( node(ivertex,iface) ) = .true.
      enddo
    enddo
    do i=1,numNodes
      if( used(i) ) write(unit,'(10x,3(es16.8,1x))') x(i), y(i), z(i)
    enddo
    deallocate( used )

  endif

end subroutine mm_write_points


!***********************************************************************
!
subroutine lu_factor( n, a, piv, info )
!
!***********************************************************************
!
!  LU factorization with partial pivoting of the dense matrix a(n,n), in place. info = 0 on success, otherwise the column of
!  the zero pivot. Row interchanges are stored in piv.
!
  implicit none

  integer, intent(in) :: n
  real(dp), dimension(n,n), intent(inout) :: a
  integer, dimension(n), intent(out) :: piv
  integer, intent(out) :: info

  integer :: j, k, p
  real(dp) :: amax
  real(dp), dimension(n) :: row

  info = 0

  do k=1,n

    ! Pivot
    p = k - 1 + maxloc( abs( a(k:n,k) ), 1 )
    amax = abs( a(p,k) )

    if( amax < 1.0e-300_dp ) then
      info = k
      return
    endif

    piv(k) = p

    if( p /= k ) then
      row = a(k,:)
      a(k,:) = a(p,:)
      a(p,:) = row
    endif

    a(k+1:n,k) = a(k+1:n,k)/a(k,k)

    do j=k+1,n
      a(k+1:n,j) = a(k+1:n,j) - a(k+1:n,k)*a(k,j)
    enddo

  enddo

end subroutine lu_factor


!***********************************************************************
!
subroutine lu_solve( n, a, piv, b )
!
!***********************************************************************
!
!  Solves a x = b with the LU factors from lu_factor, b is overwritten by the solution.
!
  implicit none

  integer, intent(in) :: n
  real(dp), dimension(n,n), intent(in) :: a
  integer, dimension(n), intent(in) :: piv
  real(dp), dimension(n), intent(inout) :: b

  integer :: j, k
  real(dp) :: t

  ! Row interchanges (all of them first, the L factor stored in a is the one of the fully permuted matrix)
  do k=1,n
    if( piv(k) /= k ) then
      t = b(k);  b(k) = b(piv(k));  b(piv(k)) = t
    endif
  enddo

  ! Forward substitution (unit lower triangular L)
  do k=1,n
    b(k+1:n) = b(k+1:n) - a(k+1:n,k)*b(k)
  enddo

  ! Back substitution (U)
  do j=n,1,-1
    b(j) = b(j)/a(j,j)
    b(1:j-1) = b(1:j-1) - a(1:j-1,j)*b(j)
  enddo

end subroutine lu_solve

end module moving_mesh
