!***********************************************************************
!
subroutine init
!
!***********************************************************************
!     Contents:
!
! Various initialisations
!     1)  Parameter Initialisation
!     2)  Field Initialisation
! Read Restart File And Set Field Values
! Initial Gradient Calculation
! Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use gradients
  use sparse_matrix, only: su,sv
  use utils, only: get_unit,i4_to_s_left
  ! use LIS_linear_solver_library
  use field_initialization
  use output

  implicit none

  ! 
  ! Local variables 
  !
  character( len = 5) :: nproc_char
  integer :: i, ijp, ijn, ijb, ib, iface, ipro
  ! integer :: output_unit
  integer :: nsw_backup
  real(dp) :: fxp, fxn, ui, vi, wi
  real(dp) :: sor_backup

!
!***********************************************************************
!

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

!
! Various initialisations
!

! 1) Parameter Initialisation

  !
  ! Pressure reference cell
  !

  iPrefProcess = 0
  pRefCell = 1


  ! Switches which define wether we look for some files in folder 0.
  if (lturb) then
    if ( TurbModel==1 .or. TurbModel==2 ) then
      solveEpsilon = .true.
      solveTKE = .true.
    elseif ( TurbModel==3 .or. TurbModel==4 ) then
      solveOmega = .true.
      solveTKE = .true.
    elseif( TurbModel==6 ) then
      solveTKE = .true.
    else
      solveOmega = .false.
      solveEpsilon = .false.
      solveTKE = .false.
    endif
  endif

  ! Set the string for second scalar equation, which appears in linear solver log, default is 'epsilon',
  ! Note, change also the script 'plotResiduals'
  if ( solveOmega ) chvarSolver(6) = 'Omega  '


  ! Reciprocal and complementary to one values of underrelaxation factors.
  do i=1,nphi
    urfr(i)=1.0_dp / urf(i)
    urfm(i)=1.0_dp - urf(i)
  enddo

  ! Initial time
  if(.not.lread) time = 0.0_dp

  ! Set to zero cumulative error in continuity
  cumulativeContErr = 0.0_dp

  ! Bulk velocity - important const_mflux flow!
  magUbar = uin


! 2)  Field Initialisation
 
  if(myid .eq. 0) then 
    write(*,'(a)') ' '
    write(*,'(a)') '  Initializing internal field and boundaries (reading 0/.. ):'
    write(*,'(a)') ' '
  endif

  ! Field initialisation vector field

  call initialize_vector_field(u,v,w,dUdxi,'U')

  ! Initialize previous time step value to current value.
  uo = u
  vo = v
  wo = w

  uoo = u
  voo = v
  woo = w

  if (bdf3) then
    uooo = u
    vooo = v
    wooo = w
  endif

  ! Field initialisation scalars

  ! 
  ! > TE Turbulent kinetic energy.
  !   
  call initialize_scalar_field(te,dTEdxi,'k')

  ! 
  ! > ED Specific turbulent kinetic energy dissipation rate, also turbulence frequency - omega
  !  
  if(solveOmega) then   
    call initialize_scalar_field(ed,dEDdxi,'omega')
  else
    call initialize_scalar_field(ed,dEDdxi,'epsilon')
  endif


  ! Density
  den = densit

  ! Effective viscosity
  vis = viscos

  visw = viscos

  ! Temperature
  if(lcal(ien)) t = tin

  ! Temperature variance
  if(lcal(ivart)) vart = vartin
  
  ! Concentration
  if(lcal(icon)) con=conin

  ! ! Reynolds stress tensor components
  ! if (lturb) then
  !   uu = 0.0_dp
  !   vv = 0.0_dp
  !   ww = 0.0_dp
  !   uv = 0.0_dp
  !   uw = 0.0_dp
  !   vw = 0.0_dp
  ! endif

  ! ! Turbulent heat fluxes
  ! if(lcal(ien).and.lbuoy) then
  !   utt = 0.0_dp
  !   vtt = 0.0_dp
  !   wtt = 0.0_dp
  ! endif

  ! ! Reynolds stress anisotropy
  ! if(lturb.and.lasm) bij = 0.0_dp

  ! ! Pressure and pressure correction
  ! p = 0.0_dp
  ! pp = p


  !
  ! > Initialize mass flow
  !

  ! Mass flow trough inner faces
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)

    fxn = facint(i)
    fxp = 1.0_dp-facint(i)

    ui = u(ijp)*fxp + u(ijn)*fxn
    vi = v(ijp)*fxp + v(ijn)*fxn
    wi = w(ijp)*fxp + w(ijn)*fxn

    flmass(i) = den(ijp)*( arx(i)*ui + ary(i)*vi + arz(i)*wi )

  enddo

   
  !
  ! Mass flow at boundaries of inner domain and buffer cells
  !
  call exchange( u )
  call exchange( v )
  call exchange( w )


  ! Initialize mass flux at faces on process boundary

  ipro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'process') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijn = iBndValueStart(ib) + i
        ipro = ipro + 1

        fxn = fpro(ipro)
        fxp = 1.0_dp-fxn

        ui = u(ijp)*fxp + u(ijn)*fxn
        vi = v(ijp)*fxp + v(ijn)*fxn
        wi = w(ijp)*fxp + w(ijn)*fxn

        flmass(iface) = den(ijp)*( arx(iface)*ui + ary(iface)*vi + arz(iface)*wi )

      enddo
    endif
  enddo


!
! Read Restart File And Set Field Values
!
  if(lread) call readfiles

!
! Initial Gradient Calculation
!

  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_grad_matrix(U,dUdxi)
  endif

  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  !
  ! Distance to the nearest wall (needed for some turbulence models).
  !

  if(myid .eq. 0) then
    write(*,*) ' '
    write(*,*) ' Calculate distance to the nearest wall:'
    write(*,*) ' '
  endif

  ! Source term
  su(1:numCells) = -Vol(1:numCells)

  ! Initialize solution
  p = 0.0_dp

  !  Coefficient array for Laplacian (we choose apu because it is numPCells long,
  !  what is required for 'mu' argument in Laplacian)
  ! apu = 1.0_dp       

  ! Laplacian operator and BCs         
  call laplacian(p) 

  sor_backup = sor(ip)
  nsw_backup = nsw(ip)

  sor(ip) = 1e-10
  nsw(ip) = 500

  ! Solve system
  ! call jacobi(p,ip)
  ! call dpcg(p,ip) 
  call iccg(p,ip) 
  ! call bicgstab(p,ip) 
  ! call solve_csr(numCells,nnz,ioffset,ja,a,su,p)

  ! Update values at constant gradient bc faces - we need these values for correct gradients

  call exchange ( p ) 

  ipro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) /= 'wall' .and. bctype(ib) /= 'proces' ) then
    ! All other boundary faces besides wall which has to be zero, and process, which is set above

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        p(ijb) = p(ijp)

      enddo

    endif
  enddo

  sor(ip) = sor_backup
  nsw(ip) = nsw_backup

  ! Gradient of solution field stored in p (gradient stored in dPdxi) :
  call grad(p,dPdxi)

  ! Wall distance computation from Poisson eq. solution stored in pp:
  wallDistance = -sqrt(  dPdxi(1,1:numCells)**2+dPdxi(2,1:numCells)**2+dPdxi(3,1:numCells)**2  ) + &
                  sqrt(  dPdxi(1,1:numCells)**2+dPdxi(2,1:numCells)**2+dPdxi(3,1:numCells)**2 + 2*p(1:numCells)  )

  ! do i = 1,numCells
  ! write(*,*) wallDistance(i)
  ! enddo

  ! Clear arrays
  su = 0.0_dp
  sv = 0.0_dp 
  p = 0.0_dp
  dPdxi = 0.0_dp


  ! ! Write wall distance field.
  ! !+-----------------------------------------------------------------------------+
  ! call get_unit( output_unit )

  ! open(unit=output_unit,file='processor'//trim(nproc_char)//'/VTK'// &
  ! &                          '/wallDistance-'//'proc'//trim(nproc_char)//'.vtu')

  ! ! Header
  ! call vtu_write_XML_header ( output_unit )
  ! ! Scalar field
  ! call vtu_write_XML_scalar_field ( output_unit, 'wallDistance', wallDistance )
  ! ! Mesh data
  ! call vtu_write_XML_meshdata ( output_unit )

  ! close( output_unit )
  ! !+-----------------------------------------------------------------------------+


end subroutine
