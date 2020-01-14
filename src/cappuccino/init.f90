! ***********************************************************************
!
subroutine init
!
!***********************************************************************
!   Contents:
!
!   Various initialisations
!       Parameter Initialisation
!       Field Initialisation
!   Read Restart File And Set Field Values
!   Initial Gradient Calculation
!   Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use gradients
  use sparse_matrix
  use utils, only: get_unit
  ! use LIS_linear_solver_library
  use field_initialization
  use output

  implicit none

  ! 
  ! Local variables 
  !

  ! character(len=5)  :: maxno
  ! character(10) :: tol
  integer :: i, ijp, ijn, iface, output_unit, ijb, ib
  integer :: nsw_backup
  real(dp) :: fxp, fxn, ui, vi, wi
  real(dp) :: sor_backup

!
!***********************************************************************
!

!
! Various initialisations
!

! Parameter Initialisation

  ! Set the string for second scalar equation, which appears in linear solver log, default is 'epsilon',
  ! Note, change also the script 'plotResiduals'
  if (lturb) then

    if ( TurbModel==3 .or. TurbModel==4 ) then

      solveOmega = .true.

    endif
    
  endif

  if ( solveOmega ) chvarSolver(6) = 'Omega  '


  ! Reciprocal values of underrelaxation factors
  do i=1,nphi
    urfr(i)=1.0_dp / urf(i)
    urfm(i)=1.0_dp - urf(i)
  enddo

  ! Initial time iz zero.
  if(.not.lread) time = 0.0_dp

  ! Set to zero cumulative error in continuity
  cumulativeContErr = 0.0_dp

  ! Bulk velocity - important const_mflux flow!
  magUbar = uin



! 1.2)  Field Initialisation
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

  !-------------------------------------------------------    
  ! Field initialisation over inner cells + boundary faces
  !-------------------------------------------------------

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
  if(lcal(icon)) con = conin

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

  ! Pressure and pressure correction
  p = 0.0_dp
  if ( simple ) pp = p
  if ( piso )   po = p
  if ( piso )   poo = p

  ! Initialize mass flow
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)

    fxn = facint(i)
    fxp = 1.0_dp-facint(i)

    ui = u(ijp)*fxp + u(ijn)*fxn
    vi = v(ijp)*fxp + v(ijn)*fxn
    wi = w(ijp)*fxp + w(ijn)*fxn

    flmass(i) = den(ijp)*(arx(i)*ui+ary(i)*vi+arz(i)*wi)

  enddo

  if (piso) then
    flmasso = flmass
    flmassoo = flmass
  endif

!
! Read Restart File And Set Field Values
!
  if(lread) call readfiles


!
! Initial Gradient Calculation
!

  ! Set matrix with geometric coefficients for Least-squares gradient reconstruction.
  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_grad_matrix(U,dUdxi)
  endif

  ! Initial velocity gradient calculation with given velocity field.
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

! 
! Initialization of residual for all variables
!
  ! do i=1,nphi
  !   rnor(i) = 1.0_dp
  !   resor(i)= 0.0_dp
  ! enddo

  ! rnor(iu)  = 1.0_dp/(flomom+small)
  ! rnor(iv)  = rnor(iu)
  ! rnor(iw)  = rnor(iu)

  ! rnor(ip)  = 1.0_dp/(flomas+small)
  ! rnor(ite) = 1.0_dp/(flowte+small)
  ! rnor(ied) = 1.0_dp/(flowed+small)

!
! Distance to the nearest wall (needed for some turbulence models) for all cells via Poisson equation.
!

  write(*,*) ' '
  write(*,*) ' Calculate distance to the nearest wall:'
  write(*,*) ' '

  ! Source term
  su(1:numCells) = -Vol(1:numCells)

  ! Initialize solution
  p = 0.0_dp

  !  Coefficient array for Laplacian
  sv = 1.0_dp       

  ! Laplacian operator and BCs         
  call laplacian(sv,p) 

  sor_backup = sor(ip)
  nsw_backup = nsw(ip)

  sor(ip) = 1e-10
  nsw(ip) = 500

  ! Solve system
  call iccg(p,ip) 
  ! call bicgstab(p,ip) 
  ! call pmgmres_ilu ( numCells, nnz, ioffset, ja, a, diag, p(1:numCells), ip, su, 100, 4, 1e-8, sor(ip) )
  ! write(maxno,'(i5)') nsw(ip)
  ! write(tol,'(es9.2)') sor(ip)
  ! ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
  ! write(options,'(a)') "-i cg -p ilu -ilu_fill 1 -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
  ! call solve_csr( numCells, nnz, ioffset, ja, a, su, p )


  ! Update values at constant gradient bc faces - we need these values for correct gradients

  do ib=1,numBoundaries

    if ( bctype(ib) /= 'wall' ) then
    ! All other boundary faces besides wall which has to be zero.

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
  wallDistance = -sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:)  ) + &
                  sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) + 2*p(1:numCells)  )

  ! Clear arrays
  su = 0.0_dp
  sv = 0.0_dp 
  p = 0.0_dp
  dPdxi = 0.0_dp


  ! Write wall distance field.
  !+-----------------------------------------------------------------------------+
  call get_unit( output_unit )

  open(unit=output_unit,file='VTK/wallDistance.vtu')

  ! Header
  call vtu_write_XML_header ( output_unit )
  ! Scalar field
  call vtu_write_XML_scalar_field ( output_unit, 'wallDistance', wallDistance )
  ! Mesh data
  call vtu_write_XML_meshdata ( output_unit )

  close( output_unit )
  !+-----------------------------------------------------------------------------+


end subroutine
