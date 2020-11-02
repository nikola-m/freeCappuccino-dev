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
  use mhd

  implicit none

  ! 
  ! Local variables 
  !

  ! character(len=5)  :: maxno
  ! character(10) :: tol
  integer :: i, ijp, ijn, iface, ijb, ib
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

  ! Parameters used in postprocessing phase
  if (lturb) then

    select case (TurbModel)

      case (1)
        solveTKE = .true.
        solveEpsilon = .true.

      case (2)
        solveTKE = .true.
        solveEpsilon = .true.

      case (3) 
        solveTKE = .true.
        solveOmega = .true.

      case (4)
        solveTKE = .true.
        solveOmega = .true.

      case (6)
        solveTKE = .true.

      case (7)
        solveTKE = .true.
        solveEpsilon = .true.

      case (8)
        solveTKE = .true.
        solveEpsilon = .true.

      case (9)
        solveTKE = .true.
        solveEpsilon = .true.

      case default
          
      end select
    
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


! 1.2)  Field Initialisation

  call initialize_vector_field(u,v,w,dUdxi,iu,'U')

  if (levm) then
  ! 
  ! > TE Turbulent kinetic energy.
  !   
  call initialize_scalar_field(te,dTEdxi,ite,'k')

  ! 
  ! > ED Specific turbulent kinetic energy dissipation rate, also turbulence frequency - omega
  !  
  if(solveOmega) then   
    call initialize_scalar_field(ed,dEDdxi,ied,'omega')
  else
    call initialize_scalar_field(ed,dEDdxi,ied,'epsilon')
  endif

  endif

  ! 
  ! > Temperature
  !
  if( lcal(ien) )   call initialize_scalar_field(t,dTdxi,ien,'T')

  ! 
  ! > Magnetic field
  ! 
  if ( lcal(iep) ) call initialize_vector_field(bmagx,bmagy,bmagz,dEpotdxi,ibmag,'B')

  ! Density
  den = densit

  ! Effective viscosity
  vis = viscos
  visw = viscos  

  ! Initialize mass flow on inner faces
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
! Distance to the nearest wall (needed for some turbulence models) for all cells via Poisson equation.
!

  if( TurbModel == 3 .or. TurbModel == 4 ) then

    write(*,*) ' '
    write(*,*) ' Calculate distance to the nearest wall:'
    write(*,*) ' '

    ! Source term
    su(1:numCells) = -Vol(1:numCells)

    ! Initialize solution
    pp = 0.0_dp

    !  Coefficient array for Laplacian
    sv = 1.0_dp       

    ! Laplacian operator and BCs         
    call laplacian(sv,pp) 

    sor_backup = sor(ip)
    nsw_backup = nsw(ip)

    sor(ip) = 1e-10
    nsw(ip) = 500

    ! Solve system
    call iccg(pp,ip) 
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

          pp(ijb) = pp(ijp)

        enddo

      endif
    enddo

    sor(ip) = sor_backup
    nsw(ip) = nsw_backup

    ! Gradient of solution field stored in p (gradient stored in dPdxi) :
    call grad(pp,dPdxi)

    ! Wall distance computation from Poisson eq. solution stored in pp:
    wallDistance = -sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) ) + &
                    sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) + 2*pp(1:numCells) )

    ! Clear arrays
    su = 0.0_dp
    sv = 0.0_dp 
    dPdxi = 0.0_dp
    pp = p

  end if

end subroutine
