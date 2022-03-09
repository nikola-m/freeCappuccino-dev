module concentration
!
! Implementation of sclar transport equation for concentration.
!
  use types
  use scalar_fluxes

  implicit none

  !
  ! Discrtetization and solution parameters - modified trough input.nml file
  !
  logical  :: calcCon = .False.                               ! To activate the solution of this field in the main function. 
  real(dp) :: urfCon  = 0.9                                   ! Under-relaxation factors.
  real(dp) :: gdsCon = 1.0                                    ! Deferred correction factor.
  character ( len=30 ) :: cSchemeCon = 'boundedLinearUpwind'  ! Convection scheme - default is second order upwind.
  character ( len=12 ) :: dSchemeCon = 'skewness' ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
  integer :: nrelaxCon = 0 ! Type of non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed. 1/0/-1.
  character ( len=12 ) :: lSolverCon = 'bicgstab'             ! Linear algebraic solver.
  integer  :: maxiterCon = 10                                 ! Max number of iterations in linear solver.
  real(dp) :: tolAbsCon = 1e-13                               ! Absolute residual level.
  real(dp) :: tolRelCon = 0.025                               ! Relative drop in residual to exit linear solver.
  real(dp) :: sigCon = 0.9_dp                                 ! Prandtl-Schmidt
  
  public :: calcCon, urfCon, gdsCon, cSchemeCon, dSchemeCon, nrelaxCon, lSolverCon, maxiterCon, tolAbsCon, tolRelCon, sigCon ! params

  public :: calculate_concentration_field

contains


subroutine calculate_concentration_field
!
! Purpose: 
!   Assemble and solve transport equation for a passive scalar field.
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use linear_solvers
  
  implicit none

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iface, iWall
  integer :: iPer, l, if, iftwin
  real(dp) :: gam, prtr, apotime
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: fimax,fimin
  real(dp) :: urfr,urfm
  real(dp) :: are,nxf,nyf,nzf,dcoef,dCn


  ! Variable specific coefficients:
  gam = gdsCon          !< first to higher order convection scheme blending parameter gamma.
  prtr = 1.0d0/sigCon   !< reciprocal value of Prandtl-Schmidt number

! Calculate gradient: 
  call grad(Con,dCondxi)

! Initialize matrix and source arrays
  a  = 0.0_dp
  su = 0.0d0
  sp = 0.0d0

  !
  ! Volume source terms 
  !

  do inp=1,numCells

    ! UNSTEADY TERM
    if( bdf .or. cn ) then
      apotime = den(inp)*vol(inp)/timestep
      su(inp) = su(inp) + apotime*Cono(inp)
      sp(inp) = sp(inp) + apotime
    elseif( bdf2 ) then
      apotime=den(inp)*vol(inp)/timestep
      su(inp) = su(inp) + apotime*( 2*Cono(inp) - 0.5_dp*Conoo(inp) )
      sp(inp) = sp(inp) + 1.5_dp*apotime
    endif

  enddo


  !
  ! Calculate terms integrated over faces
  !

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, cSchemeCon, dSchemeCon, nrelaxCon,  &
                     Con, dCondxi, prtr, cap, can, suadd )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap


    ! > Elements on main diagonal:

    ! ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources:

    su(ijp) = su(ijp) + suadd
    su(ijn) = su(ijn) - suadd 

  enddo


  !
  ! Boundary conditions
  !
  iWall = 0
  iPer = 0
  l = numInnerFaces

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         Con, dCondxi, prtr, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Con(ijb) + suadd

      end do

    elseif ( bctype(ib) == 'wall' .and. bcDFraction(9,ib)==0.0_dp ) then

      ! Wall with zdro grad

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Con(ijb) = Con(ijp)

      enddo

    elseif ( bctype(ib) == 'wall' .and. bcDFraction(9,ib)==1.0_dp ) then

      ! Prescribed concentration flux wall boundaries (that's actually non homo Neumann)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        dcoef = viscos/pranl+(vis(ijp)-viscos)/sigCon
        dcoef = dcoef*srdw(iWall)

        ! Value of the concentration gradient in normal direction is set trough 
        ! proper choice of component values. Let's project in to normal direction
        ! to recover normal gradient.

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Gradient in face-normal direction        
        dCn = dCondxi(1,ijb)*nxf+dCondxi(2,ijb)*nyf+dCondxi(3,ijb)*nzf

        ! Explicit source
        su(ijp) = su(ijp) + dcoef*dCn

      enddo

    elseif ( bctype(ib) == 'periodic') then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face


        ! face flux scalar but for periodic boundaries - it will be recognized by arguments
        call facefluxsc_periodic( ijp, ijn, &
                                  xf(if), yf(if), zf(if), &
                                  arx(if), ary(if), arz(if), &
                                  flmass(if), gam, &
                                  Con, dCondxi, prtr, cap, can, suadd )


        ! > Off-diagonal elements:

        ! l is in interval [numInnerFaces+1, numInnerFaces+numPeriodic]
        l = l + 1

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_index(l)
        a(k) = can

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_index(l)
        a(k) = cap

        ! > Elements on main diagonal:

        ! ! (icell,icell) main diagonal element
        k = diag(ijp)
        a(k) = a(k) - can

        ! ! (jcell,jcell) main diagonal element
        k = diag(ijn)
        a(k) = a(k) - cap

        ! > Sources:

        su(ijp) = su(ijp) + suadd
        su(ijn) = su(ijn) - suadd 


      end do 


    endif
    
  enddo ! Boundary conditions loop  



  ! Modify coefficients for Crank-Nicolson
  if (cn) then

    a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

    do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      k = icell_jcell_csr_index(i)
      su(ijp) = su(ijp) - a(k)*Cono(ijn)

      k = jcell_icell_csr_index(i)
      su(ijn) = su(ijn) - a(k)*Cono(ijp)

    enddo

    do ijp=1,numCells
      apotime=den(ijp)*vol(ijp)/timestep
      off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) - a(diag(ijp))
      su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*Cono(ijp)
      sp(ijp) = sp(ijp)+apotime

    enddo

  endif

  ! Underrelaxation factors
  urfr = 1.0_dp / urfCon
  urfm = 1.0_dp - urfCon

  ! Main diagonal term assembly:
  do inp = 1,numCells
    ! Main diagonal term assembly:
    ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
    ! we substract it from the sum, to eliminate it from the sum.
    off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(ijp))
    a(diag(inp)) = sp(inp) - off_diagonal_terms

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = su(inp) + urfm*a(diag(inp))*Con(inp)

  enddo

  ! Solve linear system:
  call csrsolve(lSolverCon, Con, su, resor(9), maxiterCon, tolAbsCon, tolRelCon, 'Conc')

  !
  ! Update field values at boundaries
  !
  call updateBoundary( Con )

! Report range of scalar values and clip if negative
  fimin = minval( Con(1:numCells))
  fimax = maxval( Con(1:numCells))
  
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= Concentration <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) Con(1:numCells) = max( Con(1:numCells),small)

end subroutine


end module