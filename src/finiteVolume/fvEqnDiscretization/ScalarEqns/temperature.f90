module temperature
!
! Implementation of sclar transport equation for temperature.
!
  use types
  use scalar_fluxes, only: facefluxsc

  implicit none

  !
  ! Discrtetization and solution parameters - modified trough input.nml file
  !
  logical  :: calcT = .False.                        ! To activate the solution of this field in the main function. 
  real(dp) :: urfT  = 0.9_dp                         ! Under-relaxation factors.
  real(dp) :: gdsT = 1.0_dp                          ! Deferred correction factor.
  character ( len=30 ) :: cSchemeT = 'linearUpwind'  ! Convection scheme - default is second order upwind.
  character ( len=12 ) :: dSchemeT = 'skewness'      ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
  integer :: nrelaxT = 0                             ! Face-normal gradient non-orthogonal correction relaxation parameter.
  character ( len=12 ) :: lSolverT = 'bicgstab'      ! Linear algebraic solver.
  integer  :: maxiterT = 10                          ! Max number of iterations in linear solver.
  real(dp) :: tolAbsT = 1e-13                        ! Absolute residual level.
  real(dp) :: tolRelT = 0.02_dp                      ! Relative drop in residual to exit linear solver.
  real(dp) :: sigt = 1.0_dp                          ! sigma_t


  private 

  public :: calculate_temperature
  public :: calcT, urfT, gdsT, cSchemeT, dSchemeT, nrelaxT, lSolverT, maxiterT, tolAbsT, tolRelT, sigT ! params


contains


subroutine calculate_temperature
!
! Purpose:
!   Assemble and solve transport equation for temperature scalar field.
!

  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use linear_solvers
  use scalar_fluxes, only: facefluxsc  

  implicit none

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iWall
  integer :: iPer, l, if, iftwin
  real(dp) :: gam, prl, prtr, apotime
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: dcoef
  real(dp) :: fimax,fimin
  real(dp) :: are,nxf,nyf,nzf,dTn
  real(dp) :: urfr, urfm


! Variable specific coefficients:
  gam = gdsT
  prl = 1.0d0/pranl ! <-inverse Prandlt number of the working fluid.
  prtr = 1./sigt

! Calculate gradient: 
  call grad(T,dTdxi)

! Initialize matrix and source arrays
  a  = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! > Volume source terms
!
  do inp=1,numCells


    ! Unsteady Term
    if (ltransient) then

      apotime = den(inp)*vol(inp)/timestep

      if( bdf .or. cn ) then

        su(inp) = su(inp) + apotime*To(inp)
        sp(inp) = sp(inp) + apotime

      elseif( bdf2 ) then

        su(inp) = su(inp) + apotime*( 2*To(inp) - 0.5_dp*Too(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime

      endif

    endif


  enddo

  if( lturb ) then
    ! We need the freshest utt,vtt,wtt - turbulent heat fluxes
    call calcheatflux 
    call heatflux_source
  end if

!
! > Calculate terms integrated over faces
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    ! Diffusion coefficient
    dcoef = viscos*prl+( ( vis(ijp) + (vis(ijn)-vis(ijp))*facint(i) ) - viscos )*prtr

    call facefluxsc( i, ijp, ijn, &
                     xf(i), yf(i), zf(i), &
                     arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, cSchemeT, &
                     T, dTdxi, dcoef, cap, can, suadd )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap


    ! > Explicit part of convection and difussion fluxes

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

        dcoef = viscos*prl+(vis(ijp)-viscos)/sigt

        call facefluxsc( ijp, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         dTdxi, dcoef, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*T(ijb) + suadd

      end do


    elseif ( bctype(ib) == 'wall' .and. bcDFraction(7,ib)==1.0_dp ) then

      ! Isothermal wall boundaries (that's Dirichlet on temperature)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        dcoef = viscos/pranl+(vis(ijp)-viscos)/sigt
        dcoef = dcoef*srdw(iWall)

        sp(ijp) = sp(ijp) + dcoef
        su(ijp) = su(ijp) + dcoef*T(ijb)

      enddo

    elseif ( bctype(ib) == 'wall' .and. bcDFraction(7,ib)==0.0_dp ) then

      ! Adiabatic wall boundaries (that's actually zero grad on temperature)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        T(ijb) = T(ijp)

      enddo

    elseif ( bctype(ib) == 'wallQFlux') then

      ! Prescribed heat flux wall boundaries (that's actually non homo Neumann on temperature)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        dcoef = viscos/pranl+(vis(ijp)-viscos)/sigt
        dcoef = dcoef*srdw(iWall)

        ! Value of the temperature gradient in normal direction is set trough 
        ! proper choice of component values. Let's project in to normal direction
        ! to recover normal gradient.

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Gradient in face-normal direction        
        dTn = dTdxi(1,ijb)*nxf+dTdxi(2,ijb)*nyf+dTdxi(3,ijb)*nzf

        ! Explicit source
        su(ijp) = su(ijp) + dcoef*dTn

      enddo

    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face


        dcoef = viscos*prl + ( ( vis(ijp) + (vis(ijn)-vis(ijp))*facint(i) ) - viscos )*prtr

        ! face flux scalar but for periodic boundaries - it will be recognized by arguments
        call facefluxsc( i, ijp, ijn, &
                         xf(if), yf(if), zf(if), &
                         arx(if), ary(if), arz(if), &
                         flmass(if), gam, &
                         T, dTdxi, dcoef, cap, can, suadd )


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
      su(ijp) = su(ijp) - a(k)*To(ijn)

      k = jcell_icell_csr_index(i)
      su(ijn) = su(ijn) - a(k)*To(ijp)
    enddo
    
    do ijp=1,numCells
      apotime=den(ijp)*vol(ijp)/timestep
      off_diagonal_terms = sum( a( ia(ijp) : ia(ijp+1)-1 ) ) - a(diag(ijp))
      su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*To(ijp)
      sp(ijp) = sp(ijp)+apotime
    enddo

  endif

  ! Underrelaxation factors
  urfr = 1.0_dp / urfT
  urfm = 1.0_dp - urfT

  ! Main diagonal term assembly:
  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
    ! we substract it from the sum, to eliminate it from the sum.
    off_diagonal_terms  = sum( a(ia(inp) : ia(inp+1)-1) ) - a(diag(inp))
    a(diag(inp)) = sp(inp) - off_diagonal_terms

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = su(inp) + urfm*a(diag(inp))*T(inp)

  enddo

  ! Solve linear system:
  call csrsolve(lSolverT, T, su, resor(7), maxiterT, tolAbsT, tolRelT, 'Temp')

  
  ! Update field values at boundaries
  call updateBoundary( T )


! Report range of scalar values and clip if negative
  fimin = minval( T(1:numCells) )
  fimax = maxval( T(1:numCells) )
  
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= Temperature <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) T(1:numCells) = max( T(1:numCells), small ) 

end subroutine


end module temperature