module velocity

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing discretization and solution of momentum equations.
!
!  Discussion:
!
!    This one is called for laminar and turbulent flow.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types

  implicit none

  !
  ! Discrretization and solution parameters - modified trough input.nml file
  !
  logical  :: calcU = .True.                         ! To activate the solution of this field in the main function. 
  real(dp) :: urfU(3) = (/ 0.7, 0.7, 0.7 /)          ! Under-relaxation factors.
  real(dp) :: gdsU = 1.0                             ! Deferred correction factor.
  character ( len=30 ) :: cSchemeU = 'linearUpwind'  ! Convection scheme - default is second order upwind.
  character ( len=12 ) :: dSchemeU = 'skewness'      ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
  integer :: nrelaxU = 0  ! Type of non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed (-1/0/1).
  character ( len=12 ) :: lSolverU = 'bicgstab'  ! Linear algebraic solver.
  integer  :: maxiterU = 5                       ! Max number of iterations in linear solver.
  real(dp) :: tolAbsU = 1e-13                    ! Absolute residual level.
  real(dp) :: tolRelU = 0.025                    ! Relative drop in residual to exit linear solver.

  ! !...or =>
  ! type( FieldEquation ) :: Vel
  ! ! Now access data as: Vel%urf, Vel%cScheme, Vel%maxiter,...etc.
  
  private 

  public :: calcuvw, updateVelocityAtBoundary, calc_wall_shear
  public :: calcU, urfU, gdsU, cSchemeU, dSchemeU, nrelaxU, lSolverU, maxiterU, tolAbsU, tolRelU ! params



contains



subroutine calcuvw
!
!  Purpose: 
!
!    Assembles and solves discrete momentum equations.
!
!  Discussion:
!
!    Made for seqential solution method such as SIMPLE or PISO.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    --
!
!  Author:
!
!    Nikola Mirkov
!    Email: largeddysimulation@gmail.com/nikolamirkov@yahoo.com
!
!  Reference:
!
!    Ferziger & Perić, 
!    Mirkov et al JCP 2015
!
!  Parameters:
!
!    None.
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients, only: grad
  use linear_solvers
  use nablap
  use mhd
  use fieldManipulation, only: explDiv
  
  implicit none


!
! Local variables
!
  integer :: i, k, l, inp, ijp, ijn, ijb, iface
  integer :: ib,if,iftwin
  integer :: iWall, iSym, iPer      
  real(dp) :: apotime, heat
  real(dp) :: sup, svp, swp
  real(dp) :: sum_off_diagonal_terms
  real(dp) :: cp, cb        ! Temp. for eq. coeficients
  real(dp) :: cap, can      ! Temp. for eq. coeficients
  real(dp) :: cf            ! Diffusion coefficient
  real(dp) :: nxf, nyf, nzf ! Boundary face normals
  real(dp) :: are,arer      ! Face area, reciprocal area
  real(dp) :: dpb           ! distance cell center to face center
  real(dp) :: vsol          ! Diffusion coefficient (II)
  real(dp) :: Upb, Vpb, Wpb ! Velocity difference
  real(dp) :: viss
  real(dp) :: urfr,urfm

  ! logical :: ScndOrderWallBC_Model = .false.
  ! real(dp) :: FdUi,FdVi,FdWi! Diffusive flux auxiliary
  ! real(dp) :: Utp, Vtp, Wtp
  ! real(dp) :: Vnp,fdne

  ! Initialize sources
  su = 0.0_dp
  sv = 0.0_dp
  sw = 0.0_dp

  ! For u  sp => spu; for v  sp => spv; for w  sp => sp 
  spu = 0.0_dp
  spv = 0.0_dp
  sp  = 0.0_dp


  if ( piso .and. bdf ) then

    ! If you want to use midpoint timestepping method to improve piso to 2nd order,
    ! extrapolate using Adams-Bashfort method extrapolation to time interval midpoint, t^n+1/2
    ! flmass = 1.5_dp*flmasso - 0.5_dp*flmassoo
    ! p = 1.5_dp*po - 0.5_dp*poo
    ! u = 1.5_dp*uo - 0.5_dp*uoo
    ! v = 1.5_dp*vo - 0.5_dp*voo
    ! w = 1.5_dp*wo - 0.5_dp*woo

  ! For consistent 2nd order PISO algorithm.
  ! Estimate mass flux, velocity components and pressure for the current timestep 
  ! using 2nd order Adams-Bashfort extrapolation from previous two time steps.
  elseif ( piso .and. bdf2 ) then 

    flmass = 2*flmasso - flmassoo
    p = 2*po - poo
    u = 2*uo - uoo
    v = 2*vo - voo
    w = 2*wo - woo

  elseif ( piso .and. bdf3 ) then 

    ! Extrapolation from three time levels to n+1 time level
    flmass = 3*flmasso - 8*flmassoo + 6*flmassooo
    p = 3*po - 8*poo + 6*pooo
    u = 3*uo - 8*uoo + 6*uooo
    v = 3*vo - 8*voo + 6*vooo
    w = 3*wo - 8*woo + 6*wooo

  endif

  ! Velocity gradients: 
  call updateVelocityAtBoundary
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  ! Update values of Lorentz force vector components in MHD case.
  if( calcEpot)  call calculate_Lorentz_force

  !----------------------------------------------------------------------------
  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  !----------------------------------------------------------------------------

  do inp=1,numCells

    !
    ! Constant mass flow forcing - used only on U velocity component
    !
    if(const_mflux) then
      su(inp) = su(inp) + gradPcmf*vol(inp)
    endif

    !
    ! Buoyancy source terms
    !
    if( lbuoy ) then

      if( boussinesq ) then
        ! Boussinesq-ova aproximacija:
        heat = beta*densit*(t(inp)-tref)*vol(inp)
      else
        heat = (densit-den(inp))*vol(inp)
      endif

      su(inp) = su(inp) - gravx*heat
      sv(inp) = sv(inp) - gravy*heat
      sw(inp) = sw(inp) - gravz*heat

    endif

    !
    ! MHD: Lorentz force source terms
    !
    if( calcEpot )  then

      ! Add Lorentz force volume source terms.
      su(inp) = su(inp) + sigma*florx(inp)*vol(inp)
      sv(inp) = sv(inp) + sigma*flory(inp)*vol(inp)
      sw(inp) = sw(inp) + sigma*florz(inp)*vol(inp)  

    endif


    !
    ! Unsteady term
    !
    if( ltransient ) then

      if( bdf .or. cn ) then
      !
      ! Backward differentiation formula of 1st order.
      !
        apotime = den(inp)*vol(inp)/timestep

        ! RHS vector contribution
        su(inp) = su(inp) + apotime*uo(inp)
        sv(inp) = sv(inp) + apotime*vo(inp)
        sw(inp) = sw(inp) + apotime*wo(inp)

        ! Matrix diagonal element contribution
        spu(inp) = spu(inp) + apotime
        spv(inp) = spv(inp) + apotime
        sp(inp)  = sp(inp)  + apotime

      elseif( bdf2 ) then
      !
      ! Three Level Implicit Time Integration (BDF2) - 2nd order.
      !
        apotime = den(inp)*vol(inp)/timestep

        ! RHS vector contribution
        su(inp) = su(inp) + apotime*( 2*uo(inp) - 0.5_dp*uoo(inp) )
        sv(inp) = sv(inp) + apotime*( 2*vo(inp) - 0.5_dp*voo(inp) ) 
        sw(inp) = sw(inp) + apotime*( 2*wo(inp) - 0.5_dp*woo(inp) )

        ! Matrix diagonal element contribution
        spu(inp) = spu(inp) + 1.5_dp*apotime
        spv(inp) = spv(inp) + 1.5_dp*apotime
        sp(inp)  = sp(inp)  + 1.5_dp*apotime


      elseif( bdf3 ) then

        apotime = den(inp)*vol(inp)/timestep

        ! RHS vector contribution
        su(inp) = su(inp) + apotime*( 3*uo(inp) - 1.5_dp*uoo(inp) + 1./3.0_dp*uooo(inp) )
        sv(inp) = sv(inp) + apotime*( 3*vo(inp) - 1.5_dp*voo(inp) + 1./3.0_dp*vooo(inp) ) 
        sw(inp) = sw(inp) + apotime*( 3*wo(inp) - 1.5_dp*woo(inp) + 1./3.0_dp*wooo(inp) )

        ! Matrix diagonal element contribution
        spu(inp) = spu(inp) + 11./6.0_dp*apotime
        spv(inp) = spv(inp) + 11./6.0_dp*apotime
        sp(inp)  = sp(inp)  + 11./6.0_dp*apotime

      endif

    endif ! unsteady term

  end do

  !
  ! > Fluxes trough faces:
  !

  ! Inner faces
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxuvw(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), flmass(i), facint(i), gdsU, &
      cap, can, sup, svp, swp)

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! > Sources: 

    su(ijp) = su(ijp) + sup
    sv(ijp) = sv(ijp) + svp
    sw(ijp) = sw(ijp) + swp

    su(ijn) = su(ijn) - sup
    sv(ijn) = sv(ijn) - svp
    sw(ijn) = sw(ijn) - swp

  end do 


  iSym = 0
  iWall = 0
  iPer = 0
  l = numInnerFaces
      
  ! Implement boundary conditions
 
  do ib=1,numBoundaries

    if (  bctype(ib) == 'inlet' .or. &
          bctype(ib) == 'outlet' .or. &
          bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxuvw_bnd( ijp, ijb, &
                              xf(iface), yf(iface), zf(iface), &
                              arx(iface), ary(iface), arz(iface), &
                              flmass(iface), &
                              cp, cb, sup, svp, swp)

        spu(ijp) = spu(ijp) - cb
        spv(ijp) = spv(ijp) - cb
        sp(ijp)  = sp(ijp)  - cb

        su(ijp) = su(ijp) - cb*u(ijb) + sup
        sv(ijp) = sv(ijp) - cb*v(ijb) + svp
        sw(ijp) = sw(ijp) - cb*w(ijb) + swp

      end do

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i


        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
        arer = 1.0_dp/are

        ! Face normals
        nxf = arx(iface)*arer
        nyf = ary(iface)*arer
        nzf = arz(iface)*arer

        ! dpb - Dist from cc of owner cell to cf @boundary,
        ! cannot expect dpb to be normal to boundary face in general...
        ! ->    ->
        ! dpb . nf
        dpb = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        cf = 2*vis(inp)*are/dpb

        su(ijp) = su(ijp)-cf*nxf*(nyf*v(ijp)+nzf*w(ijp))
        sv(ijp) = sv(ijp)-cf*nyf*(nxf*u(ijp)+nzf*w(ijp))
        sw(ijp) = sw(ijp)-cf*nzf*(nxf*u(ijp)+nyf*v(ijp))

        spu(ijp) = spu(ijp)+cf*nxf**2
        spv(ijp) = spv(ijp)+cf*nyf**2
        sp(ijp)  = sp(ijp) +cf*nzf**2


      end do


    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face


        call facefluxuvw_periodic(ijp, ijn, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), flmass(if), gdsU, &
          cap, can, sup, svp, swp)

        ! > Off-diagonal elements:

        ! l is in interval [numInnerFaces+1, numInnerFaces+numPeriodic]
        l = l + 1

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_index(l)
        a(k) = can

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_index(l)
        a(k) = cap

        ! > Sources: 

        su(ijp) = su(ijp) + sup
        sv(ijp) = sv(ijp) + svp
        sw(ijp) = sw(ijp) + swp

        su(ijn) = su(ijn) - sup
        sv(ijn) = sv(ijn) - svp
        sw(ijn) = sw(ijn) - swp


      end do 


    elseif ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        ! viss = viscos ! viskoznost interpolirana na boundary face
        ! if(lturb.and.ypl(iWall).gt.ctrans) viss=visw(iWall)
        viss=max(viscos,visw(iWall))

        ! Face area
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
        arer = 1.0_dp/are

        ! Face normals
        nxf = arx(iface)*arer
        nyf = ary(iface)*arer
        nzf = arz(iface)*arer

        ! dpb - Dist from cc of owner cell to cf @boundary,
        ! cannot expect dpb to be normal to boundary face in general...
        ! ->    ->
        ! dpb . nf
        dpb = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        ! Diffusion coef. 
        vsol = viss*are/dpb

        ! Velocity difference vector components
        upb = u(ijp)-u(ijb)
        vpb = v(ijp)-v(ijb)
        wpb = w(ijp)-w(ijb)

        ! Explicit contribution to main diagonal
        spu(ijp) = spu(ijp) + vsol*(1.-nxf**2)
        spv(ijp) = spv(ijp) + vsol*(1.-nyf**2)
        sp(ijp)  = sp(ijp)  + vsol*(1.-nzf**2)

        su(ijp) = su(ijp) + vsol*( u(ijb)*(1.-nxf**2) + vpb*nyf*nxf         + wpb*nzf*nxf )
        sv(ijp) = sv(ijp) + vsol*( upb*nxf*nyf        + v(ijb)*(1.-nyf**2)  + wpb*nzf*nyf )
        sw(ijp) = sw(ijp) + vsol*( upb*nxf*nzf        + vpb*nyf*nzf         + w(ijb)*(1.-nzf**2) )

      enddo

      ! *** Another way to implement wall bc - based on the code of Gabriel Usera's (Uruguay) caffa3dmbri code.

      ! do i=1,nfaces(ib)

      !   iface = startFace(ib) + i
      !   ijp = owner(iface)
      !   ijb = iBndValueStart(ib) + i

      !   viss = viscos ! viskoznost interolirana na boundary face
      !   if(lturb.and.ypl(i).gt.ctrans) viss=visw(i)
      
      !   cf=viss*srdw(i) ! cf v.s. vsol -> cf is calculated using normal distance in srdw!

      !   ! Face area 
      !   are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

      !   ! Face normals
      !   nxf = arx(iface)/are
      !   nyf = ary(iface)/are
      !   nzf = arz(iface)/are

      !   ! Dist from cc of owner cell to cf @boundary, cannot expect dpb to be normal to boundary face in general
      !   dpb = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

      !   ! Diffusion coef. 
      !   vsol = viss*are/dpb

      !   ! Velocity difference vector components
      !   upb = u(ijp)-u(ijb)
      !   vpb = v(ijp)-v(ijb)
      !   wpb = w(ijp)-w(ijb)

      !   ! Velocity difference vector projected to wall face normal.
      !   vnp = upb*nxf+vpb*nyf+wpb*nzf

      !   ! Velocity difference in tangential direction.
      !   utp = upb-vnp*nxf
      !   vtp = vpb-vnp*nyf
      !   wtp = wpb-vnp*nzf

      !   if (ScndOrderWallBC_Model) then

      !     ! Eksplicitna difuzija
      !     FdUi=viss*((dUdxi(1,ijp)+dUdxi(1,ijp))*nxf+(dUdxi(2,ijp)+dVdxi(1,ijp))*nyf+(dUdxi(3,ijp)+dWdxi(1,ijp))*nzf)
      !     FdVi=viss*((dVdxi(1,ijp)+dUdxi(2,ijp))*nxf+(dVdxi(2,ijp)+dVdxi(2,ijp))*nyf+(dVdxi(3,ijp)+dWdxi(2,ijp))*nzf)
      !     FdWi=viss*((dWdxi(1,ijp)+dUdxi(3,ijp))*nxf+(dWdxi(2,ijp)+dVdxi(3,ijp))*nyf+(dWdxi(3,ijp)+dWdxi(3,ijp))*nzf)
      !     ! Projektujes eksplicitnu difuziju na nomalu
      !     FdNe = FdUi*nxf + FdVi*nyf + FdWi*nzf
      !     ! oduzmes od eksplicitne difuzije onu komponentu duz normale
      !     FdUi = FdUi-FdNe*nxf
      !     FdVi = FdVi-FdNe*nyf
      !     FdWi = FdWi-FdNe*nzf

      !     spu(ijp) = spu(ijp) + vsol
      !     spv(ijp) = spv(ijp) + vsol
      !     sp(ijp)  = sp(ijp)  + vsol

      !     Su(ijp) = Su(ijp)+Vsol*U(ijp)-(2*cf*Utp+FdUi*Are)
      !     Sv(ijp) = Sv(ijp)+Vsol*V(ijp)-(2*cf*Vtp+FdVi*Are)
      !     Sw(ijp) = Sw(ijp)+Vsol*W(ijp)-(2*cf*Wtp+FdWi*Are)

      !   else

      !     spu(ijp) = spu(ijp) + vsol
      !     spv(ijp) = spv(ijp) + vsol
      !     sp(ijp)  = sp(ijp)  + vsol

      !     su(ijp) = su(ijp) + vsol*u(ijp) - cf*utp
      !     sv(ijp) = sv(ijp) + vsol*v(ijp) - cf*vtp
      !     sw(ijp) = sw(ijp) + vsol*w(ijp) - cf*wtp

      !   endif

      ! enddo


    endif 

  enddo


  ! For PISO, we need to save only EXPLICIT TERMS OF CONVECTION AND DIFFUSION.
  if (piso) then
    rU = su
    rV = sv
    rW = sw
  endif

  !
  ! If flow is compressible add 2/3*mu*div(U) to pressure field
  !
  if (compressible) p(1:numCells) = p(1:numCells) + 2./3.*vis(1:numCells)*explDiv( U, V, W )


  ! Contribution to source of the -grad(p) term
  if (CN) then

    call surfaceIntegratePressureCrankNicolson

    ! Modify coefficients for Crank-Nicolson (doesn't affect the main diagonal because is still zero).
    a = 0.5_dp*a 
    
  else

    call surfaceIntegratePressure

  endif



!
! Assemble and solve system for U component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if ( CN ) then

    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_index(i)
        su(ijp) = su(ijp) - a(k)*uo(ijn)

        k = jcell_icell_csr_index(i)
        su(ijn) = su(ijn) - a(k)*uo(ijp)
    enddo

    do ijp=1,numCells
        apotime = den(ijp)*vol(ijp)/timestep
        sum_off_diagonal_terms = sum( a(ioffset(ijp) : ioffset(ijp+1)-1) ) - a(diag(ijp))
        su(ijp) = su(ijp) + (apotime + sum_off_diagonal_terms)*uo(ijp)
        spu(ijp) = spu(ijp) + apotime
    enddo

  endif

  urfr = 1.0_dp/urfU(1)
  urfm = 1.0_dp-urfU(1)

  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
    ! we substract it from the sum, to eliminate it from the sum.
    ! We could also write sum( a(ioffset(inp)) : a(ioffset(inp+1)-1) ) because all diagonal terms are zero.
    sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp)) 
    
    a(diag(inp)) = spu(inp) - sum_off_diagonal_terms

    apu(inp) = 1./(a(diag(inp))+small)

    ! Alternative way, loop instead of array slice:
    ! a(diag(inp)) = spu(inp) 
    ! do k = ioffset(inp),ioffset(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo


    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = su(inp) + urfm*a(diag(inp))*u(inp)


  enddo

  ! Solve linear system of equations
  call csrsolve(lSolverU, u, su, resor(1), maxiterU, tolAbsU, tolRelU, 'U')

!
! Assemble and solve system for V component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if( CN ) then

    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_index(i)
        sv(ijp) = sv(ijp) - a(k)*vo(ijn)

        k = jcell_icell_csr_index(i)
        sv(ijn) = sv(ijn) - a(k)*vo(ijp)
    enddo

    do ijp=1,numCells
        apotime=den(ijp)*vol(ijp)/timestep
        sum_off_diagonal_terms = sum( a(ioffset(ijp) : ioffset(ijp+1)-1) ) - a(diag(ijp))
        sv(ijp) = sv(ijp) + (apotime + sum_off_diagonal_terms)*vo(ijp)
        spv(ijp) = spv(ijp)+apotime
    enddo

  endif

  ! Clear main diagonal and source vector
  do inp=1,numCells
    a(diag(inp)) = 0.0_dp
    su(inp) = 0.0_dp
  enddo

  ! Under-relaxation factors
  urfr = 1.0_dp/urfU(2)
  urfm = 1.0_dp-urfU(2)

  do inp = 1,numCells

    ! Main diagonal term assembly:
    sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp))

    a(diag(inp)) = spv(inp) - sum_off_diagonal_terms

    apv(inp) = 1./(a(diag(inp))+small)

    ! Alternative way, loop instead of array slice:
    ! a(diag(inp)) = spv(inp) 
    ! do k = ioffset(inp),ioffset(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = sv(inp) + urfm*a(diag(inp))*v(inp)

  enddo

  ! Solve linear system of equations
  call csrsolve(lSolverU, v, su, resor(2), maxiterU, tolAbsU, tolRelU, 'V') 

!
! Assemble and solve system for W component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if( CN ) then

    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_index(i)
        sw(ijp) = sw(ijp) - a(k)*wo(ijn)

        k = jcell_icell_csr_index(i)
        sw(ijn) = sw(ijn) - a(k)*wo(ijp)
    enddo

    do ijp=1,numCells
        apotime = den(ijp)*vol(ijp)/timestep
        sum_off_diagonal_terms = sum( a(ioffset(ijp) : ioffset(ijp+1)-1) ) - a(diag(ijp))
        sw(ijp) = sw(ijp) + (apotime + sum_off_diagonal_terms)*wo(ijp)
        sp(ijp) = sp(ijp) + apotime
    enddo

  endif 

  ! Clear main diagonal and source vector
  do inp=1,numCells
    a(diag(inp)) = 0.0_dp
    su(inp) = 0.0_dp
  enddo
  
  ! Set under-relaxation factors
  urfr = 1.0_dp/urfU(3)
  urfm = 1.0_dp-urfU(3)

  do inp = 1,numCells

    ! Main diagonal term assembly:
    sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp)) 
   
    a(diag(inp)) = sp(inp) - sum_off_diagonal_terms

    apw(inp) = 1./(a(diag(inp))+small)

    ! Alternative way, loop instead of array slice:
    ! a(diag(inp)) = sp(inp) 
    ! do k = ioffset(inp),ioffset(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = sw(inp) + urfm*a(diag(inp))*w(inp)

  enddo

  ! Solve linear system of equations
  call csrsolve(lSolverU, w, su, resor(3), maxiterU, tolAbsU, tolRelU, 'W')

end subroutine



subroutine facefluxuvw(ijp, ijn, xf, yf, zf, arx, ary, arz, flomass, lambda, gam, cap, can, sup, svp, swp)
!
!  Puprose: 
!   Face fluxes of momentum eq for inner faces.
!
! 
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables
  use gradients, only: sngrad
  use interpolation, only: face_value

  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

  ! Local variables
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: duxi,duyi,duzi, &
              dvxi,dvyi,dvzi, &
              dwxi,dwyi,dwzi
  real(dp) :: duxii,dvxii,dwxii, &
              duyii,dvyii,dwyii, &
              duzii,dvzii,dwzii
  real(dp) :: de, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
!----------------------------------------------------------------------


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)



  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijp)*fxp+vis(ijn)*fxn

  ! Difusion coefficient
  ! de = game*are/dpn
  de = game*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)

  ! > Equation coefficients - implicit diffusion and convection
  ce = min(flomass,zero) 
  cp = max(flomass,zero)

  can = -de + ce
  cap = -de - cp



  ! > Explicit diffusion: 


  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              u, dudxi, nrelaxU, dSchemeU, duxi, duyi, duzi,  &
              duxii, duyii, duzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              v, dvdxi, nrelaxU, dSchemeU, dvxi, dvyi, dvzi, &
              dvxii, dvyii, dvzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              w, dwdxi, nrelaxU, dSchemeU, dwxi, dwyi, dwzi, &
              dwxii, dwyii, dwzii)

!---------------------------------------------------------------------------------------
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector:
!     su = su + (fdue-fdui)
!     sv = sv + (fdve-fdvi)
!     sw = sw + (fdwe-fdwi)
!---------------------------------------------------------------------------------------

  ! Explicit diffussion: 
  fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
  fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
  fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )

  ! Implicit diffussion:
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)



  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds = cp*u(ijp)+ce*u(ijn)
  fvuds = cp*v(ijp)+ce*v(ijn)
  fwuds = cp*w(ijp)+ce*w(ijn)


! EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
! Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_f[kg/s] * Phi_f[m/s]$

  if( flomass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    ue = face_value(ijp, ijn, xf, yf, zf, fxp, u, dUdxi, cSchemeU)
    ve = face_value(ijp, ijn, xf, yf, zf, fxp, v, dVdxi, cSchemeU)
    we = face_value(ijp, ijn, xf, yf, zf, fxp, w, dWdxi, cSchemeU)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ue = face_value(ijn, ijp, xf, yf, zf, fxn, u, dUdxi, cSchemeU)
    ve = face_value(ijn, ijp, xf, yf, zf, fxn, v, dVdxi, cSchemeU)
    we = face_value(ijn, ijp, xf, yf, zf, fxn, w, dWdxi, cSchemeU)
  endif

  fuhigh = flomass*ue
  fvhigh = flomass*ve
  fwhigh = flomass*we



! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = -gam*(fuhigh-fuuds)+fdue-fdui
  svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  swp = -gam*(fwhigh-fwuds)+fdwe-fdwi

end subroutine



subroutine facefluxuvw_bnd(ijp, ijb, xf, yf, zf, arx, ary, arz, flomass, cap, can, sup, svp, swp)
!
!
! Facefluxuvw routine used for inlet and outlet boundaries. 
! Constant gradient is assumed between cell center and boundary face center.
!

  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables

  implicit none

  integer, intent(in) :: ijp, ijb
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ixi1,ixi2,ixi3
  real(dp) :: dpn,costheta,costn
  real(dp) :: xi,yi,zi

  real(dp) :: duxi,duyi,duzi, &
              dvxi,dvyi,dvzi, &
              dwxi,dwyi,dwzi

  real(dp) :: duxii,dvxii,dwxii, &
              duyii,dvyii,dwyii, &
              duzii,dvzii,dwzii

  real(dp) :: d2x,d2y,d2z,d1x,d1y,d1z

  real(dp) :: de, vole, game
  real(dp) :: fxp,fxn
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
!----------------------------------------------------------------------

  ! [NOTE]: Leaving things commented out so people can see the difference 
  !         between this and the version for iner faces.

  ! > Geometry:

  ! Face interpolation factor:
  ! fxn=lambda 
  ! fxp=1.0_dp-lambda
  fxn=1.0_dp
  fxp=0.0_dp

  ! Distance vector between cell centers:
  ! xpn=xc(ijb)-xc(ijp)
  ! ypn=yc(ijb)-yc(ijp)
  ! zpn=zc(ijb)-zc(ijp)
  xpn=xf-xc(ijp)
  ypn=yf-yc(ijp)
  zpn=zf-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! Minimal correction: nrelax = +1 :
  !costn = costheta
  ! Orthogonal correction: nrelax =  0 : 
  costn = 1.0_dp
  ! Over-relaxed approach: nrelax = -1 :
  !costn = 1./costheta
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! dpn * sf
  vole=xpn*arx+ypn*ary+zpn*arz

  ! Overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn

  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn


  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijb) !<- it's (vis(ijp)*fxp+vis(ijb)*fxn) with fxn = 1.0

  ! Difusion coefficient
  de = game*are/dpn


  ! Equation coefficients - implicit diffusion and convection
  can = -de + min(flomass,zero)
  cap = -de - max(flomass,zero)


  ! > Face velocity components and Explicit diffusion: 

  ! Coordinates of point j'
  xi = xf
  yi = yf
  zi = zf


  ! interpolate gradients defined at cv centers to faces

  ![NOTE]: Commented out version is for inner faces..
  !        fxn = 1.0, so we should take values of gradient from
  !        face center (ijb index, b as 'boundary'), 
  !        but since we have constant gradient for
  !        inlet and outlet, for which this routine is called,
  !        we take adjecent cell center value of gradient instead.

  ! duxi = dUdxi(1,ijp)*fxp+dUdxi(1,ijb)*fxn
  ! duyi = dUdxi(2,ijp)*fxp+dUdxi(2,ijb)*fxn
  ! duzi = dUdxi(3,ijp)*fxp+dUdxi(3,ijb)*fxn
  duxi = dUdxi(1,ijp)
  duyi = dUdxi(2,ijp)
  duzi = dUdxi(3,ijp) !...because constant gradient

  ! du/dx_i interpolated at cell face:
  duxii = duxi*d1x + arx/vole*( u(ijb)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
  duyii = duyi*d1y + ary/vole*( u(ijb)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
  duzii = duzi*d1z + arz/vole*( u(ijb)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 


  ! dvxi = dVdxi(1,ijp)*fxp+dVdxi(1,ijb)*fxn
  ! dvyi = dVdxi(2,ijp)*fxp+dVdxi(2,ijb)*fxn
  ! dvzi = dVdxi(3,ijp)*fxp+dVdxi(3,ijb)*fxn
  dvxi = dVdxi(1,ijp)
  dvyi = dVdxi(2,ijp)
  dvzi = dVdxi(3,ijp) !...because constant gradient

  ! dv/dx_i interpolated at cell face:
  dvxii = dvxi*d1x + arx/vole*( v(ijb)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
  dvyii = dvyi*d1y + ary/vole*( v(ijb)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
  dvzii = dvzi*d1z + arz/vole*( v(ijb)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 


  ! dwxi = dWdxi(1,ijp)*fxp+dWdxi(1,ijb)*fxn
  ! dwyi = dWdxi(2,ijp)*fxp+dWdxi(2,ijb)*fxn
  ! dwzi = dWdxi(3,ijp)*fxp+dWdxi(3,ijb)*fxn
  dwxi = dWdxi(1,ijp)
  dwyi = dWdxi(2,ijp)
  dwzi = dWdxi(3,ijp) !...because constant gradient

  ! dw/dx_i interpolated at cell face:
  dwxii = dwxi*d1x + arx/vole*( w(ijb)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
  dwyii = dwyi*d1y + ary/vole*( w(ijb)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
  dwzii = dwzi*d1z + arz/vole*( w(ijb)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! explicit diffussion 
  fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
  fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
  fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )

  ! implicit diffussion 
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  ! > Explicit convection: 
  ! - None -


  ! Explicit part of diffusion fluxes and sources due to deffered correction.

  ![NOTE]: Deferred correction blending coefficient (gam) is taken to be zero
  !        so we have only explicit diffusion below.
  ! sup = -gam*(fuhigh-fuuds)+fdue-fdui
  ! svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  ! swp = -gam*(fwhigh-fwuds)+fdwe-fdwi
  sup = fdue-fdui
  svp = fdve-fdvi
  swp = fdwe-fdwi 

end subroutine



subroutine facefluxuvw_periodic(ijp, ijn, xf, yf, zf, arx, ary, arz, flomass, gam, cap, can, sup, svp, swp)
!
!  Puprose: 
!   Face fluxes of momentum eq for periodic boundary faces.
!
!  Discussion:
!    ijp and ijn are paired cells on different sides of the domain.
!
! 
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables
  use interpolation, only: face_value_cds

  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

  ! Local variables
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: duxi,duyi,duzi,dvxi,dvyi,dvzi,dwxi,dwyi,dwzi
  real(dp) :: de, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
!----------------------------------------------------------------------


  ! > Geometry:

  ! Face interpolation factor
  fxn = 0.5_dp ! Assumption for periodic boundaries
  fxp = fxn

  ! Distance vector between cell centers
  xpn = 2*( xf-xc(ijp) )
  ypn = 2*( yf-yc(ijp) )
  zpn = 2*( zf-zc(ijp) )

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)


  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijp)*fxp+vis(ijn)*fxn

  ! Difusion coefficient
  ! de = game*are/dpn
  de = game*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)

  ! > Equation coefficients - implicit diffusion and convection
  ce = min(flomass,zero) 
  cp = max(flomass,zero)

  can = -de + ce
  cap = -de - cp



  ! > Explicit diffusion: 
  !-----------------------------------------------------------------------------
  !     We calculate explicit and implicit diffusion fde and fdi,
  !     later we put their difference (fde-fdi) to rhs vector:
  !     su = su + (fdue-fdui)
  !     sv = sv + (fdve-fdvi)
  !     sw = sw + (fdwe-fdwi)
  !-----------------------------------------------------------------------------

  ! Interpolate gradients defined at CV centers to faces
  dUxi = dUdxi(1,ijp)*fxp+dUdxi(1,ijn)*fxn
  dUyi = dUdxi(2,ijp)*fxp+dUdxi(2,ijn)*fxn
  dUzi = dUdxi(3,ijp)*fxp+dUdxi(3,ijn)*fxn

  dVxi = dVdxi(1,ijp)*fxp+dVdxi(1,ijn)*fxn
  dVyi = dVdxi(2,ijp)*fxp+dVdxi(2,ijn)*fxn
  dVzi = dVdxi(3,ijp)*fxp+dVdxi(3,ijn)*fxn

  dWxi = dWdxi(1,ijp)*fxp+dWdxi(1,ijn)*fxn
  dWyi = dWdxi(2,ijp)*fxp+dWdxi(2,ijn)*fxn
  dWzi = dWdxi(3,ijp)*fxp+dWdxi(3,ijn)*fxn

  ! Explicit diffussion: 
  fdue = game*( (duxi+duxi)*arx + (duyi+dvxi)*ary + (duzi+dwxi)*arz )
  fdve = game*( (duyi+dvxi)*arx + (dvyi+dvyi)*ary + (dvzi+dwyi)*arz )
  fdwe = game*( (duzi+dwxi)*arx + (dwyi+dvzi)*ary + (dwzi+dwzi)*arz )

  ! Implicit diffussion:
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)


  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds = cp*u(ijp)+ce*u(ijn)
  fvuds = cp*v(ijp)+ce*v(ijn)
  fwuds = cp*w(ijp)+ce*w(ijn)


! EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
! Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_f[kg/s] * Phi_f[m/s]$

  if( flomass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    ue = face_value_cds(ijp, ijn, fxp, u) 
    ve = face_value_cds(ijp, ijn, fxp, v) 
    we = face_value_cds(ijp, ijn, fxp, w) 
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ue = face_value_cds(ijn, ijp, fxn, u) 
    ve = face_value_cds(ijn, ijp, fxn, v) 
    we = face_value_cds(ijn, ijp, fxn, w) 
  endif

  fuhigh = flomass*ue
  fvhigh = flomass*ve
  fwhigh = flomass*we



! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = -gam*(fuhigh-fuuds)+fdue-fdui
  svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  swp = -gam*(fwhigh-fwuds)+fdwe-fdwi

end subroutine




subroutine updateVelocityAtBoundary
!  
!******************************************************************************
!
!     Updates values at symmetry boundaries 
! 
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

!
!     Local variables
!
  integer :: i,ijp,ijb,ib,iface,iper
  real(dp) :: Unmag,flowo

  ! Update velocity components along outlet boundaries
  ! and correct mass flux to satisfy global mass conservation

  flowo=0.0_dp

  iper = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'empty' .or. bctype(ib) == 'periodic' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        U(ijb) = U(ijp)
        V(ijb) = V(ijp)
        W(ijb) = W(ijp)

      enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! Project velocity vector to face normal direction:
        Unmag = u(ijp)*arx(iface)+v(ijp)*ary(iface)+w(ijp)*arz(iface)

        U(ijb) = U(ijp)-Unmag*arx(iface)
        V(ijb) = V(ijp)-Unmag*ary(iface)
        W(ijb) = W(ijp)-Unmag*arz(iface)

      end do

    ! elseif (  bctype(ib) == 'periodic' ) then

    !   iPer = iPer + 1

    !   ! Faces trough periodic boundaries, Taiwo first
    !   do i=1,nfaces(ib)

    !     iface = startFace(ib) + i
    !     ijp = owner(iface)
    !     ijb = iBndValueStart(ib) + i

    !     iface = startFaceTwin(iPer) + i
    !     ijn = owner(iface)

    !     U(ijb) = half*( U(ijp)+U(ijn) )
    !     V(ijb) = half*( U(ijp)+U(ijn) )
    !     W(ijb) = half*( U(ijp)+U(ijn) )

    !     ! Now find where is the twin in field array
    !     ijbt = numCells + ( startFaceTwin(iPer) - numInnerFaces ) + i
        
    !     ! Twin takes the same values
    !     U(ijbt) = U(ijb)
    !     V(ijbt) = V(ijb)
    !     W(ijbt) = W(ijb)

    !   enddo


    endif 

  enddo

end subroutine

subroutine calc_wall_shear
!
! Purpose: 
! Returns wall shear stress - tau and non-dimensional wall distance of the first cell layer - y+.
! Usually needed for post processing.
!
  use types
  use parameters
  use geometry
  use variables

  implicit none 

  integer :: iWall, ib, i, iface, ijp, ijb
  real(dp) :: are 
  real(dp) :: nxf,nyf,nzf 
  real(dp) :: upb,vpb,wpb,vsol,fshearx,fsheary,fshearz

  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Velocity difference vector components
        upb = u(ijp)-u(ijb)
        vpb = v(ijp)-v(ijb)
        wpb = w(ijp)-w(ijb)

        vsol = max(viscos,visw(iWall))*srdw(iWall)

        ! Shear forces at wall in x, y and z direction.
        fshearx = vsol * ( (u(ijb)-u(ijp))*(1.-nxf**2) + vpb*nyf*nxf                  + wpb*nzf*nxf )
        fsheary = vsol * ( upb*nxf*nyf                 + (v(ijb)-v(ijp))*(1.-nyf**2)  + wpb*nzf*nyf )
        fshearz = vsol * ( upb*nxf*nzf                 + vpb*nyf*nzf                  + (w(ijb)-w(ijp))*(1.-nzf**2) )

        tau(iWall) = sqrt(fshearx**2+fsheary**2+fshearz**2)/are

        ! ypl(iWall) = den(ijb)*sqrt( Tau(iWall) / den(ijb) )*dnw(iWall)/viscos ! reduce to =>
          ypl(iWall) =          sqrt( Tau(iWall) * den(ijb) )*dnw(iWall)/viscos


      enddo

      ! ! Let's print the skin friction coefficient - Cf
      ! ! Below we print Cf_wall v.s. Rex ( this is hardcoded for turb flat plate case)
      ! write(*,'(a)') ' '  
      ! write(*,'(a)') '  Skin friction coefficient: '//trim( bcname(ib) )  
      ! write(*,'(a)') ' '  

      ! iwall = iwall - nfaces(ib)
      ! do i=1,nfaces(ib)
      !   iface = startFace(ib) + i
      !   iWall = iWall + 1       
      !   write(*,*) rhoref*xf(iFace)*Uref/viscos, tau(iWall) / (0.5*rhoref*uref**2)
      ! enddo



    endif 

  enddo


end subroutine


end module


