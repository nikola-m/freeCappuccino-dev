!***********************************************************************
!
subroutine calcuvw
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use hcoef
  use gradients, only: grad
  use faceflux_velocity
  use fieldManipulation, only: calcPressDiv
  use mhd

  implicit none
!
!***********************************************************************

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, iface
  ! integer :: istage
  real(dp) :: urfrs, urfms, apotime, heat
  real(dp) :: sup, svp, swp
  real(dp) :: sum_off_diagonal_terms

  ! logical :: ScndOrderWallBC_Model
  integer :: ib
  integer :: iWall, iSym
  integer :: ijb            ! Boundary field value indexes
  real(dp) :: cp, cb        ! Temp. for eq. coeficients
  real(dp) :: cap, can      ! Temp. for eq. coeficients
  real(dp) :: vsi           ! Interpolated dynamic viscosity
  real(dp) :: cf            ! Diffusion coefficient
  real(dp) :: nxf, nyf, nzf ! Boundary face normals
  real(dp) :: Are           ! Face area
  real(dp) :: dpb           ! distance cell center to face center
  real(dp) :: vsol          ! Diffusion coefficient (II)
  real(dp) :: fdne          ! Diffusive flux auxiliary
  real(dp) :: Upb, Vpb, Wpb ! Velocity difference
  real(dp) :: viss
  ! real(dp) :: FdUi,FdVi,FdWi! Diffusive flux auxiliary
  ! real(dp) :: Utp, Vtp, Wtp
  ! real(dp) :: Vnp

  
  ! ScndOrderWallBC_Model = .true.

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
    ! extrapolate using Adams-Abshfort method extrapolation to time interval midpoint, t^n+1/2
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
  if(lcal(iep)) call calculate_Lorentz_force

  !----------------------------------------------------------------------------
  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  !----------------------------------------------------------------------------

  do inp=1,numCells

    !
    ! Constant mass flow forcing - used only on U velocity component
    !
    if(const_mflux) su(inp) = su(inp) + gradPcmf*vol(inp)

    !
    ! Buoyancy source terms
    !
    if(lcal(ien).and.lbuoy) then

      if(boussinesq) then
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
    if(lcal(iep)) then

      ! Add Lorentz force volume source terms.
      su(inp) = su(inp) + sigma*florx(inp)*vol(inp)
      sv(inp) = sv(inp) + sigma*flory(inp)*vol(inp)
      sw(inp) = sw(inp) + sigma*florz(inp)*vol(inp)  

    endif


    !
    ! Unsteady term
    !
    if(ltransient) then

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

  !----------------------------------------------------------------------------
  ! Calculate Reynols stresses explicitly and additional asm terms:
  !----------------------------------------------------------------------------
  ! if(lturb) then
  !   call calcstress
  !   if (lasm) call Additional_algebraic_stress_terms
  ! end if




  !----------------------------------------------------------------------------
  ! CALCULATE TERMS INTEGRATED OVER FACES
  !----------------------------------------------------------------------------

  !
  ! > Fluxes trough faces:
  !

  ! Inner faces
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxuvw(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), flmass(i), facint(i), gds(iu), &
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

      
  ! Implement boundary conditions
 
  iSym = 0
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxuvw(ijp, ijb, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), flmass(iface), &
                         cp, cb, sup, svp, swp)

        spu(ijp) = spu(ijp) - cb
        spv(ijp) = spv(ijp) - cb
        sp(ijp)  = sp(ijp)  - cb

        su(ijp) = su(ijp) - cb*u(ijb) + sup
        sv(ijp) = sv(ijp) - cb*v(ijb) + svp
        sw(ijp) = sw(ijp) - cb*w(ijb) + swp

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxuvw(ijp, ijb, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), flmass(iface), &
          cp, cb, sup, svp, swp)  

        spu(ijp) = spu(ijp) - cb
        spv(ijp) = spv(ijp) - cb
        sp(ijp)  = sp(ijp)  - cb

        su(ijp) = su(ijp) - cb*u(ijb) + sup
        sv(ijp) = sv(ijp) - cb*v(ijb) + svp
        sw(ijp) = sw(ijp) - cb*w(ijb) + swp

      enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iSym = iSym + 1

        ! Diffusion coef.
        vsi = vis(ijb)
        cf = vsi*srds(iSym) ! cf v.s. vsol -> cf is calculated using normal distance in srds!

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Dist from cc of owner cell to cf @boundary, cannot expect dpb to be normal to boundary face in general
        dpb = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

        ! Diffusion coef. 
        vsol = vsi*are/dpb

        ! Velocity difference vector components
        upb = u(ijp)-u(ijb)
        vpb = v(ijp)-v(ijb)
        wpb = w(ijp)-w(ijb)

        fdne = 2*cf*( upb*nxf + vpb*nyf + wpb*nzf )

        spu(ijp) = spu(ijp) + vsol
        spv(ijp) = spv(ijp) + vsol
        sp(ijp)  = sp(ijp)  + vsol

        su(ijp) = su(ijp)+vsol*u(ijp)-fdne*nxf
        sv(ijp) = sv(ijp)+vsol*v(ijp)-fdne*nyf
        sw(ijp) = sw(ijp)+vsol*w(ijp)-fdne*nzf

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

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

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

      ! *** Another way to implement wall bc

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



  !----------------------------------------------------------------------------
  ! PRESSURE DIVERGENCE CONTRIBUTION TO SOURCE (using Gauss rule and linear interpolation):
  !----------------------------------------------------------------------------
  call calcPressDiv
  !----------------------------------------------------------------------------




  ! Modify coefficients for Crank-Nicolson
  if (cn) then
    a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.
  endif


!
!.....Assemble and solve system for U component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if (cn) then

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

  urfrs=urfr(iu)
  urfms=urfm(iu)

  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
    ! we substract it from the sum, to eliminate it from the sum.
    ! We could also write sum( a(ioffset(inp)) : a(ioffset(inp+1)-1) ) because all diagonal terms are zero.
    sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp)) 
    a(diag(inp)) = spu(inp) - sum_off_diagonal_terms

    ! a(diag(inp)) = spu(inp) 
    ! do k = ioffset(inp),ioffset(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*u(inp)

    apu(inp) = 1./(a(diag(inp))+small)

  enddo

  ! Solve fvm equations
  call bicgstab(u,iu)

!
!.....Assemble and solve system for V component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if(cn) then

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

  urfrs=urfr(iv)
  urfms=urfm(iv)
  
  do inp=1,numCells
    a(diag(inp)) = 0.0_dp
    su(inp) = 0.0_dp
  enddo

  do inp = 1,numCells

    ! Main diagonal term assembly:
    sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp))
    a(diag(inp)) = spv(inp) - sum_off_diagonal_terms

    ! a(diag(inp)) = spv(inp) 
    ! do k = ioffset(inp),ioffset(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = sv(inp) + urfms*a(diag(inp))*v(inp)

    apv(inp) = 1./(a(diag(inp))+small)

  enddo

  ! Solve fvm equations
  call bicgstab(v,iv)
 
!
!.....Assemble and solve system for W component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if(cn) then

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

  urfrs=urfr(iw)
  urfms=urfm(iw)

  do inp=1,numCells
    a(diag(inp)) = 0.0_dp
    su(inp) = 0.0_dp
  enddo
    
  do inp = 1,numCells

    ! Main diagonal term assembly:
    sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp)) 
    a(diag(inp)) = sp(inp) - sum_off_diagonal_terms

    ! a(diag(inp)) = sp(inp) 
    ! do k = ioffset(inp),ioffset(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = sw(inp) + urfms*a(diag(inp))*w(inp)

    apw(inp) = 1./(a(diag(inp))+small)


  enddo

  ! Solve fvm equations
  call bicgstab(w,iw)

end subroutine