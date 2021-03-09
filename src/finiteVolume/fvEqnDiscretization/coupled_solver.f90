module coupled_solver
!
! Purpose:
!  Coupled pressure based solver for all speeds.
!
! Description:
!  Assemples and solves block-coupled system for pressure and velocity.
!  
! Author:
!  Nikola Mirkov 07/2020.
!
use types
use parameters
use geometry
use sparse_matrix
use variables
use gradients
use LIS_linear_solver_library

implicit none

  real(prec), parameter :: RVOZD = 287.058  ! J/(kg K)
  real(prec), parameter :: _33 = 1./3.0_dp

! Ide u sparse_matrix:
real(dp), dimension(:), allocatable :: auu,avv,aww,aup,avp,awp,apu,apv,apw,app



public

contains


subroutine calcuvwp

  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use hcoef
  use gradients, only: grad
  use faceflux_velocity
  use mhd ! maybe later we do calcuvwpEp, we could add Electric potential to large system?

  implicit none


! Local variables
  integer :: i, k, inp, ijp, ijn, iface, if
  ! logical :: ScndOrderWallBC_Model
  integer :: ib
  integer :: iWall, iSym
  integer :: ijb            ! Boundary field value indexes
  real(dp) :: urfrs, urfms, apotime, heat
  real(dp) :: sup, svp, swp
  real(dp) :: sum_off_diagonal_terms
  real(dp) :: const_11_6
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
  real(prec) :: C_rho
  ! real(dp) :: FdUi,FdVi,FdWi! Diffusive flux auxiliary
  ! real(dp) :: Utp, Vtp, Wtp
  ! real(dp) :: Vnp

  real(dp) :: fxn, fxp
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn,dene
  real(dp) :: nxx,nyy,nzz
  real(dp) :: smdpn,sfdpnr
  real(dp) :: xpp,ypp,zpp,xep,yep,zep
  real(dp) :: dpe
  real(dp) :: dpex,dpey,dpez
  real(dp) :: dpxi,dpyi,dpzi
  real(dp) :: Dpu,Dpv,Dpw,dpecorr
  
  ! ScndOrderWallBC_Model = .true.

  ! Initialize sources - note: sp is rhs vector here in coupled solver.
  su = 0.0_dp
  sv = 0.0_dp
  sw = 0.0_dp
  sp = 0.0_dp


  ! Update velocity at boudary and update velocity gradients: 
  call updateVelocityAtBoundary

  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

!@# Maybe Lorentz force in next iteration? Extended coupled MHD solver.
  ! Update values of Lorentz force vector components in MHD case.
  if(lcal(iep)) call calculate_Lorentz_force

  
  !
  ! > Volume sources 
  !

  do inp=1,numCells


    ! Unsteady term
  
    if(ltransient) then

      if( bdf .or. cn ) then
      !
      ! Backward differentiation formula of 1st order.
      !
        apotime = deno(inp)*vol(inp)/timestep

        ! RHS vector contribution
        su(inp) = su(inp) + apotime*uo(inp)
        sv(inp) = sv(inp) + apotime*vo(inp)
        sw(inp) = sw(inp) + apotime*wo(inp)

        ! NAPOMENA - NIJE DEN CONST!, PROMENI GORE ISTO, NEMOZE BITI APOTIME ISTO ZA SVE TIMESTEPOVE
        sp(inp) = deno(inp)*vol(inp)/timestep

!@# Ispravi da idu na dijagonalu matrice!

        ! Matrix diagonal element contribution
        ! spu(inp) = spu(inp) + apotime
        ! spv(inp) = spv(inp) + apotime
        ! sp(inp)  = sp(inp)  + apotime

        apotime = den(inp)*vol(inp)/timestep

        k = diag( inp )

        auu(k) = auu(k) + apotime
        avv(k) = avv(k) + apotime
        aww(k) = aww(k) + apotime

        C_rho = 1./(RVOZD*T(inp))
        app(k) = C_rho*vol(inp)/timestep

      elseif( bdf2 ) then
      !
      ! Three Level Implicit Time Integration (BDF2) - 2nd order.
      !
        apot = deno(inp)*vol(inp)/timestep
        apoot = denoo(inp)*vol(inp)/timestep

        ! RHS vector contribution
        su(inp) = su(inp) + apot*2*uo(inp) - apoot*0.5_dp*uoo(inp)
        sv(inp) = sv(inp) + apot*2*vo(inp) - apoot*0.5_dp*voo(inp) 
        sw(inp) = sw(inp) + apot*2*wo(inp) - apoot*0.5_dp*woo(inp) 

        ! NAPOMENA - NIJE DEN CONST!, PROMENI GORE ISTO, NEMOZE BITI APOTIME ISTO ZA SVE TIMESTEPOVE
        sp(inp) = ( 2*deno(inp) - 0.5_dp*denoo(inp) )*vol(inp)/timestep

        ! Matrix diagonal element contribution
        ! spu(inp) = spu(inp) + 1.5_dp*apotime
        ! spv(inp) = spv(inp) + 1.5_dp*apotime
        ! sp(inp)  = sp(inp)  + 1.5_dp*apotime

        apotime = den(inp)*vol(inp)/timestep

        k = diag( inp )

        auu(k) = auu(k) + 1.5_dp*apotime
        avv(k) = avv(k) + 1.5_dp*apotime
        aww(k) = aww(k) + 1.5_dp*apotime

        C_rho = 1./(RVOZD*T(inp))
        app(inp) = C_rho*1.5_dp*vol(inp)/timestep

      elseif( bdf3 ) then
      
        apot = deno(inp)*vol(inp)/timestep
        apoot = denoo(inp)*vol(inp)/timestep
        apooot = denooo(inp)*vol(inp)/timestep

        ! RHS vector contribution
        su(inp) = su(inp) + apot*3*uo(inp) - apoot*1.5_dp*uoo(inp) + apooot*_33*uooo(inp) 
        sv(inp) = sv(inp) + apot*3*vo(inp) - apoot*1.5_dp*voo(inp) + apooot*_33*vooo(inp) 
        sw(inp) = sw(inp) + apot*3*wo(inp) - apoot*1.5_dp*woo(inp) + apooot*_33*wooo(inp)
        sp(inp) = sp(inp)  ( 3*deno(inp) - 1.5_dp*denoo(inp) + _33*denooo(inp) )*vol(inp)/timestep

        ! Matrix diagonal element contribution
        ! spu(inp) = spu(inp) + 11./6.0_dp*apotime
        ! spv(inp) = spv(inp) + 11./6.0_dp*apotime
        ! sp(inp)  = sp(inp)  + 11./6.0_dp*apotime


        k = diag( inp )
 
        const_11_6 = 11./6.0_dp
        apotime = den(inp)*vol(inp)/timestep
        auu(k) = auu(k) + const_11_6*apotime
        avv(k) = avv(k) + const_11_6*apotime
        aww(k) = aww(k) + const_11_6*apotime

        C_rho = 1./(RVOZD*T(inp))
        app(inp) = C_rho*const_11_6*vol(inp)/timestep

      endif

    endif ! unsteady term

    
    ! Buoyancy source terms for momentum eqns.
    if(lcal(ien).and.lbuoy) then

      if(boussinesq) then ! Boussinesq-ova aproximacija  
        heat = beta*densit*(t(inp)-tref)*vol(inp)
      else
        heat = (densit-den(inp))*vol(inp)
      endif

      su(inp) = su(inp) - gravx*heat
      sv(inp) = sv(inp) - gravy*heat
      sw(inp) = sw(inp) - gravz*heat
      
    endif

    
    ! MHD: Lorentz force source terms for momentum eqns.
    if(lcal(iep)) then

      ! Add Lorentz force volume source terms.
      su(inp) = su(inp) + sigma*florx(inp)*vol(inp)
      sv(inp) = sv(inp) + sigma*flory(inp)*vol(inp)
      sw(inp) = sw(inp) + sigma*florz(inp)*vol(inp)  

    endif

  end do


!  __             __
! |                 |
! | auu auv auw aup |
! |                 |
! | avu avv avw avp |
! |                 |
! | awu awv aww awp |
! |                 |
! | apu apv apw app |
!  -               -
!
! auu, etc. are numCells x numCells sparse matrix each!
!  __             __
! |                 |
! | auu  o   o  aup |
! |                 |
! |  o  avv  o  avp |
! |                 |
! |  o   o  aww awp |
! |                 |
! | apu apv apw app |
!  -               -
!

  !
  ! > Fluxes trough faces:
  !

  ! Inner faces
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    fxn=facint(i)
    fxp=1.0_dp-facint(i)

    call facefluxuvw( ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), flmass(i), facint(i), gds(iu), &
                      cap, can, sup, svp, swp )


    ! Where is (icell,jcell) matrix element in long 'aval' CSR format array, of size [1,nnz]?
    ! Answer: icell_jcell_csr_index(iface)
 
    ! > (Icell,Jcell) matrix element:
    k = icell_jcell_csr_index(i)

    auu(k) = auu(k) + can
    avv(k) = avv(k) + can
    aww(k) = aww(k) + can
 
    aup(k) = aup(k) + fxn*arx(i)
    avp(k) = avp(k) + fxn*ary(i) 
    awp(k) = awp(k) + fxn*arz(i)

    ! > (Icell,Icell) matrix element: 
    k = diag(ijp)

    auu(k) = auu(k) - can
    avv(k) = avv(k) - can
    aww(k) = aww(k) - can

    aup(k) = aup(k) + fxp*arx(i)
    avp(k) = avp(k) + fxp*ary(i)
    awp(k) = awp(k) + fxp*arz(i)

    ! > (Jcell,Icell) matrix element:
    k = jcell_icell_csr_index(i)

    ! a(k) = cap
    auu(k) = auu(k) + cap
    avv(k) = avv(k) + cap
    aww(k) = aww(k) + cap

    aup(k) = fxn*arx(i)
    avp(k) = fxn*ary(i) 
    awp(k) = fxn*arz(i)

    ! > (Jcell,Jcell) matrix element: 
    k = diag(ijn)

    auu(k) = auu(k) - cap
    avv(k) = avv(k) - cap
    aww(k) = aww(k) - cap

    aup(k) = aup(k) + fxp*arx(i)
    avp(k) = avp(k) + fxp*ary(i)
    awp(k) = awp(k) + fxp*arz(i)

    !
    ! > Sources (deffered correction, diffrerence to complete difffusion term,etc.): 
    !
 
    su(ijp) = su(ijp) + sup
    sv(ijp) = sv(ijp) + svp
    sw(ijp) = sw(ijp) + swp

    su(ijn) = su(ijn) - sup
    sv(ijn) = sv(ijn) - svp
    sw(ijn) = sw(ijn) - swp

!@# Kad obilazimo unutrasnje fejsove mpzemo odmah da damo doprinos za pritisak


    ! call facefluxmass(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), capd, cand, flmass(i))


  ! > Geometry:

  ! Face interpolation factor
  fxn = lambda 
  fxp = 1.0_dp-lambda
  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)
  ! Distance between cell centers
  ! dpn = sqrt(xpn**2+ypn**2+zpn**2)
  ! cell face area
  are = sqrt(arx**2+ary**2+arz**2)
  ! Unit vectors of the normal
  nxx = arx/are
  nyy = ary/are
  nzz = arz/are
  !  ______
  ! (Vol/Ap)f
  Dpu = (fxn*Vol(ijn)*Apu(ijn)+fxp*Vol(ijp)*Apu(ijp))
  Dpv = (fxn*Vol(ijn)*Apv(ijn)+fxp*Vol(ijp)*Apv(ijp))
  Dpw = (fxn*Vol(ijn)*Apw(ijn)+fxp*Vol(ijp)*Apw(ijp))
  ! Density at the cell face
  dene = den(ijp)*fxp+den(ijn)*fxn
  !
  ! COEFFICIENTS OF PRESSURE-CORRECTION EQUATION
  sfdpnr = 1.0d0/(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz+small)
  ! sfdpnr = 1.0d0/(arx*xpn+ary*ypn+arz*zpn)  
  !
  ! (Sf.Sf) / (dpn.Sf)
  smdpn = (arx*arx+ary*ary+arz*arz)*sfdpnr
  ! smdpn = (arx*arx+ary*ary+arz*arz)/(arx*xpn+ary*ypn+arz*zpn)
  !
  capd = -dene*Dpu*smdpn
  cand = cap


!////////////////////////////////////////////////////////////
!   RHIE-CHOW velocity interolation at face    
!         __     ______________      ______
!   Uf = (U)f + (dPdxi*(vol/ap))f - (vol/ap)f*(Pn'-Pp')/dp'n'
!   
!   Fluxmass = Densit*dot(Uf,Sf)
!////////////////////////////////////////////////////////////

  !  ______     _____    
  ! (vol/ap)f *(dPdxi)f 
  dpxi = Dpu*(fxn*dPdxi(1,ijn)+fxp*dPdxi(1,ijp))
  dpyi = Dpv*(fxn*dPdxi(2,ijn)+fxp*dPdxi(2,ijp))
  dpzi = Dpw*(fxn*dPdxi(3,ijn)+fxp*dPdxi(3,ijp))

  !.....Values at points p' and e' due to non-orthogonality. 
  xpp = xf-(xf-xc(ijp))*nxx
  ypp = yf-(yf-yc(ijp))*nyy
  zpp = zf-(zf-zc(ijp))*nzz
  xep = xf-(xf-xc(ijn))*nxx
  yep = yf-(yf-yc(ijn))*nyy
  zep = zf-(zf-zc(ijn))*nzz
  !.....Distances |P'P| and |E'E| projected ionto x,y,z-axis
  xpp = xpp-xc(ijp)
  ypp = ypp-yc(ijp)
  zpp = zpp-zc(ijp)
  xep = xep-xc(ijn)
  yep = yep-yc(ijn)
  zep = zep-zc(ijn)

  dpe = dPdxi(1,ijn)*xep+dPdxi(2,ijn)*yep+dPdxi(3,ijn)*zep - &
        dPdxi(1,ijp)*xpp+dPdxi(2,ijp)*ypp+dPdxi(3,ijp)*zpp 

  ! Pressure gradient along normal between N' and P' point which are on face normal direction.
  dpn = dpe/(xpn*nxx+ypn*nyy+zpn*nzz)

  ! Rhie-Chow Interpolation 
  ue = -Dpu * (dpn - dpxi)
  ve = -Dpv * (dpn - dpyi)
  we = -Dpw * (dpn - dpzi)

  ! rhs term
  bpc = dene*(ue*arx+ve*ary+we*arz)


    ! call facefluxsc( ijp, ijn, &
    !                  xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
    !                  flmass(i), facint(i), gam, &
    !                  fi, dFidxi, prtr, capc, canc, suadd )



    ! Convection fluxes - uds
    Crhof = 1./(RVOZD*T(inp)*den(inp))*fxp + 1./(RVOZD*T(ine)*den(ine))*fxe

    canc =  min( flmass(i)*Crhof, zero )
    capc = -max( flmass(i)*Crhof, zero )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    app(k) = cand + canc

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    app(k) = capd + capc

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    app(k) = app(k) - cand - canc

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    app(k) = app(k) - capd - capc

    ! > Sources:

    sp(ijp) = sp(ijp) - bpc - flmass(i) + suadd
    sp(ijn) = sp(ijn) + bpc + flmass(i) - suadd 

  end do ! end - flux terms

      
  ! Implement boundary conditions
 
  iSym = 0
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)
        ijb = iBndValueStart(ib) + i

        call facefluxuvw(ijp, ijb, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), flmass(if), cp, cb, sup, svp, swp)

        spu(ijp) = spu(ijp) - cb
        spv(ijp) = spv(ijp) - cb
        sp(ijp)  = sp(ijp)  - cb

        su(ijp) = su(ijp) - cb*u(ijb) + sup
        sv(ijp) = sv(ijp) - cb*v(ijb) + svp
        sw(ijp) = sw(ijp) - cb*w(ijb) + swp

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)
        ijb = iBndValueStart(ib) + i

        call facefluxuvw(ijp, ijb, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), flmass(if), cp, cb, sup, svp, swp)  

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


  ! Modify coefficients for Crank-Nicolson
  if (cn) then
    a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.
  endif


!
! Assembling and Under-relaxation 
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

    apu(inp) = 1.0/(a(diag(inp))+small)

  enddo

  ! ! Solve fvm equations
  ! call bicgstab(u,iu)

!
! Assemble and underrelax system for V component of velocity
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
  
!@# Ovo je trebalod a se resetuje dijagonala matrice A i rhs vektor Su
  ! do inp=1,numCells
  !   a(diag(inp)) = 0.0_dp
  !   su(inp) = 0.0_dp
  ! enddo

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

    apv(inp) = 1.0/(a(diag(inp))+small)

  enddo

  ! ! Solve fvm equations
  ! call bicgstab(v,iv)
 
!
! Assemble and underrelax system for W component of velocity
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

!@# Ovo je trebalod a se resetuje dijagonala matrice A i rhs vektor Su
  ! do inp=1,numCells
  !   a(diag(inp)) = 0.0_dp
  !   su(inp) = 0.0_dp
  ! enddo
    
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

  ! ! Solve fvm equations
  ! call bicgstab(w,iw)

end subroutine

end module