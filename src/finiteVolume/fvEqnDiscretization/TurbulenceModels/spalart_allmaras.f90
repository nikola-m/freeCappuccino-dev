module Spalart_Allmaras
!
! Implementation of Spalart-Allmaras turbulence model.
! Reference:  
! * Spalart, P. R. and Allmaras, S. R., "A One-Equation Turbulence Model for Aerodynamic Flows",
!   Recherche Aerospatiale, No. 1, 1994, pp. 5-21.
! * https://turbmodels.larc.nasa.gov/spalart.html (for implementation)
! ***
! * NOTE: We use TKE array to store nu_tilda!
! ***
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use TurbModelData, only: TurbModel

  implicit none


  ! Turbulence model constants
  real(dp), parameter :: Cb1 = 0.1355_dp
  real(dp), parameter :: Cb2 = 0.622_dp
  real(dp), parameter :: sigma = twothirds
  real(dp), parameter :: Cw1 = 3.2390678167757287_dp ! = Cb1/cappa**2 + (1.0_dp+Cb2)/sigma
  real(dp), parameter :: Cw2 = 0.3_dp
  real(dp), parameter :: Cw3 = 2.0_dp
  real(dp), parameter :: Cv1 = 7.1_dp
  real(dp), parameter :: Ct3 = 1.2_dp
  real(dp), parameter :: Ct4 = 0.5_dp

  real(dp), parameter :: one_sixth = 1._dp/6._dp
  logical :: notf2 = .false. ! Set to .true. for SA model without f_t2 term. Also for DES and DDES.

  private 

  public :: modify_viscosity_spalart_allmaras
  public :: modify_viscosity_inlet_spalart_allmaras


contains


subroutine modify_viscosity_spalart_allmaras()
!
! Main module routine to solve turbulence model equations and update effective viscosity.
! ***
! * NOTE: We use TKE array to store nu_tilda here!
! ***
!

  implicit none

  call calcsc(TE,dTEdxi,1) ! Assemble and solve nu_tilda eq.
  call modify_mu_eff()

end subroutine


subroutine calcsc(Fi,dFidxi,ifi)
!
! Assemble and solve transport equation for scalar field.
!

  use sparse_matrix
  use linear_solvers
  
  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, l, inp, ijp, ijn, ijb, iper, ib, iwall, iface, if, iftwin
  real(dp) :: prtr, apotime, urfrs, urfms, genp, genn
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: fimax,fimin
  real(dp) :: nu_tilda,nu,xi,fv1,fv2,ft2,r,g,fw,strain_tilda
  real(dp) :: dnutdx,dnutdy,dnutdz
  real(dp) :: urf, tolAbs, tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 

! Variable specific coefficients:
  urf = TurbModel%Scalar(ifi)%urf
  lSolver = TurbModel%Scalar(ifi)%lSolver 
  maxiter = TurbModel%Scalar(ifi)%maxiter
  tolAbs = TurbModel%Scalar(ifi)%tolAbs
  tolRel = TurbModel%Scalar(ifi)%tolRel

  prtr=1./sigma

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

  !
  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  !

  do inp=1,numCells

    nu_tilda = te(inp)   ! we store nu_tilda in tke field.
    nu = viscos/den(inp)   ! kinematic viscosity
    xi = nu_tilda / nu   ! xi parameter
    fv1 = xi**3 / (xi**3 + Cv1**3) ! parameter
    fv2 = 1.0_dp - xi / (1.0_dp + xi*fv1) ! parameter
    strain_tilda = vorticity(inp) + (nu_tilda / (cappa*wallDistance(inp))**2) * fv2 ! Shat

    ! Limiting strain - this is very important, there are other ways besides clipping implemented here.
    ! Clipping implemented here originates from a note referencing private communication with P. Spalart on turbmodels.larc.nasa.gov website
    strain_tilda = max(strain_tilda, 0.3_dp*vorticity(inp)) 

    ft2 = Ct3 * exp(-Ct4 * xi**2)
    if (notf2) ft2 = 0.0_dp ! If we want to eliminate ft2 term, other method is setting ct3 = 0.

    gen(inp) = Cb1 * (1.0_dp - ft2) * strain_tilda * nu_tilda

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! > Add production term to the rhs:
    su(inp)=genp*vol(inp)  


    ! > Cross diffusion: 
    dnutdx=dTEdxi(1,inp)
    dnutdy=dTEdxi(2,inp)
    dnutdz=dTEdxi(3,inp)

    su(inp) = su(inp) + ( Cb2/sigma * (dnutdx**2+dnutdy**2+dnutdz**2) ) * vol(inp)


    ! > Add destruction term to the lhs:

    r = min(nu_tilda/(strain_tilda*(cappa*wallDistance(inp))**2 + small), 10.0_dp) ! r parameter

    g = r + Cw2*(r**6-r) ! g parameter

    fw = g*((1.0_dp+Cw3**6)/(g**6+Cw3**6))**one_sixth ! fw parameter

    ! > destruction term
    sp(inp) = (Cw1*fw-(Cb1/cappa**2)*ft2)*(nu_tilda/wallDistance(inp))**2 * den(inp)*vol(inp)/(nu_tilda+small)

    ! > If production negative move it to the lhs:
    sp(inp) = sp(inp)-genn*vol(inp)/(nu_tilda+small)

    ! > UNSTEADY TERM
    if (ltransient) then
      if ( bdf .or. cn ) then
        apotime = den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*teo(inp)
        sp(inp) = sp(inp) + apotime
      elseif( bdf2 ) then
        apotime=den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*( 2*teo(inp) - 0.5_dp*teoo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime
      endif
    endif

  ! End of volume source terms
  enddo


!
! CALCULATE TERMS INTEGRATED OVER FACES
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), &
                     arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), TurbModel%Scalar(ifi)%gds, &
                     fi, dFidxi, prtr, cap, can, suadd )

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

        call facefluxsc_boundary( ijp, ijb, &
                                  xf(iface), yf(iface), zf(iface), &
                                  arx(iface), ary(iface), arz(iface), &
                                  flmass(iface), &
                                  Fi, dFidxi, prtr, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

      end do


    ! elseif ( bctype(ib) == 'wall') then

    !   do i=1,nfaces(ib)

    !     iface = startFace(ib) + i
    !     ijp = owner(iface)
    !     ijb = iBndValueStart(ib) + i
    !     iWall = iWall + 1

    !     call facefluxsc_boundary( ijp, ijb, &
    !                               xf(iface), yf(iface), zf(iface), &
    !                               arx(iface), ary(iface), arz(iface), &
    !                               flmass(iface), &
    !                               FI, dFidxi, prtr, cap, can, suadd )

    !     Sp(ijp) = Sp(ijp)-can

    !     Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

    !   enddo


    elseif (  bctype(ib) == 'periodic' ) then

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
                                  flmass(if), &
                                  TurbModel%Scalar(ifi)%gds, &
                                  fi, dFidxi, prtr, cap, can, suadd )


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
        su(ijp) = su(ijp) - a(k)*teo(ijn)

        k = jcell_icell_csr_index(i)
        su(ijn) = su(ijn) - a(k)*teo(ijp)
    enddo

    do ijp=1,numCells
        apotime=den(ijp)*vol(ijp)/timestep
        off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) - a(diag(ijp))
        su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*teo(ijp)
        sp(ijp) = sp(ijp)+apotime
    enddo
        
  endif

  ! Underrelaxation factors
  urfrs=1.0_dp/urf
  urfms=1.0_dp-urf

  ! Main diagonal term assembly:
  do inp = 1,numCells

    ! Main diagonal term assembly:
    a(diag(inp)) = sp(inp) 
    do k = ioffset(inp),ioffset(inp+1)-1
      if (k.eq.diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) -  a(k)
    enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)
                    
  enddo

  ! Solve linear system:
  call csrsolve(lSolver, fi, su, resor(5), maxiter, tolAbs, tolRel, 'nutilda' )


  ! Update field values at boundaries
  call updateBoundary( fi )

  ! iPer = 0

  ! do ib=1,numBoundaries

  !   if ( bctype(ib) == 'outlet' .or. &
  !       bctype(ib) == 'symmetry' .or. &
  !        bctype(ib) == 'pressure' .or. &
  !       bctype(ib) == 'empty' ) then

  !     do i=1,nfaces(ib)

  !       iface = startFace(ib) + i
  !       ijp = owner(iface)
  !       ijb = iBndValueStart(ib) + i

  !       fi(ijb)=fi(ijp)

  !     enddo

  !   elseif (  bctype(ib) == 'periodic' ) then

  !     iPer = iPer + 1

  !     ! Faces trough periodic boundaries, Taiwo first
  !     do i=1,nfaces(ib)

  !       iface = startFace(ib) + i
  !       ijp = owner(iface)
  !       ijb = iBndValueStart(ib) + i

  !       iface = startFaceTwin(iPer) + i
  !       ijn = owner(iface)

  !       fi(ijb) = half*( fi(ijp)+fi(ijn) )

  !       ! Now find where is the twin in field array
  !       ijbt = numCells + ( startFaceTwin(iPer) - numInnerFaces ) + i
        
  !       ! Twin takes the same values
  !       fi(ijbt) = fi(ijb)

  !     enddo

  !   endif

  ! enddo


! Report range of scalar values and clip if negative
  fimin = minval(fi)
  fimax = maxval(fi)
  
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= nutilda <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


subroutine modify_mu_eff()
!
! Update turbulent and effective viscosity.
!

  implicit none

  integer :: i,ib,inp
  integer :: iface, iwall,ijp,ijb
  real(dp) :: urf, visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau,viscw
  ! real(dp) :: Utauvis,Utaulog,Upl
  real(dp) :: Uplblend
  real(dp) :: nu_tilda,nu,xi,fv1

!
! Loop trough cells 
!
  urf = TurbModel%urfVis

  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    nu_tilda = te(inp)
    nu = viscos/den(inp)
    xi = nu_tilda / nu 
    fv1 = xi**3 / (xi**3 + Cv1**3)
    ! Update effective viscosity:
    vis(inp)=viscos + den(inp)*nu_tilda*fv1

    ! Underelaxation
    vis(inp)=urf*vis(inp)+(1.0_dp-urf)*visold

  enddo

  !
  ! Boundary faces 
  !

  ! Update value at every bounary face type except 'wall', which is treated below.
  call updateBoundary( vis )

  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        !
        ! Wall boundaries - update Visw and Ypl
        !

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Magnitude of a cell center velocity projected on boundary face normal
        Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

        ! Tangential velocity components 
        xtp = U(ijp)-Vnp*nxf
        ytp = V(ijp)-Vnp*nyf
        ztp = W(ijp)-Vnp*nzf

        ! Its magnitude
        Vtp = xtp*xtp+ytp*ytp+ztp*ztp

        ! Tangent direction
        xtp = xtp/vtp
        ytp = ytp/vtp
        ztp = ztp/vtp

        ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        Tau(iWall) = viscos*Ut2/dnw(iWall)
        Utau = sqrt( Tau(iWall) / den(ijb) )
        ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
        ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+
        ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos

        ! viscw = zero
        ! if(ypl(iWall) > ctrans) then
        !   viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
        ! endif

        ! Potpuno neka druga logika sad

        ! *** Automatic wall treatment ***
 
        ! Ne mogu ove dve linije dole ako nemam TKE!!
        ! utau = sqrt( viscos*Vtp/(densit*dnw(iWall)) + cmu25*te(ijp) ) ! It's actually u* in original reference...

        ! ypl(iWall) = den(ijp)*Utau*dnw(iWall)/viscos 

        ! Utauvis=ypl(iWall)
        ! Utaulog=1.0/cappa*log(Elog*ypl(iWall)) 

        ! Upl=sqrt(sqrt(Utauvis**4+Utaulog**4)) 
          
        ! viscw = den(ijp)*utau*dnw(iWall)/Upl 

        ! Blended version of shear stress - probati ovo(!?)
        ! tau(iWall) = den(ijp) * (Vtp/Uplblend)**2

        ! Varijanta 2, u originalnoj referenci...
        ! tau(iWall) = den(ijp) * Vtp*Utau/Upl

        !*** END: Automatic wall treatment ***

        ! *** Enhanced wall treatment - Reichardt blending ***

        ! ! Below is a variant where we use Reichardt blending
        ! ! for whole span of y+ values.
        ! ! Some authors say that Reichardt function for u+ approximates
        ! ! the composite u+(y+) curve, better that Kader blending function.
 
        ! utau = sqrt( viscos*Vtp/(densit*dnw(iWall)) + cmu25*te(ijp) ) ! It's actually u* in original reference...

        ! ypl(iWall) = den(ijp)*Utau*dnw(iWall)/viscos 

        Uplblend = one/cappa*log(one+cappa*ypl(iWall)) + &
                   7.8_dp*(1.-exp(-ypl(iWall)/11.0_dp)-(ypl(iWall)/11.0_dp)*exp(-ypl(iWall)/3.0_dp))
          
        viscw = den(ijp)*utau*dnw(iWall)/Uplblend  

        ! ! Blended version of shear stress - probati ovo(!?)
        ! ! tau(iWall) = den(ijp) * (Vtp/Uplblend)**2

        ! Varijanta 2, u originalnoj referenci...
        tau(iWall) = den(ijp) * Vtp*Utau/Uplblend

        !*** END: Enhanced wall treatment - Reichardt blending ***

        visw(iWall) = max(viscos,viscw)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo


end subroutine


subroutine modify_viscosity_inlet_spalart_allmaras()
!
! Update turbulent and effective viscosity at inlet.
!

  implicit none

  integer :: i,ib,ijb
  real(dp) :: nu_tilda,nu,fv1,xi

  !
  ! Boundary faces - loop over inlet boundaries
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        ijb = iBndValueStart(ib) + i

        nu_tilda = te(ijb) ! nu_tilda stored in TKE array..
        nu = viscos/den(ijb)
        xi = nu_tilda / nu 
        fv1 = xi**3 / (xi**3 + Cv1**3)

        Vis(ijb) = viscos+den(ijb)*nu_tilda*fv1

      end do

    endif

  enddo    


end subroutine



!***********************************************************************
!
subroutine facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                      fm, lambda, gam, FI, dFidxi, &
                      prtr, cap, can, suadd)
!
!***********************************************************************
!
  use interpolation

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn                       ! owner and neighbour indices
  real(dp), intent(in) :: xf,yf,zf                      ! face centroid coordinates
  real(dp), intent(in) :: arx, ary, arz                 ! area vector
  real(dp), intent(in) :: fm                            ! mass flow at face
  real(dp), intent(in) :: lambda                        ! interpolation factor
  real(dp), intent(in) :: gam                           ! deferred correction factor [0,1], 1-high order.
  real(dp), dimension(numTotal), intent(in) :: Fi       ! Scalar field in question
  real(dp), dimension(3,numCells), intent(in) :: dFidxi ! Gradient of the scalar field in question.
  real(dp), intent(in) :: prtr                          ! One over prandtl coefficient
  real(dp), intent(inout) :: cap, can, suadd            ! On return - matrix coeffs for owner and neighbour and source contribution.


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn!, xi,yi,zi,r1,r2,psie,psiw
  real(dp) :: dpn
  real(dp) :: cp,ce
  real(dp) :: fii
  real(dp) :: fdfie,fdfii,fcfie,fcfii,ffic
  real(dp) :: de, game, nu, nu_tilda
  real(dp) :: fxp,fxn
  real(dp) :: dfixi,dfiyi,dfizi
  real(dp) :: dfixii,dfiyii,dfizii
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


  ! Cell face diffussion coefficient
  nu = viscos/den(ijp) 
  nu_tilda = fi(ijp)*fxp+fi(ijn)*fxn
  game = prtr*(nu+nu_tilda) ! prtr is 1/sigma, is set in calcsc above

  ! Difusion coefficient for linear system
  ! de = game*are/dpn
  de = game*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)

  ! Convection fluxes - uds
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  cap = -de - cp
  can = -de + ce
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              Fi, dFidxi, TurbModel%Scalar(1)%nrelax, TurbModel%Scalar(1)%dScheme, dfixi, dfiyi, dfizi, &
              dfixii, dfiyii, dfizii)
  

  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)  

  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)


  !-------------------------------------------------------
  ! Explicit higher order convection; cScheme set in input
  !-------------------------------------------------------
  if( fm .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dFidxi, TurbModel%Scalar(1)%cScheme)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dFidxi, TurbModel%Scalar(1)%cScheme)
  endif

  fcfie = fm*fii

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deferred correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  ffic = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -ffic+fdfie-fdfii 

end subroutine



!***********************************************************************
!
subroutine facefluxsc_periodic(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                               flmass, gam, &
                               FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use interpolation

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass
  real(dp), intent(in) :: gam 
  ! character( len=30 ), intent(in) :: cScheme ! Convection scheme.
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: dpn
  real(dp) :: Cp,Ce
  real(dp) :: fii,fm
  real(dp) :: fdfie,fdfii,fcfie,fcfii,ffic
  real(dp) :: de, game, nu, nu_tilda
  real(dp) :: fxp,fxn
  real(dp) :: dfixi,dfiyi,dfizi
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


  ! Cell face diffussion coefficient

  nu = viscos/den(ijp) 
  nu_tilda = fi(ijp)*fxp+fi(ijn)*fxn
  game = prtr*(nu+nu_tilda) ! prtr is 1/sigma, is set in calcsc above


  ! Difusion coefficient for linear system
  ! de = game*are/dpn
  de = game*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  cap = -de - max(fm,zero)
  can = -de + min(fm,zero)
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

  ! Explicit diffusion
  fdfie = game*(dfixi*arx + dfiyi*ary + dfizi*arz)  

  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)


  !-------------------------------------------------------
  ! Explicit higher order convection, cSheme is set in input
  !-------------------------------------------------------
  if( flmass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value_cds(ijp, ijn, fxp, fi)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value_cds(ijn, ijp, fxn, fi)
  endif

  fcfie = fm*fii

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deffered correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  ffic = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -ffic+fdfie-fdfii 

end subroutine




!***********************************************************************
!
subroutine facefluxsc_boundary(ijp, ijn, xf, yf, zf, arx, ary, arz, fm, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numTotal), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
  real(dp) :: cp,ce
  real(dp) :: fdfie,fdfii
  real(dp) :: d1x,d1y,d1z,d2x,d2y,d2z
  real(dp) :: de, vole, game, nu, nu_tilda
  real(dp) :: fxp,fxn
  real(dp) :: dfixi,dfiyi,dfizi
  real(dp) :: dfixii,dfiyii,dfizii

!----------------------------------------------------------------------

  dfixi = 0.0_dp
  dfiyi = 0.0_dp
  dfizi = 0.0_dp

  ! > Geometry:

  ! Face interpolation factor
  fxn=1.0_dp
  fxp=0.0_dp

  ! Distance vector between cell center and face center
  xpn=xf-xc(ijp)
  ypn=yf-yc(ijp)
  zpn=zf-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! Cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectors n and i_xi - we need cosine
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

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Cell face diffussion coefficient
  nu = viscos/den(ijp)
  nu_tilda = fi(ijp)*fxp+fi(ijn)*fxn
  game = prtr*(nu+nu_tilda) ! prtr is 1/sigma, is set in calcsc above


  !-- Skewness correction --

  ! Overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn

  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn

  ! Interpolate gradients defined at CV centers to faces..
  ! It should be dFidxi(:,ijn), because fxn=1.0, but we write dFidxi(:,ijp) because constant gradient
  ! is applied between cell center and boundary cell face.
  dfixi = dFidxi(1,ijp)
  dfiyi = dFidxi(2,ijp)
  dfizi = dFidxi(3,ijp) 

  !.....du/dx_i interpolated at cell face:
  dfixii = dfixi*d1x + arx/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfiyii = dfiyi*d1y + ary/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfizii = dfizi*d1z + arz/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
 

  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)   

  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)


  ! Difusion coefficient
  de = game*are/dpn

  ! Convection fluxes - uds
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  cap = -de - cp
  can = -de + ce

  ! Explicit part of fluxes
  suadd = fdfie-fdfii 

end subroutine


end module Spalart_Allmaras