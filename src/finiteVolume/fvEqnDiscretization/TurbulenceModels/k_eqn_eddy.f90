module k_eqn_eddy
!
! Implementation of one-equation eddy viscosity SGS, 
! Yoshizawa et al. - Bridging between eddy-viscosity-type and second-order models using an two-scale DIA,
! 9th International symposium on turbulent shear flow, Kyoto, Japan, Vol. 3., 1993, p 23.1.1-23.1.6.
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use TurbModelData, only: TurbModel
  use scalar_fluxes, only: facefluxsc

  implicit none

  ! Turbulence model constants
  real(dp), parameter :: Ck = 0.094_dp
  real(dp), parameter :: Ceps = 1.048_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(ck))
  real(dp), parameter :: CMU75 = cmu25**3


  private 

  public :: modify_viscosity_k_eqn_eddy
  public :: modify_viscosity_inlet_k_eqn_eddy


contains


!***********************************************************************
!
subroutine modify_viscosity_k_eqn_eddy
!
!***********************************************************************
!
! Main module routine to solve turbulence model equations and update effective viscosity.
!
!***********************************************************************

  implicit none

  call calcsc(TE,dTEdxi,1) ! Assemble and solve turbulence kinetic energy eq.
  call modify_mu_eff

end subroutine



!***********************************************************************
!
subroutine calcsc(fi,dfidxi,ifi)
!
!***********************************************************************
!
! Discretization of scalar equation
!
!***********************************************************************

  use sparse_matrix
  use linear_solvers
  
  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, iftwin
  real(dp) :: prtr, apotime, urfrs, urfms
  real(dp) :: utp, vtp, wtp, utn, vtn, wtn
  real(dp) :: genp, genn
  real(dp) :: uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2
  ! real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real(dp) :: viss,viste,dcoef
  real(dp) :: fimax,fimin
  real(dp) :: delta,third,three_halfs
  real(dp) :: gam, urf, tolAbs, tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 

  third = 1./3.0_dp
  three_halfs = 3./2.0_dp

  ! Variable specific coefficients:
  gam = TurbModel%Scalar( ifi )%gds
  urf = TurbModel%Scalar( ifi )%urf
  lSolver = TurbModel%Scalar( ifi )%lSolver 
  maxiter = TurbModel%Scalar( ifi )%maxiter
  tolAbs = TurbModel%Scalar( ifi )%tolAbs
  tolRel = TurbModel%Scalar( ifi )%tolRel

  prtr = 1.0_dp

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
!

  !=========================================================
  ! STANDARD PRODUCTION
  !=========================================================
  do inp=1,numCells

    ! dudx = dudxi(1,inp)
    ! dudy = dudxi(2,inp)
    ! dudz = dudxi(3,inp)

    ! dvdx = dvdxi(1,inp)
    ! dvdy = dvdxi(2,inp)
    ! dvdz = dvdxi(3,inp)

    ! dwdx = dwdxi(1,inp)
    ! dwdy = dwdxi(2,inp)
    ! dwdz = dwdxi(3,inp)

    ! ! Minus here in fron because UU,UV,... calculated in calcstress hold -tau_ij
    ! ! So the exact production is calculated as tau_ij*dui/dxj
    ! gen(inp) = -den(inp)*( uu(inp)*dudx+uv(inp)*(dudy+dvdx)+uw(inp)*(dudz+dwdx) &
    !                                    +vv(inp)*dvdy       +vw(inp)*(dvdz+dwdy) &
    !                                                        +ww(inp)*dwdz )

    gen(inp) = 2*abs(vis(inp)-viscos)*magStrain(inp)*magStrain(inp)

  enddo

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================
  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Add production term to the rhs:
    su(inp)=genp*vol(inp)  

    ! Add destruction term to the lhs:
    delta = vol(inp)**third
    sp(inp)=den(inp) * (Ceps*te(inp)**three_halfs/delta) * vol(inp)/(te(inp)+small)
    sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)


    !=====================================
    ! VOLUME SOURCE TERMS: buoyancy
    !=====================================
    if(lbuoy) then
      
      ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
      call calcheatflux 

      if(boussinesq) then
         uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)*beta
         vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)*beta
         wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)*beta
      else
         uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)/(t(inp)+small)
         vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)/(t(inp)+small)
         wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)/(t(inp)+small)
      end if

      utp=max(uttbuoy,zero)
      vtp=max(vttbuoy,zero)
      wtp=max(wttbuoy,zero)
      utn=min(uttbuoy,zero)
      vtn=min(vttbuoy,zero)
      wtn=min(wttbuoy,zero)

      su(inp)=su(inp)+utp+vtp+wtp
      sp(inp)=sp(inp)-(utn+vtn+wtn)/(te(inp)+small)

    end if


    ! UNSTEADY TERM
    if(ltransient) then

      apotime = den(inp)*vol(inp)/timestep

      if( bdf .or. cn ) then

        su(inp) = su(inp) + apotime*teo(inp)
        sp(inp) = sp(inp) + apotime

      elseif( bdf2 ) then

        su(inp) = su(inp) + apotime*( 2*teo(inp) - 0.5_dp*teoo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime

      endif

    endif

  ! End of TKE volume source terms
  enddo



!
! CALCULATE TERMS INTEGRATED OVER FACES
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    ! Diffusion coefficient
    viste = ( vis(ijp) + (vis(ijn)-vis(ijp))*facint(i) )-viscos
    dcoef = viscos+viste*prtr

    call facefluxsc( i, ijp, ijn, &
                     xf(i), yf(i), zf(i), &
                     arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, &
                     TurbModel%Scalar( ifi )%cScheme, &
                     fi, dfidxi, dcoef, cap, can, suadd )

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
         bctype(ib) == 'pressure') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! Diffusion coefficient
        viste = vis(ijb)-viscos
        dcoef = viscos + viste*prtr 

        call facefluxsc( ijp, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         dfidxi, dcoef, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp) - can*fi(ijb) + suadd

      end do


    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic bnd pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin peridoic face

        ! Diffusion coefficient
        viste = half*(vis(ijp)+vis(ijn))-viscos
        dcoef = viscos+viste*prtr

        call facefluxsc( i, ijp, ijn, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), gam, &
                         fi, dfidxi, dcoef, cap, can, suadd )

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



    elseif ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        !
        ! > Wall boundary conditions for turbulence kinetic energy eq.
        !

        su(ijp )= su(ijp)-gen(ijp)*vol(ijp) ! take out standard production from wall ajdecent cell.
        viss=viscos
        if(ypl(i).gt.ctrans) viss=visw(iwall)
        ! viss = max(viscos,visw(iWall))

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
        Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)

        ! Tangent direction - unit vector
        xtp = xtp/vtp
        ytp = ytp/vtp
        ztp = ztp/vtp

        ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        Tau(iwall) = viss*Ut2/dnw(iwall)

        gen(ijp)=abs(tau(iwall))*cmu25*sqrt(te(ijp))/(dnw(iwall)*cappa)
        su(ijp)=su(ijp)+gen(ijp)*vol(ijp)


      enddo


    endif

  enddo


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
        off_diagonal_terms = sum( a( ia(ijp) : ia(ijp+1)-1 ) ) - a(diag(ijp))
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
    a(diag(inp)) = sp(inp) - sum( a(ia(inp) : ia(inp+1)-1) ) - a(diag(inp))
    ! a(diag(inp)) = sp(inp) 
    ! do k = ia(inp),ia(inp+1)-1
    !   if (k.eq.diag(inp)) cycle
    !   a(diag(inp)) = a(diag(inp)) -  a(k)
    ! enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)
                    
  enddo

  ! Solve linear system:
  call csrsolve(lSolver, te, su, resor(5), maxiter, tolAbs, tolRel, 'k' )


  ! Update field values at boundaries
  call updateBoundary( fi )


! Report range of scalar values and clip if negative
  fimin = minval(fi)
  fimax = maxval(fi)

  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= k <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


subroutine modify_mu_eff
!
! Update turbulent and effective viscosity.
!

  implicit none

  integer :: i,inp
  integer :: iface, ijp,ijb,ib,iwall
  real(dp) :: urf,visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau,viscw
  real(dp) :: third,delta

  third = 1./3.0_dp

  urf = TurbModel%urfVis

  !
  ! Loop trough cells 
  !

  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    delta = vol(inp)**third

    ! Update effective viscosity:
    ! \mu_{eff}=\mu+\mu_t; \mu_t = C_k * k^{0.5} * Delta
    vis(inp)=viscos + den(inp)*ck*sqrt(te(inp))*delta

    ! Underelaxation
    vis(inp)=urf*vis(inp)+(1.0_dp-urf)*visold

  enddo

  !
  ! Boundary faces 
  !

  ! Update value at every bounary face type except "wall', which is treated below.
  call updateBoundary( vis )

  ! Now 'wall' type.
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
        Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)

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
        !ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
        ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

        viscw = zero

        if(ypl(iWall) > ctrans) then
          viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
        endif

        visw(iWall) = max(viscos,viscw)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo


end subroutine


subroutine modify_viscosity_inlet_k_eqn_eddy
!
! Update turbulent and effective viscosity at inlet.
!

  implicit none

  integer :: i,ib,ijb,iface
  real(dp) :: delta,are

  ! Loop over inlet boundaries
 
  !
  ! Boundary faces 
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijb = iBndValueStart(ib) + i

        are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )

        Delta = sqrt(are)

        ! Update effective viscosity:
        ! \mu_{eff}=\mu+\mu_t; \mu_t = C_k * k^{0.5} * Delta
        Vis(ijb) = viscos + den(ijb)*ck*sqrt(te(ijb))*delta

      end do

    endif

  enddo 


end subroutine


end module k_eqn_eddy