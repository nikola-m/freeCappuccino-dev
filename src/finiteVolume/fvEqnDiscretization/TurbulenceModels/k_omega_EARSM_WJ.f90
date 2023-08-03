module k_omega_earsm_wj
!
! Implementation of k-omega Shear Stress Transport (SST) two equation turbulence model.
!
! REFERENCES:
!     * https://turbmodels.larc.nasa.gov/easmko.html
!     * Hellsten, A., "New Advanced k-omega Turbulence Model for High-Lift Aerodynamics," 
!       AIAA Journal, Vol. 43, No. 9, 2005, pp. 1857-1869, https://doi.org/10.2514/1.13754.
!     * Wallin, S. and Johansson, A. V., "An Explicit Algebraic Reynolds Stress Model for 
!       Incompressible and Compressible Turbulent Flows," J. Fluid Mechanics, Vol. 403, 2000, pp. 89-132, 
!       https://doi.org/10.1017/S0022112099007004.
!     * V. Hamalainen, Implementing an Explicit Algebraic Reynolds Stress Model 
!       Into The Three-Dimensional FINFLO Flow Solver, Report B-52, Espoo, 2001, Finland.
!    
!
!======================================================================
!     Adopting:
!     A subroutine to calculate the extra anisotropy components and
!     the effective eddy-viscosity coefficient used with the EARSM.
!     Date:      Author:     Affiliation:
!     23.11.2001 Ville H. /  Laboratory of Aerodynamics, Espoo, Finland 
!     13.12.2011 Nikola M. / Laboratory for Thermal Engineering and Energy,
!                            Institute of Nuclear Sciences "Vinca", Belgrade, Serbia  
!
!     There are two versions here: the original Wallin-Johansson and Menter et. al 2009
!======================================================================
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use TurbModelData, only: TurbModel
  use scalar_fluxes, only: facefluxsc

  implicit none

  logical :: LowRe = .false. ! Has to be set in calling routine or in main program.
  logical :: EARSM_WJ = .true.
  logical :: EARSM_M = .false.


! Often used numbers
  real(dp), parameter :: P3= 1./3.   ! One third
  real(dp), parameter :: P6= 0.5*P3  ! One sixth
  real(dp), parameter :: TT= 2.0*P3  ! Two thirds

! Model coefficient and its variants
  real(dp), parameter :: C1E = 1.8_dp
  real(dp), parameter :: C1P = 1.8_dp  ! 9.0*(C1E - 1.0)/4.0 ! C1_prime
  real(dp), parameter :: C1P3 = 0.6_dp ! C1P*P3
  real(dp), parameter :: C1PSQ = 3.24_dp ! C1P**2
  real(dp), parameter :: CT = 6.0_dp

  real(dp), parameter :: BETTAST=0.09_dp   
  real(dp), parameter :: A1=0.31_dp

  ! Turbulence coeffficients (a bit different for EARSM_M and EARSM_WJ)
  real(dp) :: SIGMK1
  real(dp) :: SIGMK2
  real(dp) :: SIGMOM1
  real(dp) :: SIGMOM2
  real(dp) :: BETAI1
  real(dp) :: BETAI2
  real(dp) :: ALPHA1
  real(dp) :: ALPHA2

!  SST-1994 coefficients:
!  ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM1)
!  ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM2)
!
!  SST-2003 coefficients: The first is higher than the original constant
!  definition by approximately 0.43%, and the second is lower by less than 0.08%. 
!  ALPHA1=5./9.0_dp; ALPHA2=0.44_dp
!

  real(dp), parameter :: C3 = 1.44_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(BETTAST))
  real(dp), parameter :: CMU75 = cmu25**3

  ! real(dp) :: TEIN ! Used on We would neeed to define it from the inlet boundary data for TKE.

  real(dp), dimension(:), allocatable :: fsst,cmueff ! Blending function and effecive Cmu coefficient

  private 

  public :: EARSM_M, EARSM_WJ, LowRe
  public :: modify_viscosity_k_omega_earsm_wj
  public :: modify_viscosity_inlet_k_omega_earsm_wj

contains



subroutine modify_viscosity_k_omega_earsm_wj
!
! Main module routine to solve turbulence model equations and update effective viscosity.
!

  implicit none

  if(.not.allocated(fsst)) then
    allocate( fsst(numTotal) )
    write(*,'(a)') "  **Allocated fsst blending function for k-omega-EARSM."
  endif

  if(.not.allocated(cmueff)) then
    allocate( cmueff(numCells) )
    write(*,'(a)') "  **Allocated cmueff for k-omega-EARSM."
  endif

  ! Turbulence model constants 

  IF (EARSM_M) THEN
      SIGMK1=0.856_dp
      SIGMK2=1.1_dp
      SIGMOM1=0.5_dp
      SIGMOM2=0.856_dp
      BETAI1=0.075
      BETAI2=0.0828
      ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMOM1)
      ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMOM2)
  ELSE ! EARSM_WJ
      SIGMK1=1.1_dp
      SIGMK2=1.1_dp
      SIGMOM1=0.53_dp
      SIGMOM2=1.0
      BETAI1=0.0747
      BETAI2=0.0828
      ALPHA1=0.518_dp
      ALPHA2=0.44
  END IF

  call calcsc(TE,dTEdxi,1) ! Assemble and solve turbulence kinetic energy eq.
  call calcsc(ED,dEDdxi,2) ! Assemble and solve specific dissipation rate (omega [1/s]) of tke eq.
  call modify_mu_eff

end subroutine



subroutine calcsc(Fi,dFidxi,ifi)
!
!  Discretization and solution of scalar equation.
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
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, prtr_ijp, prtr_ijn, apotime, urfrs, urfms
  real(dp) :: utp, vtp, wtp, utn, vtn, wtn
  real(dp) :: genp, genn
  real(dp) :: uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2
  real(dp) :: viss,viste,dcoef
  real(dp) :: fimax,fimin
  real(dp) :: wldist,domegapl,ksi                           
  real(dp) :: dtedx,dtedy,dtedz,deddx,deddy,deddz
  real(dp) :: alphasst,bettasst,domega,vist
  real(dp) :: wlog,wvis
  ! real(dp) :: W_S,Ri,F4
  real(dp) :: gam,urf,tolAbs,tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 


  ! Variable specific coefficients:
  urf = TurbModel%Scalar(ifi)%urf
  gam = TurbModel%Scalar(ifi)%gds
  lSolver = TurbModel%Scalar(ifi)%lSolver 
  maxiter = TurbModel%Scalar(ifi)%maxiter
  tolAbs = TurbModel%Scalar(ifi)%tolAbs
  tolRel = TurbModel%Scalar(ifi)%tolRel


  ! Calculate gradient: 
  call grad(fi,dfidxi)

  ! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
!

  ! TKE volume source terms
  if(ifi.eq.1) then


  !=========================================================
  ! STANDARD PRODUCTION
  ! Note: 
  !   In find_strain_rate we calculate strain rate as:
  !   S = sqrt (2*Sij*Sij).
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
    ! gen(inp) = -den(inp)*( uu(inp)*dudx+uv(inp)*(dudy+dvdx)+ &
    !                        uw(inp)*(dudz+dwdx)+vv(inp)*dvdy+ &
    !                        vw(inp)*(dvdz+dwdy)+ww(inp)*dwdz )


    gen(inp)=abs(vis(inp)-viscos)*magStrain(inp)*magStrain(inp)

    ! Production limiter
    gen(inp)=min(gen(inp),0.9_dp*den(inp)*te(inp)*ed(inp))        

  enddo


  ! VOLUME SOURCE TERMS 
  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Add production term to the rhs:
    su(inp)=genp*vol(inp)  

    ! Add destruction term to the lhs:
    sp(inp)=bettast*ed(inp)*den(inp)*vol(inp)    

    ! If production negative move to lhs
    sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)



    ! VOLUME SOURCE TERMS: buoyancy
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
    if (ltransient) then

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


!****************************************
  elseif(ifi.eq.2) then
!****************************************


! Omega volume source terms


  do inp=1,numCells

    ! Wall distance
    wldist = walldistance(inp)

    ! Gradient of turbulence kinetic energy
    dtedx=dTEdxi(1,inp)
    dtedy=dTEdxi(2,inp)
    dtedz=dTEdxi(3,inp)

    ! Gradient of turbulence kinetic energy specific dissipation rate 
    deddx=dEDdxi(1,inp)
    deddy=dEDdxi(2,inp)
    deddz=dEDdxi(3,inp)

    ! Find $d_{\Omega}^{+}$ D_omega+
    domegapl=max( wldist*wldist/(ed(inp)+small) * (dtedx*deddx+dtedy*deddy+dtedz*deddz), 1e-20) ! ORIGINAL: 200*TEIN)

    ! Find ksi
    if (EARSM_M) then
      ksi=min(  max(                                                   &
                    sqrt(te(inp))/(BETTAST*wldist*ed(inp)+small),      &
                    500*viscos/den(inp)/(wldist**2*ed(inp)+small)      &
                   ),                                                  &
                2*te(inp)/domegapl                                     &   
              )
    else ! EARSM_WJ
      ksi=min(  max(                                                   &
                    sqrt(te(inp))/(BETTAST*wldist*ed(inp)+small),      &
                    500*viscos/den(inp)/(wldist**2*ed(inp)+small)      &
                   ),                                                  &
                20*te(inp)/domegapl                                    &
              )
    endif

    ! Find the blending function f_sst:
    if (EARSM_M) then
      fsst(inp) = tanh(ksi**4)
    else ! EARSM_WJ
      fsst(inp) = tanh(1.5_dp*ksi**4)
    endif

  enddo



  ! VOLUME SOURCE TERMS 
  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)


    ! Production of dissipation
    vist = (vis(inp)-viscos)/densit

    ! Production coefficient alpha_sst
    alphasst=fsst(inp)*alpha1+(1.0_dp-fsst(inp))*alpha2               

    su(inp)=alphasst*genp*vol(inp)/(vist+small)

    ! FIND $D_omega$ CROSS DIFFUSION MODIFICATION: 

    ! Gradient of turbulence kinetic energy
    dtedx=dTEdxi(1,inp)
    dtedy=dTEdxi(2,inp)
    dtedz=dTEdxi(3,inp)

    ! Gradient of turbulence kinetic energy specific dissipation rate 
    deddx=dEDdxi(1,inp)
    deddy=dEDdxi(2,inp)
    deddz=dEDdxi(3,inp)



    if (EARSM_M) then  
      domega = 2*(1.0_dp-fsst(inp))*den(inp)*SIGMOM2/(ed(inp)+small)*(dtedx*deddx+dtedy*deddy+dtedz*deddz)
    else ! EARSM_WJ
      ! (fsst(inp)+0.4_dp*(1.0_dp-fsst(inp))) = 0.4_dp+0.6_dp*fsst(inp)
      domega = (0.4_dp+0.6_dp*fsst(inp))*den(inp)/(ed(inp)+small)*(dtedx*deddx+dtedy*deddy+dtedz*deddz)
    endif

    domega = max(domega,0.0_dp)

    su(inp)=su(inp)+domega*vol(inp)



    ! Destruction of dissipation. 

    ! Destruction coefficient beta_sst
     bettasst=fsst(inp)*betai1+(1.0_dp-fsst(inp))*betai2


    ! Add destruction term (-beta*rho*w**2) to the lhs :   
    sp(inp)=bettasst*den(inp)*ed(inp)*vol(inp) 
    !
    !..or using the destruction term that incorporates Simplified Curvature Correction:
    ! Multiply destruction by F4 Simplified Curvature Correction term b Hellsten
    ! to obtain SST-2003RC-Hellsten model
    ! W_S = Vorticity(inp)/magStrain(inp)
    ! Ri = W_S*(W_S-one)
    ! F4 = one/(one + 1.4_dp*Ri)    
    ! sp(inp)=F4*bettasst*den(inp)*ed(inp)*vol(inp)

    ! Negative value of production moved to lhs.
    sp(inp)=sp(inp)-alphasst*genn*vol(inp)/(vist*ed(inp)+small) 



    ! UNSTEADY TERM
    if (ltransient) then

      apotime = den(inp)*vol(inp)/timestep

      if( bdf .or. cn ) then

        su(inp) = su(inp) + apotime*edo(inp)
        sp(inp) = sp(inp) + apotime

      elseif( bdf2 ) then

        su(inp) = su(inp) + apotime*( 2*edo(inp) - 0.5_dp*edoo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime

      endif

    endif

  ! End of Dissipation volume source terms
  enddo

!--------------------------------------
  end if ! ite or ied conditional



!
! CALCULATE CONVECTIVE AND DIFFUSIVE FLUXES BY INTEGRATING OVER FACES
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    ! In SST model the Effective diffusivity is a field variable:
    if(ifi.eq.1) then
      prtr_ijp = fsst(ijp)*sigmk1  + (1.0_dp-fsst(ijp))*sigmk2
      prtr_ijn = fsst(ijn)*sigmk1  + (1.0_dp-fsst(ijn))*sigmk2
    else
      prtr_ijp = fsst(ijp)*sigmom1 + (1.0_dp-fsst(ijp))*sigmom2
      prtr_ijn = fsst(ijn)*sigmom1 + (1.0_dp-fsst(ijn))*sigmom2
    endif

    prtr = half*(prtr_ijp+prtr_ijn)

    ! Diffusion coefficient
    viste = ( vis(ijp) + (vis(ijn)-vis(ijp))*facint(i) )-viscos
    dcoef = viscos+viste*prtr

    call facefluxsc(  i, ijp, ijn, &
                      xf(i), yf(i), zf(i), &
                      arx(i), ary(i), arz(i), &
                      flmass(i), facint(i), gam, &
                      TurbModel%Scalar(ifi)%cScheme, &
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
         bctype(ib) == 'pressure'  ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! In SST model the Effective diffusivity is a field variable:
        if(ifi.eq.1) then
          prtr=fsst(ijp)*sigmk1 + (1.0_dp-fsst(ijp))*sigmk2
        else
          prtr=fsst(ijp)*sigmom1 + (1.0_dp-fsst(ijp))*sigmom2
        endif

        ! Diffusion coefficient
        viste = vis(ijb)-viscos
        dcoef = viscos + viste*prtr 

        call facefluxsc( ijp, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         dfidxi, prtr, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp) - can*fi(ijb) + suadd

      end do



    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face



        ! In SST model the Effective diffusivity is a field variable:
        if(ifi.eq.1) then
          prtr_ijp = fsst(ijp)*SIGMK1  + (1.0_dp-fsst(ijp))*SIGMK2
          prtr_ijn = fsst(ijn)*SIGMK1  + (1.0_dp-fsst(ijn))*SIGMK2
        else
          prtr_ijp = fsst(ijp)*SIGMOM1 + (1.0_dp-fsst(ijp))*SIGMOM2
          prtr_ijn = fsst(ijn)*SIGMOM1 + (1.0_dp-fsst(ijn))*SIGMOM2
        endif

        prtr = half*(prtr_ijp+prtr_ijn )

        ! Diffusion coefficient
        viste = half*(vis(ijp)+vis(ijn))-viscos
        dcoef = viscos+viste*prtr

        ! face flux scalar but for periodic boundaries - it will be recognized by arguments
        call facefluxsc(  i, ijp, ijn, &
                          xf(if), yf(if), zf(if), &
                          arx(if), ary(if), arz(if), &
                          flmass(if), gam, &
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

        if (ifi .eq. 1) then

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

          Tau(iWall) = viss*Ut2/dnw(iwall)

          gen(ijp)=abs(tau(iWall))*cmu25*sqrt(te(ijp))/(dnw(iwall)*cappa)

          su(ijp)=su(ijp)+gen(ijp)*vol(ijp)

        else

          ! Wall boundaries approximated with wall functions
          ! for correct values of dissipation all coefficients have
          ! to be zero, su equal the dissipation, and ap = 1

          ! Automatic wall treatment - quadratic blend of log-layer and vis sublayer value:
          wlog=sqrt(te(ijp))/(cmu25*cappa*dnw(iWall))

          wvis=6.0_dp*(viscos/den(ijp))/(betai1*dnw(iWall)**2) 

          ed(ijp) = sqrt(wvis**2+wlog**2)
          su(ijp)=ed(ijp)

          a( ia(ijp):ia(ijp+1)-1 ) = 0.0_dp
          sp(ijp) = 1.0_dp

        endif

      enddo

    endif

  enddo

  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

      if(ifi.eq.1) then

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

      else ! dissipation

        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_index(i)
            su(ijp) = su(ijp) - a(k)*edo(ijn)

            k = jcell_icell_csr_index(i)
            su(ijn) = su(ijn) - a(k)*edo(ijp)
        enddo

        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ia(ijp) : ia(ijp+1)-1 ) ) - a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*edo(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo

      endif

  endif

  ! Underrelaxation factors
  urfrs=1.0_dp/urf
  urfms=1.0_dp-urf

  ! Main diagonal term assembly:
  do inp = 1,numCells

    ! Main diagonal term assembly:
    a(diag(inp)) = sp(inp) 
    do k = ia(inp),ia(inp+1)-1
      if (k.eq.diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) -  a(k)
    enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)
                    
  enddo

  

  ! Solve linear system:
  if (ifi.eq.1) then
    call csrsolve(lSolver, te, su, resor(5), maxiter, tolAbs, tolRel, 'k' )
  else
    call csrsolve(lSolver, ed, su, resor(6), maxiter, tolAbs, tolRel, 'Omega' )
  endif

  ! Update field values at boundaries
  call updateBoundary( fi )


  ! Report range of scalar values and clip if negative
  fimin = minval(fi)
  fimax = maxval(fi)

  if ( ifi == 1 ) then 
    write(*,'(2x,es11.4,a,es11.4)') fimin,' <= k <= ',fimax
  else
    write(*,'(2x,es11.4,a,es11.4)') fimin,' <= Omega <= ',fimax
  endif

  ! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


!***********************************************************************
!
subroutine modify_mu_eff
!
!***********************************************************************
!
! Update turbulent and effective viscosity.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  implicit none

  integer :: i,inp
  integer :: ib,iface,ijp,ijb,iwall
  real(dp) :: visold,urf
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Utau,viscw
  real(dp) :: Utauvis,Utaulog,Upl
  ! real(dp) :: fimax,fimin

  REAL(DP) :: DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ
  REAL(DP) :: TSTUR,TSVIS,TTS,HTTS                     ! Turbulence timescale related
  REAL(DP) :: S11,S12,S13,S22,S23,S33,S21,S32,W12,W13,W23 ! Components of strain and vorticity tensors
  REAL(DP) :: S11SQ,S22SQ,S33SQ,S12SQ,S13SQ,S23SQ         ! Strain rate tensor components squared
  REAL(DP) :: W12SQ,W13SQ,W23SQ                           ! Vorticity tensor components squared
  REAL(DP) :: SII,WII,WIIP3,SWWIV,SWWIVTT,SSWWV
  REAL(DP) :: TERM3C11,TERM3C12,TERM3C13,TERM3C22,TERM3C23
  REAL(DP) :: TERM4C11,TERM4C12,TERM4C13,TERM4C22,TERM4C23
  REAL(DP) :: TERM6C11,TERM6C12,TERM6C13,TERM6C22,TERM6C23
  REAL(DP) :: TERM9C11,TERM9C12,TERM9C13,TERM9C22,TERM9C23
  REAL(DP) :: P1,P2,SQP2,PM,PMP,SQPM,FACOS,FNC,DE,FII1,FII2,FN,FNSQ,Q,PQ,PQ1,Q1
  REAL(DP) :: BETA1,BETA3,BETA4,BETA6,BETA9
  ! REAL(DP) :: TERMLR,TERMLR11,TERMLR12,TERMLR13,TERMLR22,TERMLR23
  ! REAL(DP) :: beta1eq


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Start of legacy code for Explicit Algebraic Reynolds Stress model %%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      DO INP=1,numCells

!=====Velocity Gradient Tensor Invariants Computation=============================

!     Turbulence Timescale
      TSTUR = 1./(BETTAST*ED(INP))                         
      TSVIS = CT*SQRT( VISCOS /( DENSIT*BETTAST*ED(INP)*TE(INP) )  )
      TTS = MAX(TSTUR,TSVIS) ! Limiter for turbulence time-scale

!     Half of the time scale
      HTTS = 0.5*TTS

!     Velocity Gradients:
      dudx = dudxi(1,inp)
      dudy = dudxi(2,inp)
      dudz = dudxi(3,inp)

      dvdx = dvdxi(1,inp)
      dvdy = dvdxi(2,inp)
      dvdz = dvdxi(3,inp)

      dwdx = dwdxi(1,inp)
      dwdy = dwdxi(2,inp)
      dwdz = dwdxi(3,inp)


!.....FIND SCALAR INVARIANT OF THE STRAIN RATE AND VORTICITY TENSOR

!     Strain-rate and vorticity tensor components
      S11=TTS*dUdx
      S22=TTS*dVdy
      S33=TTS*dWdz

      S12=HTTS*(dUdy + dVdx)
      S13=HTTS*(dUdz + dWdx)
      S23=HTTS*(dVdz + dWdy)
 
      S21=S12     
      S32=S23

      W12=HTTS*(dUdy - dVdx)
      W13=HTTS*(dUdz - dWdx)
      W23=HTTS*(dVdz - dWdy)

!     Squares of strain-rate and vorticity tensor components
      S11SQ = S11*S11
      S22SQ = S22*S22
      S33SQ = S33*S33

      S12SQ = S12*S12
      S13SQ = S13*S13
      S23SQ = S23*S23  

      W12SQ = W12*W12
      W13SQ = W13*W13
      W23SQ = W23*W23

!     Second invariants of the strain rate and vorticity tensors
      SII = S11SQ + S22SQ + S33SQ + 2.0*(S12SQ+S13SQ+S23SQ)
      WII =-2.0*(W12SQ + W13SQ + W23SQ)
      WIIP3 = P3*WII ! One third of the invariant

!     Third invanriant of the strain rate and vorticity tensors:
      SWWIV =-S11*(W12SQ + W13SQ) - S22*(W12SQ + W23SQ)     &
            - S33*(W13SQ + W23SQ)                           &
            + 2.0*(-S12*W13*W23 + S13*W12*W23 - S23*W12*W13)
      SWWIVTT= TT*SWWIV ! Two thirds of the invariant

!     Fourth invariant of the strain rate and vorticity tensors
      SSWWV = 2.0*(-(S12*S13 + S22*S23 + S23*S33)*W12*W13    &
                  +(S11*S13 + S12*S23 + S13*S33)*W12*W23     &
                  -(S11*S12 + S12*S22 + S13*S23)*W13*W23)    &
        - (S11SQ + S12SQ + S13SQ)*(W12SQ + W13SQ)            &
        - (S12SQ + S22SQ + S23SQ)*(W12SQ + W23SQ)            &
        - (S13SQ + S23SQ + S33SQ)*(W13SQ + W23SQ)

!     Tensor component terms for beta 3
      TERM3C11 = - W12SQ - W13SQ - WIIP3
      TERM3C22 = - W12SQ - W23SQ - WIIP3
      TERM3C12 = - W13*W23
      TERM3C13 =   W12*W23
      TERM3C23 = - W12*W13

!     Tensor component terms for beta 4
      TERM4C11 =-2.0*(S12*W12 + S13*W13)
      TERM4C22 = 2.0*(S12*W12 - S23*W23)
      TERM4C12 = (S11-S22)*W12       - S23*W13       - S13*W23
      TERM4C13 =     - S23*W12 + (S11-S33)*W13       + S12*W23
      TERM4C23 =       S13*W12       + S12*W13 + (S22-S33)*W23

!     tensor component terms for beta 6
      TERM6C11 = -2.0*((S12*W13 - S13*W12)*W23       &
                     + S11*(W12SQ + W13SQ)) - SWWIVTT- WII*S11
      TERM6C22 = -2.0*((S23*W12 + S12*W23)*W13       &
                     + S22*(W12SQ + W23SQ)) - SWWIVTT- WII*S22
      TERM6C12 = -S12*(2.0*W12SQ + W13SQ + W23SQ)    &
               - (S13*W13-S23*W23)*W12 - (S11+S22)*W13*W23- WII*S12
      TERM6C13 = -S13*(W12SQ + 2.0*W13SQ + W23SQ)    &
               - (S12*W12+S23*W23)*W13 + (S11+S33)*W12*W23- WII*S13
      TERM6C23 = -S23*(W12SQ + W13SQ + 2.0*W23SQ)    &
               + (S12*W12-S13*W13)*W23 - (S22+S33)*W12*W13- WII*S23

!     Tensor component terms for beta 9
      TERM9C11 =-2.0*(( S12*W12 + S13*W13 - S23*W23)*W12SQ  &
                     +( S12*W12 + S13*W13 + S23*W23)*W13SQ  &
                     +( S22-S33)*W12*W13*W23)
      TERM9C22 =-2.0*((-S12*W12 - S13*W13 + S23*W23)*W12SQ  &
                     +(-S12*W12 + S13*W13 + S23*W23)*W23SQ  &
                     +(-S11+S33)*W12*W13*W23)
      TERM9C12 = ((S11-S22)*W12 - 2.0*(S13*W23+S23*W13))*W12SQ  &
               + ((S11-S33)*W12 - 2.0* S13*W23)         *W13SQ  &
               + ((S33-S22)*W12 - 2.0* S23*W13)         *W23SQ
      TERM9C13 = ((S11-S22)*W13 + 2.0* S12*W23)         *W12SQ  &
               + ((S11-S33)*W13 + 2.0*(S12*W23-S23*W12))*W13SQ  &
               + ((S22-S33)*W13 - 2.0* S23*W12)         *W23SQ
      TERM9C23 = ((S22-S11)*W23 + 2.0* S12*W13)         *W12SQ  &
               + ((S11-S33)*W23 + 2.0* S13*W12)         *W13SQ  &
               + ((S22-S33)*W23 + 2.0*(S12*W13+S13*W12))*W23SQ

!     Solution of the third degree equation for N_c
!.....Correction to C1p (A3') appearing in Anti's thesis.......
!      beta1eq = -6./5.* (4.05/(16.4025 -2.*WII))              !
!      C1P = C1P + 9./4. * 2.2 * max(1+beta1eq*SII,0.)         !
!      C1PSQ = C1P**2                                          !
!      C1P3 = C1P*P3                                           !
!.............................................................!
      P1     = (C1PSQ/27.0 + 0.45*SII - TT*WII)*C1P
      P2     = P1**2 - (C1PSQ*P3**2 + 0.9*SII + TT*WII)**3
      IF (P2 .GE. 0.0) THEN
        SQP2 = SQRT(P2)
        PM = P1 - SQP2
        PMP = ABS(PM)**P3
        FNC = C1P3 + (P1+SQP2)**P3 + SIGN(PMP,PM)
      ELSE
        PM = P1**2 - P2
        SQPM = SQRT(PM)
        FACOS = P3*ACOS(P1/SQPM)
        FNC = C1P3 + 2.0*(PM**P6)*COS(FACOS)
      END IF

!.....Improvement of the approximation of the N...........................
!     Nonlinear EARSMko2005 Model with better approximation for 3D Flows !
      DE = 20.0*(FNC**4)*(FNC - 0.5*C1P)     &
        - WII*(10.0*FNC + 15.0*C1P)*FNC**2  &
        + 10.0*C1P*WII**2
      FII1 = SWWIV**2
      FII2 = SSWWV - 0.5*SII*WII
      FN = FNC + 162.0*(FII1 + FII2*FNC**2)/DE
!........................................................................!

!     The denominator of the betas
      FNSQ = FN**2
      Q = 5.0*P6*(FNSQ - 2.0*WII)*(2.0*FNSQ - WII)     !%$ <<< Wallin&Johansson A.Hellsten

      IF (EARSM_M) THEN
        Q = (1./1.245)*(FNSQ - 2.0*WII)               !%$ <<< Menter 2009. 4% increase
        Q1=Q*P6*(2.0*FNSQ - WII)                      !
        PQ1=1.0/Q1                                    !
      ENDIF

      PQ = 1.0/Q


!     Coefficients (betas) (Menter 2009)
      IF (EARSM_M) THEN
        BETA1 = - PQ*FN
        BETA3 = - PQ1*2.0*SWWIV/FN
        BETA4 = - PQ
        BETA6 = - PQ1*FN  
        BETA9 =   PQ1  
      ELSE
!     Coefficients (betas) (WJ Original)
        BETA1 = - PQ*FN*(2.0*FNSQ - 7.0*WII)  !%$ <<< Wallin&Johansson A.Hellsten
        BETA3 = - PQ*12.0*SWWIV/FN            !%$ <<<
        BETA4 = - PQ*2.0*(FNSQ - 2.0*WII)     !%$ <<<
        BETA6 = - PQ*6.0*FN                   !%$ <<<
        BETA9 =   PQ*6.0                      !%$ <<<
      ENDIF

! !=====Low Re version===================================================
!       IF (LowRe) THEN
! !.....FIRST VECTOR-DISTANCE FROM GHOST CELL CENTER TO THIS CELL'S CENTER  
!       FAC= -1.                                                 
!       DXF=FAC*(XC(INB)-XC(INP))                                
!       DYF=FAC*(YC(INB)-YC(INP))    
!       DZF=FAC*(ZC(INB)-ZC(INP))
! !.....SECOND X THIRD DEFINE  ALWAYS CF AREA
!       ARE=DSQRT(AR1X(INB)**2+AR1Y(INB)**2+AR1Z(INB)**2)
! !.....COMPONENTS OF THE UNIT NORMAL VECTOR
!       ARER=1./ARE
!       ALF=AR1X(INB)*ARER
!       BET=AR1Y(INB)*ARER
!       GAM=AR1Z(INB)*ARER
! !.....NORMAL DISTANCE FROM SCALAR PRODUCT
!       DELN=(DXF*ALF+DYF*BET+DZF*GAM)

! !     Reynolds No. - Rey
!       REYLR = Densit*sqrt(TE(inp))*DELN/Viscos             

! !     factor for Low-Re(LR) correction - f1
!       FACLR = 1.-exp(-0.092*sqrt(REYLR) - 1.2e-4*REYLR**2) 

! !     405.*1.8**2/(216.*1.8-160.)
!       SIIEQ = 5.735139863                                  

! !     Additional term in anisotropy tensor expression 
!       TERMLR = (1.-FACLR**2)*1.4/max(SII,SIIEQ)
!       TERMLR11 = TERMLR * ((S11SQ+S12SQ+S13SQ)     - SII*P3)
!       TERMLR22 = TERMLR * ((S21*S12+S22SQ+S23*S32) - SII*P3)
!       TERMLR12 = TERMLR * (S11*S12+S12*S22+S13*S32)
!       TERMLR13 = TERMLR * (S11*S13+S12*S23+S13*S33)
!       TERMLR23 = TERMLR * (S21*S13+S22*S23+S23*S33)

! !     Extra anisotropy components.
!       AIJ(1,INP) = (TERMLR11                                         &           
!            + FACLR**2*BETA3*TERM3C11                                 & 
!            +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C11  &
!            + FACLR*   BETA6*TERM6C11                                 &
!            + FACLR**2*BETA9*TERM9C11)
!       AIJ(2,INP) = (TERMLR12                                         &
!            + FACLR**2*BETA3*TERM3C12                                 &
!            +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C12  &
!            + FACLR*   BETA6*TERM6C12                                 &
!            + FACLR**2*BETA9*TERM9C12)
!       AIJ(3,INP) = (TERMLR13                                         &
!            + FACLR**2*BETA3*TERM3C13                                 &
!            +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C13  &
!            + FACLR*   BETA6*TERM6C13                                 &
!            + FACLR**2*BETA9*TERM9C13)
!       AIJ(4,INP) = (TERMLR22                                         &
!            + FACLR**2*BETA3*TERM3C22                                 &
!            +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C22  &
!            + FACLR*   BETA6*TERM6C22                                 &
!            + FACLR**2*BETA9*TERM9C22)
!       AIJ(5,INP) = (TERMLR23                                         &
!            + FACLR**2*BETA3*TERM3C23                                 &
!            +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C23  &
!            + FACLR*   BETA6*TERM6C23                                 &
!            + FACLR**2*BETA9*TERM9C23)
! !=====End Of Low Re Version============================================
!       ELSE

      ! High-Re form:
      ! Extra anisotropy components.
      AIJ(1,INP) = BETA3*TERM3C11 + BETA4*TERM4C11 + BETA6*TERM6C11 + BETA9*TERM9C11
      AIJ(2,INP) = BETA3*TERM3C12 + BETA4*TERM4C12 + BETA6*TERM6C12 + BETA9*TERM9C12
      AIJ(3,INP) = BETA3*TERM3C13 + BETA4*TERM4C13 + BETA6*TERM6C13 + BETA9*TERM9C13
      AIJ(4,INP) = BETA3*TERM3C22 + BETA4*TERM4C22 + BETA6*TERM6C22 + BETA9*TERM9C22
      AIJ(5,INP) = BETA3*TERM3C23 + BETA4*TERM4C23 + BETA6*TERM6C23 + BETA9*TERM9C23
      ! It seems that aij(3,3)= -aij(1,1)-aij(2,2), 
      ! see p. 41 in V.Hamalainen's report (see references in module header):
      AIJ(6,INP) = -aij(1,inp)-aij(4,inp)

!       ENDIF ! LowRe

!     Effective C_mu:
      CMUEFF(INP) = -0.5*(BETA1 + WII*BETA6) 

!     Low Re correction for Effective Cmu:
!       IF (LowRe) CMUEFF(INP) = CMUEFF(INP) * FACLR

!     Use Limiter that Tom Gatski & Chris Rumsey use:
!      CMUEFF(INP) = max(CMUEFF(INP), 0.0005)


      END DO ! CELL LOOP

!=====END OF Velocity Gradient Tensor Invariants computation===========



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% End of legacy code for Explicit Algebraic Reynolds Stress model %%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  urf = TurbModel%urfVis


  ! Loop trough cells 
  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    ! Update effective viscosity:
    if(EARSM_M) then
      vis(inp)=viscos+den(inp)*te(inp)/(ed(inp)+small)
    else ! EARSM_WJ
      vis(inp)=viscos+den(inp)*(cmueff(inp)/BETTAST)*te(inp)/(ed(inp)+small)
    endif

    ! Underelaxation
    vis(inp)=visold + urf*(vis(inp)-visold)

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

        ! Magnitude of tangential velocity component
        Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)

        ! ! Tangent direction
        ! xtp = xtp/vtp
        ! ytp = ytp/vtp
        ! ztp = ztp/vtp

        ! ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        ! Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        ! Tau(iWall) = viscos*Ut2/dnw(iWall)
        ! Utau = sqrt( Tau(iWall) / den(ijb) )
        ! ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        ! Ima i ova varijanta...ovo je tehnicki receno ystar iliti y* a ne y+
        ! ypl(iWall) = den(ijp)*cmu25*sqrt(te(ijp))*dnw(iWall)/viscos


        ! *** Automatic wall treatment ***
 
        utau = sqrt( viscos*Vtp/(densit*dnw(iWall)) + cmu25*te(ijp) ) ! It's actually u* in original reference...

        ypl(iWall) = den(ijp)*Utau*dnw(iWall)/viscos 

        Utauvis=ypl(iWall)
        Utaulog=1.0/cappa*log(Elog*ypl(iWall)) 

        Upl=sqrt(sqrt(Utauvis**4+Utaulog**4)) 
          
        viscw = den(ijp)*utau*dnw(iWall)/Upl 

        ! Blended version of shear stress - probati ovo(!?)
        ! tau(iWall) = den(ijp) * (Vtp/Upl)**2

        ! Varijanta 2, u originalnoj referenci...
        tau(iWall) = den(ijp) * Vtp*Utau/Upl

        !*** END: Automatic wall treatment ***

        ! ! *** Enhanced wall treatment - Reichardt blending ***

        ! ! Below is a variant where we use Reichardt blending
        ! ! for whole span of y+ values.
        ! ! Some authors say that Reichardt function for u+ approximates
        ! ! the composite u+(y+) curve, better that Kader blending function.
 
        ! utau = sqrt( viscos*Vtp/(densit*dnw(iWall)) + cmu25*te(ijp) ) ! It's actually u* in original reference...

        ! ypl(iWall) = den(ijp)*Utau*dnw(iWall)/viscos 

        ! Uplblend = one/cappa*log(one+cappa*ypl(iWall)) + &
        !            7.8_dp*(1.-exp(-ypl(iWall)/11.0_dp)-(ypl(iWall)/11.0_dp)*exp(-ypl(iWall)/3.0_dp))
          
        ! viscw = den(ijp)*utau*dnw(iWall)/Uplblend  

        ! ! Blended version of shear stress - probati ovo(!?)
        ! ! tau(iWall) = den(ijp) * (Vtp/Uplblend)**2

        ! ! Varijanta 2, u originalnoj referenci...
        ! tau(iWall) = den(ijp) * Vtp*Utau/Uplblend

        ! !*** END: Enhanced wall treatment - Reichardt blending ***
        
        ! viscw = zero
        ! if(ypl(iWall) > ctrans) then
        !   viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
        ! endif
        ! Wall shear stress
        ! tau(iwall) = cappa*den(ijp)*Vtp*cmu25*sqrt(te(ijp))/log(Elog*ypl(iWall))

        visw(iWall) = max(viscos,viscw)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  ! fimin = minval(vis/viscos)
  ! fimax = maxval(vis/viscos)

  ! write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= Viscosity ratio <= ',fimax


end subroutine


subroutine modify_viscosity_inlet_k_omega_earsm_wj
!
! Update turbulent and effective viscosity at inlet.
!
  use types
  use parameters
  use geometry, only: numBoundaries,nfaces,iBndValueStart
  use variables

  implicit none

  integer :: i,ib,ijb

  !
  ! Boundary faces 
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        ijb = iBndValueStart(ib) + i

        Vis(ijb) = viscos+den(ijb)*te(ijb)/(ed(ijb)+small)

      end do

    endif

  enddo 

end subroutine


end module