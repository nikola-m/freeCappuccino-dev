module k_omega_SST
!
! Implementation of k-omega Shear Stress Transport (SST) two equation turbulence model.
!
! REFERENCES:
!     * ANSYS FLUENT Theory Guide
!     * Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications",
!       AIAA Journal, Vol. 32, No. 8, August 1994, pp. 1598-1605. 
!     * Menter, F. R., Kuntz, M., and Langtry, R., "Ten Years of Industrial Experience with the SST Turbulence Model",
!       Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano, and M. Tummers, Begell House, Inc., 2003, pp. 625 - 632. 
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none

  logical :: LowRe = .false. ! Has to be set in calling routine or in main program.

  ! Turbulence model constants 
  real(dp), parameter :: BETTAST=0.09_dp   
  real(dp), parameter :: SIGMK1=0.85_dp
  real(dp), parameter :: SIGMK2=1.0_dp
  real(dp), parameter :: SIGMOM1=0.5_dp
  real(dp), parameter :: SIGMOM2=0.856_dp
  real(dp), parameter :: BETAI1=0.075_dp
  real(dp), parameter :: BETAI2=0.0828_dp
  real(dp), parameter :: A1=0.31_dp

! SST-1994 coefficients
!  ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM1)
!  ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM2)

! SST-2003 coefficients. The first is higher than the original constant
! definition by approximately 0.43%, and the second is lower by less than 0.08%. 
  real(dp), parameter :: ALPHA1=5./9.0_dp
  real(dp), parameter :: ALPHA2=0.44_dp

  real(dp), parameter :: C3 = 1.44_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(BETTAST))
  real(dp), parameter :: CMU75 = cmu25**3

  real(dp), dimension(:), allocatable :: fsst

  private 

  public :: LowRe
  public :: modify_viscosity_k_omega_sst
  public :: modify_viscosity_inlet_k_omega_sst

contains


!***********************************************************************
!
subroutine modify_viscosity_k_omega_sst()
!
!***********************************************************************
!
! Main module routine to solve turbulence model equations and update effective viscosity.
!
!***********************************************************************
  use types
  use parameters
  use variables
  use gradients
  implicit none

  if(.not.allocated(fsst)) then

    allocate(fsst(numTotal))

    write(*,'(a)') ' '
    write(*,'(a)') '  **Allocated SST blending function.'
    write(*,'(a)') ' '
    
  endif

  call calcsc(TE,dTEdxi,ite) ! Assemble and solve turbulence kinetic energy eq.
  call calcsc(ED,dEDdxi,ied) ! Assemble and solve specific dissipation rate (omega [1/s]) of tke eq.
  call modify_mu_eff()

end subroutine


!***********************************************************************
!
subroutine calcsc(Fi,dFidxi,ifi)
!
!***********************************************************************
!
!  Discretization of scalar equation
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use title_mod

  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  real(dp) :: gam, prtr, prtr_ijp, prtr_ijn, apotime, const, urfrs, urfms
  real(dp) :: utp, vtp, wtp, utn, vtn, wtn
  real(dp) :: genp, genn
  real(dp) :: uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: magStrainSq
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2
  ! real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real(dp) :: viss
  real(dp) :: fimax,fimin
  real(dp) :: wldist,domegapl,ksi,tmp                            
  real(dp) :: dtedx,dtedy,dtedz,deddx,deddy,deddz
  real(dp) :: alphast,alphasst,bettasst,domega,vist
  real(dp) :: wlog,wvis
  ! real(dp) :: W_S,Ri,F4
 

! Variable specific coefficients:
  gam=gds(ifi)


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
  if(ifi.eq.ite) then


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


    magStrainSq=magStrain(inp)*magStrain(inp)
    gen(inp)=abs(vis(inp)-viscos)*magStrainSq

    ! PRODUCTION LIMITER FOR SST AND SAS MODELS:
    ! 10*bettainf=10*0.09=0.9

    ! High-Re version...............................................................
      gen(inp)=min(gen(inp),0.9_dp*den(inp)*te(inp)*ed(inp))        

      if (LowRe) then
    ! Low-Re version of Wilcox and SST k-omega......................................
        tmp=10*bettast*(4./15.0_dp+(den(inp)*te(inp)/(8.0_dp*viscos*ed(inp)))**4) & !
                      /(1.0_dp    +(den(inp)*te(inp)/(8.0_dp*viscos*ed(inp)))**4)   !  
        gen(inp)=min(gen(inp),tmp*den(inp)*te(inp)*ed(inp))                         !
    !...............................................................................!
      end if 

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

    !======================================================================
    !     Note there is possibility to add a source term to eliminate 
    !     non-physical decay of turbulence variables in the freestream
    !     for external aerodynamic problems
    !     Reference:
    !     Spalart, P. R. and Rumsey, C. L., "Effective Inflow Conditions for 
    !     Turbulence Models in Aerodynamic Calculations," AIAA Journal,
    !     Vol. 45, No. 10, 2007, pp. 2544 - 2553.
    !======================================================================
    ! ADD SUSTAIN TERMS (ONLY FOR SST!):
    ! su(inp)=su(inp)+bettast*tein*edin*den(inp)*vol(inp)

    ! Add destruction term to the lhs:

    !.....High-Re version.....................................................
    sp(inp)=bettast*ed(inp)*den(inp)*vol(inp)    

    if(LowRe) then                                        
    !.....Low-Re version of Wilcox and SST k-omega.............................
        tmp = BETTAST*(4./15.0_dp+(den(inp)*te(inp)/(8*viscos*ed(inp)))**4) & !
                    /(1.0_dp    +(den(inp)*te(inp)/(8*viscos*ed(inp)))**4)    !
        sp(inp)=tmp*ed(inp)*den(inp)*vol(inp)                                 !           
    !.........................................................................!
    endif

    ! If gen negative move to lhs
    sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)


    !
    !=====================================
    ! VOLUME SOURCE TERMS: buoyancy
    !=====================================
      if(lcal(ien).and.lbuoy) then
        
        ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
        call calcheatflux 

        if(boussinesq) then
           uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)*beta
           vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)*beta
           wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)*beta
        else
           uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)/(t(inp)+273.15_dp)
           vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)/(t(inp)+273.15_dp)
           wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)/(t(inp)+273.15_dp)
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

      !
      !=====================================
      ! UNSTEADY TERM
      !=====================================
      if (ltransient) then
        if( bdf .or. cn ) then
          apotime = den(inp)*vol(inp)/timestep
          su(inp) = su(inp) + apotime*teo(inp)
          sp(inp) = sp(inp) + apotime
        elseif( bdf2 ) then
          apotime=den(inp)*vol(inp)/timestep
          su(inp) = su(inp) + apotime*( 2*teo(inp) - 0.5_dp*teoo(inp) )
          sp(inp) = sp(inp) + 1.5_dp*apotime
        endif
      endif

  ! End of TKE volume source terms
  enddo

!****************************************
  elseif(ifi.eq.ied) then
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

    ! Find $d_{\omega}^{+}$ d_omega+
    domegapl=max(2*SIGMOM2*den(inp)/(ed(inp)) * (dtedx*deddx+dtedy*deddy+dtedz*deddz),1e-20)

    ! Find ksi
    ksi=min(  max(                                                   &
                  sqrt(te(inp))/(BETTAST*wldist*ed(inp)+small),      &
                  500.0_dp*viscos/den(inp)/(wldist**2*ed(inp)+small) &
                 ),                                                  &
              4.0_dp*den(inp)*te(inp)*SIGMOM2/(domegapl*wldist**2)   &
            )

    ! Find the SST model blending function f_sst:
    fsst(inp) = tanh(ksi**4)

  enddo

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================
  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)


    ! Production of dissipation
    vist = (vis(inp)-viscos)/densit

    ! Production coefficient alpha_sst
    !.....High-Re version...............................................
      alphasst=fsst(inp)*alpha1+(1.0_dp-fsst(inp))*alpha2               !<

      If(LowRe) then
     ! Low-Re version of SST k-omega......................................
      alphast=(0.024_dp+(densit*te(inp))/(6.0_dp*viscos*ed(inp)))  &    !             
            /(1.0_dp+(densit*te(inp))/(6.0_dp*viscos*ed(inp)))          !
      tmp=alpha1/alphast*                                   &           !                                     
             (1./9.0_dp+ (densit*te(inp))/(2.95_dp*viscos*ed(inp))) &   !
            /(1.0_dp  + (densit*te(inp))/(2.95_dp*viscos*ed(inp)))      !
      alphasst=fsst(inp)*tmp + (1.0_dp-fsst(inp))*alpha2                !<
      !.................................................................!
      endif

    su(inp)=alphasst*genp*vol(inp)/(vist+small)

    ! FIND D_omega CROSS DIFFUSION MODIFICATION: 

      ! Gradient of turbulence kinetic energy
      dtedx=dTEdxi(1,inp)
      dtedy=dTEdxi(2,inp)
      dtedz=dTEdxi(3,inp)

      ! Gradient of turbulence kinetic energy specific dissipation rate 
      deddx=dEDdxi(1,inp)
      deddy=dEDdxi(2,inp)
      deddz=dEDdxi(3,inp)

      domega = 2*(1.0_dp-fsst(inp))*den(inp)*SIGMOM2/(ed(inp)+small)*(dtedx*deddx+dtedy*deddy+dtedz*deddz)
      domega = max(domega,0.0_dp)

    su(inp)=su(inp)+domega*vol(inp)



    ! Destruction of dissipation. 

    ! Destruction coefficient beta_sst
     bettasst=fsst(inp)*betai1+(1.0_dp-fsst(inp))*betai2

    ! ADD SUSTAIN TERMS
    ! su(inp)=su(inp)+bettasst*edin*edin*den(inp)*vol(inp)


    ! Add destruction term (-beta*rho*w**2) to the lhs :   
    sp(inp)=bettasst*den(inp)*ed(inp)*vol(inp) 
    !..or using the destruction term that incorporates Simplified Curvature Correction:
    ! Multiply destruction by F4 Simplified Curvature Correction term b Hellsten
    ! to obtain SST-2003RC-Hellsten model
    ! W_S = Vorticity(inp)/magStrain(inp)
    ! Ri = W_S*(W_S-one)
    ! F4 = one/(one + 1.4_dp*Ri)    
    ! sp(inp)=F4*bettasst*den(inp)*ed(inp)*vol(inp)

    ! Negative value of production moved to lhs.
    sp(inp)=sp(inp)-alphasst*genn*vol(inp)/(vist*ed(inp)+small) 

    !
    !=====================================
    ! VOLUME SOURCE TERMS: Buoyancy
    !=====================================
      if(lcal(ien).and.lbuoy) then
        const=c3*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)

        if(boussinesq) then
           uttbuoy=-gravx*utt(inp)*const*beta
           vttbuoy=-gravy*vtt(inp)*const*beta
           wttbuoy=-gravz*wtt(inp)*const*beta
        else
           uttbuoy=-gravx*utt(inp)*const/(t(inp)+273.15)
           vttbuoy=-gravy*vtt(inp)*const/(t(inp)+273.15)
           wttbuoy=-gravz*wtt(inp)*const/(t(inp)+273.15)
        end if

        utp=max(uttbuoy,zero)
        vtp=max(vttbuoy,zero)
        wtp=max(wttbuoy,zero)
        utn=min(uttbuoy,zero)
        vtn=min(vttbuoy,zero)
        wtn=min(wttbuoy,zero)

        su(inp)=su(inp)+utp+vtp+wtp
        sp(inp)=sp(inp)-(utn+vtn+wtn)/(ed(inp)+small)
      end if

      !
      !=====================================
      ! UNSTEADY TERM
      !=====================================
      if (ltransient) then
        if( bdf .or. cn ) then
          apotime = den(inp)*vol(inp)/timestep
          su(inp) = su(inp) + apotime*edo(inp)
          sp(inp) = sp(inp) + apotime
        elseif( bdf2 ) then
          apotime=den(inp)*vol(inp)/timestep
          su(inp) = su(inp) + apotime*( 2*edo(inp) - 0.5_dp*edoo(inp) )
          sp(inp) = sp(inp) + 1.5_dp*apotime
        endif
      endif

  ! End of Epsilon volume source terms
  enddo
!--------------------------------------
  end if

!
! CALCULATE TERMS INTEGRATED OVER FACES
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    ! In SST model the Effective diffusivity is a field variable:
    if(ifi.eq.ite) then
      prtr_ijp = fsst(ijp)*SIGMK1  + (1.0_dp-fsst(ijp))*SIGMK2
      prtr_ijn = fsst(ijn)*SIGMK1  + (1.0_dp-fsst(ijn))*SIGMK2
    else
      prtr_ijp = fsst(ijp)*SIGMOM1 + (1.0_dp-fsst(ijp))*SIGMOM2
      prtr_ijn = fsst(ijn)*SIGMOM1 + (1.0_dp-fsst(ijn))*SIGMOM2
    endif

    call facefluxsc(  ijp, ijn, &
                      xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                      flmass(i), facint(i), gam, &
                      fi, dFidxi, prtr_ijp, cap, can, suadd )                      
                      ! fi, dFidxi, prtr_ijp, prtr_ijn, cap, can, suadd )

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

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! In SST model the Effective diffusivity is a field variable:
        if(ifi.eq.ite) then
          prtr=fsst(ijp)*SIGMK1  + (1.0_dp-fsst(ijp))*SIGMK2
        else
          prtr=fsst(ijp)*SIGMOM1 + (1.0_dp-fsst(ijp))*SIGMOM2
        endif

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         Fi, dFidxi, prtr, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp) - can*fi(ijb) + suadd

      end do


    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! In SST model the Effective diffusivity is a field variable:
        if(ifi.eq.ite) then
          prtr=fsst(ijp)*SIGMK1  + (1.0_dp-fsst(ijp))*SIGMK2
        else
          prtr=fsst(ijp)*SIGMOM1 + (1.0_dp-fsst(ijp))*SIGMOM2
        endif

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         FI, dFidxi, prtr, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp) - can*fi(ijb) + suadd

      end do

    elseif ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        if (ifi .eq. ite) then

        !
        ! > Wall boundary conditions for turbulence kinetic energy eq.
        !

          su(ijp )= su(ijp)-gen(ijp)*vol(ijp) ! take out standard production from wall ajdecent cell.

          ! viss=viscos
          ! if(ypl(i).gt.ctrans) viss=visw(iwall)
          viss = max(viscos,visw(iWall))

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

          a( ioffset(ijp):ioffset(ijp+1)-1 ) = 0.0_dp
          sp(ijp) = 1.0_dp

        endif

      enddo

    endif

  enddo

  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

      if(ifi.eq.ite) then

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

      else ! ifi.eq.ied

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
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) - a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*edo(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo

      endif

  endif

  ! Underrelaxation factors
  urfrs=urfr(ifi)
  urfms=urfm(ifi)

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
  call bicgstab(fi,ifi)

  !
  ! Update symmetry and outlet boundaries
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'outlet' .or. bctype(ib) == 'symmetry' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        fi(ijb)=fi(ijp)

      enddo

    endif

  enddo

! Report range of scalar values and clip if negative
  fimin = minval(fi)
  fimax = maxval(fi)

  write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= ',chvar(ifi),' <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


!***********************************************************************
!
subroutine modify_mu_eff()
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
  real(dp) :: visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Utau,viscw
  real(dp) :: wldist,etha,f2_sst,alphast
  real(dp) :: Utauvis,Utaulog,Upl
  ! real(dp) :: fimax,fimin


  ! Loop trough cells 
  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    ! Update effective viscosity:

    ! Wall distance
    wldist = walldistance(inp)

    ! find etha:
    etha=max(2*sqrt(te(inp))/(bettast*wldist*ed(inp)), &
             (500*viscos/den(inp))/(wldist**2*ed(inp))) 

    ! find f2: 
    f2_sst = tanh(etha*etha)

    vis(inp)=viscos+den(inp)*a1*te(inp)/(max(a1*ed(inp), magStrain(inp)*f2_sst))

    ! Low-re version..........................................................
    if (LowRe) then                                                          !
    ! Let's find alpha*                                                      !                                         
      alphast=(0.024_dp+(densit*te(inp))/(6*viscos*ed(inp)))   &             !           
             /(1.0_dp+(densit*te(inp))/(6*viscos*ed(inp)))                   !
      vis(inp)=viscos+den(inp)*te(inp)/(ed(inp)+small)               &       !  
            *1.0_dp/max(1.0_dp/alphast, magStrain(inp)*f2_sst/(a1*ed(inp)))  !                                                   
    ! End of low-re version..................................................!
    end if


    ! Underelaxation
    vis(inp)=urf(ivis)*vis(inp)+(1.0_dp-urf(ivis))*visold

  enddo


  !
  ! Boundary faces 
  !

  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then

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
        ! tau(iWall) = den(ijp) * (Vtp/Uplblend)**2

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

        visw(iWall) = max(viscos,viscw)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  ! fimin = minval(vis/viscos)
  ! fimax = maxval(vis/viscos)

  ! write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= Viscosity ratio <= ',fimax


end subroutine modify_mu_eff


subroutine modify_viscosity_inlet_k_omega_sst()
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