module k_epsilon_rlzb_2lewt
!
! Implementation of Realizable k-epsilon two equation turbulence model (Shih et al. 1994).
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none

  ! Turbulence model constants
  real(dp), parameter :: CMU = 0.09_dp
  real(dp)            :: C1 ! Not a constant here.
  real(dp), parameter :: C2 = 1.90_dp
  real(dp), parameter :: C3 = 1.44_dp
  real(dp), parameter :: sigma_k = 1.0_dp
  real(dp), parameter :: sigma_epsilon = 1.2_dp
  real(dp), parameter :: a0 = 4.04 ! Original version 4.0, Fluent uses calibrated value of 4.04.

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(cmu))
  real(dp), parameter :: CMU75 = cmu25**3

  ! Coefficients for 2-layer approach
  ! real(dp), parameter :: Reyst = 75.0 ! Also a value of 200.0 is reported in Fluent documenttion  
  real(dp), parameter :: Reyst = 200.0 ! Also a value of 200.0 is reported in Fluent documenttion  
  real(dp), parameter :: Clst = cappa/cmu75
  real(dp), parameter :: Aeps = 2*Clst
  real(dp), parameter :: Amu = 70.0
  ! real(dp), parameter :: Ablend = 4.896499054173959  ! For Rey*=75 and alpha=0.15*Rey*
  real(dp), parameter :: Ablend = 13.057330811130559 ! For Rey*=200 and alpha=0.15*Rey*

  private 


  public :: modify_viscosity_k_epsilon_rlzb_2lewt
  public :: modify_viscosity_inlet_k_epsilon_rlzb_2lewt

contains


!***********************************************************************
!
subroutine modify_viscosity_k_epsilon_rlzb_2lewt()
!
!***********************************************************************
!
! Main module routine to solve turbulence model equations and update 
! effective viscosity.
!
!***********************************************************************
!
  use types
  use parameters
  use variables
  use gradients

  implicit none
!
!***********************************************************************
!
  call calcsc(TE,dTEdxi,ite) ! Assemble and solve turbulence kinetic energy eq.
  call calcsc(ED,dEDdxi,ied) ! Assemble and solve dissipation rate of tke eq.
  call modify_mu_eff()

end subroutine



!***********************************************************************
!
subroutine modify_viscosity_inlet_k_epsilon_rlzb_2lewt()
!
!***********************************************************************
!
! Update effective viscosity at inlet
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
!
  call modify_mu_eff_inlet()

end subroutine



!***********************************************************************
!
subroutine calcsc(Fi,dFidxi,ifi)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use title_mod

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  real(dp) :: gam, prtr, apotime, const, urfrs, urfms, &
              utp, vtp, wtp, utn, vtn, wtn, &
              genp, genn, &
              uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: etarlzb 
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real(dp) :: viss
  real(dp) :: fimax,fimin
  real(dp) :: k12,Rey,lambeps,leps,eps2l


! Variable specific coefficients:
  gam=gds(ifi)

  if(ifi.eq.ite) prtr=1.0_dp/sigma_k
  if(ifi.eq.ied) prtr=1.0_dp/sigma_epsilon

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================

  ! TKE volume source terms
  if(ifi.eq.ite) then

  !=========================================================
  ! STANDARD PRODUCTION
  !=========================================================

  do inp=1,numCells

    dudx = dudxi(1,inp)
    dudy = dudxi(2,inp)
    dudz = dudxi(3,inp)

    dvdx = dvdxi(1,inp)
    dvdy = dvdxi(2,inp)
    dvdz = dvdxi(3,inp)

    dwdx = dwdxi(1,inp)
    dwdy = dwdxi(2,inp)
    dwdz = dwdxi(3,inp)

    ! Minus here in fron because UU,UV,... calculated in calcstress hold -tau_ij
    ! So the exact production is calculated as tau_ij*dui/dxj
    gen(inp) = -den(inp)*( uu(inp)*dudx+uv(inp)*(dudy+dvdx)+ &
                           uw(inp)*(dudz+dwdx)+vv(inp)*dvdy+ &
                           vw(inp)*(dvdz+dwdy)+ww(inp)*dwdz )

    ! magStrainSq=magStrain(inp)*magStrain(inp)
    ! gen(inp)=abs(vis(inp)-viscos)*magStrainSq

  enddo


  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Add production term to the rhs:
    su(inp)=genp*vol(inp)  

    ! Add destruction term to the lhs:
    sp(inp)=ed(inp)*den(inp)*vol(inp)/(te(inp)+small)
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
           uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)/(t(inp)+273.15)
           vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)/(t(inp)+273.15)
           wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)/(t(inp)+273.15)
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
      if( bdf .or. cn ) then
        apotime = den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*teo(inp)
        sp(inp) = sp(inp) + apotime
      elseif( bdf2 ) then
        apotime=den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*( 2*teo(inp) - 0.5_dp*teoo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime
      endif

  ! End of TKE volume source terms
  enddo

!****************************************
  elseif(ifi.eq.ied) then
!****************************************

  ! Epsilon volume source terms

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================
  do inp=1,numCells

    genp=max(magStrain(inp),zero)
    genn=min(magStrain(inp),zero)

    ! Production of dissipation
    etarlzb = magStrain(inp)*te(inp)/(ed(inp)+small)
    c1 = max(0.43,etarlzb/(etarlzb+5.))
    su(inp)=c1*genp*ed(inp)*vol(inp)

    ! Destruction of dissipation
    sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/( te(inp)+sqrt(viscos/densit*ed(inp)) )

    ! Negative value of production moved to lhs.
    sp(inp) = sp(inp) - c1*genn*ed(inp)*vol(inp)

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
        else ! if(boussinesq.eq.0)
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
    !.....UNSTEADY TERM
    !=====================================
      if( bdf .or. cn ) then
        apotime = den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*edo(inp)
        sp(inp) = sp(inp) + apotime
      elseif( bdf2 ) then
        apotime=den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*( 2*edo(inp) - 0.5_dp*edoo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime
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

    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, &
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

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         Fi, dFidxi, prtr, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         FI, dFidxi, prtr, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

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

          ! viss=viscos
          ! if(ypl(iWall).gt.ctrans) viss=visw(iWall)
          viss=max(viscos,visw(iWall))

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

          Tau(iWall) = viss*Ut2/dnw(iWall)

          ! Production of TKE in wall adjecent cell
          ! First substract the standard production from source term
          su(ijp)=su(ijp)-gen(ijp)*vol(ijp)
          ! Calculate production for wall adjecent cell
          gen(ijp)=abs(tau(iWall))*cmu25*sqrt(te(ijp))/(dnw(iWall)*cappa)
          ! Add this production to source vector
          su(ijp)=su(ijp)+gen(ijp)*vol(ijp)

        ! else
        ! !
        ! ! > Wall boundary conditions for dissipation rate of turbulence kinetic energy eq.
        ! !

        !   ! Wall boundaries approximated with wall functions
        !   ! for correct values of dissipation all coefficients have
        !   ! to be zero, su equal the dissipation, and diagonal element a(diag(ijp)) = 1

        !   a( ioffset(ijp):ioffset(ijp+1)-1 ) = 0.0_dp
        !   sp(ijp) = 1.0_dp

        !   ed(ijp)=cmu75*te(ijp)**1.5/(cappa*dnw(iWall))
        !   su(ijp)=ed(ijp)

        endif

      enddo

    endif
    
  enddo ! Boundary conditions loop  



  if(ifi .eq.ied) then

  ! The two-layer apprach - enhanced wall treatment
  ! Penalty approach to set epsilon in 'near wall region' based on Rey
  ! This will overwite what is set for the wall adjecent layer in boundary conditions loop above
  ! and expand the region where epsilon is set to a specifi value - the whole region where Rey < Rey*

    ! Loop over cells
    do inp=1,numCells

      ! sqrt(tke)
      k12 = sqrt(abs(te(inp)))

      !
      ! Re number based on wall distance
      !
      Rey = den(inp)*wallDistance(inp)*k12/viscos

      if (Rey < Reyst) then

        !
        ! **Now perform Jonger blending (check parameters in the header of the module and change in needed)
        !

        ! Blending factor
        lambeps = 0.5_dp*( 1 + tanh( (Rey - Reyst)/Ablend ) )


        ! While here we will also produce new value for epsilon based on blending between two values
        leps = wallDistance(inp)*Clst*( 1.0-exp(-Rey/Aeps) ) ! the Chen-Patel length scale 

        ! Inner layer value for issipation
        eps2l = k12**3/leps 

        !
        ! **Now perform Jonger blending for dissipation
        !  
        ed(inp) = ( lambeps*ed(inp) + (1.0-lambeps)*eps2l )

        !
        ! Penalty formulation - not exactly but close
        !
        a( ioffset(inp):ioffset(inp+1)-1 ) = 0.0_dp
        sp(inp) = 1.0_dp
        su(inp)=ed(inp)

      endif

    enddo

  endif


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
  fimin = minval(fi(1:numCells))
  fimax = maxval(fi(1:numCells))
  
  write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= ',chvar(ifi),' <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi(1:numCells) = max(fi(1:numCells),small)

end subroutine calcsc


!***********************************************************************
!
subroutine modify_mu_eff()
!
! Update turbulent and effective viscosity.
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use variables
  implicit none
!
!***********************************************************************
!
  integer :: i,ib,inp
  integer :: iface, ijp,ijb,iWall,ivisc
  real(dp) :: visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: viscw
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real(dp) :: s11,s12,s13,s21,s22,s23,s31,s32,s33,w12,w13,w23
  real(dp) :: wrlzb,ffi,ass,ust,cmur,stild
  real(dp) :: k12,Rey,lmu,mut2l,mut,lambeps
  real(dp) :: Uplblend,Utau
  
  ivisc = 0

  !
  ! Loop trough cells 
  !
  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    ! sqrt(tke)
    k12 = sqrt(abs(te(inp)))

    !
    ! Re number based on wall distance
    !
    Rey = den(inp)*wallDistance(inp)*k12/viscos

    if (Rey < Reyst) ivisc = ivisc+1

    !
    ! Inner layer turbulent viscosity based on Wolfstein model
    !
    lmu = wallDistance(inp)*Clst*( 1.0-exp(-Rey/Amu) ) ! the length scale 
    mut2l = den(inp)*cmu*lmu*k12

    !
    ! Outer layer viscosity
    !
    dudx = dudxi(1,inp)
    dudy = dudxi(2,inp)
    dudz = dudxi(3,inp)

    dvdx = dvdxi(1,inp)
    dvdy = dvdxi(2,inp)
    dvdz = dvdxi(3,inp)

    dwdx = dwdxi(1,inp)
    dwdy = dwdxi(2,inp)
    dwdz = dwdxi(3,inp)

    ! Find strain rate tensor
    ! [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
    s11=dudx
    s12=0.5*(dudy+dvdx)
    s13=0.5*(dudz+dwdx)
    s22=dvdy
    s23=0.5*(dvdz+dwdy) 
    s33=dwdz
    s21=s12
    s31=s13
    s32=s23

    ! Find antisymmetric part of velocity gradient tensor
    ! [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
    w12=0.5*(dudy - dvdx)
    w13=0.5*(dudz - dwdx)
    w23=0.5*(dvdz - dwdy)

    ! FIND Stilda = sqrt (Sij*Sij) NOTE THERE'S NO FACTOR OF 2!!! like in S=sqrt (2*Sij*Sij)
    stild=sqrt(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2))

    ! W = Sij*Sjk*Ski/S
    wrlzb = (s11*s11*s11+s11*s12*s21+s11*s13*s31+ &
             s12*s21*s11+s12*s22*s21+s12*s23*s31+ &
             s13*s31*s11+s13*s32*s21+s13*s33*s31+ &
             s21*s11*s12+s21*s12*s22+s21*s13*s32+ &
             s22*s21*s12+s22*s22*s22+s22*s23*s32+ &
             s23*s31*s12+s23*s32*s22+s23*s33*s32+ &
             s31*s11*s13+s31*s12*s23+s31*s13*s33+ &
             s32*s21*s13+s32*s22*s23+s32*s23*s33+ &
             s33*s31*s13+s33*s32*s23+s33*s33*s33)/(stild**3)

    ffi = 1./3. * acos(max(-1.0d0, min( sqrt(6.0d0)*wrlzb ,1.0d0))) ! we make sure the argument of arccos lies in interval [-1,1]
    ass = dsqrt(6.0d0)*cos(ffi)
    ust = dsqrt(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2 + w12**2+w13**2+w23**2))

    cmur = 1.0d0/(a0 + ass*ust*te(inp)/(ed(inp) + small))

    mut = den(inp)*cmur*te(inp)**2/(ed(inp)+small)

    !
    ! **Now perform Jonger blending (check parameters in the header of the module and change in needed)
    !

    ! Blending factor
    lambeps = 0.5_dp*( 1 + tanh( (Rey - Reyst)/Ablend ) )

    ! Update effective viscosity:
    vis(inp)=viscos + ( lambeps*mut + (1.0-lambeps)*mut2l )

    ! Underelaxation
    vis(inp)=urf(ivis)*vis(inp)+(1.0_dp-urf(ivis))*visold

  enddo

  write(*,'(a,i0,a)') "  Two-layer model: Rey < Rey* in ",ivisc," cells."

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

        ! Its magnitude
        Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)

        ! ! Tangent direction
        ! xtp = xtp/vtp
        ! ytp = ytp/vtp
        ! ztp = ztp/vtp

        ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        ! Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        ! Tau(iWall) = viscos*Ut2/dnw(iWall)
        ! Utau = sqrt( Tau(iWall) / den(ijb) )
        ! ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        ! Ima i ova varijanta..ovo je tehnicki receno ystar iliti y* a ne y+
        ! ypl(iWall) = den(ijp)*cmu25*sqrt(te(ijp))*dnw(iWall)/viscos

        viscw = zero


        ! *** Enhanced wall treatment - Reichardt blending ***

        ! Below is a variant where we use Reichardt blending
        ! for whole span of y+ values.
        ! Some authors say that Reichardt function for u+ approximates
        ! the composite u+(y+) curve, better that Kader blending function.
 
        utau = sqrt( viscos*Vtp/(densit*dnw(iWall)) + cmu25*te(ijp) ) ! It's actually u* in original reference...

        ypl(iWall) = den(ijp)*Utau*dnw(iWall)/viscos 

        Uplblend = one/cappa*log(one+cappa*ypl(iWall)) + &
                   7.8_dp*(1.-exp(-ypl(iWall)/11.0_dp)-(ypl(iWall)/11.0_dp)*exp(-ypl(iWall)/3.0_dp))
          
        viscw = den(ijp)*utau*dnw(iWall)/Uplblend  

        ! Blended version of shear stress - probati ovo(!?)
        ! tau(iWall) = den(ijp) * (Vtp/Uplblend)**2

        ! Varijanta 2, u originalnoj referenci...
        tau(iWall) = den(ijp) * Vtp*Utau/Uplblend

        !*** END: Enhanced wall treatment - Reichardt blending ***


        ! if(ypl(iWall) > ctrans) then

        !   viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))

        !   ! Shear stress for log region
        !   tau(iwall) = cappa*den(ijp)*Vtp*cmu25*sqrt(te(ijp))/log(Elog*ypl(iWall))

        ! elseif( ypl(iWall) > 3.0 ) then

        !   !
        !   ! Enhanced wall treatment
        !   !
        !   Upvisc = ypl(iWall)
        !   Uplog  = log(Elog*ypl(iWall))/cappa
        !   Gmblend = -0.01_dp*ypl(iWall)**4/(1.+5*ypl(iWall))        
        !   Uplblend = exp(Gmblend)*Upvisc + exp(1./Gmblend)*Uplog          
        !   viscw = ypl(iWall)*viscos/Uplblend  

        !   ! Blended version of shear stress
        !   tau(iwall) = den(ijp) * (Vtp/Uplblend)**2

        ! else

        !   ! Shear stress for viscous region
        !   tau(iWall) = viscos*Ut2/dnw(iWall)
        !   Utau = sqrt( Tau(iWall) / den(ijp) )
        !   ypl(iWall) = den(ijp)*Utau*dnw(iWall)/viscos

        ! endif

        ! Set visw used in calcuvw wall bc treatment
        visw(iWall) = max(viscos,viscw)

        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo


end subroutine modify_mu_eff



!***********************************************************************
!
subroutine modify_mu_eff_inlet()
!
! Update effective viscosity for standard k-epsilon:
! \mu_{eff}=\mu+\mu_t; \mu_t = C_\mu * \frac{k^2}{\epsilon} 
! 
! NOTE: Although this is realizable k-eps we will use simple standard
! k-epsilon mu_t expression at inlet. 
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numBoundaries,nfaces,iBndValueStart
  use variables

  implicit none
!
!***********************************************************************
!
  integer :: i,ib,ijb

  !
  ! Boundary faces 
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        ijb = iBndValueStart(ib) + i

        Vis(ijb) = viscos+den(ijb)*te(ijb)**2*cmu/(ed(ijb)+small)

      end do

    endif

  enddo    

end subroutine modify_mu_eff_inlet


end module k_epsilon_rlzb_2lewt