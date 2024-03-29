module k_epsilon_rlzb
!
! Implementation of Realizable k-epsilon two equation turbulence model (Shih et al. 1994).
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


  private 


  public :: modify_viscosity_k_epsilon_rlzb
  public :: modify_viscosity_inlet_k_epsilon_rlzb

contains


subroutine modify_viscosity_k_epsilon_rlzb
!
! Purpose: 
!   Main module routine to solve turbulence model equations and update effective viscosity.
!
!

  call calcsc_tke     ! Assemble and solve turbulence kinetic energy eqn.
  call calcsc_epsilon ! Assemble and solve dissipation rate of tke eqn.
  call modify_mu_eff  ! Update viscosity

end subroutine


subroutine calcsc_tke
!
! Purpose:
!  Assemble and solve turbulence kinetic energy equation
!

  use sparse_matrix
  use linear_solvers
  
  implicit none

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, apotime, urfrs, urfms, &
              utp, vtp, wtp, utn, vtn, wtn, &
              genp, genn, &
              uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2
  real(dp) :: viss,viste,dcoef
  real(dp) :: fimax,fimin
  real(dp) :: gam, urf, tolAbs, tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 

! Variable specific coefficients:
  gam = TurbModel%Scalar(1)%gds
  urf = TurbModel%Scalar(1)%urf
  lSolver = TurbModel%Scalar(1)%lSolver 
  maxiter = TurbModel%Scalar(1)%maxiter
  tolAbs = TurbModel%Scalar(1)%tolAbs
  tolRel = TurbModel%Scalar(1)%tolRel

  prtr=1.0_dp/sigma_k

! Calculate gradient: 
  call grad(te,dTedxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

  !
  ! > VOLUME SOURCE TERMS 
  !

  !=========================================================
  ! STANDARD PRODUCTION
  !=========================================================

  do inp=1,numCells
    gen(inp)=abs(vis(inp)-viscos)*magStrain(inp)*magStrain(inp)
  enddo


  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Add production term to the rhs:
    su(inp)=genp*vol(inp)  

    ! Add destruction term to the lhs:
    sp(inp)=ed(inp)*den(inp)*vol(inp)/(te(inp)+small)
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



!
! > CALCULATE TERMS INTEGRATED OVER FACES
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
                     TurbModel%Scalar(1)%cScheme, &
                     TE, dTEdxi, dcoef, cap, can, suadd )

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

        ! Diffusion coefficient
        viste = vis(ijb)-viscos
        dcoef = viscos + viste*prtr

        call facefluxsc( ijp, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         dtedxi, dcoef, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*te(ijb) + suadd

      end do

    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face

        ! Diffusion coefficient
        viste = half*(vis(ijp)+vis(ijn))-viscos
        dcoef = viscos+viste*prtr

        ! face flux scalar but for periodic boundaries - it will be recognized by arguments
        call facefluxsc(  i, ijp, ijn, &
                          xf(if), yf(if), zf(if), &
                          arx(if), ary(if), arz(if), &
                          flmass(if), gam, &
                          te, dTedxi, dcoef, cap, can, suadd )


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

      enddo

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

    ! Main diagonal term assembly: a_p = sp-sum( a_nb )
    a(diag(inp)) = sp(inp) 
    ! now add to main diagonal the negative sum of off-diagonal elements:
    do k = ia(inp),ia(inp+1)-1
      if (k.eq.diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) -  a(k)
    enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*te(inp)
                    
  enddo

  ! Solve linear system:
  call csrsolve(lSolver, te, su, resor(5), maxiter, tolAbs, tolRel, 'k' )


  ! Update field values at boundaries
  call updateBoundary( te )


! Report range of scalar values and clip if negative
  fimin = minval( te(1:numCells) )
  fimax = maxval( te(1:numCells) )
   
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= k <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) te(1:numCells) = max( te(1:numCells),small )

end subroutine



subroutine calcsc_epsilon
!
! Assemble and solve turbulence dissipation rate equation.
!

  use sparse_matrix
  use linear_solvers
  
  implicit none

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, apotime, const, urfrs, urfms, &
              utp, vtp, wtp, utn, vtn, wtn, &
              genp, genn, &
              uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: viste,dcoef
  real(dp) :: etarlzb 
  real(dp) :: off_diagonal_terms
  real(dp) :: fimax,fimin
  real(dp) :: gam, urf, tolAbs, tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 

! Variable specific coefficients:
  gam = TurbModel%Scalar(2)%gds
  urf = TurbModel%Scalar(2)%urf
  lSolver = TurbModel%Scalar(2)%lSolver 
  maxiter = TurbModel%Scalar(2)%maxiter
  tolAbs = TurbModel%Scalar(2)%tolAbs
  tolRel = TurbModel%Scalar(2)%tolRel

  prtr=1.0_dp/sigma_epsilon

! Calculate gradient: 
  call grad(ED,dEDdxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! > VOLUME SOURCE TERMS 
!
  do inp=1,numCells

    genp=max(magStrain(inp),zero)
    genn=min(magStrain(inp),zero)

    ! Production of dissipation
    etarlzb = magStrain(inp)*te(inp)/(ed(inp)+small)
    c1 = max(0.43,etarlzb/(etarlzb+5.0))
    su(inp)=c1*genp*ed(inp)*vol(inp)

    ! Destruction of dissipation
    sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/( te(inp)+sqrt(viscos/densit*ed(inp))+small )

    ! Negative value of production moved to lhs.
    sp(inp) = sp(inp) - c1*genn*ed(inp)*vol(inp)

  
    !=====================================
    ! VOLUME SOURCE TERMS: Buoyancy
    !=====================================
    if(lbuoy) then

      const=c3*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)

      if(boussinesq) then
         uttbuoy=-gravx*utt(inp)*const*beta
         vttbuoy=-gravy*vtt(inp)*const*beta
         wttbuoy=-gravz*wtt(inp)*const*beta
      else ! if(boussinesq.eq.0)
         uttbuoy=-gravx*utt(inp)*const/(t(inp)+small)
         vttbuoy=-gravy*vtt(inp)*const/(t(inp)+small)
         wttbuoy=-gravz*wtt(inp)*const/(t(inp)+small)
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

  ! End of Epsilon volume source terms
  enddo

!
! > CALCULATE TERMS INTEGRATED OVER FACES
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
                     TurbModel%Scalar(2)%cScheme, &
                     ed, deddxi, dcoef, cap, can, suadd )

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

        ! Diffusion coefficient
        viste = vis(ijb)-viscos
        dcoef = viscos + viste*prtr

        call facefluxsc( ijp, &
                         xf(iface), yf(iface), zf(iface), &
                         arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         deddxi, dcoef, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*ed(ijb) + suadd

      end do


    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face

        ! Diffusion coefficient
        viste = half*(vis(ijp)+vis(ijn))-viscos
        dcoef = viscos+viste*prtr

        ! face flux scalar but for periodic boundaries - it will be recognized by arguments
        call facefluxsc(  i, ijp, ijn, &
                          xf(if), yf(if), zf(if), &
                          arx(if), ary(if), arz(if), &
                          flmass(if), gam, &
                          Ed, dEddxi, dcoef, cap, can, suadd )


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
        ! > Wall boundary conditions for dissipation rate of turbulence kinetic energy eq.
        !

        ! Wall boundaries approximated with wall functions
        ! for correct values of dissipation all coefficients have
        ! to be zero, su equal the dissipation, and diagonal element a(diag(ijp)) = 1

        a( ia(ijp):ia(ijp+1)-1 ) = 0.0_dp
        sp(ijp) = 1.0_dp

        ed(ijp)=cmu75*te(ijp)**1.5/(cappa*dnw(iWall))
        su(ijp)=ed(ijp)

      enddo

    endif
    
  enddo ! Boundary conditions loop  

  ! Modify coefficients for Crank-Nicolson
  if (cn) then

    a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

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
    su(inp) = su(inp) + urfms*a(diag(inp))*ed(inp)
                    
  enddo

  ! Solve linear system:
  call csrsolve(lSolver, ed, su, resor(6), maxiter, tolAbs, tolRel, 'epsilon' )



  ! Update field values at boundaries
  call updateBoundary( ed )


! Report range of scalar values and clip if negative
  fimin = minval( ed(1:numCells) )
  fimax = maxval( ed(1:numCells) )
  
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= epsilon <= ',fimax


! These field values cannot be negative
  if(fimin.lt.0.0_dp) ed(1:numCells) = max(ed(1:numCells),small)

end subroutine



subroutine modify_mu_eff
!
! Update turbulent and effective viscosity.
!

  use types
  use parameters
  use geometry
  use variables

  implicit none

  integer :: i,ib,inp
  integer :: iface, ijp,ijb,iWall
  real(dp) :: urf,visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp,Ut2
  real(dp) :: viscw
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real(dp) :: s11,s12,s13,s21,s22,s23,s31,s32,s33,w12,w13,w23
  real(dp) :: wrlzb,ffi,ass,ust,cmur,vist,stild

  urf = TurbModel%urfVis

  !
  ! Loop trough cells 
  !
  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

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

    vist = den(inp)*cmur*te(inp)**2/(ed(inp)+small)

    ! Update effective viscosity:
    vis(inp)=viscos + vist

    ! Underelaxation
    vis(inp)=urf*vis(inp)+(1.0_dp-urf)*visold

  enddo


  !
  ! Boundary faces 
  !
  ! Update value at every bounary face type except 'wall', which is treated below.
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

        ! Tau(iWall) = viscos*Ut2/dnw(iWall)
        ! Utau = sqrt( Tau(iWall) / den(ijb) )
        ! ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        ! Ima i ova varijanta...ovo je tehnicki receno ystar iliti y* a ne y+
        ypl(iWall) = den(ijp)*cmu25*sqrt(te(ijp))*dnw(iWall)/viscos

        ! Wall shear stress
        tau(iwall) = cappa*den(ijp)*Vtp*cmu25*sqrt(te(ijp))/log(Elog*ypl(iWall))

        viscw = zero

        ! Standard wall function - first cell in log region
        if( ypl(iWall) > ctrans ) then
          viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
        endif

        visw(iWall) = max(viscos,viscw)

        vis(ijb) = visw(iWall)



      enddo

    endif 

  enddo


end subroutine


subroutine modify_viscosity_inlet_k_epsilon_rlzb
!
! Update effective viscosity for standard k-epsilon.
! 
! NOTE: Although this is realizable k-eps we will use simple standard
! k-epsilon mu_t expression at inlet. 
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

        Vis(ijb) = viscos+den(ijb)*te(ijb)**2*cmu/(ed(ijb)+small)

      end do

    endif

  enddo    

end subroutine


end module k_epsilon_rlzb