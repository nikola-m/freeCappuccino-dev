module k_epsilon_rng
!
! Implementation of k-epsilon Renormalisation group (RNG) two equation turbulence model.
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
  real(dp), parameter :: CMU = 0.0845_dp
  real(dp), parameter :: C1 = 1.42_dp
  real(dp), parameter :: C2 = 1.68_dp
  real(dp), parameter :: C3 = 1.44_dp
  real(dp), parameter :: sigma_k = 0.7194_dp
  real(dp), parameter :: sigma_epsilon = 0.7194_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(cmu))
  real(dp), parameter :: CMU75 = cmu25**3


  private 


  public :: modify_viscosity_k_epsilon_rng
  public :: modify_viscosity_inlet_k_epsilon_rng

contains



subroutine modify_viscosity_k_epsilon_rng()
!
! Main module routine to solve turbulence model equations and update effective viscosity.
!

  implicit none

  call calcsc(TE,dTEdxi,1) ! Assemble and solve turbulence kinetic energy eq.
  call calcsc(ED,dEDdxi,2) ! Assemble and solve specific dissipation rate (omega [1/s]) of tke eq.
  call modify_mu_eff

end subroutine


!***********************************************************************
!
subroutine calcsc(fi,dfidxi,ifi)
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
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, apotime, const, urfrs, urfms
  real(dp) :: utp, vtp, wtp, utn, vtn, wtn
  real(dp) :: genp, genn
  real(dp) :: uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2
  real(dp) :: viss,viste,dcoef
  real(dp) :: fimax,fimin
  real(dp) :: reta,etarng
  real(dp) :: gam, urf, tolAbs, tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 

! Variable specific coefficients:
  gam = TurbModel%Scalar(ifi)%gds
  urf = TurbModel%Scalar(ifi)%urf
  lSolver = TurbModel%Scalar(ifi)%lSolver 
  maxiter = TurbModel%Scalar(ifi)%maxiter
  tolAbs = TurbModel%Scalar(ifi)%tolAbs
  tolRel = TurbModel%Scalar(ifi)%tolRel

  if(ifi.eq.1) then
    prtr=1.0_dp/sigma_k
  else
    prtr=1.0_dp/sigma_epsilon
  endif

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
  !   In 'find_strain_rate' we calculate strain rate as:
  !   S = sqrt (2*Sij*Sij),
  !   for RNG model we need S = sqrt (Sij*Sij).
  !=========================================================

  ! sqrt2r = 1./sqrt(2.0_dp)

  do inp=1,numCells

    ! magStrain(inp)  = magStrain(inp) / sqrt(2.0_dp); when squared ==>>

    gen(inp)=abs(vis(inp)-viscos)*half*magStrain(inp)*magStrain(inp)

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

!****************************************
  elseif(ifi.eq.2) then
!****************************************

  ! Epsilon volume source terms

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================
  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Production of dissipation
    su(inp)=c1*genp*ed(inp)*vol(inp)/(te(inp)+small)

    ! Destruction of dissipation
    sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)

    ! Negative value of production moved to lhs.
    sp(inp)=sp(inp)-c1*genn*vol(inp)/(te(inp)+small) 

    ! Additional term of RNG model
    etarng = magStrain(inp)*te(inp)/(ed(inp)+small)
    reta = cmu*etarng**3*(1.0_dp-etarng/4.38_dp)/(1.0_dp+0.012_dp*etarng**3)
    if (etarng.le.4.38_dp) then
      sp(inp)=sp(inp)+reta*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
    else
      su(inp)=su(inp)-reta*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
    end if

    !
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
!--------------------------------------
  end if

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
                         dfidxi, dcoef, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*fi(ijb) + suadd

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

        if ( ifi .eq. 1 ) then
        !
        ! > Wall boundary conditions for turbulence kinetic energy eq.
        !

          viss=viscos
          if(ypl(iWall).gt.ctrans) viss=visw(iWall)

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

        else
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

        endif

      enddo

    endif
    
  enddo ! Boundary conditions loop 

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
    call csrsolve(lSolver, ed, su, resor(6), maxiter, tolAbs, tolRel, 'epsilon' )
  endif


  ! Update field values at boundaries
  call updateBoundary( fi )


! Report range of scalar values and clip if negative
  fimin = minval(fi)
  fimax = maxval(fi)
  
  if ( ifi.eq.1 ) then 
    write(6,'(2x,es11.4,a,es11.4)') fimin,' <= k <= ',fimax
  else
    write(6,'(2x,es11.4,a,es11.4)') fimin,' <= epsilon <= ',fimax
  endif

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


!***********************************************************************
!
subroutine modify_mu_eff
!
! Update turbulent and effective viscosity.
!
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables

  implicit none
!
!***********************************************************************
!
  integer :: i,ib,inp
  integer :: iface, ijp,ijb,iwall
  real(dp) :: urf,visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau,viscw

  urf = TurbModel%urfVis

  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    ! Update effective viscosity:
    vis(inp)=viscos + den(inp)*cmu*te(inp)**2/(ed(inp)+small)

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
        ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
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



end subroutine modify_mu_eff


!***********************************************************************
!
subroutine modify_viscosity_inlet_k_epsilon_rng()
!
! Update turbulent and effective viscosity at inlet.
!
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
!
  integer :: i,ib,ijb

  !
  ! Loop over inlet boundaries
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


end module k_epsilon_rng