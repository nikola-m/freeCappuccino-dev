module k_epsilon_zeta_f
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of the k-epsilon-zeta-f2004 turbulence model.
!
!  Discussion:
!
!  Reference:
!    [1] K. Hanjalic, M. Popovac, M. Hadziabdic, A Robust Near-Wall Elliptic-Relaxation Eddy-Viscosity Turbulence Model for CFD,
!        Int. J. Heat and Fluid Flow, Vol. 25, 2004, pp.1047-1051.
!    [2] M. Popovac, K. Hanjalic, Compound Wall Treatment for RANS Computation of Complex Turbulent Flows and Heat Transfer,
!        Flow, Turbulcence and Combustion, Vol. 78, 2007, pp.177-202.
!    [3] Langley Research Center-Turbulence Modeling Resource, turbmodels.larc.nasa.gov
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    17th February 2023
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
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
  real(dp), parameter :: CMU = 0.22_dp
  real(dp), parameter :: SIXCMU = sqrt(6.0)*0.22_dp  
  real(dp), parameter :: C1 = 1.4_dp
  real(dp)            :: CEPS1 ! = 1.4*(1.+(0.012/zeta))
  real(dp), parameter :: CEPS2 = 1.9_dp
  real(dp), parameter :: C2P = 0.65_dp
  real(dp), parameter :: CT = 6.0_dp
  real(dp), parameter :: CL = 0.36_dp
  real(dp), parameter :: CETHA = 85.0_dp

  real(dp), parameter :: sigma_k = 1.0_dp
  real(dp), parameter :: sigma_epsilon = 1.3_dp
  real(dp), parameter :: sigma_zeta = 1.2_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(cmu))
  real(dp), parameter :: CMU75 = cmu25**3

  real(dp), dimension(:), allocatable :: frlx ! eliptic relaxation scalar field
  real(dp), dimension(:), allocatable :: zeta, zetao, zetaoo ! zeta scalar variable, zeta = v^2/k
  real(dp), dimension(:,:), allocatable :: dZetadxi


  private 


  public :: modify_viscosity_k_epsilon_zeta_f
  public :: modify_viscosity_inlet_k_epsilon_zeta_f

contains



subroutine modify_viscosity_k_epsilon_zeta_f
!
! Purpose: 
!   Main module routine to solve turbulence model equations and update effective viscosity.
!
!
  implicit none


  if(.not.allocated(frlx)) then

    allocate(frlx(numTotal))

    write(*,'(a)') ' '
    write(*,'(a)') '  **Allocated eliptic relaxation function for k-epsilon-zeta-f2004 model.'
    write(*,'(a)') ' '
   
  endif

  if(.not.allocated(zeta)) then

    allocate(zeta(numTotal), zetao(numTotal), zetaoo(numTotal), dZetadxi(3,numTotal))

    write(*,'(a)') ' '
    write(*,'(a)') '  **Allocated zeta field for k-epsilon-zeta-f2004 model.'
    write(*,'(a)') ' '
   
  endif

  call calcsc_tke     ! Assemble and solve turbulence kinetic energy eq.
  call calcsc_epsilon ! Assemble and solve dissipation rate of tke eq.
  call calcsc_zeta    ! Assemble and solve zeta equation.
  call calcsc_frelax  ! Assemble and solve elliptic relaxation eqn.
  call modify_mu_eff  ! Update efective viscosity

end subroutine


subroutine calcsc_frelax
!
! Purpose: 
!   Assemble and solve eliptic relaxation equation.
!
!
  use sparse_matrix
  use linear_solvers

  implicit none

  integer :: i,inp,iwall,iper,l,ijp,ijb,iface,ib
  real(dp) :: nu,lenscale,timescale,dummy
  character( len=12 ) :: lSolver 
  integer :: maxiter
  real(dp) :: tolAbs, tolRel

  do inp=1,numCells

    nu = viscos/densit

    timescale = max(   min( te(inp)/(ed(inp)+small), 0.6_dp/( SIXCMU*magStrain(inp)*zeta(inp)+small ) ), &
                       CT*sqrt( nu/(ed(inp)+small) ) )

    ! Source term
    su(inp) = -( C1-1.0_dp+C2P*gen(inp)/( densit*ed(inp)+small ) ) * ( zeta(inp)-twothirds ) / timescale

    ! Coefficient array for Laplacian - we will use existing array sv.

    lenscale = CL*max( &
                        min( te(inp)**0.75_dp/(ed(inp)+small) , sqrt(te(inp))/( SIXCMU*magStrain(inp)*zeta(inp)+small ) ), &
                        CETHA*sqrt(sqrt( nu**3/(ed(inp)+small) )) &
                      )
    sv(inp) = lenscale*lenscale

  enddo

  ! Initialize solution
  frlx = 0.0_dp
     

  ! Laplacian operator and BCs         
  call laplacian(sv,frlx) 

  ! Modify to Helmholtz eq.
  a(diag(1:numCells)) = a(diag(1:numCells))+1.0_dp



  !
  ! Boundary conditions
  !

  iWall = 0
  iPer = 0
  l = numInnerFaces

  do ib=1,numBoundaries


    if ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        a( ia(ijp):ia(ijp+1)-1 ) = 0.0_dp
        sp(ijp) = 1.0_dp

        ! first wall cell fixed value of epsilon here:
        frlx(ijp) = -2*viscos/densit*zeta(ijp)/(dnw(iWall)*dnw(iWall)) 
        su(ijp) = frlx(ijp)

      enddo

    endif
    
  enddo ! Boundary conditions loop  


! Variable specific coefficients.
! (Note: frelax is forth scalar in this model)
  lSolver = TurbModel%Scalar(4)%lSolver 
  maxiter = TurbModel%Scalar(4)%maxiter
  tolAbs = TurbModel%Scalar(4)%tolAbs
  tolRel = TurbModel%Scalar(4)%tolRel

  ! Solve system
  call csrsolve(lSolver, frlx, su, dummy, maxiter, tolAbs, tolRel, 'frlx' )

end subroutine



subroutine calcsc_tke
!
! Purpose:
!  Assemble and solve turbulence kinetic energy equation.
!
  use sparse_matrix
  use linear_solvers
  
  implicit none

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, apotime, urfrs, urfms
  real(dp) :: genp, genn
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,vtp,xtp,ytp,ztp,ut2
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


  ! CALCULATE STANDARD PRODUCTION
  gen(1:numCells) = abs(vis(1:numCells)-viscos)*magStrain(1:numCells)*magStrain(1:numCells)


  !
  ! > Add Volume Source Terms:
  !

  ! ! Maybe vectorised like this:

  ! su(1:numCells) = max(gen(1:numCells),zero)*vol(1:numCells) &
  ! + den(1:numCells)*vol(1:numCells)/timestep*teo(1:numCells)

  ! sp(1:numCells) = ( ( ed(1:numCells)*den(1:numCells)-min(gen(1:numCells),zero) )/(te(1:numCells)+small) &
  !                    + den(1:numCells)/timestep )*vol(1:numCells)

  ! elseif( bdf2 ) then

  !   su(1:numCells) = su(1:numCells) + apotime*( 2*teo(1:numCells) - half*teoo(1:numCells) )
  !   sp(1:numCells) = sp(1:numCells) + 1.5_dp*apotime

  ! endif



  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Add production term to the rhs:
    su(inp)=genp*vol(inp)  

    ! Add destruction term to the lhs:
    sp(inp)=ed(inp)*den(inp)*vol(inp)/(te(inp)+small)
    sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)



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
  ! > Calculate Convective And Diffusive Fluxes
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
                     te, dtedxi, dcoef, cap, can, suadd )

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
                         dtedxi, dcoef, cap, can, suadd ) 

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
                          te, dtedxi, dcoef, cap, can, suadd )


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
        ! First subtract the standard production from source term
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

    ! Main diagonal term assembly:
    a(diag(inp)) = sp(inp) 
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

end subroutine calcsc_tke



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
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, apotime, urfrs, urfms
  real(dp) :: genp, genn
  real(dp) :: viste,dcoef
  real(dp) :: cap, can, suadd
  real(dp) :: timescale
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


    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Timescale
    timescale = max(   min( te(inp)/(ed(inp)+small), 0.6_dp/(SIXCMU*magStrain(inp)*zeta(inp)) ), &
                       CT*sqrt( viscos/densit/(ed(inp)+small) ) )

    CEPS1  = 1.4*( 1. + ( 0.012/(zeta(inp)+small) ) )

    ! Production of dissipation
    su(inp) = CEPS1*genp*vol(inp)/(timescale+small)

    ! Destruction of dissipation
    sp(inp) = CEPS2*den(inp)*ed(inp)*vol(inp)/(timescale+small)

    ! Negative value of production moved to lhs.
    sp(inp) = sp(inp)-CEPS1*genn*vol(inp)/(ed(inp)*timescale+small) 


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
                     ED, dEDdxi, dcoef, cap, can, suadd )

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
                         dEDdxi, dcoef, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*ED(ijb) + suadd

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

        ! For correct values of dissipation all coefficients have
        ! to be zero, su equal the dissipation, and diagonal element a(diag(ijp)) = 1

        a( ia(ijp):ia(ijp+1)-1 ) = 0.0_dp
        sp(ijp) = 1.0_dp

        ! first wall cell fixed value of epsilon here:
        ed(ijp) = 2*viscos/densit*te(ijp)/(dnw(iWall)*dnw(iWall)) 
        su(ijp) = ed(ijp)

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



subroutine calcsc_zeta
!
! Purpose:
!  Assemble and solve equation for zeta scalar field (zeta = v^2/k).
!
  use sparse_matrix
  use linear_solvers
  
  implicit none

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, ijb, ib, iface, iwall
  integer :: iper, l, if, iftwin
  real(dp) :: prtr, apotime, urfrs, urfms
  real(dp) :: dummy
  real(dp) :: viste,dcoef
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: fimax,fimin
  real(dp) :: gam, urf, tolAbs, tolRel
  integer :: maxiter
  character( len=12 ) :: lSolver 

! Variable specific coefficients:
  gam = TurbModel%Scalar(3)%gds
  urf = TurbModel%Scalar(3)%urf
  lSolver = TurbModel%Scalar(3)%lSolver 
  maxiter = TurbModel%Scalar(3)%maxiter
  tolAbs = TurbModel%Scalar(3)%tolAbs
  tolRel = TurbModel%Scalar(3)%tolRel

  prtr=1.0_dp/sigma_zeta


  ! Calculate gradient: 
  call grad(zeta,dZetadxi)
 
  ! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp


  !
  ! > Add Volume Source Terms:
  !

  do inp=1,numCells

    ! Add production term to the rhs:
    su(inp) = den(inp)*frlx(inp)*vol(inp)  

    ! Add destruction term to the lhs, divide by solution variable (zeta) and change sign to plus:
    sp(inp) = gen(inp)*vol(inp)/(te(inp)+small) 


    ! UNSTEADY TERM
    if (ltransient) then

      apotime = den(inp)*vol(inp)/timestep
      
      if( bdf .or. cn ) then

        su(inp) = su(inp) + apotime*zetao(inp)
        sp(inp) = sp(inp) + apotime

      elseif( bdf2 ) then

        su(inp) = su(inp) + apotime*( 2*zetao(inp) - 0.5_dp*zetaoo(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime

      endif

    endif

  ! End of volume source terms
  enddo


  !
  ! > Calculate Convective And Diffusive Fluxes
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
                     TurbModel%Scalar(3)%cScheme, &
                     Zeta, dZetadxi, dcoef, cap, can, suadd )

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
                         dZetadxi, dcoef, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*zeta(ijb) + suadd

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
                          zeta, dZetadxi, dcoef, cap, can, suadd )


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

        ! For correct values of dissipation all coefficients have
        ! to be zero, su equal the dissipation, and diagonal element a(diag(ijp)) = 1

        a( ia(ijp):ia(ijp+1)-1 ) = 0.0_dp
        sp(ijp) = 1.0_dp

        ! first wall cell fixed value of epsilon here:
        zeta(ijp) = 2*viscos/densit*te(ijp)/(dnw(iWall)*dnw(iWall)) 
        su(ijp) = ed(ijp)

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
      su(ijp) = su(ijp) - a(k)*zetao(ijn)

      k = jcell_icell_csr_index(i)
      su(ijn) = su(ijn) - a(k)*zetao(ijp)
    enddo

    do ijp=1,numCells
      apotime=den(ijp)*vol(ijp)/timestep
      off_diagonal_terms = sum( a( ia(ijp) : ia(ijp+1)-1 ) ) - a(diag(ijp))
      su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*zetao(ijp)
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
    su(inp) = su(inp) + urfms*a(diag(inp))*zeta(inp)
                    
  enddo


  ! Solve linear system:
  call csrsolve(lSolver, zeta, su, dummy, maxiter, tolAbs, tolRel, 'zeta' )


  ! Update field values at boundaries
  call updateBoundary( zeta )


! Report range of scalar values and clip if negative
  fimin = minval( zeta(1:numCells) )
  fimax = maxval( zeta(1:numCells) )
   
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= zeta <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) zeta(1:numCells) = max( zeta(1:numCells),small )

end subroutine



subroutine modify_mu_eff
!
! Purpose: 
!   Update turbulent and effective viscosity.
!

  implicit none

  integer :: i,inp
  integer :: iface,ijp,ijb,ib,iwall
  real(dp) :: urf,visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau,viscw
  real(dp) :: timescale


  urf = TurbModel%urfVis

!
! Loop trough cells 
!
  do inp=1,numCells

    ! Store old value
    visold=vis(inp)

    ! Timescale
    timescale = max(   min( te(inp)/(ed(inp)+small), 0.6_dp/( SIXCMU*magStrain(inp)*zeta(inp)+small ) ), &
                       CT*sqrt( viscos/densit/(ed(inp)+small) ) )
    
    ! Update effective viscosity:
    vis(inp) = viscos + den(inp)*cmu*zeta(inp)*te(inp)*timescale

    ! Underelaxation
    vis(inp) = urf*vis(inp)+(1.0_dp-urf)*visold

  enddo

  !
  ! > Boundary faces 
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

        ! Tau(iWall) = viscos*Ut2/dnw(iWall)
        ! Utau = sqrt( Tau(iWall) / den(ijb) )
        ! ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
        ypl(iWall) = den(ijp)*cmu25*sqrt(te(ijp))*dnw(iWall)/viscos
        ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

        viscw = zero
     
        ! Standard wall function - first cell in log region
        if( ypl(iWall) > ctrans ) then

          viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
          tau(iwall) = cappa*den(ijp)*Vtp*cmu25*sqrt(te(ijp))/log(Elog*ypl(iWall))

        else

          ! Viscous sublayer
          Tau(iWall) = viscos*Ut2/dnw(iWall)
          Utau = sqrt( Tau(iWall) / den(ijb) )
          ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        endif

        visw(iWall) = max(viscos,viscw)

        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

end subroutine modify_mu_eff



subroutine modify_viscosity_inlet_k_epsilon_zeta_f
!
! Update turbulent and effective viscosity at inlet.
!

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


end module k_epsilon_zeta_f