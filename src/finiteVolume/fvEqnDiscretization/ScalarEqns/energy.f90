module energy
!
! Purpose:
! Implementation of scalar transport equation for mean (in the sense of Reynolds averaging) 
! 1) Total Enthalpy or
! 2) Total Energy,
! 
! Discussion:
!   The approach taken is that of Murthy and Mathur as documented in papers [1] and [2].
!   In this approach the equations for Total Energy or Total Enthalpy are solved for temperature
!   and on the right hand side, as source terms, we have unsteady and convective terms of difference
!   between the total energy written in terms of temperature and explicit total energy,
!   or between the total enthalpy written in terms of temperature and explicit total enthaply.
!   Look into the papers [1] and [2] for clarification.
!   As I understand this is how it is implemented in Fluent.
!   Also we consulted the web page where Total Energy equation is discussed and some of the ideas are 
!   also present here.
!
!
!  References:
!    [1] J.Y. Murthy and S.R. Mathur, A Conservative Numerical Scheme for the Energy Equation, 
!        J. of Heat Transfer (Nov. 1998) Vol. 120, p. 1081.
!    [2] S.R. Mathur and J.Y. Murthy, Unstructured Finite Volume Methods for Multi-Mode Heat Transfer,
!        Chapter Two in Advances of Numerical Heat Transfer, Vol.2 pp.37-70.
!    [3] 'Energy Equation in OpenFOAM', https://doc.cfd.direct/openfoam/energy-equation/
!
!
!  Author: 
!    Nikola Mirkov (largeddysimulation@gmail.com)
!
!
  use types

  implicit none

  !
  ! Discrtetization and solution parameters - modified trough input.nml file
  !
  logical  :: calcEn = .False.                      ! To activate the solution of this field in the main function. 
  real(dp) :: urfEn  = 0.9_dp                       ! Under-relaxation factors.
  real(dp) :: gdsEn = 1.0_dp                        ! Deferred correction factor.
  character ( len=30 ) :: cSchemeEn = 'linearUpwind'! Convection scheme - default is second order upwind.
  character ( len=12 ) :: dSchemeEn = 'skewness'    ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
  integer :: nrelaxEn = 0                           ! Type of non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed. 1/0/-1
  character ( len=12 ) :: lSolverEn = 'bicgstab'    ! Linear algebraic solver.
  integer  :: maxiterEn = 10                        ! Max number of iterations in linear solver.
  real(dp) :: tolAbsEn = 1e-13                      ! Absolute residual level.
  real(dp) :: tolRelEn = 0.025                      ! Relative drop in residual to exit linear solver.
  real(dp) :: sigtEn = 0.85_dp                       ! sigma_t
  logical :: solveTotalEnergy = .True.              !@ What we solve here - default is total energy, change it 
  logical :: solveInternalEnergy = .False.          !@ in input.nml file.
  logical :: solveEnthalpy = .False.                !@
  logical :: addViscDiss = .False.                  ! Add viscous dissipation term? T/F

  real(dp), parameter :: c_p = 1006.0               ! Heat capacity at constant pressure [J/kgK] ! From input?


  ! type, extends(fieldEqn) :: energyEqn
  !   real(dp) :: sigtEn = 0.9_dp                       ! sigma_t
  !   logical :: solveTotalEnergy = .True.              !@ What we solve here - default is total energy, change it 
  !   logical :: solveInternalEnergy = .False.          !@ in input.nml file.
  !   logical :: solveEnthalpy = .False.                !@
  !   logical :: addViscDiss = .False.                  ! Add viscous dissipation term? T/F
  ! end type energyEqn

  ! type(energyEqn) :: Energy


  private 

  public :: calc_energy

  public :: calcEn, urfEn, gdsEn, cSchemeEn, dSchemeEn, nrelaxEn, lSolverEn, maxiterEn, tolAbsEn, tolRelEn, sigtEn, & ! params
            solveEnthalpy,solveInternalEnergy,solveTotalEnergy, &
            addViscDiss
contains



subroutine calc_energy
!
! Purpose:
!   Ansemble and solve total energy/internal energy/enthalpy equation.
!
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use fieldManipulation, only: explDiv,fieldInterpolate
  use linear_solvers
  
  implicit none

!
! Local variables
!
  integer ::  i, k, ib, inp, ijp, ijn, ijb, iface, iwall
  integer :: iPer, l, if, iftwin
  real(dp) :: gam, prtr, apotime
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: coef,dcoef
  real(dp) :: fimax,fimin
  real(dp) :: Ke, Keo, Keoo
  ! real(dp) :: dpdt
  real(dp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
  real(dp) :: Psi, Phi
  real(dp) :: are,nxf,nyf,nzf,dTn
  real(dp) :: urfr,urfm
  real(dp) :: viscDiss
  real(dp), dimension(:), allocatable :: DivRhoUK,DivUp
  real(dp), dimension(:), allocatable :: up,vp,wp
  real(dp), dimension(:), allocatable :: deni

! Variable specific coefficients:
  gam = gdsEn
  prtr = 1.0_dp/sigtEn

! Calculate gradient: 
  call grad(En,dEndxi)
  call grad(T,dTdxi)

  ! Velocity gradients update (is it necessary here?) 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! > Volume source terms - vectors
!

  ! We need to add material derivative of mechnical energy D/Dt(rho*K), K=1/2|u|**2
  ! For the sake of efficiency, explicit divergence is done as field fieldManipulation
  ! operation, and explicit time derivative is given below inside loop over cells.
  ! Finally we get DdtRhoK + DivRhoUK
  if ( solveEnthalpy .or. solveInternalEnergy ) then

    ! Face density - allocate array
    allocate( deni(numTotal), DivRhoUK(numCells) )

    ! Interpolate density to faces
    deni = fieldInterpolate( den )

    ! Explicit divergence of mechanical energy
    DivRhoUK = explDiv( flmass/deni, 0.5*U**2, 0.5*V**2, 0.5*W**2 ) !<- K define

    deallocate( deni )

  endif

  if ( solveTotalEnergy .or. solveInternalEnergy ) then

    ! Prebaci gore, a ove neka idu kao private promenljive modula, da ih ne bi stalno alocirao
    allocate( up(numTotal), vp(numTotal), wp(numTotal), DivUp(numCells) )

    !The thing below is a new vector field now, and we're going to find its divergence, explicitly.
    ! Component-wise multiplication first:
    up = u*p 
    vp = v*p
    wp = w*p

    !                                       _
    ! This holds a vector of sources of div(u*p)
    !
    DivUp = explDiv(up,vp,wp)

    deallocate(up,vp,wp)

  endif


    ! ! > Convective term correction
    ! ! NAPOMENA!!:
    ! !  Za ovaj modul uvek stavljamo da je upwind schema 1st upwind i tako cemo racunati
    ! !  i konvekciju ispod za Temp, ali eksplicitna konvekcija za En (totalna energija ili totalna entalpija)
    ! !  mora da bude high order.
    ! !  Ako ljudi stave drugaciju upwind schemu overridovacemo u 1st Upwind. Da pojednostavimo.
    ! su = su + explDiv( flmass, (c_p*T - En) ) ! PROVERITI DA LI OVAKO MOZE



!
! > Volume source terms - loop over cells
!
  do inp=1,numCells


    !
    ! > Viscous heating term (compressible flows)
    !
    if (addViscDiss) then

      dudx = dUdxi(1,inp)
      dudy = dUdxi(2,inp)
      dudz = dUdxi(3,inp)
      
      dvdx = dVdxi(1,inp)
      dvdy = dVdxi(2,inp)
      dvdz = dVdxi(3,inp)

      dwdx = dWdxi(1,inp)
      dwdy = dWdxi(2,inp)
      dwdz = dWdxi(3,inp)

      Psi = (dudx+dvdy+dwdz)**2

      Phi = 2*(dudx**2+dvdy**2+dwdz**2)+ (dudy+dvdx)**2+(dudz+dwdx)**2+(dwdy+dvdz)**2
      
      viscDiss = vis(inp)*(Phi-2./3.*Psi)*Vol(inp)

    endif


    
    !
    ! > Source terms caused by writing the energy variable (total energy or total enthaply) in terms of temperature.
    !

    ! > Unsteady term correction
    if (ltransient) then

      apotime = den(inp)*vol(inp)/timestep

      if( bdf .or. cn ) then
        
        su(inp) = su(inp) + apotime*( c_p*( T(inp) - To(inp) ) - ( En(inp)-Eno(inp) ) )

      elseif( bdf2 ) then

        su(inp) = su(inp) + apotime*( &
                                      ( c_p*( 1.5_dp*T(inp) - 2*To(inp) + 0.5_dp*Too(inp) ) ) - &
                                      ( 1.5_dp*En(inp) - 2*Eno(inp) + 0.5_dp*Enoo(inp) ) &
                                    )
      endif

    endif


    !
    ! > Source terms
    !

    if (ltransient) then

      if ( solveEnthalpy .or. solveInternalEnergy ) then
        ! > ddt(rho*K)

        if( bdf .or. cn ) then
          Ke  = 0.5_dp*(u(inp)**2+v(inp)**2+w(inp)**2)        ! 1/2 * |u|^2
          Keo = 0.5_dp*(uo(inp)**2+vo(inp)**2+wo(inp)**2)     ! 1/2 * |uo|^2

          su(inp) = su(inp) + den(inp)*vol(inp) * ( Ke - Keo )

        elseif( bdf2 ) then
          Ke   = 0.5_dp*(u(inp)**2+v(inp)**2+w(inp)**2)       ! 1/2 * |u|^2
          Keo  = 0.5_dp*(uo(inp)**2+vo(inp)**2+wo(inp)**2)    ! 1/2 * |uo|^2
          Keoo = 0.5_dp*(uoo(inp)**2+voo(inp)**2+woo(inp)**2) ! 1/2 * |uoo|^2

          su(inp) = su(inp) + den(inp)*vol(inp) * ( 1.5_dp*Ke - 2*Keo + 0.5_dp*Keoo )

        endif
      endif

      if ( solveEnthalpy ) then
        ! > dpdt - the time derivative of pressure - find explicitly as source

        if( bdf .or. cn ) then
          su(inp) = su(inp) + vol(inp) * ( p(inp) - po(inp) )

        elseif( bdf2 ) then
          su(inp) = su(inp) + vol(inp) * ( 1.5_dp*p(inp) - 2*po(inp) + 0.5_dp*poo(inp) )

        endif
      endif

    endif ! Unsteady


    ! Unsteady Term for Temperature - we are writting the equation in terms of Temperature
    if (ltransient) then

      if( bdf .or. cn ) then
        apotime = den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*To(inp)
        sp(inp) = sp(inp) + apotime

      elseif( bdf2 ) then
        apotime=den(inp)*vol(inp)/timestep
        su(inp) = su(inp) + apotime*( 2*To(inp) - 0.5_dp*Too(inp) )
        sp(inp) = sp(inp) + 1.5_dp*apotime

      endif

    endif


  enddo


!
! > Calculate terms integrated over faces
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, &
                     T, dTdxi, prtr, cap, can, suadd ) ! Promenio u Temp iz En

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
    ! sp(ijp) = sp(ijp) - can

    ! ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap
    !sp(ijn) = sp(ijn) - cap

    ! > Explicit part of convection and difussion fluxes

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
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         T, dTdxi, prtr, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*T(ijb) + suadd

      end do


    elseif ( bctype(ib) == 'wall' .and. bcDFraction(7,ib)==1.0_dp ) then

      ! Isothermal wall boundaries (that's Dirichlet on temperature)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        dcoef = viscos/pranl+(vis(ijp)-viscos)/sigtEn
        coef=dcoef*srdw(iWall)
        !# are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
        !# coef = dcoef*are/dnw(iWall)
        sp(ijp) = sp(ijp) + coef
        su(ijp) = su(ijp) + coef*T(ijb)

      enddo

    elseif ( bctype(ib) == 'wall' .and. bcDFraction(7,ib)==0.0_dp ) then

      ! Adiabatic wall boundaries (that's actually zero grad on temperature)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        T(ijb)=T(ijp)

      enddo

    elseif ( bctype(ib) == 'wallQFlux') then

      ! Prescribed heat flux wall boundaries (that's actually non homo Neumann on temperature)

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        dcoef = viscos/pranl+(vis(ijp)-viscos)/sigtEn
        coef=dcoef*srdw(iWall)

        ! Value of the temprature gradient in normal direction is set trough 
        ! proper choice of component values. Let's project in to normal direction
        ! to recover normal gradient.

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Gradient in face-normal direction        
        dTn = dTdxi(1,ijb)*nxf+dTdxi(2,ijb)*nyf+dTdxi(3,ijb)*nzf

        ! Explicit source
        su(ijp) = su(ijp) + coef*dTn

      enddo

    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! count periodic boundary pairs

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i ! Where does the face data for twin start, looking at periodic boundary pair with index iPer.
        ijn = owner(iftwin)              ! Owner cell of the twin periodic face


        ! face flux scalar but for periodic boundaries - it will be recognized by arguments
        call facefluxsc_periodic( ijp, ijn, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), &
                                  flmass(if), gam, &
                                  T, dTdxi, prtr, cap, can, suadd )


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
          su(ijp) = su(ijp) - a(k)*To(ijn)

          k = jcell_icell_csr_index(i)
          su(ijn) = su(ijn) - a(k)*To(ijp)
      enddo
      
      do ijp=1,numCells
          apotime=den(ijp)*vol(ijp)/timestep
          off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) - a(diag(ijp))
          su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*To(ijp)
          sp(ijp) = sp(ijp)+apotime
      enddo

  endif

  ! Underrelaxation factors
  urfr = 1.0/urfEn
  urfm = 1.0_dp-urfEn

  ! Main diagonal term assembly:
  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! Sum all off-diagonal coefs in a row of a sparse matrix.
    ! we substract the diagonal element from the sum, to eliminate it from the sum.
    off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp))
    a(diag(inp)) = sp(inp) - off_diagonal_terms

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfr
    su(inp) = su(inp) + urfm*a(diag(inp))*T(inp)

  enddo

  ! Solve linear system:
  call csrsolve(lSolverEn, T, su, resor(7), maxiterEn, tolAbsEn, tolRelEn, 'Energy')


  ! Update field values at boundaries
  call updateBoundary(T)


! Report range of scalar values and clip if negative
  fimin = minval( T(1:numCells) )
  fimax = maxval( T(1:numCells) )
  
  write(6,'(2x,es11.4,a,es11.4)') fimin,' <= Temperature <= ',fimax

  ! These field values cannot be negative
  if(fimin.lt.0.0_dp) T(1:numCells) = max( T(1:numCells), small )  


  !
  ! > Recalculate energy variable (total energy, ethalpy or internal energy)
  !
  call correct_energy


end subroutine


!***********************************************************************
!
subroutine facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                      fm, lambda, gam, FI, dFidxi, &
                      prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters, only: zero, viscos
  use geometry, only: numTotal, numCells, xc,yc,zc
  use variables, only: vis
  use interpolation
  use gradients

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn ! owner and neighbour indices
  real(dp), intent(in) :: xf,yf,zf ! face centroid coordinates
  real(dp), intent(in) :: arx, ary, arz ! area vector
  real(dp), intent(in) :: fm ! mass flow at face
  real(dp), intent(in) :: lambda ! interpoaltion factor
  real(dp), intent(in) :: gam  ! deferred correction factor [0,1], 1-high order.
  real(dp), dimension(numTotal), intent(in) :: Fi ! Scalar field in question
  real(dp), dimension(3,numCells), intent(in) :: dFidxi ! Gradient of the scalar field in question.
  real(dp), intent(in) :: prtr ! One over prandtl coefficient
  real(dp), intent(inout) :: cap, can, suadd ! On return - matrix coeffs for owner and neighbour and source contribution.


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn!, xi,yi,zi,r1,r2,psie,psiw
  real(dp) :: dpn
  real(dp) :: cp,cn
  real(dp) :: fii
  real(dp) :: fdfie,fdfii,fcfie,fcfii,ffic
  real(dp) :: de, game, viste
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
  viste = (vis(ijp)-viscos)*fxp+(vis(ijn)-viscos)*fxn
  game = viscos*prtr+viste/sigtEn


  ! Difusion coefficient for linear system
  ! de = game*are/dpn
  de = game*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)

  ! Convection fluxes - uds
  cn = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  cap = -de - cp
  can = -de + cn
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              Fi, dFidxi, nrelaxEn, dSchemeEn, dfixi, dfiyi, dfizi, &
              dfixii, dfiyii, dfizii)
  

  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)  

  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)


  !-------------------------------------------------------
  ! Explicit higher order convection
  !-------------------------------------------------------
  if( fm .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dFidxi, cSchemeEn )
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dFidxi, cSchemeEn )
  endif

  fcfie = fm*fii

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = cn*fi(ijn)+cp*fi(ijp)

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
subroutine facefluxsc_periodic(ijp, ijn, &
                               xf, yf, zf, arx, ary, arz, &
                               flmass, gam, &
                               FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use variables, only: vis
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
  real(dp) :: de, game, viste
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
  viste = (vis(ijp)-viscos)*fxp+(vis(ijn)-viscos)*fxn
  game = viscos*prtr+viste/sigtEn


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
  use types
  use parameters, only: zero,viscos
  use geometry, only: numTotal, numCells, xc,yc,zc
  use variables, only: vis

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
  real(dp) :: cp,cn
  real(dp) :: fdfie,fdfii
  real(dp) :: d1x,d1y,d1z,d2x,d2y,d2z
  real(dp) :: de, vole, game, viste
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

  ! Turbulent viscosity
  viste = vis(ijn)-viscos

  ! Cell face diffussion coefficient for TEMPERATURE
  game = viscos*prtr + viste/sigtEn


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
  cn = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  cap = -de - cp
  can = -de + cn

  ! Explicit part of fluxes
  suadd = fdfie-fdfii 

end subroutine


subroutine correct_energy
!
! Purpose:
!  Recalculate energy variable (total energy, ethalpy or internal energy).
!
! Author:
!  Nikola Mirkov
!
  ! use variables, only: u,v,w,p,T,En,den

  implicit none

  if ( solveEnthalpy ) then

    ! En = 

  else if ( solveInternalEnergy ) then

    ! En = 

  elseif ( solveTotalEnergy .or. solveInternalEnergy ) then

    ! En = 

  else
    write(*,*) "Error: choose one energy variable!"
    stop
  endif


end subroutine


end module
