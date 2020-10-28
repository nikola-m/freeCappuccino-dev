module interpolation
!
! Module which defines function for interpolation of variables to cell faces.
! Various approaches are implemented
!
  use types
  use geometry, only: numTotal,numCells,xc,yc,zc

  implicit none

  ! Choosing discretization scheme applied for all variables at the moment.
  logical :: lcds
  logical :: lcdsc
  logical :: lkappa ! Hardcoded kappa scheme as in Ansys Fluent
  logical :: lsmart
  logical :: lavl
  logical :: lmuscl
  logical :: lumist
  logical :: lkoren
  logical :: lcharm
  logical :: lvanleer
  logical :: lospre
  logical :: lminmod
  logical :: lcentral
  logical :: llinearUpwind
  logical :: lcdsBounded
  logical :: lludsBounded
  logical :: lludsBounded02
  logical :: lluds, lfromm, lcui, lquick ! Various kappa schemes
  logical :: lspl_13, lspl_max_12, lspl_max_13 ! Various symmetric piecewise linear schemes

  logical :: flux_limiter

  public

  private :: lcds, lcdsc, lluds, lkappa, lminmod, lsmart, lavl, lmuscl, lumist, lkoren, lcharm, lvanleer, &
             lospre, lcentral, llinearUpwind, lcdsBounded, lludsBounded, lludsBounded02, &
             lfromm, lcui, lquick,lspl_13, lspl_max_12, lspl_max_13
  private :: flux_limiter

  contains


!***********************************************************************
!
subroutine set_convective_scheme(name)
!
! Set convective scheme according to given name.
!
!***********************************************************************
!
  implicit none

  character(len=25), intent(in) :: name

  ! Set initial false for everything
  lcds = .false.
  lcdsc = .false.
  lluds = .false.
  lkappa = .false.
  lsmart = .false.
  lavl = .false.
  lmuscl = .false.
  lumist = .false.
  lkoren = .false.
  lcharm = .false.
  lospre = .false.
  lminmod = .false.
  lcentral = .false.
  llinearUpwind = .false.
  lcdsBounded = .false.
  lludsBounded = .false.

  flux_limiter = .false.

  select case( trim(name) ) 

    case ('cds')
      lcds = .true.
    case ('cdscorr')
      lcdsc = .true.
    case ('luds')
      lluds = .true.
    case('kappa')
      lkappa = .true.
    case ('smart')
      lsmart = .true.
    case ('avl-smart')
      lavl = .true.
    case ('muscl')
      lmuscl = .true.
    case ('umist')
      lumist = .true.
    case ('koren')
      lkoren = .true.
    case ('charm')
      lcharm = .true.
    case ('ospre')
      lospre = .true.
    case ('minmod')
      lminmod = .true.
    case ('central')
      lcentral = .true.
    case ('linearUpwind')
      llinearUpwind = .true.
    case ('boundedLinearUpwind')
      lludsBounded = .true.
    case ('boundedLinearUpwind02')
      lludsBounded02 = .true.
    case ('boundedCentral')
      lcdsBounded = .true.
    case default
      write(*,'(a)') "Using default convective scheme - 2nd order upwind."
      llinearUpwind = .true.

  end select
  
  ! Set value for flux_limiter logical
  if(lluds.or.lsmart.or.lavl.or.lmuscl.or.lumist.or.lkoren.or.lcharm.or.lospre.or. &
     lminmod.or.lludsBounded.or.lludsBounded02.or.lcdsBounded) then
    flux_limiter = .true.
  endif

end subroutine



!***********************************************************************
!
function face_value(ijp,ijn,xf,yf,zf,lambda,u,dUdxi) result(ue)
!
!***********************************************************************
!
  implicit none

  ! Result
  real(dp) :: ue

  ! Input
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf,lambda
  real(dp), dimension(numTotal) :: u
  real(dp), dimension(3,numCells) :: dUdxi


  if (lcds) then 
    ue = face_value_cds(ijp,ijn, lambda, u)  

  elseif (lcdsc) then 
    ue = face_value_cds_corrected(ijp, ijn, xf, yf, zf, lambda, u, dUdxi)  

  elseif (lcentral) then 
    ue = face_value_central(ijp, ijn, xf, yf, zf, u, dUdxi)

  elseif (llinearUpwind) then 
    ue = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)

  elseif (lkappa) then
    ue = face_value_kappa(ijp, ijn, xf, yf, zf, u, dUdxi)

  elseif (flux_limiter) then
    ue = face_value_flux_limiter(ijp, ijn, u, dUdxi)

  else
    ue = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)
 

  endif 

end function

!***********************************************************************
!
function face_value_w_option(ijp,ijn,xf,yf,zf,lambda,u,dUdxi,scheme) result(ue)
!
!***********************************************************************
!
  implicit none

  ! Result
  real(dp) :: ue

  ! Input
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf,lambda
  real(dp), dimension(numTotal) :: u
  real(dp), dimension(3,numCells) :: dUdxi
  character(len=10), intent(in) :: scheme


  if (scheme == 'cds') then 
    ue = face_value_cds(ijp,ijn, lambda, u)  

  elseif (scheme == 'cdscorr') then 
    ue = face_value_cds_corrected(ijp, ijn, xf, yf, zf, lambda, u, dUdxi)  

  elseif (scheme == 'central') then 
    ue = face_value_central(ijp, ijn, xf, yf, zf, u, dUdxi)

  elseif (scheme == 'linearUpwind') then 
    ue = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)

  elseif (scheme == 'kappa') then
    ue = face_value_kappa(ijp, ijn, xf, yf, zf, u, dUdxi)

  else
    ue = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi) 

  endif 

end function

!***********************************************************************
!
  function face_value_cds(ijp, ijn, lambda, fi) result(face_value)
!
!***********************************************************************
!
! Calculates face value using values of variables at neighbour cell-centers.
!
!***********************************************************************

  implicit none

  !     Result
  real(dp) :: face_value

  !      Input
  integer :: ijp, ijn
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: fi

  !     Locals
  real(dp) :: fxn,fxp

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  !            |________uj'___________|
  face_value = fi(ijp)*fxp+fi(ijn)*fxn

  end function



!***********************************************************************
!
  function face_value_cds_corrected(ijp,ijn, xf, yf, zf, lambda, fi, dfidxi) result(face_value)
!
!***********************************************************************
!
!     Calculates face value using values of variables and their gradients
!     at centers of adjecent cells..
!
!***********************************************************************

  implicit none

  !     Result
  real(dp) :: face_value

  !      Input
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf, lambda
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

  !     Locals
  real(dp) :: fxn,fxp,xi,yi,zi,dfixi,dfiyi,dfizi

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Coordinates of point j'
  xi = xc(ijp)*fxp+xc(ijn)*fxn
  yi = yc(ijp)*fxp+yc(ijn)*fxn
  zi = zc(ijp)*fxp+zc(ijn)*fxn

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

  !            |________uj'___________|_________________ucorr____________________|
  face_value = fi(ijp)*fxp+fi(ijn)*fxn+(dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi))

  end function



!***********************************************************************
!
  function face_value_central(inp,inn, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells..
!
!***********************************************************************

  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: inp, inn
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: gradfi

  ! Locals
  real(dp) :: gradfidr

  gradfidr=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) &
          +gradfi(1,inn)*(xf-xc(inn))+gradfi(2,inn)*(yf-yc(inn))+gradfi(3,inn)*(zf-zc(inn))

  face_value = 0.5_dp*( fi(inp) + fi(inn) + gradfidr)

  end function


!***********************************************************************
!
  function face_value_2nd_upwind(inp, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells..
!    Corresponds to unlimited second order upwind scheme as 
!    used in ANSYS FLUENT.
!
!***********************************************************************
!

  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: inp
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: gradfi

  ! Locals
  real(dp) :: gradfidr

  gradfidr = gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp))

  face_value = fi(inp) + gradfidr

  end function


!***********************************************************************
!
  function face_value_kappa(inp,inn, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells..
!    Corresponds to MUSCL scheme as used in ANSYS FLUENT.
!
!***********************************************************************
!

  implicit none

  !.....Result
  real(dp) :: face_value

  !     Input
  integer :: inp, inn
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: gradfi

  !     Locals
  real(dp) :: gradfidr_2nd_upwind,gradfidr_central,face_value_2nd_upwind,face_value_central
  real(dp) :: theta

  ! theta = 0.125_dp ! Fluent theta = 1/8
  ! theta = 2./3.    ! CUI
  theta = 0.5_dp   ! Fromm


  gradfidr_2nd_upwind=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) 

  gradfidr_central=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) &
                  +gradfi(1,inn)*(xf-xc(inn))+gradfi(2,inn)*(yf-yc(inn))+gradfi(3,inn)*(zf-zc(inn))

  face_value_2nd_upwind = ( fi(inp) + gradfidr_2nd_upwind )
  face_value_central = 0.5_dp*( fi(inp) + fi(inn) + gradfidr_central)

  face_value = theta*face_value_central + (1.0_dp-theta)*face_value_2nd_upwind
  
  end function


!***********************************************************************
!
  function face_value_flux_limiter(ijp, ijn, u, dUdxi) result(face_value)
!
!***********************************************************************
!
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells.
!    Flux limited versionof BOUNDED HIGH-ORDER CONVECTIVE SCHEMES.
!    Reference paper is Waterson & Deconinck JCP 224 (2007) pp. 182-207
!    Also Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003., 
!    for definition of 'r' ratio expression on unstructured meshes.
!
!***********************************************************************
!
  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: ijn, ijp
  real(dp), dimension(numTotal) :: u
  real(dp), dimension(3,numCells) :: dUdxi

  ! Locals
  real(dp) :: r,psi,xpn,ypn,zpn

  psi = 1.0_dp ! Initial

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)


  ! Gradient ratio expression taken from Darwish-Moukalled 'TVD schemes for unstructured grids' paper.
  r = (2*dUdxi(1,ijp)*xpn + 2*dUdxi(2,ijp)*ypn + 2*dUdxi(3,ijp)*zpn)/(u(ijn)-u(ijp)+1e-30) - 1.0_dp


  if(lsmart) then
    psi = max(0., min(2*r, 0.75*r+0.25, 4.0))

  elseif(lavl) then
    psi = max(0., min(1.5*r, 0.75*r+0.25, 2.5))

  elseif(lmuscl) then
    psi = max(0., min(2*r, 0.5*r+0.5, 2.0))
 
  elseif(lumist) then
    psi = max(0., min(2*r, 0.75*r+0.25, 0.25*r+0.75, 2.0))

  elseif(lkoren) then
    psi = max(0., min(2*r, 2./3._dp*r+1./3.0_dp, 2.0))

  elseif(lcharm) then
    psi = max(0.,(r+abs(r))*(3*r+1.0)/(2*(r+1.0)**2)) ! Charm

  elseif(lvanleer) then
    psi = max(0., min( (r+abs(r))/(r+1.0), 2.0 )) ! Harmonic - Van Leer

  elseif(lospre) then
    psi = max(0.,3*r*(r+1.0)/(2*(r**2+r+1.0) ))

  elseif(lminmod) then
    psi = max(0., min(r, 1.0))

  elseif(lludsBounded) then
    psi = max(0., min(2*r, 1.0)) ! <-See Waterson-Deconninck JCP paper how to get this from Chakravarty-Osher

  elseif(lludsBounded02) then
    psi = max(0., min(10*r, 1.0)) ! <- this is the same as limitedLinear 0.2 in OpenFOAM.

  elseif(lcdsBounded) then
    psi = max(0., min(r, 2.0)) ! <-See Waterson-Deconninck JCP paper how to get this from Chakravarty-Osher

  elseif(lluds) then
    psi = 1.0_dp ! unbounded 2nd order upwind (luds) scheme, also a kappa scheme with kappa = -1

  elseif(lfromm) then ! Fromm scheme - unique symmetric kappa scheme, kappa=0
    psi = 0.5_dp*r+0.5_dp

  elseif(lcui) then ! Unique third order kappa scheme, kappa=1/3
    psi = 2./3.*r+1./3.

  elseif(lquick) then ! Quick scheme, as unique quadratic kappa scheme, kappa=1/2
    psi = 3./4.*r + 1./4.

  elseif(lspl_13) then ! Symmetric piecewise linear scheme, see Waterson and Deconinck JCP 224 (2007)
    psi = max(0., min(2*r, 1./3.*r+2./3., 2./3.*r+1./3., 2.0))

  ! else(lspl_max_12) then
  !   psi = 

  ! else(lspl_max_13) then
  !   psi = 

  ! else(lsuperbee) then
  !   psi = 

  ! else(lcopla) then
  !   psi = 

  ! else(lcubista) then
  !   psi = 

  ! else(lwaceb) then
  !   psi = 

  ! else(lstoic) then
  !   psi = 

  ! else(lvonos) then
  !   psi = 

  ! else(lhoab) then
  !   psi = 

  ! else(lvanalbada) then
  !   psi = 

  ! else() then
  !   psi = 
  ! else() then
  !   psi = 

  ! else() then
  !   psi = 

  ! else() then
  !   psi = 

  ! else() then
  !   psi = 

  ! else() then
  !   psi = 

  ! else() then
  !   psi = 

  end if

  face_value = u(ijp) + 0.5_dp*psi*( u(ijn) - u(ijp) ) 


  end function


end module
