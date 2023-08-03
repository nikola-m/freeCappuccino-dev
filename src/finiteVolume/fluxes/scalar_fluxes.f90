module scalar_fluxes
!
! Implementation of functions for obtaining discretized convection and diffusion
! fluxes for transport equations of salar fields.
! To the outside world we show only the interface function 'facefluxsc', 
! wisely checking the arguments, module decides what function to call.
!
  use types
  use parameters
  use geometry, only: numTotal,numCells,xc,yc,zc,Df
  use gradients, only: sngrad
  use interpolation, only: face_value, face_value_cds

  implicit none


  interface facefluxsc
    module procedure facefluxsc
    module procedure facefluxsc_boundary
    module procedure facefluxsc_periodic
  end interface

 
  public :: facefluxsc


contains


!***********************************************************************
!
subroutine facefluxsc(i, ijp, ijn, xf, yf, zf, arx, ary, arz, &
                      fm, lambda, gam, cScheme,  &
                      fi, dfidxi, dcoef, cap, can, suadd)
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: i, ijp, ijn    ! Face ID, owner cell, neighbor cell
  real(dp), intent(in) :: xf,yf,zf      ! Face centroid coordinates
  real(dp), intent(in) :: arx, ary, arz ! Face area vector components
  real(dp), intent(in) :: fm            ! Mass flow
  real(dp), intent(in) :: lambda        ! Face interpolation factor
  real(dp), intent(in) :: gam           ! gamma - deferred correction coef
  character( len=30 ), intent(in) :: cScheme            ! Convection scheme - name
  real(dp), dimension(numTotal), intent(in) :: fi       ! Scalar field in question
  real(dp), dimension(3,numCells), intent(in) :: dfidxi ! Gradient of the scalar field in question
  real(dp), intent(in) :: dcoef                         ! Diffusion coefficient
  real(dp), intent(inout) :: cap, can, suadd            ! Matrix coefs and source term


! Local variables
  real(dp) :: xpn,ypn,zpn
  real(dp) :: Cp,Ce
  real(dp) :: fii
  real(dp) :: fdfie
  real(dp) :: fcfie,fcfii
  real(dp) :: de
  real(dp) :: fxp,fxn
  real(dp) :: dfixi,dfiyi,dfizi
!----------------------------------------------------------------------

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)


  ! Difusion coefficient for linear system
  de = dcoef*Df(i)

  ! Convection fluxes - uds
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  can = -de + ce
  cap = -de - cp


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  ! Explicit diffusion - ONLY nonorthogonal correction to the implicit diffussion:

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dfidxi(1,ijp)*fxp+dfidxi(1,ijn)*fxn
  dfiyi = dfidxi(2,ijp)*fxp+dfidxi(2,ijn)*fxn
  dfizi = dfidxi(3,ijp)*fxp+dfidxi(3,ijn)*fxn

  ! Create non-orthogonal correction to gradient only
  dfixi = dfixi*(arx-Df(i)*xpn)
  dfiyi = dfiyi*(ary-Df(i)*ypn)
  dfizi = dfizi*(arz-Df(i)*zpn)

  fdfie = dcoef*( dfixi + dfiyi + dfizi )



  !-------------------------------------------------------
  ! Explicit higher order convection, cSheme is set in input
  !-------------------------------------------------------
  if( fm .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dfidxi, cScheme)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dfidxi, cScheme)
  endif

  fcfie = fm*fii

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deffered correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  fcfie = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -fcfie+fdfie

end subroutine



!***********************************************************************
!
subroutine facefluxsc_periodic(i, ijp, ijn, xf, yf, zf, arx, ary, arz, &
                               fm, gam, &
                               fi, dfidxi, dcoef, cap, can, suadd)
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: i, ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: dcoef
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: xpn,ypn,zpn
  real(dp) :: Cp,Ce
  real(dp) :: fii
  real(dp) :: fdfie,fcfie,fcfii
  real(dp) :: de
  real(dp) :: fxp,fxn
  real(dp) :: dfixi,dfiyi,dfizi
!----------------------------------------------------------------------

  ! > Geometry:

  ! Face interpolation factor
  fxn = half ! Assumption for periodic boundaries
  fxp = fxn

  ! Distance vector between cell centers
  xpn = 2*( xf-xc(ijp) )
  ypn = 2*( yf-yc(ijp) )
  zpn = 2*( zf-zc(ijp) )


  ! Difusion coefficient for linear system
  de = dcoef * Df(i)

  ! Convection fluxes - uds
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  can = -de + ce
  cap = -de - cp


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  ! Explicit diffussion - ONLY nonorthogonal correction to the implicit diffussion:

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dfidxi(1,ijp)*fxp+dfidxi(1,ijn)*fxn
  dfiyi = dfidxi(2,ijp)*fxp+dfidxi(2,ijn)*fxn
  dfizi = dfidxi(3,ijp)*fxp+dfidxi(3,ijn)*fxn

  ! Create non-orthogonal correction to gradient only
  dfixi = dfixi*(arx-Df(i)*xpn)
  dfiyi = dfiyi*(ary-Df(i)*ypn)
  dfizi = dfizi*(arz-Df(i)*zpn)

  fdfie = dcoef*( dfixi + dfiyi + dfizi )

  !-------------------------------------------------------
  ! Explicit higher order convection, cSheme is set in input usually, but here is set to CDS
  !-------------------------------------------------------
  if( fm .ge. zero ) then 
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
  fcfie = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -fcfie+fdfie

end subroutine


!***********************************************************************
!
subroutine facefluxsc_boundary(ijp, &
                               xf, yf, zf, &
                               arx, ary, arz, &
                               fm, &
                               dfidxi, dcoef, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm 
  real(dp), dimension(3,numTotal), intent(in) :: dfidxi
  real(dp), intent(in) :: dcoef
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: de, Dfi
  real(dp) :: dfixi,dfiyi,dfizi

!----------------------------------------------------------------------

  ! > Geometry:

  ! Distance vector between cell center and face center
  xpn=xf-xc(ijp)
  ypn=yf-yc(ijp)
  zpn=zf-zc(ijp)

  ! (sf*sf)/(dpn * sf)
  Dfi = (arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz) 

  ! Difusion coefficient
  de = dcoef * Dfi

  ! Convection fluxes - uds
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  can = -de + ce
  cap = -de - cp


  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------

  ! NOTE: No need for deferred correction for convection - value of the scalar at boundary face is known.

  ! Explicit diffussion - ONLY nonorthogonal correction to the implicit diffussion:

  ! Interpolate gradients defined at CV centers to faces
  ! dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  ! dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  ! dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn
  ! It should be dfidxi(:,ijn), because fxn=1.0, but we write dfidxi(:,ijp) because constant gradient
  ! is applied between cell center and boundary cell face.
  dfixi = dfidxi(1,ijp)
  dfiyi = dfidxi(2,ijp)
  dfizi = dfidxi(3,ijp) 

  ! Create non-orthogonal correction to gradient only
  dfixi = dfixi*(arx-Dfi*xpn)
  dfiyi = dfiyi*(ary-Dfi*ypn)
  dfizi = dfizi*(arz-Dfi*zpn)

  suadd = dcoef*( dfixi + dfiyi + dfizi )
  

end subroutine




!***********************************************************************
!
subroutine facefluxscalar(i, ijp, ijn, xf, yf, zf, arx, ary, arz, &
                          fm, lambda, gam, cScheme,  &
                          fi, dfidxi, diffcoef, cap, can, suadd)
!
!***********************************************************************
!
! Purpose: In this routine diffusion coefficient is passed as an argument
!          in the form of an array, as opposed in facefluxsc. 
!          I will still figure out what is the best way. Nikola 
!
!
!***********************************************************************
! 

  implicit none


  integer, intent(in) :: i, ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  character( len=30 ), intent(in) :: cScheme 
  real(dp), dimension(numTotal), intent(in) :: fi
  real(dp), dimension(3,numCells), intent(in) :: dfidxi
  real(dp), dimension(numCells), intent(in) :: diffcoef
  real(dp), intent(inout) :: cap, can, suadd


  ! Local variables
  real(dp) :: xpn,ypn,zpn
  real(dp) :: Cp,Ce
  real(dp) :: fii
  real(dp) :: fdfie
  real(dp) :: fcfie,fcfii
  real(dp) :: de, game
  real(dp) :: fxp,fxn
  real(dp) :: dfixi,dfiyi,dfizi
!----------------------------------------------------------------------

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)


  ! Cell face diffussion coefficient
  game = diffcoef(ijp)*fxp+diffcoef(ijn)*fxn

  ! Difusion coefficient for linear system
  de = game * Df(i)


  ! Convection fluxes - uds
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  can = -de + ce
  cap = -de - cp



  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  ! Explicit diffusion - ONLY nonorthogonal correction to the implicit diffussion:

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dfidxi(1,ijp)*fxp+dfidxi(1,ijn)*fxn
  dfiyi = dfidxi(2,ijp)*fxp+dfidxi(2,ijn)*fxn
  dfizi = dfidxi(3,ijp)*fxp+dfidxi(3,ijn)*fxn

  ! Create non-orthogonal correction to gradient only
  dfixi = dfixi*(arx-Df(i)*xpn)
  dfiyi = dfiyi*(ary-Df(i)*ypn)
  dfizi = dfizi*(arz-Df(i)*zpn)

  fdfie = game*( dfixi + dfiyi + dfizi )



  !-------------------------------------------------------
  ! Explicit higher order convection, cSheme is set in input
  !-------------------------------------------------------
  if( fm .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dfidxi, cScheme)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dfidxi, cScheme)
  endif

  fcfie = fm*fii

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deffered correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  fcfie = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -fcfie+fdfie

end subroutine


end module