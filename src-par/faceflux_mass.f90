module faceflux_mass
!
! Implementation of common functions for obtaining discretized fluxes for
! Navier-Stokes equation.
! To the outside world we show only the interface function 'facefluxmass', 
! wisely checking the arguments, module decides what function to call.
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,vol
  use variables, only: den,U,V,W,dUdxi,dVdxi,dWdxi,p,dpdxi
  use sparse_matrix, only: apu,apv,apw
  use interpolation

  implicit none


  private 

  public :: facefluxmass,facefluxmass2,facefluxmass_piso,fluxmc


contains


!***********************************************************************
!
subroutine facefluxmass(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, fluxmass)
!
!***********************************************************************
!
  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in):: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: fluxmass

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn,dene
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ui,vi,wi,ue,ve,we
  real(dp) :: smdpn,sfdpnr
  real(dp) :: xpp,ypp,zpp,xep,yep,zep
  real(dp) :: dpe
  real(dp) :: dpex,dpey,dpez
  real(dp) :: dpxi,dpyi,dpzi
  real(dp) :: Dpu,Dpv,Dpw,dpecorr


  ! > Geometry:

  ! Face interpolation factor
  fxn = lambda 
  fxp = 1.0_dp-lambda

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)

  ! Distance between cell centers
  dpn = sqrt(xpn**2+ypn**2+zpn**2)

  ! cell face area
  are = sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx = arx/are
  nyy = ary/are
  nzz = arz/are

  !  ______
  ! (Vol/Ap)_f
  Dpu = (fxn*Vol(ijn)*Apu(ijn)+fxp*Vol(ijp)*Apu(ijp))
  Dpv = (fxn*Vol(ijn)*Apv(ijn)+fxp*Vol(ijp)*Apv(ijp))
  Dpw = (fxn*Vol(ijn)*Apw(ijn)+fxp*Vol(ijp)*Apw(ijp))

  ! Density at the cell face
  dene = den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE-CORRECTION EQUATION
  ! sfdpnr = 1./(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz)
  sfdpnr = 1./(arx*xpn+ary*ypn+arz*zpn)  
  smdpn = (arx*arx+ary*ary+arz*arz)*sfdpnr
  cap = -dene*Dpu*smdpn
  can = cap

!
! CELL FACE PRESSURE GRADIENTS AND VELOCITIES
!

!////////////////////////////////////////////////////////
!     RHIE-CHOW velocity interolation at face
!
!   Uf = UI+DPDXI-API*Sf*(Pn-Pp)
!         _
!   UI-> (U)f -> second order interpolation at face
!            _______________ 
!   DPDXI-> (dPdx*Vol*(1/ap))f -> second order interpolation at cell face
!          ______
!   API*Sf*(Pn-Pp) -> (1/ap)f*Sf*(Pn-Pp) cell face coefficient 1/Ap x Area_f x (p1-p2)
!  
!   Finally:     
!         __     _______________     ______
!   Uf = (U)f + (dPdx*Vol*(1/ap))f - (1/ap)f*Sf*(Pn-Pp)
!
!   and:
!   Fluxmass = Densit*dot(Uf,Sf)
!/////////////////////////////////////////////////////////


  ! UI-> (U)f -> second order interpolation at face
  !+Interpolate velocities to face center:+++++++++++++++++++++++++

  ! ui = face_value_cds( ijp, ijn, lambda, u )
  ! vi = face_value_cds( ijp, ijn, lambda, v )
  ! wi = face_value_cds( ijp, ijn, lambda, w )
  
  ! ui = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, u, dUdxi )
  ! vi = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, v, dVdxi )
  ! wi = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, w, dWdxi )

  ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )

  
  !+END: Interpolate velocities to face center:+++++++++++++++++++++++++


  ! DPDXI-> (dPdx*Vol*(1/ap))f -> second order interpolation at cell face
  !+Interpolate pressure gradients to cell face center++++++++++++++++++
  dpxi = Dpu*(fxn*dPdxi(1,ijn)+fxp*dPdxi(1,ijp))*xpn*nxx
  dpyi = Dpv*(fxn*dPdxi(2,ijn)+fxp*dPdxi(2,ijp))*ypn*nyy
  dpzi = Dpw*(fxn*dPdxi(3,ijn)+fxp*dPdxi(3,ijp))*zpn*nzz
  !+END: Interpolate pressure gradients to cell face center+++++++++++++


  ! (1/ap)f*Sf*(Pn-Pp)
  !+Pressure deriv. along normal+++++++++++++++++++++++++++++++++++++++++ 

  !.....Values at points p' and e' due to non-orthogonality. 
  xpp = xf-(xf-xc(ijp))*nxx
  ypp = yf-(yf-yc(ijp))*nyy
  zpp = zf-(zf-zc(ijp))*nzz
  xep = xf-(xf-xc(ijn))*nxx
  yep = yf-(yf-yc(ijn))*nyy
  zep = zf-(zf-zc(ijn))*nzz
  !.....Distances |P'P| and |E'E| projected ionto x,y,z-axis
  xpp = xpp-xc(ijp)
  ypp = ypp-yc(ijp)
  zpp = zpp-zc(ijp)
  xep = xep-xc(ijn)
  yep = yep-yc(ijn)
  zep = zep-zc(ijn)

  dpe = (p(ijn)-p(ijp))
  dpecorr = ( dPdxi(1,ijn)*xep+dPdxi(2,ijn)*yep+dPdxi(3,ijn)*zep - & !<<--Correction
              dPdxi(1,ijp)*xpp+dPdxi(2,ijp)*ypp+dPdxi(3,ijp)*zpp  )  !<<|
  dpe = dpe+dpecorr

  dpex = Dpu * dpe*sfdpnr * arx
  dpey = Dpv * dpe*sfdpnr * ary
  dpez = Dpw * dpe*sfdpnr * arz
  !+END: Pressure deriv. along normal++++++++++++++++++++++++++++++++++++

  ! Rhie-Chow Interpolation 
  ue = ui - dpex + dpxi
  ve = vi - dpey + dpyi
  we = wi - dpez + dpzi

  ! MASS FLUX via Rhie-Chow Interpolation
  fluxmass = dene*(ue*arx+ve*ary+we*arz)

end subroutine


!***********************************************************************
!
subroutine facefluxmass2(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, fluxmass)
!
!***********************************************************************
!
  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: fluxmass

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn,dene
  real(dp) :: ui,vi,wi
  real(dp) :: dpxi,dpyi,dpzi
  real(dp) :: Kj ! notation from Muzaferija&Gosman JCP paper


  ! > Geometry:

  ! Face interpolation factor
  fxn = lambda 
  fxp = 1.0_dp-lambda

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)

  ! Distance between cell centers
  dpn = sqrt(xpn**2+ypn**2+zpn**2)

  ! cell face area
  are = sqrt(arx**2+ary**2+arz**2)


  ! density at the cell face
  dene = den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE-CORRECTION EQUATION
  Kj = vol(ijp)*apu(ijp)*fxp + vol(ijn)*apu(ijn)*fxn
  cap = -dene*Kj*are/dpn
  can = cap


!
! CELL FACE PRESSURE GRADIENTS AND VELOCITIES
!

!////////////////////////////////////////////////////////
!     RHIE-CHOW velocity interolation at face
!
!   Uf = UI+DPDXI-API*Sf*(Pn-Pp)
!         _
!   UI-> (U)f -> second order interpolation at face
!            _______________ 
!   DPDXI-> (dPdx*Vol*(1/ap))f -> second order interpolation at cell face
!          ______
!   API*Sf*(Pn-Pp) -> (1/ap)f*Sf*(Pn-Pp) cell face coefficient 1/Ap x Area_f x (p1-p2)
!  
!   Finally:     
!         __     _______________     ______
!   Uf = (U)f + (dPdx*Vol*(1/ap))f - (1/ap)f*Sf*(Pn-Pp)
!
!   and:
!   Fluxmass = Densit*dot(Uf,Sf)
!/////////////////////////////////////////////////////////


  ! UI-> (U)f -> second order interpolation at face
  ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )

  ! ui = face_value_cds( ijp, ijn, lambda, u )
  ! vi = face_value_cds( ijp, ijn, lambda, v )
  ! wi = face_value_cds( ijp, ijn, lambda, w )

  dpxi = ( dPdxi(1,ijn)*fxp + dPdxi(1,ijp)*fxn ) * xpn
  dpyi = ( dPdxi(2,ijn)*fxp + dPdxi(2,ijp)*fxn ) * ypn
  dpzi = ( dPdxi(3,ijn)*fxp + dPdxi(3,ijp)*fxn ) * zpn

  ! MASS FLUX via Rhie-Chow Interpolation
  fluxmass = dene*(ui*arx+vi*ary+wi*arz) + cap*(p(ijn)-p(ijp)-dpxi-dpyi-dpzi)

end subroutine


!***********************************************************************
!
subroutine facefluxmass_piso(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, flmass)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,vol
  use sparse_matrix, only: apu
  use variables, only: den,U,V,W,dUdxi,dVdxi,dWdxi
  use gradients
  use interpolation

  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: flmass

  ! Local variables
  real(dp) :: fxn, fxp
  ! real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn,dene
  real(dp) :: ui,vi,wi
  real(dp) :: Kj

  ! > Geometry:

  ! Face interpolation factor
  fxn = lambda 
  fxp = 1.0_dp-lambda

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  ! dpn = sqrt(xpn**2+ypn**2+zpn**2) 

  ! cell face area
  ! are = sqrt(arx**2+ary**2+arz**2)

  ! density at the cell face
  dene = den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE EQUATION
  Kj = vol(ijp)*apu(ijp)*fxp + vol(ijn)*apu(ijn)*fxn
  ! cap = - dene*Kj*are/dpn
  cap = -dene*Kj*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)
  can = cap


  ! Interpolate velocities to face center:
  ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )


  ! MASS FLUX
  ! Calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
  flmass = dene*(ui*arx+vi*ary+wi*arz)

end subroutine



!***********************************************************************
!
subroutine fluxmc(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, fmcor)
!
!***********************************************************************
!
!   This routine calculates mass flux correction in the
!   second pressure-correction step which accounts for the
!   effects of non-orthogonality.
!            ___________                             ->       ->         ->      ->
!   FMCOR = (rho*Vol/Apu)f * |Sf| / |d_{P'N'}| * [ grad(p)_N.r_{NN'} - grad(p)_P.r_{PP'} ]
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,vol
  use variables
  use sparse_matrix, only: apu

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: fmcor
  !
  ! Local variables
  !
  real(dp) :: fxn,fxp
  real(dp) :: xpn,ypn,zpn
  real(dp) :: rapr
  real(dp) :: are
  real(dp) :: nxx,nyy,nzz
  real(dp) :: xep,yep,zep,xpp,ypp,zpp
  real(dp) :: dppnnr

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Distance from P' to N'-reciprocal value
  dppnnr = 1.0_dp/((xpn*nxx)+(ypn*nyy)+(zpn*nzz))

  ! Values at points p' and e' due to non-orthogonality. 
  xpp=xf-(xf-xc(ijp))*nxx
  ypp=yf-(yf-yc(ijp))*nyy
  zpp=zf-(zf-zc(ijp))*nzz

  xep=xf-(xf-xc(ijn))*nxx
  yep=yf-(yf-yc(ijn))*nyy
  zep=zf-(zf-zc(ijn))*nzz

  ! Distances |P'P| and |E'E| projected ionto x,y,z-axis
  xpp=xpp-xc(ijp)
  ypp=ypp-yc(ijp)
  zpp=zpp-zc(ijp)

  xep=xep-xc(ijn)
  yep=yep-yc(ijn)
  zep=zep-zc(ijn)

  ! APU==1./AP x density - interpolated            
  rapr = (apu(ijp)*den(ijp)*vol(ijp)*fxp+apu(ijn)*den(ijn)*vol(ijn)*fxn)


  ! Mass flux correction for the second p'-equation (source term)
  fmcor = rapr*are*((dPdxi(1,ijn)*xep-dPdxi(1,ijp)*xpp)   & 
                   +(dPdxi(2,ijn)*yep-dPdxi(2,ijp)*ypp)   &
                   +(dPdxi(3,ijn)*zep-dPdxi(3,ijp)*zpp))*dppnnr 

end subroutine



end module