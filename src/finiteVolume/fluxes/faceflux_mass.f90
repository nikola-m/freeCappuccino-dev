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
  use variables!, only: den,U,V,W,dUdxi,dVdxi,dWdxi,p,dpdxi
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
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: fluxmass

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: are
  ! real(dp) :: dpn
  real(dp) :: xpn,ypn,zpn,dene
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ui,vi,wi
  ! real(dp) :: ue,ve,we
  real(dp) :: smdpn
  ! real(dp) :: sfdpnr
  real(dp) :: xpp,ypp,zpp,xep,yep,zep
  real(dp) :: dpe
  ! real(dp) :: dpx,dpy,dpz
  real(dp) :: dpxi,dpyi,dpzi
  real(dp) :: Dpu!,Dpv,Dpw


  ! > Geometry:

  ! Face interpolation factor
  fxn = lambda 
  fxp = 1.0_dp-lambda

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)

  ! Distance between cell centers
  ! dpn = sqrt(xpn**2+ypn**2+zpn**2)

  ! cell face area
  are = sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx = arx/are
  nyy = ary/are
  nzz = arz/are

  !  ______
  ! (Vol/Ap)_f
  Dpu = (fxn*Vol(ijn)*Apu(ijn)+fxp*Vol(ijp)*Apu(ijp))
  ! Dpv = (fxn*Vol(ijn)*Apv(ijn)+fxp*Vol(ijp)*Apv(ijp))
  ! Dpw = (fxn*Vol(ijn)*Apw(ijn)+fxp*Vol(ijp)*Apw(ijp))

  ! Density at the cell face
  dene = den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE-CORRECTION EQUATION

  ! (Sf.Sf) / (dpn.Sf)  
  !# smdpn = (arx*arx+ary*ary+arz*arz)/(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz+small)
  smdpn = (arx*arx+ary*ary+arz*arz)/(arx*xpn+ary*ypn+arz*zpn)

  ! |Sf| / (dpn.nf)
  !smdpn = are/(xpn*nxx+ypn*nyy+zpn*nzz)

  cap = -dene*Dpu*smdpn
  can = cap



!////////////////////////////////////////////////////////////
!   RHIE-CHOW velocity interolation at face    
!         __     ______________      ______
!   Uf = (U)f + (dPdxi*(vol/ap))f - (vol/ap)f*(Pn'-Pp')/dp'n'
!   
!   Fluxmass = Densit*dot(Uf,Sf)
!////////////////////////////////////////////////////////////


  ! UI-> (U)f -> second order interpolation at face
  !+Interpolate velocities to face center:

  ! ui = face_value_cds( ijp, ijn, lambda, u )
  ! vi = face_value_cds( ijp, ijn, lambda, v )
  ! wi = face_value_cds( ijp, ijn, lambda, w )
  
  ! ui = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, u, dUdxi )
  ! vi = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, v, dVdxi )
  ! wi = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, w, dWdxi )

  ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )

  ! Try this - like in Fluent Theory Guide
  ! ui = ( u(ijp)/Apu(ijp)+u(ijn)/Apu(ijn) ) / ( 1./Apu(ijp) + 1./Apu(ijn) )
  ! vi = ( v(ijp)/Apv(ijp)+v(ijn)/Apv(ijn) ) / ( 1./Apv(ijp) + 1./Apv(ijn) )
  ! wi = ( w(ijp)/Apw(ijp)+w(ijn)/Apw(ijn) ) / ( 1./Apw(ijp) + 1./Apw(ijn) )

  !+Interpolate pressure gradients to cell face center
  dpxi = ( dPdxi(1,ijn)*fxp + dPdxi(1,ijp)*fxn ) * xpn
  dpyi = ( dPdxi(2,ijn)*fxp + dPdxi(2,ijp)*fxn ) * ypn
  dpzi = ( dPdxi(3,ijn)*fxp + dPdxi(3,ijp)*fxn ) * zpn

  !+Face centered pressure gradient along cell-face normal

  ! Values at points p' and e' due to non-orthogonality. 
  xpp = xf-(xf-xc(ijp))*nxx
  ypp = yf-(yf-yc(ijp))*nyy
  zpp = zf-(zf-zc(ijp))*nzz
  xep = xf-(xf-xc(ijn))*nxx
  yep = yf-(yf-yc(ijn))*nyy
  zep = zf-(zf-zc(ijn))*nzz

  ! Distances |P'P| and |E'E| projected ionto x,y,z-axis
  xpp = xpp-xc(ijp)
  ypp = ypp-yc(ijp)
  zpp = zpp-zc(ijp)
  xep = xep-xc(ijn)
  yep = yep-yc(ijn)
  zep = zep-zc(ijn)

  dpe = p(ijn) + dPdxi(1,ijn)*xep+dPdxi(2,ijn)*yep+dPdxi(3,ijn)*zep - &
        p(ijp) + dPdxi(1,ijp)*xpp+dPdxi(2,ijp)*ypp+dPdxi(3,ijp)*zpp 

  ! Pressure gradient along normal between N' and P' point which are on face normal direction.
  ! dpn = dpe/(xpn*nxx+ypn*nyy+zpn*nzz+small)

  ! Rhie-Chow Interpolation 
  ! ue = ui - Dpu * (dpn - dpxi)
  ! ve = vi - Dpv * (dpn - dpyi)
  ! we = wi - Dpw * (dpn - dpzi)

  ! MASS FLUX via Rhie-Chow Interpolation
  ! fluxmass = dene*(ue*arx+ve*ary+we*arz)
  fluxmass = dene*(ui*arx+vi*ary+wi*arz) + cap*(dpe-dpxi-dpyi-dpzi)  

end subroutine


!***********************************************************************
!
subroutine facefluxmass2(ijp, ijn, arx, ary, arz, lambda, cap, can, fluxmass)
!
!***********************************************************************
!
! You will have to call it like this:
! subroutine facefluxmass2(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, fluxmass)
! with additional arguments xf, yf, zf, compare to the case above,
! if you are using
! ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
! ...
! to find velocity interpolation value at face.
! Now we use weighted interpolation as described in Fluent Theory Guide
! for testing purposes.
!
!***********************************************************************
  implicit none

  integer, intent(in) :: ijp, ijn
  ! real(dp), intent(in) :: xf,yf,zf
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
  ! cap = -dene*Kj*are/dpn
  !...or..(maybe better)
  cap = -dene*Kj*(arx*arx+ary*ary+arz*arz)/(arx*xpn+ary*ypn+arz*zpn) 
  can = cap


!
! CELL FACE PRESSURE GRADIENTS AND VELOCITIES
!

!////////////////////////////////////////////////////////
!     RHIE-CHOW velocity interolation at face
!     Fluxmass = Densit*dot(Uf,Sf)
!/////////////////////////////////////////////////////////


  ! UI-> (U)f -> second order interpolation at face
  ! ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  ! vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  ! wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )

  ! ! Try this - like in Fluent Theory Guide
  ui = ( u(ijp)/Apu(ijp)+u(ijn)/Apu(ijn) ) / ( 1./Apu(ijp) + 1./Apu(ijn) )
  vi = ( v(ijp)/Apv(ijp)+v(ijn)/Apv(ijn) ) / ( 1./Apv(ijp) + 1./Apv(ijn) )
  wi = ( w(ijp)/Apw(ijp)+w(ijn)/Apw(ijn) ) / ( 1./Apw(ijp) + 1./Apw(ijn) )

  dpxi = ( dPdxi(1,ijn)*fxp + dPdxi(1,ijp)*fxn ) * xpn
  dpyi = ( dPdxi(2,ijn)*fxp + dPdxi(2,ijp)*fxn ) * ypn
  dpzi = ( dPdxi(3,ijn)*fxp + dPdxi(3,ijp)*fxn ) * zpn

  ! MASS FLUX via Rhie-Chow Interpolation
  fluxmass = dene*(ui*arx+vi*ary+wi*arz) + cap*(p(ijn)-p(ijp)-dpxi-dpyi-dpzi)

end subroutine



!***********************************************************************
!
subroutine facefluxmass_piso(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, flmass)!, fmo, fmoo,fmooo)
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
  ! real(dp), intent(in) :: fmo,fmoo,fmooo

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: are,dpn
  real(dp) :: nx,ny,nz
  real(dp) :: xpn,ypn,zpn,dene
  real(dp) :: ui,vi,wi
  ! real(dp) :: ufo, vfo, wfo
  real(dp) :: Kj


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

  ! Face unit normal
  nx = arx/are
  ny = ary/are
  nz = arz/are

  ! density at the cell face
  dene=den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE EQUATION
  Kj = vol(ijp)*apu(ijp)*fxp + vol(ijn)*apu(ijn)*fxn
  ! cap = - dene*Kj*are/dpn
  cap = -dene*Kj*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)
  can = cap

  ! Interpolate velocities (H/Ap) to face center:

  ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )

  ! Try this - like in Fluent Theory Guide
  ! ui = ( u(ijp)/Apu(ijp)+u(ijn)/Apu(ijn) ) / ( 1./Apu(ijp) + 1./Apu(ijn) )
  ! vi = ( v(ijp)/Apv(ijp)+v(ijn)/Apv(ijn) ) / ( 1./Apv(ijp) + 1./Apv(ijn) )
  ! wi = ( w(ijp)/Apw(ijp)+w(ijn)/Apw(ijn) ) / ( 1./Apw(ijp) + 1./Apw(ijn) )

  ! MASS FLUX
  ! Calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
  flmass = dene*(ui*arx+vi*ary+wi*arz) 

  !
  ! Create the difference of face velocity derived from previous time-step mass flows, satisfying  continuity equation, i.e. flmaso/(den*are) and flmassoo(den*are),
  ! and projection to face normal of face interpolated value of velocity from previous timesteps, i.e. dot_protuct((Uo)f, nf) , and dot_product((Uoo)f, nf),
  ! and correct mass flow trough face with it.
  ! This is needed FOR TIME-STEPPING CONSISTENT INTERPOLATION.
  !

  ! den*Vol/(ap*dt)
  Kj = dene*Kj/timestep

  ! if (bdf) then
    
  !   ! From timestep n-1
  !   ufo = face_value_cds( ijp, ijn, lambda, uo ) 
  !   vfo = face_value_cds( ijp, ijn, lambda, vo ) 
  !   wfo = face_value_cds( ijp, ijn, lambda, wo )  

  !   ui = Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nx )
  !   vi = Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * ny )
  !   wi = Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nz )

  !   flmass = flmass + (ui*arx+vi*ary+wi*arz)

  ! endif

  ! if (bdf2) then

  !   ! From timestep n-1
  !   ufo = face_value_cds( ijp, ijn, lambda, uo ) 
  !   vfo = face_value_cds( ijp, ijn, lambda, vo ) 
  !   wfo = face_value_cds( ijp, ijn, lambda, wo )  

  !   ui = 2*Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nx )
  !   vi = 2*Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * ny )
  !   wi = 2*Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nz )

  !   flmass = flmass + (ui*arx+vi*ary+wi*arz)

  !   ! From timestep n-2
  !   ufo = face_value_cds( ijp, ijn, lambda, uoo ) 
  !   vfo = face_value_cds( ijp, ijn, lambda, voo ) 
  !   wfo = face_value_cds( ijp, ijn, lambda, woo )  

  !   ui = -0.5_dp * Kj * ( ( fmoo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nx )
  !   vi = -0.5_dp * Kj * ( ( fmoo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * ny )
  !   wi = -0.5_dp * Kj * ( ( fmoo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nz )

  !   flmass = flmass + (ui*arx+vi*ary+wi*arz)

  ! endif 


  ! if (bdf3) then

  !   ! From timestep n-1
  !   ufo = face_value_cds( ijp, ijn, lambda, uo ) 
  !   vfo = face_value_cds( ijp, ijn, lambda, vo ) 
  !   wfo = face_value_cds( ijp, ijn, lambda, wo )  

  !   ui = 3*Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nx )
  !   vi = 3*Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * ny )
  !   wi = 3*Kj * ( ( fmo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nz )

  !   flmass = flmass + (ui*arx+vi*ary+wi*arz)

  !   ! From timestep n-2
  !   ufo = face_value_cds( ijp, ijn, lambda, uoo ) 
  !   vfo = face_value_cds( ijp, ijn, lambda, voo ) 
  !   wfo = face_value_cds( ijp, ijn, lambda, woo )  

  !   ui = -1.5_dp * Kj * ( ( fmoo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nx )
  !   vi = -1.5_dp * Kj * ( ( fmoo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * ny )
  !   wi = -1.5_dp * Kj * ( ( fmoo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nz )

  !   flmass = flmass + (ui*arx+vi*ary+wi*arz)


  !   ! From timestep n-3
  !   ufo = face_value_cds( ijp, ijn, lambda, uooo ) 
  !   vfo = face_value_cds( ijp, ijn, lambda, vooo ) 
  !   wfo = face_value_cds( ijp, ijn, lambda, wooo )  

  !   ui = 1./3.0_dp * Kj * ( ( fmooo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nx )
  !   vi = 1./3.0_dp * Kj * ( ( fmooo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * ny )
  !   wi = 1./3.0_dp * Kj * ( ( fmooo / (dene*are) - (ufo*nx+vfo*ny+wfo*nz) ) * nz )

  !   flmass = flmass + (ui*arx+vi*ary+wi*arz)

  ! endif

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

  ! (den*Vol/ap)f x Sf/(dpn.nf)         
  rapr = (apu(ijp)*den(ijp)*vol(ijp)*fxp+apu(ijn)*den(ijn)*vol(ijn)*fxn) * are/(xpn*nxx+ypn*nyy+zpn*nzz)


  ! Mass flux correction for the second p'-equation (source term)
  fmcor = rapr*((dPdxi(1,ijn)*xep-dPdxi(1,ijp)*xpp)+&
                (dPdxi(2,ijn)*yep-dPdxi(2,ijp)*ypp)+&
                (dPdxi(3,ijn)*zep-dPdxi(3,ijp)*zpp))

end subroutine


end module