module interpolation
!
! Module which defines function for interpolation of variables to cell faces.
! Various approaches are implemented
!
  use types
  use parameters
  use geometry, only: numTotal,numCells,numPCells,xc,yc,zc

  implicit none

  public

  contains



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
  real(dp), dimension(3,numPCells) :: dUdxi


  if (lcds) then 
    ue = face_value_cds(ijp,ijn, lambda, u)  

  elseif (lcdsc) then 
    ue = face_value_cds_corrected(ijp, ijn, xf, yf, zf, lambda, u, dUdxi)  

  elseif (lcds_flnt) then 
    ue = face_value_central(ijp, ijn, xf, yf, zf, u, dUdxi)

  elseif (l2nd_flnt) then 
    ue = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)

  elseif (lmuscl_flnt) then
    ue = face_value_muscl(ijp, ijn, xf, yf, zf, u, dUdxi)

  elseif (lcubic) then
    ue = face_value_cubic(ijp, ijn, lambda, u, dUdxi)

  elseif (flux_limiter) then
    ue = face_value_2nd_upwind_flux_limiter(ijp, ijn, lambda, u, dUdxi)

  else
    ue = face_value_muscl(ijp, ijn, xf, yf, zf, u, dUdxi)
 

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
!     at neighbours cell-centers.
!
!***********************************************************************

  implicit none

  !     Result
  real(dp) :: face_value

  !      Input
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf, lambda
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numPCells) :: dfidxi

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
!    at neighbours cell-centers.
!
!***********************************************************************

  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: inp, inn
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numPCells) :: gradfi

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
!    at neighbours cell-centers.
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
  real(dp), dimension(3,numPCells) :: gradfi

  ! Locals
  real(dp) :: gradfidr

  gradfidr = gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp))

  face_value = fi(inp) + gradfidr
  end function


!***********************************************************************
!
  function face_value_muscl(inp,inn, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!    Calculates face value using values of variables and their gradients
!    at neighbours cell-centers.
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
  real(dp), dimension(3,numPCells) :: gradfi

  !     Locals
  real(dp) ::  phi_p, phi_n
  real(dp) :: xcp,ycp,zcp
  real(dp) :: xcn,ycn,zcn
  real(dp) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
  real(dp) :: gradfi_n_x,gradfi_n_y,gradfi_n_z
  real(dp) :: gradfidr_2nd_upwind,gradfidr_central,face_value_2nd_upwind,face_value_central
  real(dp) :: theta

  ! theta = 1/8
  theta = 0.125_dp

  !.....Values at cell center's of neighbouring cells:
  phi_p = fi(inp)

  phi_n = fi(inn)

  xcp = xc(inp)
  ycp = yc(inp)
  zcp = zc(inp)

  xcn = xc(inn)
  ycn = yc(inn)
  zcn = zc(inn)

  gradfi_p_x = gradfi(1,inp)
  gradfi_p_y = gradfi(2,inp)
  gradfi_p_z = gradfi(3,inp)

  gradfi_n_x = gradfi(1,inn)
  gradfi_n_y = gradfi(2,inn)
  gradfi_n_z = gradfi(3,inn)


  ! gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
  gradfidr_2nd_upwind=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) 
  gradfidr_central=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
                  +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

  face_value_2nd_upwind = ( phi_p + gradfidr_2nd_upwind )
  face_value_central = 0.5_dp*( phi_p + phi_n + gradfidr_central)

  face_value = theta*face_value_central + (1.0_dp-theta)*face_value_2nd_upwind
  
  end function


!***********************************************************************
!
  function face_value_cubic(inp,inn, lambda, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!    CUBIC interpolation function.
!    Calculates face value using values of variables and their gradients
!    at neighbours cell-centers fitted into a cubic polynomial. 
!    Implementation corresponds to one found tin OpenFOAM, but the details
!    are taken from derivation in an cfd-online thread:
!    https://www.cfd-online.com/Forums/openfoam-programming-development/184766-factors-cubic-inerpolation-scheme-openfoam.html
!
!***********************************************************************
!

  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: inp, inn
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numPCells) :: gradfi

  ! Locals
  real(dp) :: central,corr1,corr2,lambda_,fP,fN,gx,gy,gz,dx,dy,dz,dpn

  ! Face interpolation factor
  lambda_=1.0_dp-lambda
  
  ! The CUBIC interpolation has follwong for according to derivation
  !\phi_f = \phi(1-\lambda) = \lambda \phi_P + (1-\lambda)\phi_N + \lambda (1-\lambda)(1-2\lambda)(\phi_N -\phi_P)  - \lambda^2(\lambda-1) \nabla \phi_P - \lambda (1-\lambda)^2 \nabla \phi_N
  ! Here \lambda is our lambda_=1-lambda.
  !

  central  = lambda_ * fi(inp) + (1.0-lambda_)*fi(inn) 

  corr1    = lambda_*(1-lambda_)*(1-2*lambda_)*(fi(inn)-fi(inp)) 

  fP = lambda_**2*(lambda_-1.0)
  fN = lambda_*(1.0-lambda_)**2 
  gx = fP * gradfi(1,inp) + fN * gradfi(1,inn)
  gy = fP * gradfi(2,inp) + fN * gradfi(2,inn)
  gz = fP * gradfi(3,inp) + fN * gradfi(3,inn)


  ! Project to face normal direction by doting with face area vector and then mulitpy with distance dpn
  !  ->    ->
  ! Grad * Sf/|Sf| *
  !aresq = arx * arx + ary * ary + arz * arz
  !corr2 = ( gx * arx + gy * ary + gz * arz ) / are * dpn  ! Maybe we don't need dpn factor?

  dx = xc(inn) - xc(inp)
  dy = yc(inn) - yc(inp)
  dz = zc(inn) - zc(inp)
  dpn = sqrt(dx*dx + dy*dy +dz*dz)

  corr2 = (gx*dx + gy*dy + gz*dz)/dpn

  ! Finally cubic scheme, central scheme + high order corrections
  face_value = central + corr1 - corr2

  end function


!***********************************************************************
!
  function face_value_2nd_upwind_flux_limiter(ijp, ijn, lambda, u, dUdxi) result(face_value)
!
!***********************************************************************
!
!    Calculates face value using values of variables and their gradients
!    at neighbours cell-centers.
!    Flux limited versionof BOUNDED HIGH-ORDER CONVECTIVE SCHEMES.
!    Reference paper is Waterson & Deconinck JCP 224 (2007) pp. 182-207
!    Also Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003., 
!    for definition of 'r' ratio expression on unstructured meshes.
!
!***********************************************************************
!
  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: ijn, ijp
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: u
  real(dp), dimension(3,numPCells) :: dUdxi

  ! Locals
  real(dp) :: r,psi,xpn,ypn,zpn,fxp

  ! Face interpolation factor
  fxp = 1.0_dp-lambda

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)

  ! Gradient ratio expression taken from Darwish-Moukalled TVD schemes paper
  r = (2*dUdxi(1,ijp)*xpn + 2*dUdxi(2,ijp)*ypn + 2*dUdxi(3,ijp)*zpn)/(u(ijn)-u(ijp)) - 1.0_dp


  if(lsmart) then
    psi = max(0., min(2.*r, 0.75*r+0.25, 4.))

  elseif(lavl) then
    psi = max(0., min(1.5*r, 0.75*r+0.25, 2.5)) 

  elseif(lmuscl) then
    psi = max(0., min(2.*r, 0.5*r+0.5, 2.))
 
  elseif(lumist) then
    psi = max(0., min(2.*r, 0.75*r+0.25, 0.25*r+0.75, 2.))

  elseif(lkoren) then
    psi = max(0., min(2.*r, 2./3._dp*r+1./3._dp, 2.))

  elseif(lcharm) then
    psi = max(0., min(2*r,(r+abs(r))*(3*r+1.0)/(2*(r+1.0)**2), 2.0))

  elseif(lospre) then
    psi = max(0., min(2*r,1.5*r*(r+1.0)/(r**2+r+1.0), 2.0))

  else
    ! psi = 1.0_dp   ! psi for 2nd order upwind (luds) scheme
    psi = max(0., min(r, 4.0)) ! boundedCDS - just temporary here until I implement it
    ! psi = 2./3.*r+1./3. ! Unique third order kappa scheme, kappa=1/3
  end if

  face_value = u(ijp) + fxp*psi*(u(ijn)-u(ijp))

  end function


end module
