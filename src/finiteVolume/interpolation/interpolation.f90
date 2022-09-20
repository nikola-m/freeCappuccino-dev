module interpolation

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!     Module which defines functions for interpolation of variables from cell centers to cell faces.
!
!  Discussion:
!
!    We implement various interpolation schemes, see below.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use geometry, only: numTotal,numCells,xc,yc,zc

  implicit none

  public


contains


function face_value(ijp,ijn,xf,yf,zf,lambda,u,dUdxi,scheme) result(vf)
!
!  Purpose: 
!
!   Computes interpolation value at face centroid of the given field 'u'.
!
!  Discussion:
!
!    This function is called within face loop. Therefore function parameters
!    are set within this loop.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    --
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Reference:
!
!    See below important references for flux limiter approach.
!
!  Parameters:
!
!    ijp               - Input, integer, index of the owner cell (of the face).
!    ijn               - Input, integer, index of the neighbor cell (of the face).
!    xf,yf,zf          - Input, double, face centroid coordinates.
!    lambda            - Input, double, face centroid interpolation distance weighted factor.
!    u(numTotal)       - Input, double, scalar field values at cell centers.
!    dUdxi(3,numTotal) - Input, double, gradient vector of the scalar field at cell centers.
!    scheme            - Input, character string, interpolation scheme employed. 
!    vf                - Output, double, intepolated value of the field at face centroid.
!
  implicit none

!
! Result
!
  real(dp) :: vf

!
! Parameters
!
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf,lambda
  real(dp), dimension(numTotal) :: u
  real(dp), dimension(3,numCells) :: dUdxi
  character(len=30), intent(in) :: scheme


  select case( trim(scheme) ) 

    case ( 'cds' )
      vf = face_value_cds(ijp,ijn, lambda, u) 

    ! case ( 'cdscorr' )
    !   vf = face_value_cds_corrected(ijp, ijn, xf, yf, zf, lambda, u, dUdxi)

    case ( 'central' )
      vf = face_value_central(ijp, ijn, xf, yf, zf, u, dUdxi)

    ! case ('harmonic')
    !   vf = face_value_harmonic(ijp,ijn, lambda, u) 

    case ('linearUpwind')
      vf = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)

    case('kappa')
      vf = face_value_kappa(ijp, ijn, xf, yf, zf, u, dUdxi)

    ! case('cubic')
    !   vf = face_value_cubic(ijp, ijn, lambda, u, dUdxi)

    case default 
      ! Defaults to flux limiter schemes to speed-up the search.
      vf = face_value_flux_limiter(ijp, ijn, lambda, u, dUdxi, scheme )

  end select


end function


function face_value_cds(ijp, ijn, lambda, fi) result(vf)
!
!  Purpose: 
!    Calculates face value using central differencing scheme.
!
!  Discussion:
!    Actually the face value is obtained using distance weighted interpolation,
!    which is CDS equivalent.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    --
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Parameters:
!    ijp               - Input, integer, index of the owner cell (of the face).
!    ijn               - Input, integer, index of the neighbor cell (of the face).
!    lambda            - Input, double, face centroid interpolation distance weighted factor.
!    fi(numTotal)      - Input, double, scalar field values at cell centers.
!    vf                - Output, double, intepolated value of the field at face centroid.
!
  implicit none
!
! Result
!
  real(dp) :: vf
!
! Parameters
!
  integer :: ijp, ijn
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: fi
!
! Local variables
!
  real(dp) :: fxn,fxp

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  vf = fi(ijp)*fxp+fi(ijn)*fxn

end function


function face_value_cds_corrected(ijp,ijn, xf, yf, zf, lambda, fi, dfidxi) result(vf)
!
!  Purpose: 
!     Calculates face value using values of variables and their gradients
!     at centers of adjecent cells.
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Reference:
!    Demirdzic & Muzaferija papers.
!
!  Parameters:
!
!    ijp               - Input, integer, index of the owner cell (of the face).
!    ijn               - Input, integer, index of the neighbor cell (of the face).
!    xf,yf,zf          - Input, double, face centroid coordinates.
!    lambda            - Input, double, face centroid interpolation distance weighted factor.
!    fi(numTotal)      - Input, double, scalar field values at cell centers.
!    dfidxi(3,numTotal)- Input, double, gradient vector of the scalar field at cell centers.
!    vf                - Output, double, intepolated value of the field at face centroid.
!

  implicit none

  ! Result
  real(dp) :: vf

  ! Input
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf, lambda
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

  ! Locals
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

  !    |________uj'___________|_________________ucorr____________________|
  vf = fi(ijp)*fxp+fi(ijn)*fxn+(dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi))

end function


function face_value_central(inp,inn, xf, yf, zf, fi, gradfi) result(vf)
!
!  Purpose: 
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells.
!
!  Discussion:
!    This is two sided averaged extrapolation. 
!    This considered cetral scheme in Fluent.
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Reference:
!    ANSYS Fluent Theory guide. Section on numerics.
!
!  Parameters:
!
!    inp               - Input, integer, index of the owner cell (of the face).
!    inn               - Input, integer, index of the neighbor cell (of the face).
!    xf,yf,zf          - Input, double, face centroid coordinates.
!    fi(numTotal)      - Input, double, scalar field values at cell centers.
!    gradfi(3,numTotal)- Input, double, gradient vector of the scalar field at cell centers.
!    vf                - Output, double, intepolated value of the field at face centroid.
!

  implicit none

  ! Result
  real(dp) :: vf

  ! Input
  integer :: inp, inn
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: gradfi

  ! Locals
  real(dp) :: gradfidr

  gradfidr=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) &
          +gradfi(1,inn)*(xf-xc(inn))+gradfi(2,inn)*(yf-yc(inn))+gradfi(3,inn)*(zf-zc(inn))

  vf = 0.5_dp*( fi(inp) + fi(inn) + gradfidr)

end function


function face_value_harmonic(ijp, ijn, lambda, fi) result(vf)
!
!  Purpose: 
!    Calculates face value using HARMONIC differencing scheme.
!
!  Discussion:
!    This sceheme is proposed in S.Patankar's book and recommended
!    for diffusion problems when diffusion coefs have greater spatial variation.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    --
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Parameters:
!    ijp               - Input, integer, index of the owner cell (of the face).
!    ijn               - Input, integer, index of the neighbor cell (of the face).
!    lambda            - Input, double, face centroid interpolation distance weighted factor.
!    fi(numTotal)      - Input, double, scalar field values at cell centers.
!    vf                - Output, double, intepolated value of the field at face centroid.
!
  implicit none
!
! Result
!
  real(dp) :: vf
!
! Parameters
!
  integer :: ijp, ijn
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: fi
!
! Local variables
!
  real(dp) :: fxn,fxp

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

    vf = ( fi(ijp)*fi(ijn) ) / ( fi(ijp)*fxp+fi(ijn)*fxn + 1e-30 ) 

end function

function face_value_2nd_upwind(inp, xf, yf, zf, fi, gradfi) result(vf)
!
!  Purpose:
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells.
!
!  Discussion:
!    Corresponds to unlimited second order upwind scheme as 
!    used in ANSYS FLUENT.
!
!  Reference:
!    ANSYS Fluent Theory guide. Section on numerics.
!
!    inp               - Input, integer, index of the owner cell (of the face).
!    xf,yf,zf          - Input, double, face centroid coordinates.
!    fi(numTotal)      - Input, double, scalar field values at cell centers.
!    gradfi(3,numTotal)- Input, double, gradient vector of the scalar field at cell centers.
!    vf                - Output, double, intepolated value of the field at face centroid.
!

  implicit none

  ! Result
  real(dp) :: vf

  ! Input
  integer :: inp
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: gradfi

  ! Locals
  real(dp) :: gradfidr

  gradfidr = gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp))

  vf = fi(inp) + gradfidr

end function


function face_value_kappa(inp,inn, xf, yf, zf, fi, gradfi) result(vf)
!
!  Purpose:
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells.
!
!  Discussion:
!    Kappa scheme is obtained as a parametrized blend of central and linear upwind schemes.
!    By changing the parameter you get QUCK, Fromm, CUI..
!    
!    It corresponds to MUSCL scheme as used in ANSYS FLUENT if kappa=1/8.
!    I don't know how is that MUSCL scheme?!
!
!  Reference:
!    Waterson & Deconinck JCP 224 (2007) pp. 182-207
!    ANSYS Fluent Theory guide. Section on numerics.
!
!  Parameters:
!    inp               - Input, integer, index of the owner cell (of the face).
!    inn               - Input, integer, index of the neighbor cell (of the face).
!    xf,yf,zf          - Input, double, face centroid coordinates.
!    fi(numTotal)      - Input, double, scalar field values at cell centers.
!    gradfi(3,numTotal)- Input, double, gradient vector of the scalar field at cell centers.
!    vf                - Output, double, intepolated value of the field at face centroid.
!

  implicit none

  ! Result
  real(dp) :: vf

  ! Input
  integer :: inp, inn
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: gradfi

  ! Locals
  real(dp) :: gradfidr_2nd_upwind,gradfidr_central,vf_2nd_upwind,vf_central
  real(dp) :: theta

  ! theta = 0.125_dp ! Fluent theta = 1/8
  theta = 2./3.      ! CUI
  ! theta = 0.5_dp   ! Fromm


  gradfidr_2nd_upwind=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) 

  gradfidr_central=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) &
                  +gradfi(1,inn)*(xf-xc(inn))+gradfi(2,inn)*(yf-yc(inn))+gradfi(3,inn)*(zf-zc(inn))

  vf_2nd_upwind = ( fi(inp) + gradfidr_2nd_upwind )
  vf_central = 0.5_dp*( fi(inp) + fi(inn) + gradfidr_central)

  vf = theta*vf_central + (1.0_dp-theta)*vf_2nd_upwind

  
end function



! function face_value_cubic(inp,inn, lambda, fi, gradfi) result(vf)
! !
! !  Purpose:
! !    Calculates face value using values of variables and their gradients
! !    at centers of adjecent cells.
! !
! !  Discussion:
! !    Cubic scheme from an online reference (see below)
! !
! !  Reference:
! !    https://www.cfd-online.com/Forums/openfoam-programming-development/184766-factors-cubic-inerpolation-scheme-openfoam.html
! !
! !  Parameters:
! !    inp               - Input, integer, index of the owner cell (of the face).
! !    inn               - Input, integer, index of the neighbor cell (of the face).
! !    lambda            - Input, double, face interpolation factor.
! !    fi(numTotal)      - Input, double, scalar field values at cell centers.
! !    gradfi(3,numTotal)- Input, double, gradient vector of the scalar field at cell centers.
! !    vf                - Output, double, intepolated value of the field at face centroid.
! !

!   implicit none

!   ! Result
!   real(dp) :: vf

!   ! Input
!   integer :: inp, inn
!   real(dp) :: lambda
!   real(dp), dimension(numTotal) :: fi
!   real(dp), dimension(3,numCells) :: gradfi

!   ! Local variables
!   real(dp) :: xpn,ypn,zpn,dpn,gradP,gradN

!   ! Face interpolation factor
!   ! fxn=lambda 
!   ! fxp=1.0_dp-lambda

!   lambda = 1.-lambda

!   ! Distance vector between cell centers
!   xpn=xc(inn)-xc(inp)
!   ypn=yc(inn)-yc(inp)
!   zpn=zc(inn)-zc(inp)

!   ! Distance from P to neighbor N
!   dpn=sqrt(xpn**2+ypn**2+zpn**2)  

!   ! Normal on the line connecting cell centers
!   xpn = xpn/dpn
!   ypn = ypn/dpn
!   zpn = zpn/dpn

!   ! Projektuj gradijente u P i u N na pravac dpn da bi dobio skalarne vrednosti i onda to ubaci u jedn.
!   gradP = gradfi(1,inp)*xpn+gradfi(2,inp)*ypn+gradfi(3,inp)*zpn
!   gradN = gradfi(1,inn)*xpn+gradfi(2,inn)*ypn+gradfi(3,inn)*zpn

!   ! Cubic interpolation formula
!   ! vf = fxp*fi(inp) + fxn*fi(inn) + fxp*fxn*(1.-2*fxp)*(fi(inn)-fi(inp)) - fxp**2*fxn*gradP - fxp*fxn*(1+fxp)*gradN

!   vf = lambda*fi(inp)+(1.-lambda)*fi(inn) &
!      +(lambda*(1. - lambda*(3. - 2*lambda)))*(fi(inn)-fi(inp)) &
!      -( sqrt(lambda)*(lambda - 1.)) *gradP &
!      -( sqrt(1. - lambda)*lambda )*gradN
    
! end function


function face_value_flux_limiter(ijp, ijn, lambda, u, dUdxi, scheme) result(vf)
!
!  Purpose: 
!
!   Computes interpolation value at face centroid of the given field 'u'.
!
!  Discussion:
!
!    Calculates face value using values of variables and their gradients
!    at centers of adjecent cells.
!    Flux limited version of BOUNDED HIGH-ORDER CONVECTIVE SCHEMES.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    --
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Reference:
!
!    Reference paper is Waterson & Deconinck JCP 224 (2007) pp. 182-207
!    Also Darwis h& Moukalled, TVD schemes for unstructured grids, IJHMT, 2003., 
!    for definition of 'r' ratio on unstructured meshes.
!
!  Parameters:
!
!    ijp               - Input, integer, index of the owner cell (of the face).
!    ijn               - Input, integer, index of the neighbor cell (of the face).
!    lambda            - Input, double, face centroid interpolation distance weighted factor.
!    u(numTotal)       - Input, double, scalar field values at cell centers.
!    dUdxi(3,numTotal) - Input, double, gradient vector of the scalar field at cell centers.
!    scheme            - Input, character string, interpolation scheme employed.   
!    vf                - Output, double, intepolated value of the field at face centroid.
!
  implicit none

  ! Result
  real(dp) :: vf

  ! Input
  integer :: ijn, ijp
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: u
  real(dp), dimension(3,numCells) :: dUdxi
  character(len=30), intent(in) :: scheme

  ! Locals
  real(dp) :: r,psi,xpn,ypn,zpn,fxp

  ! Face interpolation factor
  fxp = 1.0_dp-lambda

  psi = 1.0_dp ! Initial

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)


  ! Gradient ratio expression taken from Darwish-Moukalled 'TVD schemes for unstructured grids' paper.
  r = (2*dUdxi(1,ijp)*xpn + 2*dUdxi(2,ijp)*ypn + 2*dUdxi(3,ijp)*zpn)/(u(ijn)-u(ijp)+1e-30) - 1.0_dp


  select case( trim(scheme) ) 

  case ( 'muscl' )
    psi = max(0., min(2*r, 0.5*r+0.5, 2.0))

  case ( 'umist' ) 
    psi = max(0., min(2*r, 0.75*r+0.25, 0.25*r+0.75, 2.0))

  case ( 'koren' )
    psi = max(0., min(2*r, 2./3._dp*r+1./3.0_dp, 2.0))

  case ( 'smart' )
    psi = max(0., min(2*r, 0.75*r+0.25, 4.0))

  case ( 'avl-smart' ) 
    psi = max(0., min(1.5*r, 0.75*r+0.25, 2.5))

  case ( 'charm' )
    psi = max(0.,(r+abs(r))*(3*r+1.0)/(2*(r+1.0)**2)) ! Charm

  case ( 'vanleer' ) 
    psi = max(0., min( (r+abs(r))/(r+1.0), 2.0 )) ! Harmonic - Van Leer

  case ( 'ospre' )
    psi = max(0.,3*r*(r+1.0)/(2*(r**2+r+1.0) ))

  case ( 'minmod' ) 
    psi = max(0., min(r, 1.0))

  case ( 'boundedLinearUpwind' )
    psi = max(0., min(2*r, 1.0)) ! <-See Waterson-Deconninck JCP paper how to get this from Chakravarty-Osher

  case ( 'boundedLinearUpwind02' )
    psi = max(0., min(10*r, 1.0)) ! <- this is the same as limitedLinear 0.2 in OpenFOAM.

  case ('boundedCentral' ) 
    psi = max(0., min(r, 4.0)) ! <-See Waterson-Deconninck JCP paper how to get this from Chakravarty-Osher

  case ('linearUpwindFL') 
    psi = 1.0_dp ! unbounded 2nd order upwind (luds) scheme, also a kappa scheme with kappa = -1

  case ( 'fromm' )  ! Fromm scheme - unique symmetric kappa scheme, kappa=0
    psi = 0.5_dp*r+0.5_dp

  case ( 'cui' ) ! Unique third order kappa scheme, kappa=1/3
    psi = 2./3.*r+1./3.

  case ( 'quick' ) ! Quick scheme, as unique quadratic kappa scheme, kappa=1/2
    psi = 3./4.*r + 1./4.

  case ( 'spl13' ) ! Symmetric piecewise linear scheme, see Waterson and Deconinck JCP 224 (2007)
    psi = max(0., min(2*r, 1./3.*r+2./3., 2./3.*r+1./3., 2.0))

  ! case ( 'splmax12' )
  !   psi = 

  ! case ( 'splmax13') 
  !   psi = 

  ! case ( 'superbee' ) 
  !   psi = 

  ! case ( 'copla' ) 
  !   psi = 

  ! case ( 'cubista' ) 
  !   psi = 

  ! case ( 'waceb' ) 
  !   psi = 

  ! case ( 'stoic' ) 
  !   psi = 

  ! case ( 'vonos' ) 
  !   psi = 

  ! case ( 'hoab' ) 
  !   psi = 

  ! case( 'vanalbada' )
  !   psi = 

    case default 

      write(*,'(a)') ' '
      write(*,'(a)') 'Fatal error: non-existing interpolation scheme!'
      stop

  end select

  vf = u(ijp) + fxp * psi*(u(ijn)-u(ijp))


end function


end module
