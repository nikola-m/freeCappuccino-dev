module fviDivergence
!
! Purpose:
!  Module for IMPLICIT operations on discrete tensor fields: The Divergence oeprator.
!
! Description:
!  Module contains procedures for IMPLICIT manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!  That means we are expecting discrete volume/surface scalar/vector/tensor fields.
!  Included operations are:
!  fviDiv          (vol<Scalar/Vector/Tensor>Field -> type(fvEquation)/type(fvVectorEquation))
!
!  Author: Nikola Mirkov
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!

use types
use geometry
use tensorFields
use fvxInterpolation
use gradients

implicit none

interface fviDiv
  module procedure fvi_div_volScalarField
  module procedure fvi_div_volVectorField
  module procedure fvi_div_surfaceVectorField
end interface

public

contains



function fvi_div_volScalarField(fm,phi) result(fvEq)
!
! Description:
!    Implicit FVM discretization of divergence operator.
! Usage:
!    [type(fvEquation)] fvEq = fviDiv( [type(surfaceScalarField)] fm, [type(volumeScalarField)] phi )
!

  implicit none

  type(surfaceScalarField), intent(in) :: fm
  type(volScalarField),     intent(in) :: phi
!
! > Result
!
  type(fvVectorEquation) :: fvEq

!
! > Local
!
  integer :: i, k, ijp, ijn, ijb, iface
  real(dp) :: cap, can
  real(dp) :: are,dpw
  real(dp) :: gam
  real(dp) :: suadd


  ! Initialize matrix array and rhs vector
  fvEq%a = 0.0_dp
  fvEq%su = 0.0_dp

  ! > Assemble Laplacian system matrix

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxdiv(ijp, ijn, xf(i), yf(i), zf(i), fm%mag(i), facint(i), gam, phi, dPhidxi, cap, can, suadd)

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    fvEq%a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    fvEq%a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    fvEq%a(k) = fvEq%a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    fvEq%a(k) = fvEq%a(k) - cap

    ! > Sources:

    fvEq%su(ijp) = fvEq%su(ijp) + suadd
    fvEq%su(ijn) = fvEq%su(ijn) - suadd 

  end do


  ! > Boundary conditions

  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' ) then
      do i=1,nfaces(ib)
        if  = startFace(ib) + i
        ijp = owner(if)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, &
                         xf(if), yf(if), zf(if), &
                         fm(if), 0.0, 1.0, &
                         phi, dPhidxi, cap, can, suadd)

        fvEq%sp(ijp) = fvEq%sp(ijp) - can
        fvEq%su(ijp) = fvEq%su(ijp) - can*phi(ijb) + suadd

      enddo

    !if ( bctype(ib) == '') then

      ! If boundary is fixed value type
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
 
        fvEq%su(ijp) = fvEq%su(ijp) - fm(iface)*phi(ijb)

        ! can = min(fm,zero)
        ! fvEq%sp(ijp) = fvEq%sp(ijp) - can
        ! fvEq%su(ijp) = fvEq%su(ijp) - can*fi(ijb) + suadd

      enddo

    !if ( bctype(ib) == '') then

      ! If boundary is fixed gradient type
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! First determine boundary value by extrapolation
        dpw = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )
        phi(ijb) = phi(ijp) + gradVal*dpw

        ! Move to RHS vector
        fvEq%su(ijp) = fvEq%su(ijp) - fm(iface)*phi(ijb)

      enddo

    endif

  enddo 


end subroutine


function fvi_div_volVectorField(fm,U) result(fvEq)
!
! Description:
!    Implicit FVM discretization of divergence operator.
! Usage:
!    [type(fvVectorEquation)] fvEq = fviDiv( [type(surfaceScalarField)] fm, [type(volumeVectorField)] U )
!

  implicit none

  type(surfaceScalarField), intent(in) :: fm
  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(fvVectorEquation) :: fvEq

!
! > Local
!
  integer :: i, k, ijp, ijn, ijb, iface
  real(dp) :: cap, can
  real(dp) :: are,dpw
  real(dp) :: gam
  real(dp) :: suadd


  ! Initialize matrix array
  fvEq%a = 0.0_dp

  ! > Assemble Laplacian system matrix

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

!
! [NOTE] Sta sa gds(iu) kako to iskoristiti. Mozda to treba da bude svojstvo FvEqn? jer tu je odredjeno koliko
!        je implicit koliko explicit...
!

    call facefluxdivuvw(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), fm(i), facint(i), gds(iu), cap, can, sup, svp, swp)

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    fvEq%a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    fvEq%a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    fvEq%a(k) = fvEq%a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    fvEq%a(k) = fvEq%a(k) - cap

    ! > Sources: 

    fvEq%su(ijp) = fvEq%su(ijp) + sup
    fvEq%sv(ijp) = fvEq%sv(ijp) + svp
    fvEq%sw(ijp) = fvEq%sw(ijp) + swp

    fvEq%su(ijn) = fvEq%su(ijn) - sup
    fvEq%sv(ijn) = fvEq%sv(ijn) - svp
    fvEq%sw(ijn) = fvEq%sw(ijn) - swp

  end do


  ! > Boundary conditions

  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' ) then
      do i=1,nfaces(ib)
        if  = startFace(ib) + i
        ijp = owner(if)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, xf(if), yf(if), zf(if), fm(if), 0.0, 1.0, fi, dFidxi, cap, can, suadd)

        fvEq%sp(ijp) = fvEq%sp(ijp) - can
        fvEq%su(ijp) = fvEq%su(ijp) - can*fi(ijb) + suadd

      enddo

    endif

  enddo 


end subroutine
  

!***********************************************************************
!
subroutine facefluxdivuvw(ijp, ijn, xf, yf, zf, arx, ary, arz, fm, lambda, gam, cap, can, sup, svp, swp)
!
!***********************************************************************
!
! Face fluxes of divergence for inner faces.
!
!***********************************************************************
! 
  implicit none


  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

  ! Local variables
  real(dp) :: fxp,fxn  
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: fuuds,fvuds,fwuds
  real(dp) :: fuhigh,fvhigh,fwhigh
  real(dp) :: ui, vi, wi
!----------------------------------------------------------------------


  ! > Geometry:

  ! Face interpolation factor
  fxn = lambda 
  fxp = 1.0_dp-lambda

  ! Distance vector between cell centers
  xpn = xc(ijn)-xc(ijp)
  ypn = yc(ijn)-yc(ijp)
  zpn = zc(ijn)-zc(ijp)

  ! > Equation coefficients:

  ! > Equation coefficients - implicit convection
  ce = min(fm,zero) 
  cp = max(fm,zero)

  cap = - cp
  can =   ce


  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds = cp*u(ijp)+ce*u(ijn)
  fvuds = cp*v(ijp)+ce*v(ijn)
  fwuds = cp*w(ijp)+ce*w(ijn)


! EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
! Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_f[kg/s] * Phi_f[m/s]$

  if( fm .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    ui = face_value(ijp, ijn, xf, yf, zf, fxp, u, dUdxi)
    vi = face_value(ijp, ijn, xf, yf, zf, fxp, v, dVdxi)
    wi = face_value(ijp, ijn, xf, yf, zf, fxp, w, dWdxi)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ui = face_value(ijn, ijp, xf, yf, zf, fxn, u, dUdxi)
    vi = face_value(ijn, ijp, xf, yf, zf, fxn, v, dVdxi)
    wi = face_value(ijn, ijp, xf, yf, zf, fxn, w, dWdxi)
  endif

  fuhigh = fm*ui
  fvhigh = fm*vi
  fwhigh = fm*wi


! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = -gam*(fuhigh-fuuds)
  svp = -gam*(fvhigh-fvuds)
  swp = -gam*(fwhigh-fwuds)

end subroutine



!***********************************************************************
!
subroutine facefluxdiv(ijp, ijn, xf, yf, zf, fm, lambda, gam, fi, dFidxi, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use interpolation, only: face_value

  implicit none
!
!***********************************************************************
! 

! Arguments
  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: fm
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: fi
  real(dp), dimension(3,numTotal), intent(in) :: dFidxi
  real(dp), intent(inout) :: cap, can, suadd

! Local variables
  real(dp) :: fxp,fxn
  real(dp) :: cp,ce
  real(dp) :: fii
  real(dp) :: fcfie,fcfii
  real(dp), parameter :: zero = 0.0

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! > Convection fluxes - upwind scheme
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! > Matrix coefficients - implicit part
  cap = - cp
  can =   ce

  ! > Explicit higher order approximation of face value

  if( fm .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dFidxi)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dFidxi)
  endif

  ! > Explicit, hight order, approximation of convection
  fcfie = fm*fii

  ! > Explicit, first order, approximation of convection
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  ! > Deffered correction source for convection = gamablend*(high-low)
  suadd = -gam*(fcfie-fcfii)


end subroutine

end module


