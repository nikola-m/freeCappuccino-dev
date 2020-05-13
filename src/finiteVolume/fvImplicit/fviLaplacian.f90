module fviLaplacian
!
! Purpose:
!  Module for IMPLICIT operations on discrete tensor fields: The Laplacian operator.
!
! Description:
!  Module contains procedures for IMPLICIT manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!  That means we are expecting discrete volume/surface scalar/vector/tensor fields.
!  Included operations are:
!  fviLaplacian    (vol<Scalar/Vector/Tensor>Field -> type(fvEquation)/type(fvVectorEquation))
!
!  Author: Nikola Mirkov
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!

use types
use geometry
use tensorFields
use fvEquation
! use fvxInterpolation
! use gradients

implicit none

interface fviLaplacian
  module procedure fvi_Lapl_volScalarField
  ! module procedure fvi_Lapl_volVectorField
  ! module procedure fvi_Lapl_surfaceVectorField
end interface

public

contains

function fvi_Lapl_volScalarField(mu,phi) result(fvEq)
!
! Description:
!    Implicit FVM discretization of Laplacian operator.
! Usage:
!    [type(fvVectorEquation)] fvEq = fviLaplacian( [type(volumeScalarField)] phi )
!
! System of linear equations is written as:
! $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!

  implicit none

  type(volumeScalarField), intent(in) :: phi
  type(volumeScalarField), intent(in) :: mu

!
! > Result
!
  type(fvVectorEquation) :: fvEq

  !
  ! Local variables
  !

  integer :: i, k, ijp, ijn
  integer :: ib,ijb,iface
  real(dp) :: cap, can
  real(dp) :: are,dpw,dcoef


  ! Initialize matrix array
  fvEq%a = 0.0_dp

  ! > Assemble Laplacian system matrix

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxlaplacian(ijp, ijn, arx(i), ary(i), arz(i), facint(i), mu%mag, cap, can)

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

  end do


!.....Modify matrix coefficients to reflect presence of Boundary Conditions in PDE problem.

  do ib=1,numBoundaries

    !if ( bctype(ib) == '') then

      ! If boundary is fixed value type
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
        dpw = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

        dcoef   = mu%mag(ijp)*are/dpw
        fvEq%a( diag(ijp) ) = fvEq%a( diag(ijp) ) - dcoef  
        fvEq%su(ijp) = fvEq%su(ijp) - dcoef*phi(ijb)

      enddo

    !if ( bctype(ib) == '') then

      ! If boundary is fixed gradient type
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        fvEq%su(ijp) = fvEq%su(ijp) - gradVal*are

      enddo


    !endif
    
  enddo


end subroutine



!***********************************************************************
!
subroutine facefluxuvw(ijp, ijn, xf, yf, zf, arx, ary, arz, flomass, lambda, gam, cap, can, sup, svp, swp)
!
!***********************************************************************
!
! Face fluxes of velocity for inner faces.
!
!***********************************************************************
! 
  implicit none


  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

  ! Local variables
  integer  :: nrelax
  character(len=12) :: approach
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: duxi,duyi,duzi, &
              dvxi,dvyi,dvzi, &
              dwxi,dwyi,dwzi
  real(dp) :: duxii,dvxii,dwxii, &
              duyii,dvyii,dwyii, &
              duzii,dvzii,dwzii
  real(dp) :: de, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
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


  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijp)*fxp+vis(ijn)*fxn

  ! Difusion coefficient
  ! de = game*are/dpn
  de = game*(arx*arx+ary*ary+arz*arz)/(xpn*arx+ypn*ary+zpn*arz)

  can = -de
  cap = -de


  ! > Explicit diffusion: 

  nrelax = 0
  approach  = 'skewness'

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              u, dudxi, nrelax, approach, duxi, duyi, duzi,  &
              duxii, duyii, duzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              v, dvdxi, nrelax, approach, dvxi, dvyi, dvzi, &
              dvxii, dvyii, dvzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              w, dwdxi, nrelax, approach, dwxi, dwyi, dwzi, &
              dwxii, dwyii, dwzii)

!---------------------------------------------------------------------------------------
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector:
!     su = su + (fdue-fdui)
!     sv = sv + (fdve-fdvi)
!     sw = sw + (fdwe-fdwi)
!---------------------------------------------------------------------------------------

  ! Explicit diffussion: 
  fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
  fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
  fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )

  ! Implicit diffussion:
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)


! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = fdue-fdui
  svp = fdve-fdvi
  swp = fdwe-fdwi

end subroutine



!***********************************************************************
!
subroutine faceFluxLaplacian(ijp, ijn, arx, ary, arz, lambda, mu, cap, can)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numCells,xc,yc,zc

  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numCells), intent(in) :: mu
  real(dp), intent(inout) :: cap, can

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: xpn,ypn,zpn
  real(dp) :: smdpn

  ! real(dp) :: are
  ! real(dp) :: dpn

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  !             _    _
  ! First way: |Sf|/|dpn|
  !
  ! distance from p to p_j
  ! dpn = sqrt(xpn**2+ypn**2+zpn**2)
  ! cell face area
  ! are=sqrt(arx**2+ary**2+arz**2)
  ! smdpn = are/dpn

  !             _    _     _    _
  ! Second way: Sf . Sf / (Sf . dpn)
  smdpn = (arx*arx+ary*ary+arz*arz)/(arx*xpn+ary*ypn+arz*zpn)

  ! Coefficients of discretized Laplace equation
  cap = (fxp*mu(ijp)+fxn*mu(ijn))*smdpn
  can = cap

end subroutine

end module
