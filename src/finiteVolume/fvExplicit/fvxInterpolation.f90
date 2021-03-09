module fvxInterpolation
!
! Purpose:
!  Module for explicit operations on discrete tensor fields - Interpolation operator.
!
! Description:
!  Module contains procedures for explicit manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!  That means we are expecting discrete volume/surface scalar/vector/tensor fields.
!  Included operations are:
!  fvxInterpolation (vol<Scalar/Vector/Tensor>Field -> surface<Scalar/Vector/Tensor>Field)
!
!  Author: Nikola Mirkov
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!  
 use tensorFields
 use geometry
 use fvxGradient, only: grad,grad_scalar_field_w_option

 implicit none

  ! The inerpolations that we do for e.g. divergence is hardcoded here
  ! Recommended schemes are: 'central', 'cds', 'cdscorr'
  character(len=10), parameter :: interpolation_scheme = 'central' 

  interface fvxInterpolate
    module procedure fvxInterpolateScalar
    module procedure fvxInterpolateVector
  end interface

 public fvxInterpolate

 contains

function fvxInterpolateScalar(phi) result(psi)
!
! Description:
!  Creates surfaceScalarField from volScalarField by linear interpolation.
! Usage:
!   [type(surfaceScalarField)] psi = fvxInterpolate ( [type(volScalarField)] phi )
!

  implicit none

  type(volScalarField), intent(in) :: phi
!
! > Result
!
  type(surfaceScalarField) :: psi

!
! > Local
!
  type(volVectorField) :: dU

  integer :: i,ijp,ijn,ijb,iface

!+-----------------------------------------------------------------------------+

  psi = new_surfaceScalarField( numFaces )

  psi%field_name = 'interpolated_field'

    ! Calculate cell-centered gradient
    dU = Grad( phi )

  ! Interpolate to face

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      psi%mag (i) = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), phi%mag(i), dU%x, dU%y, dU%z )
    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijb = numCells + i
      psi%mag(iface) = phi%mag(ijb)
    enddo

end function fvxInterpolateScalar



function fvxInterpolateVector(U) result(psi)
!
! Description:
!  Creates surfaceVectorField from volVectorField by linear interpolation.
! Usage:
!   [type(surfaceVectorField)] psi = fvxInterpolate ( [type(volVectorField)] U )
!

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(surfaceVectorField) :: psi
!
! > Locals
!  
  type(volVectorField) :: dU

  integer :: i,ijp,ijn,ijb,iface

!+-----------------------------------------------------------------------------+

  psi = new_surfaceVectorField( numFaces )

  psi%field_name = 'interpolated_field'

  ! Calculate cell-centered gradient
  dU = new_volVectorField(numTotal)

  ! Calculate gradient
    call grad_scalar_field_w_option( U%x, dU%x, dU%y, dU%z, 'gauss', 'no-limit' )

  ! Interpolate to face

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    psi%x (i) = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), U%x, dU%x, dU%y, dU%z )
  enddo

  ! Contribution from boundaries
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijb = numCells + i
    psi%x (iface) = U%x(ijb)
  enddo


  ! Calculate cell-centered gradient
    call grad_scalar_field_w_option( U%y, dU%x, dU%y, dU%z, 'gauss', 'no-limit' )

  ! Interpolate to face

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    psi%y (i) = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), U%y, dU%x, dU%y, dU%z )
  enddo

  ! Contribution from boundaries
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijb = numCells + i
    psi%y(iface) = U%y(ijb)
  enddo

  ! Calculate cell-centered gradient
    call grad_scalar_field_w_option( U%z, dU%x, dU%y, dU%z, 'gauss', 'no-limit' )

  ! Interpolate to face

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    psi%z (i) = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), U%z, dU%x, dU%y, dU%z )
  enddo

  ! Contribution from boundaries
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijb = numCells + i
    psi%z(iface) = U%z(ijb)
  enddo

end function fvxInterpolateVector


!***********************************************************************
!
function face_value(ijp,ijn,xf,yf,zf,lambda,u,dUdx,dUdy,dUdz) result(ue)
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
  real(dp), dimension(numCells) :: dUdx,dUdy,dUdz

  if (interpolation_scheme == 'cdscorr') then 
    ue = face_value_cds_corrected(ijp, ijn, xf, yf, zf, lambda, u, dUdx,dUdy,dUdz)  

  elseif (interpolation_scheme == 'central') then 
    ue = face_value_central(ijp, ijn, xf, yf, zf, u, dUdx,dUdy,dUdz)

  elseif (interpolation_scheme == 'cds') then 
    ue = face_value_cds(ijp,ijn, lambda, u) 

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

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: ijp, ijn
  real(dp) :: lambda
  real(dp), dimension(numTotal) :: fi

  ! Locals
  real(dp) :: fxn,fxp

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  face_value = fi(ijp)*fxp+fi(ijn)*fxn

  end function



!***********************************************************************
!
  function face_value_cds_corrected(ijp,ijn, xf, yf, zf, lambda, fi, dfidx,dfidy,dfidz) result(face_value)
!
!***********************************************************************
!
! Calculates face value using values of variables and their gradients
! at centers of adjecent cells..
!
!***********************************************************************

  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: ijp, ijn
  real(dp) :: xf, yf, zf, lambda
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(numCells) :: dfidx,dfidy,dfidz

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
  dfixi = dFidx(ijp)*fxp+dFidx(ijn)*fxn
  dfiyi = dFidy(ijp)*fxp+dFidy(ijn)*fxn
  dfizi = dFidz(ijp)*fxp+dFidz(ijn)*fxn

  !            |________uj'___________|_________________ucorr____________________|
  face_value = fi(ijp)*fxp+fi(ijn)*fxn+(dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi))

  end function


!***********************************************************************
!
  function face_value_central(inp,inn, xf, yf, zf, fi, dfidx,dfidy,dfidz) result(face_value)
!
!***********************************************************************
!
! Calculates face value using values of variables and their gradients
! at centers of adjecent cells..
!
!***********************************************************************

  implicit none

  ! Result
  real(dp) :: face_value

  ! Input
  integer :: inp, inn
  real(dp) :: xf, yf, zf
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(numCells) :: dfidx,dfidy,dfidz

  ! Locals
  real(dp) :: gradfidr

  gradfidr = dfidx(inp)*(xf-xc(inp))+dfidy(inp)*(yf-yc(inp))+dfidz(inp)*(zf-zc(inp)) &
           + dfidx(inn)*(xf-xc(inn))+dfidy(inn)*(yf-yc(inn))+dfidz(inn)*(zf-zc(inn))

  face_value = 0.5_dp*( fi(inp) + fi(inn) + gradfidr)

  end function



end module