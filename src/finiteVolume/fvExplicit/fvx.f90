module fvx
!
! Purpose:
!  Module for explicit operations on discrete tensor fields.
!
! Description:
!  Module contains procedures for explicit manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!  That means we are expecting discrete volume/surface scalar/vector/tensor fields.
!  Included operations are:
!  fvxIterpolation (vol<Scalar/Vector/Tensor>Field -> surface<Scalar/Vector/Tensor>Field)
!  fvxGrad         (vol<Scalar/Vector>Field ->vol<Vector/Tensor>Field)
!  fvxDiv          (vol<Vector>Field ->vol<Scalar>Field) or (surface<Vector>Field ->vol<Scalar>Field)
!
!  Author: Nikola Mirkov
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!

use types
use geometry
use tensor_fields
use interpolation
use gradients
use fieldManipulation, only: explDiv

implicit none



interface fvxGrad
  module procedure fvx_grad_volScalarField
  module procedure fvx_grad_VolVectorField
end interface

interface fvxInterpolate
  module procedure fvcInterpolateScalar
  module procedure fvcInterpolateVector
end interface
  
interface fvxDiv
  module procedure fvx_div_volVectorField
  module procedure fvx_div_surfaceVectorField
end interface

interface fvxSurfaceSum
  module procedure fvcSurfaceSumGivenVolField
  module procedure fvcSurfaceSumGivenSurfaceField
end interface

interface fvxSurfaceIntegrate
  module procedure fvcSurfaceIntegrateGivenVolField
  module procedure fvcSurfaceIntegrateGivenSurfaceField
end interface

interface fvxAverage
  module procedure fvcAverageGivenVolField
  module procedure fvcAverageGivenSurfaceField
end interface


public

contains


function fvx_div_volVectorField(U) result(phi)
!
! Description:
!  Divergence of a volume vector field.
! Usage:
!    [type(volumeScalarField)] phi = fvxDiv( [type(volumeVectorField)] U )
!
  use fieldManipulation

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(volScalarField) :: phi

!+-----------------------------------------------------------------------------+

  phi = new_volScalarField( numCells )

  phi%field_name = 'divergence_field'

  !phi % mag = explDiv( U%x, U%y, U%z )


    ! Calculate cell-centered gradient
    !call updateBoundary(u,'boundary_region_name', zeroGrad/value/nGrad)
    call grad(u,dUdxi)
    call grad(v,dVdxi)
    call grad(w,dWdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u,v,w, dUdxi,dvdxi,dWdxi, dfxe)
        ! Value of the variable at cell-face center
  uf = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, u, du, scheme )
  vf = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, v, dv, scheme )
  wf = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, w, dw, scheme )

  ! (interpolated mid-face value)x(area)
  dfxe = uf*sx + vf*sy + wf*sz

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp)+dfxe
      div(ijn) = div(ijn)-dfxe

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb),v(ijb),w(ijb), dfxe)
      !
    dfx = uf*sx + vf*sy + wf*sz
      ! I'm putting minus here, because the face normal is allways facing out
      ! and if we get positive alignement with that vector we will have a positive contribution 
      ! to divergence, but the thing is in fact opposite since divergence 
      ! grows if something goes into the volume      
      div(ijp) = div(ijp) - dfxe
    enddo

end function fvx_div_volVectorField


function fvx_div_surfaceVectorField(U) result(phi)
!
! Description:
!  Divergence of a surface vector field.
! Usage:
!    [type(volumeScalarField)] phi = fvxDiv( [type(surfaceVectorField)] U )
!
  implicit none

  type(surfaceVectorField), intent(in) :: U
!
! > Result
!
  type(volScalarField) :: phi

!
! > Locals
!
  integer :: i,iface
  integer :: icell,jcell

  real(dp) :: sfuf

!
! > Initialize
!
  phi = new_volScalarField( numCells )

  phi%field_name = 'divergence_field'

  phi % mag = 0.0_dp

!
! > Inner faces contribution
!
  inner_face_loop: do iface=1,numInnerFaces 
!+-----------------------------------------------------------------------------+

  icell = owner(iface)
  jcell = neighbour(iface)

  ! Dot face area vector and interpolated vector to face, Sf.(u)_f
  sfuf = u%x (iface)*arx(iface) + u%y (iface)*ary(iface) + u%z (iface)*arz(iface)

  ! Contribution to owner and neighbour cell
  phi%mag (icell) = phi%mag (icell) + sfuf
  phi%mag (jcell) = phi%mag (jcell) - sfuf

!+-----------------------------------------------------------------------------+
  enddo inner_face_loop

!
! > Boundary faces contribution
!

  boundary_face_loop: do i=1,numBoundaryFaces
!+-----------------------------------------------------------------------------+  	
  iface = numInnerFaces + i
  icell = owner(iface)

  ! Dot face area vector and interpolated vector to face, Sf.(u)_f
  sfuf = u%x (iface)*arx(iface) + u%y (iface)*ary(iface) + u%z (iface)*arz(iface)

  ! Contribution to owner 
  ! NOTE:
  ! I'm putting minus here, because the face normal is allways facing out
  ! and if we get positive alignement with that vector we will have a positive contribution 
  ! to divergence, but the thing is in fact opposite since divergence 
  ! grows if something goes into the volume   
  phi%mag (icell) = phi%mag (icell) - sfuf

!+-----------------------------------------------------------------------------+
  enddo boundary_face_loop


end function fvx_div_surfaceVectorField




end module