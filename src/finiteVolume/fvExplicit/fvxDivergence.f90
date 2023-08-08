module fvxDivergence
!
! Purpose:
!  Module for explicit operations on discrete tensor fields.
!
! Description:
!  Module contains procedures for explicit manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!
!  Included operations are:
!  fvxDiv  
!        
!  Input -> Output:
!  (vol<Vector>Field     -> vol<Scalar>Field) or 
!  (surface<Vector>Field -> vol<Scalar>Field)
!
!  Author: Nikola Mirkov
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!

use types
use geometry
use tensorFields
use fvxInterpolation, only: face_value, face_value_cds
use fvxGradient, only : grad

implicit none

interface fvxDiv
  module procedure fvx_div_volVectorField
  ! module procedure fvx_div_volVectorField2
  module procedure fvx_div_surfaceVectorField
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

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(volScalarField) :: phi

!
! > Local
!
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: uf,vf,wf,dfxe


  phi = new_volScalarField( numCells )

  phi%field_name = 'divergence_field'

  phi%mag = 0.0


  !
  ! > Inner faces contribution
  !

  inner_face_loop: do iface=1,numInnerFaces 

    ijp = owner(i)
    ijn = neighbour(i)

    ! Value of the variable at cell-face center
    uf = face_value_cds( ijp, ijn, facint(i), U%x)
    vf = face_value_cds( ijp, ijn, facint(i), U%y)
    wf = face_value_cds( ijp, ijn, facint(i), U%z)
        
    ! (interpolated mid-face value)x(area)
    dfxe = uf*arx(i) + vf*ary(i) + wf*arz(i)

    ! Accumulate contribution at cell center and neighbour.
    phi%mag(ijp) = phi%mag(ijp) + dfxe
    phi%mag(ijn) = phi%mag(ijn) - dfxe


  enddo inner_face_loop


  !
  ! > Boundary faces contribution
  !

  boundary_face_loop: do i=1,numBoundaryFaces

    iface = numInnerFaces + i
    ijp = owner(iface)
    ijb = numCells + i

    dfxe = U%x(ijb)*arx(iface)+ U%y(ijb)*ary(iface) + U%z(ijb)*arz(iface) 

    phi%mag(ijp) = phi%mag(ijp) + dfxe

  enddo boundary_face_loop

end function fvx_div_volVectorField



function fvx_div_volVectorField2(U) result(phi)
!
! Description:
!  Divergence of a volume vector field.
! Usage:
!    [type(volumeScalarField)] phi = fvxDiv( [type(volumeVectorField)] U )
!

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(volScalarField) :: phi

!
! > Local
!
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: uf,vf,wf,dfxe

  type(volTensorField) :: D


  phi = new_volScalarField( numCells )

  phi%field_name = 'divergence_field'

  phi%mag = 0.0

  D = Grad( U )

  !
  ! > Inner faces contribution
  !

  inner_face_loop: do iface=1,numInnerFaces 

    ijp = owner(i)
    ijn = neighbour(i)

    ! Value of the variable at cell-face center
    uf = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), U%x, D%xx, D%xy, D%xz )
    vf = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), U%y, D%yx, D%yy, D%yz )
    wf = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), U%z, D%zx, D%zy, D%zz )
        
    ! (interpolated mid-face value)x(area)
    dfxe = uf*arx(i) + vf*ary(i) + wf*arz(i)

    ! Accumulate contribution at cell center and neighbour.
    phi%mag(ijp) = phi%mag(ijp) + dfxe
    phi%mag(ijn) = phi%mag(ijn) - dfxe


  enddo inner_face_loop


  !
  ! > Boundary faces contribution
  !

  boundary_face_loop: do i=1,numBoundaryFaces

    iface = numInnerFaces + i
    ijp = owner(iface)
    ijb = numCells + i

    dfxe = U%x(ijb)*arx(iface)+ U%y(ijb)*ary(iface) + U%z(ijb)*arz(iface) 

    phi%mag(ijp) = phi%mag(ijp) + dfxe

  enddo boundary_face_loop

end function fvx_div_volVectorField2


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
  integer :: i,ijb,iface
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

  icell = owner(iface)
  jcell = neighbour(iface)

  ! Dot face area vector and interpolated vector to face, Sf.(u)_f
  sfuf = u%x (iface)*arx(iface) + u%y (iface)*ary(iface) + u%z (iface)*arz(iface)

  ! Contribution to owner and neighbour cell
  phi%mag (icell) = phi%mag (icell) + sfuf
  phi%mag (jcell) = phi%mag (jcell) - sfuf

  enddo inner_face_loop

!
! > Boundary faces contribution
!

  boundary_face_loop: do i=1,numBoundaryFaces
   
  iface = numInnerFaces + i
  icell = owner(iface)
  ijb = numCells + i

  ! Dot face area vector and interpolated vector to face, Sf.(u)_f
  sfuf = u%x (ijb)*arx(iface) + u%y (ijb)*ary(iface) + u%z (ijb)*arz(iface)

  phi%mag (icell) = phi%mag (icell) + sfuf

  enddo boundary_face_loop


end function fvx_div_surfaceVectorField


end module