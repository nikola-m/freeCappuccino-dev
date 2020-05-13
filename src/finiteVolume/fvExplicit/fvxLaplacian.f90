module fvxLaplacian
!
! Purpose:
!  Module for explicit operations on discrete tensor fields  - explicit Laplacian.
!
! Description:
!  Module contains procedures for explicit manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!  That means we are expecting discrete volume/surface scalar/vector/tensor fields.
!
!  fvxLaplacian          (vol<Scalar/Vector>Field ->vol<Scalar>Field) or (surface<Vector>Field ->vol<Scalar>Field)
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

interface fvxLaplacian
  module procedure fvx_Lapl_volScalarField
  module procedure fvx_Lapl_volVectorField
  ! module procedure fvx_Lapl_surfaceVectorField
end interface

public

contains


function fvx_Lapl_volScalarField(phi) result(Lapl)
!
! Description:
!  Laplacian of a volume scalar field.
! Usage:
!    [type(volumeScalarField)] Lapl = fvxDiv( fvxGrad ( [type(volumeScalarField)] phi ) )
!

  implicit none

  type(volScalarField), intent(in) :: phi
!
! > Result
!
  type(volScalarField) :: Lapl

!+-----------------------------------------------------------------------------+

  Lapl = new_volScalarField( numCells )

  Lapl%field_name = 'Laplacian_field'

  Lapl%mag = 0.0

  Lapl = fvxDiv( fvxGrad( phi ) )

end function


function fvx_Lapl_volVectorField(U) result(Lapl)
!
! Description:
!  Laplacian of a volume vector field.
! Usage:
!    [type(volumeVectorField)] Lapl = fvxDiv( fvxGrad ( [type(volumeVectorField)] U ) )
!

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(volVectorField) :: Lapl

!+-----------------------------------------------------------------------------+

  Lapl = new_volVectorField( numCells )

  Lapl%field_name = 'Laplacian_field'

  Lapl%x = 0.0
  Lapl%y = 0.0
  Lapl%z = 0.0

  Lapl = fvxDiv( fvxGrad( phi ) )

end function


end module