module fvxDdt
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

implicit none


interface fvxDdt
  module procedure fvx_ddt_volVectorField
  module procedure fvx_ddt_volScalarField
end interface

public

contains

subroutine fvx_ddt_volVectorField
!  
!******************************************************************************
!
!     Explicit discretisation of time derivative.
!
!******************************************************************************
!
! Description:
!  Ddt operator of a volume vector field.
! Usage:
!    [type(volumeVectorField)] DUDt = fvxDdt( [type(volumeVectorField)] U )
!
  use fieldManipulation

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(volVectorField) :: DUDt

!+-----------------------------------------------------------------------------+

  DUDt = new_volVectorField( numCells )

  DUDt%field_name = 'DUDt_field'

  do inp=1,numCells

    if (bdf) then 

      apotime=den(inp)*vol(inp)/timestep

      DUDt%x(inp) = apotime*( u(inp) - uo(inp) )
      DUDt%y(inp) = apotime*( v(inp) - vo(inp) )
      DUDt%z(inp) = apotime*( w(inp) - wo(inp) )

    elseif (bdf2) then 

      apotime=den(inp)*vol(inp)/timestep

      DUDt%x(inp) = apotime*( 1.5*u(inp) - 2*uo(inp) + 0.5*uoo(inp) )
      DUDt%y(inp) = apotime*( 1.5*v(inp) - 2*vo(inp) + 0.5*voo(inp) )
      DUDt%z(inp) = apotime*( 1.5*w(inp) - 2*wo(inp) + 0.5*woo(inp) )

    endif

  end do


end subroutine


subroutine fvx_ddt_volScalarField
!  
!******************************************************************************
!
!   Explicit discretisation of time derivative.
!
!******************************************************************************
!
! Description:
!  Ddt operator of a volume scalar field.
! Usage:
!    [type(volumeScalarField)] DphiDt = fvxDdt( [type(volumeScalarField)] phi )
!
  use fieldManipulation

  implicit none

  type(volScalarField), intent(in) :: phi
!
! > Result
!
  type(volScalarField) :: DphiDt

!+-----------------------------------------------------------------------------+

  DphiDt = new_volScalarField( numCells )

  DphiDt%field_name = 'DUDt_field'

  do inp=1,numCells

    if (bdf) then 

      apotime=den(inp)*vol(inp)/timestep

      DphiDt%mag(inp) = apotime*( phi(inp) - phio(inp) )

    elseif (bdf2) then 

      apotime=den(inp)*vol(inp)/timestep

      DphiDt%mag(inp) = apotime*( 1.5*phi(inp) - 2*phio(inp) + 0.5*phioo(inp) )

    endif

  end do


end subroutine



subroutine fvx_d2dt2_scalar_field
!  
!******************************************************************************
!
!   Description: 
!     Second order time differentiation for scalar fields.
!
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix

  implicit none

  real(dp) :: apotime

  ! spu(1:numCells) = den(1:numCells)*vol(1:numCells)/timestep**2
  ! su(1:numCells) = den(1:numCells)*vol(1:numCells)/timestep**2 * (2*uo(1:numCells) - uoo(1:numCells))


  do inp=1,numCells

    if(bdf) then

      apotime=den(inp)*vol(inp)/timestep**2

      su(inp) = su(inp) + apotime*(2*uo(inp) - uoo(inp))

      spu(inp) = spu(inp) + apotime

    endif

  end do


end subroutine