module fviDdt
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

implicit none

  real(dp), parameter :: small = 1e-30
  real(dp), parameter :: zero = 0.0_dp

  interface fviDdt
    module procedure fvi_ddt_volScalarField
    module procedure fvi_ddt_volVectorField
    module procedure fvi_ddt_surfaceVectorField
  end interface

  interface fviD2dt2
    module procedure fvi_d2dt2_scalar_field
  end interface

public :: fviDdt,fviD2dt2

contains

function fvi_ddt_vector_field(den,U) result(fvEqn)
!  
!******************************************************************************
!
! Finds source and matrix coefficients representing FVM discretization 
! of time derivative operator: \partial / \partial t.
!
! System of linear equations is written as:
! $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!

  implicit none

  type(volVectorField), intent(in), optional :: den
  type(volVectorField), intent(in) :: U

!
! > Result
!
  type(fvVectorEquation) :: fvEqn

  integer :: icell 
  real(dp) :: apotime

!
! > Shift in time
! 
  ! fvEqn%xoo(:) = fvEqn%xo(:)
  ! fvEqn%yoo(:) = fvEqn%yo(:)
  ! fvEqn%zoo(:) = fvEqn%zo(:)

  ! fvEqn%xo(:) = U%x(:)
  ! fvEqn%yo(:) = U%y(:)
  ! fvEqn%zo(:) = U%z(:)


  do icell =1,numCells

    if( present(den) ) then
      apotime = den( icell )*vol( icell )/timestep
    else
      apotime = vol( icell )/timestep
    endif

    if( bdf .or. cn ) then
    !
    ! Backward differentiation formula of 1st order or Crank-Nicolson.
    !

      ! RHS vector contribution
      fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*uo( icell )
      fvEqn%sv( icell ) = fvEqn%sv( icell ) + apotime*vo( icell )
      fvEqn%sw( icell ) = fvEqn%sw( icell ) + apotime*wo( icell )

      ! Matrix diagonal element contribution
      fvEqn%spu( icell ) = fvEqn%spu( icell ) + apotime
      fvEqn%spv( icell ) = fvEqn%spv( icell ) + apotime
      fvEqn%spw( icell ) = fvEqn%spw( icell ) + apotime

    elseif( bdf2 ) then
    !
    ! Backward Differentiation (BDF2) - 2nd order.
    !

      ! RHS vector contribution
      fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*( 2*uo( icell ) - 0.5_dp*uoo( icell ) )
      fvEqn%sv( icell ) = fvEqn%sv( icell ) + apotime*( 2*vo( icell ) - 0.5_dp*voo( icell ) ) 
      fvEqn%sw( icell ) = fvEqn%sw( icell ) + apotime*( 2*wo( icell ) - 0.5_dp*woo( icell ) )

      ! Matrix diagonal element contribution
      fvEqn%spu( icell ) = fvEqn%spu( icell ) + 1.5_dp*apotime
      fvEqn%spv( icell ) = fvEqn%spv( icell ) + 1.5_dp*apotime
      fvEqn%spw( icell ) = fvEqn%spw( icell ) + 1.5_dp*apotime


    elseif( bdf3 ) then
    !
    ! Backward Differentiation (BDF3) - 3nd order.
    !

      ! RHS vector contribution
      fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*( 3*uo( icell ) - 1.5_dp*uoo( icell ) + 1./3.0_dp*uooo( icell ) )
      fvEqn%sv( icell ) = fvEqn%sv( icell ) + apotime*( 3*vo( icell ) - 1.5_dp*voo( icell ) + 1./3.0_dp*vooo( icell ) ) 
      fvEqn%sw( icell ) = fvEqn%sw( icell ) + apotime*( 3*wo( icell ) - 1.5_dp*woo( icell ) + 1./3.0_dp*wooo( icell ) )

      ! Matrix diagonal element contribution
      fvEqn%spu( icell ) = fvEqn%spu( icell ) + 11./6.0_dp*apotime
      fvEqn%spv( icell ) = fvEqn%spv( icell ) + 11./6.0_dp*apotime
      fvEqn%spw( icell ) = fvEqn%spw( icell ) + 11./6.0_dp*apotime

    endif

  end do

end function


function fvi_ddt_scalar_field(den, psi) result(fvEqn)
!  
!******************************************************************************
!
!     Finds source and matrix coefficients representing FVM discretization 
!     of time derivative operator: \partial / \partial t.
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
  use parameters

  implicit none

  type(volScalarField), intent(in), optional :: den
  type(volScalarField), intent(in) :: U

!
! > Result
!
  type(fvEquation) :: fvEqn

!
! > Local
! 
  integer ::  icell 
  real(dp) :: apotime, sut

!
! > Shift in time
! 
  ! fvEqn%oo(:) = fvEqn%o(:)
  ! fvEqn%o(:) = phi%mag(:)

  do icell = 1,numCells

    if( present(den) ) then
      apotime = den( icell )*vol( icell )/timestep
    else
      apotime = vol( icell )/timestep
    endif

    if( bdf .or. cn ) then
    !
    ! Backward differentiation formula of 1st order.
    !

      ! RHS vector contribution
      fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*uo( icell )

      ! Matrix diagonal element contribution
      fvEqn%sp( icell )  = fvEqn%sp( icell ) + apotime

    elseif( bdf2 ) then
    !
    ! Backward Differentiation (BDF2) - 2nd order.
    !

      ! RHS vector contribution
      fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*( 2*uo( icell ) - 0.5_dp*uoo( icell ) )

      ! Matrix diagonal element contribution
      fvEqn%sp( icell )  = fvEqn%sp( icell )  + 1.5_dp*apotime


    elseif( bdf3 ) then
    !
    ! Backward Differentiation (BDF3) - 3nd order.
    !

      ! RHS vector contribution
      fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*( 3*uo( icell ) - 1.5_dp*uoo( icell ) + 1./3.0_dp*uooo( icell ) )

      ! Matrix diagonal element contribution
      fvEqn%sp( icell ) = fvEqn%sp( icell )  + 11./6.0_dp*apotime

    endif

  end do


end function



function fvi_d2dt2_scalar_field( den, U ) result(fvEqn)
!  
!******************************************************************************
!
!   Description: 
!     Second order time differentiation for scalar fields.
!
!     Finds source and matrix coefficients representing FVM discretization 
!     of scnd. order time derivative operator: \partial^2 / \partial t^2.
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!

  use parameters

  implicit none

  type(volScalarField), intent(in), optional :: den
  type(volScalarField), intent(in) :: U
!
! > Result
!
  type(fvEquation) :: fvEqn
!
! > Locals
!  
  integer :: icell
  real(dp) :: apotime
!
! > Shift in time
! 
  ! fvEqn%oo(:) = fvEqn%o(:)
  ! fvEqn%o(:) = phi%mag(:)

  do icell = 1,numCells

    if( present(den) ) then
      apotime = den( icell )*vol( icell )/timestep**2
    else
      apotime = vol( icell )/timestep**2
    endif

    fvEqn%su( icell ) = fvEqn%su( icell ) + apotime*(2*uo( icell ) - uoo( icell ))

    fvEqn%sp( icell ) = fvEqn%sp( icell ) + apotime

  end do


end function