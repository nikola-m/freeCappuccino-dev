! 00204     fvm.diag() += mesh.V()*max(susp.field(), scalar(0));
! 00205 
! 00206     fvm.source() -= mesh.V()*min(susp.field(), scalar(0))

module fviSources
!
!  Purpose:
!    Module contanins functions for discretisation of source terms.
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


  interface fviSu
    module procedure fvi_su_scalar_field
    module procedure fvi_su_vector_field
  end interface

  interface fviSp
    module procedure fvi_sp_scalar_field
    module procedure fvi_sp_vector_field
  end interface

  interface fviSuSp
    module procedure fvi_susp_scalar_field
    module procedure fvi_susp_vector_field
  end interface

public

contains


function fvi_su_scalar_field(su,psi) result(fvEqn)
!
! Description:
!    Explicit discretization of the source term.
!
! Usage:
!   

  implicit none


!
! > Locals
!
  integer :: icell

!
! > Cell loop 
!
  cell_loop: do icell=1,numCells 

     fvEqn%su( icell ) = fvEqn%su( icell ) - su%mag(icell) * Vol(icell)
 
  enddo cell_loop

end function fvi_su_scalar_field



function fvi_su_vector_field(su,U) result(fvEqn)
!
! Description:
!    Explicit discretization of the source term for the vector field U.
!
! Usage:
!   

  implicit none

  type(volVectorField), intent(in) :: su
  type(volVectorField), intent(in), optional :: U

!
! > Result
!
  type(fvVectorEquation) :: fvEqn
!
! > Locals
!
  integer :: icell

!
! > Cell loop 
!
  cell_loop: do icell=1,numCells 

    fvEqn%su( icell ) = fvEqn%su( icell ) - su%x(icell) * Vol(icell)
    fvEqn%sv( icell ) = fvEqn%sv( icell ) - su%y(icell) * Vol(icell)
    fvEqn%sw( icell ) = fvEqn%sw( icell ) - su%z(icell) * Vol(icell)
 
  enddo cell_loop

end function fvi_su_vector_field




function fvi_sp_scalar_field(sp,psi) result(fvEqn)
!
! Description:
!    Implicit discretization of the source term.
!
! Usage:
!   

  implicit none

  type(volScalarField), intent(in) :: sp
  type(volScalarField), intent(in) :: psi

!
! > Result
!
  type(fvEquation) :: fvEqn
!
! > Locals
!
  integer :: icell

!
! > Cell loop 
!
  cell_loop: do icell=1,numCells 

     fvEqn%a( diag(icell) ) = fvEqn%a( diag(icell) ) + sp%mag(icell) / ( psi%mag(icell) + small ) * Vol(icell)
     ! ...or...
     ! fvEqn%sp( icell ) = fvEqn%sp( icell ) + sp%mag(icell) / ( psi%mag(icell) + small ) * Vol(icell)
 
  enddo cell_loop

end function fvi_sp_scalar_field





function fvi_sp_vector_field(sp,U) result(fvEqn)
!
! Description:
!    Implicit discretization of the source term for the vector field U.
!
! Usage:
!   

  implicit none

  type(volVectorField), intent(in) :: sp
  type(volVectorField), intent(in) :: U

!
! > Result
!
  type(fvVectorEquation) :: fvEqn
!
! > Locals
!
  integer :: icell

!
! > Cell loop 
!
  cell_loop: do icell=1,numCells 

    fvEqn%spu(icell) = fvEqn%spu(icell) + sp%x(icell) / ( U%x(icell) + small ) * Vol(icell)
    fvEqn%spv(icell) = fvEqn%spv(icell) + sp%y(icell) / ( U%y(icell) + small ) * Vol(icell)
    fvEqn%spw(icell) = fvEqn%spw(icell) + sp%z(icell) / ( U%z(icell) + small ) * Vol(icell)
 
  enddo cell_loop

end function fvi_sp_vector_field


function fvi_susp_scalar_field(susp,psi) result(fvEqn)
!
! Description:
!    Implicit-Explicit discretization of the source term.
!
!    Comment:
!    We don't allow source term to spoil diagonal dominance
!    of the linear system matrix, so if it is negative,
!    we automatically move it to rhs, i.e. to source of fvEquation.
!
! Usage:
!   

  implicit none

  type(volScalarField), intent(in) :: susp
  type(volScalarField), intent(in) :: psi

!
! > Result
!
  type(fvEquation) :: fvEqn
!
! > Locals
!
  integer :: icell

!
! > Cell loop 
!
  cell_loop: do icell=1,numCells 

    fvEqn%a( diag(icell) ) = fvEqn%a( diag(icell) ) + max( susp%mag(icell), zero) / ( psi%mag(icell) + small ) * Vol(icell)

    fvEqn%su( icell ) = fvEqn%su( icell ) - min( susp%mag(icell), zero ) * Vol(icell)

  enddo cell_loop

end function fvi_susp_scalar_field





function fvi_susp_vector_field(susp,U) result(fvEqn)
!
! Description:
!    Implicit-Explicit discretization of the source term.
!
!    Comment:
!    We don't allow source term to spoil diagonal dominance
!    of the linear system matrix, so if it is negative,
!    we automatically move it to rhs, i.e. to source of fvVectorEquation.
!
! Usage:
!   

  implicit none

  type(volVectorField), intent(in) :: susp
  type(volVectorField), intent(in) :: U

!
! > Result
!
  type(fvVectorEquation) :: fvEqn
!
! > Locals
!
  integer :: icell

!
! > Cell loop 
!
  cell_loop: do icell=1,numCells 

    fvEqn%spu(icell) = fvEqn%spu(icell) + max( susp%x(icell), zero) / ( U%x(icell) + small ) * Vol(icell)
    fvEqn%spv(icell) = fvEqn%spv(icell) + max( susp%y(icell), zero) / ( U%y(icell) + small ) * Vol(icell)
    fvEqn%spw(icell) = fvEqn%spw(icell) + max( susp%z(icell), zero) / ( U%z(icell) + small ) * Vol(icell)

    fvEqn%su( icell ) = fvEqn%su( icell ) - min( susp%x(icell), zero ) * Vol(icell)
    fvEqn%sv( icell ) = fvEqn%sv( icell ) - min( susp%y(icell), zero ) * Vol(icell)
    fvEqn%sw( icell ) = fvEqn%sw( icell ) - min( susp%z(icell), zero ) * Vol(icell)

 
  enddo cell_loop

end function fvi_susp_vector_field

end module