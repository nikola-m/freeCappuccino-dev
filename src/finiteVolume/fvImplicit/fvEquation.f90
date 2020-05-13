module fvEquation
!
! Definition of Finite Volume Method linear equation object, and oparations over it.
!
use types
use geometry
use sparse_matrix
use tensorFields
use linear_solver, only: spsolve

implicit none

!
! The fvEquation derived data type (all objects from fvImplicit belong to this type.)
!
type, extends(csrMatrix) :: fvEquation

  real(dp), dimension(:), allocatable :: su   ! right hand side vector (vector of sources)
  real(dp), dimension(:), allocatable :: sp   ! vector of sources that goes into main matrix diagonal (vector of sources)
  real(dp), dimension(:), allocatable :: res  ! Residual vector for linear solvers

  ! Linear solution control parameters
  character( len = 20 ) :: solver ! Name of designated solver like 'iccg', etc.
  integer :: itr_max              ! Max no. of iterations
  integer :: tol_abs              ! Absolute tolerance level to be reached before exiting iterations.
  integer :: tol_rel              ! Relative tolerance level to be reached before exiting iterations. 
  real(dp) :: resor               ! Residual norm - e.g. L1 norm of initial residual

end type fvEquation



type, extends(csrMatrix) :: fvVectorEquation

  ! Solution vector component fields
  real(dp), dimension(:), allocatable :: u
  real(dp), dimension(:), allocatable :: v
  real(dp), dimension(:), allocatable :: w

  ! Past time values of solution vector field 
  real(dp), dimension(:), allocatable :: uo
  real(dp), dimension(:), allocatable :: vo
  real(dp), dimension(:), allocatable :: wo

  real(dp), dimension(:), allocatable :: uoo
  real(dp), dimension(:), allocatable :: voo
  real(dp), dimension(:), allocatable :: woo

  ! Source terms that go into main diagonal
  real(dp), dimension(:), allocatable :: spu
  real(dp), dimension(:), allocatable :: spv
  real(dp), dimension(:), allocatable :: spw

  ! Sources
  real(dp), dimension(:), allocatable :: su
  real(dp), dimension(:), allocatable :: sv
  real(dp), dimension(:), allocatable :: sw

  real(dp), dimension(:), allocatable :: res ! Residual vector

  ! Linear solution control parameters
  character( len = 20 ) :: solver ! Name of designated solver like 'iccg', etc.
  integer :: itr_max              ! Max no. of iterations
  integer :: tol_abs              ! Absolute tolerance level to be reached before exiting iterations.
  integer :: tol_rel              ! Relative tolerance level to be reached before exiting iterations. 

end type fvVectorEquation


interface operator(==)
   module procedure subtract_source_from_fvEquation  
   module procedure subtract_volVectorFieldSource_from_fvVectorEquation
   module procedure subtract_fvEquations
   module procedure subtract_fvVectorEquations
end interface


! Overload summation to be able to add fvEquation and volScalarField
! This enables to have fvi... routines which usually return fvEquation on rhs of == sign
interface operator(+)
   module procedure add_source_to_fvEquation  
   module procedure add_volVectorFieldSource_to_fvVectorEquation
   module procedure add_fvEquations
   module procedure add_fvVectorEquations
end interface


interface operator(-)
   module procedure subtract_source_from_fvEquation  
   module procedure subtract_volVectorFieldSource_from_fvVectorEquation
   module procedure subtract_fvEquations
   module procedure subtract_fvVectorEquations
end interface

public

contains


function new_fvEquation( ) result(fvEqn)
    implicit none
    integer :: i
    type(fvEquation) :: fvEqn

    allocate(fvEqn % ioffset ( numCells+1 ))
    allocate(fvEqn % diag ( numCells ))
    allocate(fvEqn % ja ( nnz ))
    allocate(fvEqn % a ( nnz ))

    allocate(fvEqn % su ( numCells ))
    allocate(fvEqn % sp ( numCells ))
    allocate(fvEqn % res ( numCells ))

    allocate(fvEqn % x ( numTotal ))
    allocate(fvEqn % o ( numTotal ))
    allocate(fvEqn % oo ( numTotal ))

    do i=1,numCells+1
      fvEqn % ioffset(i) = ioffset(i)
    enddo
    do i=1,numCells+1
      fvEqn % diag(i) = diag(i)
    enddo    
    do i=1,nnz
      fvEqn % ja(i) = ja(i)
    enddo

end function new_fvEquation

function new_fvVectorEquation( ) result(fvEqn)
    implicit none
    integer :: i
    type(fvVectorEquation) :: fvEqn

    allocate(fvEqn % ioffset ( numCells+1 ))
    allocate(fvEqn % ja ( nnz ))
    allocate(fvEqn % a ( nnz ))

    allocate(fvEqn % su ( numCells ))
    allocate(fvEqn % sv ( numCells ))
    allocate(fvEqn % sw ( numCells ))

    allocate(fvEqn % spu ( numCells ))
    allocate(fvEqn % spv ( numCells ))
    allocate(fvEqn % spw ( numCells ))

    allocate(fvEqn % res ( numCells ))

    allocate(fvEqn % u ( numTotal ))
    allocate(fvEqn % v ( numTotal ))
    allocate(fvEqn % w ( numTotal ))

    allocate(fvEqn % uo ( numTotal ))
    allocate(fvEqn % vo ( numTotal ))
    allocate(fvEqn % wo ( numTotal ))

    allocate(fvEqn % uoo ( numTotal ))
    allocate(fvEqn % voo ( numTotal ))
    allocate(fvEqn % woo ( numTotal ))

    do i=1,numCells+1
      fvEqn % ioffset(i) = ioffset(i)
    enddo
    do i=1,nnz
      fvEqn % ja(i) = ja(i)
    enddo

end function new_fvVectorEquation



function add_source_to_fvEquation(fvEqnIn,source) result( fvEqnOut )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqnIn
  type(volScalarField), intent(in) :: source
!
! > Result
!
  type(fvEquation) :: fvEqnOut

  integer :: i

  fvEqnOut = new_fvEquation()

  do i=1,numCells
      fvEqnOut % su(i) = fvEqnIn % su(i) + source % mag(i)
  enddo

end function add_source_to_fvEquation



function add_volVectorFieldSource_to_fvVectorEquation(fvEqnIn,vecSource) result( fvEqnOut )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqnIn
  type(volVectorField), intent(in) :: vecSource
!
! > Result
!
  type(fvVectorEquation) :: fvEqnOut

  integer :: i

  fvEqnOut = new_fvVectorEquation()

  do i=1,numCells
      fvEqnOut % su(i) = fvEqnIn % su(i) + vecSource % x(i)
      fvEqnOut % sv(i) = fvEqnIn % sv(i) + vecSource % y(i)
      fvEqnOut % sw(i) = fvEqnIn % sw(i) + vecSource % z(i)
  enddo

end function add_volVectorFieldSource_to_fvVectorEquation



function subtract_source_from_fvEquation(fvEqnIn,source) result( fvEqnOut )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqnIn
  type(volScalarField), intent(in) :: source
!
! > Result
!
  type(fvEquation) :: fvEqnOut

  integer :: i

  fvEqnOut = new_fvEquation()

  do i=1,numCells
      fvEqnOut%source( i ) = fvEqnIn%source( i ) - source%mag( i )
  enddo

end function subtract_source_from_fvEquation




function subtract_volVectorFieldSource_from_fvVectorEquation(fvEqn,vecSource) result( fvEqn_new )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqn
  type(volVectorField), intent(in) :: vecSource
!
! > Result
!
  type(fvVectorEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvVectorEquation()

  do i=1,numCells
      fvEqn_new % su(i) = fvEqn % su(i) - vecSource % x(i)
      fvEqn_new % sv(i) = fvEqn % sv(i) - vecSource % y(i)
      fvEqn_new % sw(i) = fvEqn % sw(i) - vecSource % z(i)
  enddo

end function subtract_volVectorFieldSource_from_fvVectorEquation


function add_fvEquations(fvEqn1,fvEqn2) result( fvEqn_new )
!
! Adds two objects of type(fvEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqn1
  type(fvEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvEquation()

  do i=1,numCells
      fvEqn_new % source(i) = fvEqn1 % source(i) + fvEqn2 % source(i)
  enddo
  do i=1,nnz
      fvEqn_new % coef(i) = fvEqn1 % coef(i) + fvEqn2 % coef(i)
  enddo

end function add_fvEquations


function add_fvVectorEquations(fvEqn1,fvEqn2) result( fvEqn_new )
!
! Adds two objects of type(fvVectorEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqn1
  type(fvVectorEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvVectorEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvVectorEquation()

  do i=1,numCells
      fvEqn_new % su(i) = fvEqn1 % su(i) + fvEqn2 % su(i)
      fvEqn_new % sv(i) = fvEqn1 % sv(i) + fvEqn2 % sv(i)
      fvEqn_new % sw(i) = fvEqn1 % sw(i) + fvEqn2 % sw(i)

      fvEqn_new % spu(i) = fvEqn1 % spu(i) + fvEqn2 % spu(i)
      fvEqn_new % spv(i) = fvEqn1 % spv(i) + fvEqn2 % spv(i)
      fvEqn_new % spw(i) = fvEqn1 % spw(i) + fvEqn2 % spw(i)
  enddo
  do i=1,nnz
      fvEqn_new % coef(i) = fvEqn1 % coef(i) + fvEqn2 % coef(i)
  enddo

end function add_fvVectorEquations


function subtract_fvEquations(fvEqn1,fvEqn2) result( fvEqn_new )
!
! Substracts two objects of type(fvEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqn1
  type(fvEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvEquation()

  do i=1,numCells
      fvEqn_new % source(i) = fvEqn1 % source(i) - fvEqn2 % source(i)
  enddo
  do i=1,nnz
      fvEqn_new % coef(i) = fvEqn1 % coef(i) - fvEqn2 % coef(i)
  enddo

end function subtract_fvEquations


function subtract_fvVectorEquations(fvEqn1,fvEqn2) result( fvEqn )
!
! Substracts two objects of type(fvVectorEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqn1
  type(fvVectorEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvVectorEquation) :: fvEqn

  integer :: i

  fvEqn = new_fvVectorEquation()

  do i=1,numCells
    fvEqn%su(i) = fvEqn1%su(i) - fvEqn2%su(i)
    fvEqn%sv(i) = fvEqn1%sv(i) - fvEqn2%sv(i)
    fvEqn%sw(i) = fvEqn1%sw(i) - fvEqn2%sw(i)

    fvEqn%spu(i) = fvEqn1%spu(i) - fvEqn2%spu(i)
    fvEqn%spv(i) = fvEqn1%spv(i) - fvEqn2%spv(i)
    fvEqn%spw(i) = fvEqn1%spw(i) - fvEqn2%spw(i)
  enddo
  do i=1,nnz
      fvEqn%coef(i) = fvEqn1%coef(i) - fvEqn2%coef(i)
  enddo

end function subtract_fvVectorEquations


subroutine solve_fvEqn( fvEqn )
!
! Purpose:
!   Solve linear system of equations resulting from finite volume discrretization.
!
! Description:
!   A version of the linear solver routine that uses finite volume eqaution derived type - fvEqn as an input
!   
  implicit none
!
! Parameters
!
  type fvEquation :: fvEqn 

  call spsolve(fvEqn%solver, fvEqn%phi%mag, fvEqn%su, fvEqn%itr_max, fvEqn%tol_abs, fvEqn%tol_rel, fvEqn%phi%field_name)

end subroutine

end module fvEquation
