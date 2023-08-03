module pressure

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of algorithms for computation of pressure 
!    within specific pressure-velocity coupling approach.
!
!  Discussion:
!
!    The subroutines from this module are called after momentum predictor step which is common for both SIMPLE and PISO approach.
!    The approaches differ in corrector part. Implementation of SIMPLE here uses pressure-correction as a solution variable,
!    whereas PISO solves for pressure directly.
!    We use SIMPLE for stationary simulations, whereas for non-stationary simulations we use both SIMPLE and PISO.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types

  implicit none

  !
  ! Discrtetization and solution parameters - modified trough input.nml file
  !
  logical  :: calcP = .True.                     ! To activate the solution of this field in the main function. 
  real(dp) :: urfP = 0.2_dp                      ! Under-relaxation factor.
  character ( len = 60 ) :: lSolverP = 'iccg'      ! Linear algebraic solver.
  integer  :: maxiterP = 30                      ! Max number of iterations in linear solver.
  real(dp) :: tolAbsP = 1e-13                    ! Absolute residual level.
  real(dp) :: tolRelP = 0.025_dp                 ! Relative drop in residual to exit linear solver.

  private 

  public :: calcP, urfP, lSolverP, maxiterP, tolAbsP, tolRelP ! params

  public :: calcp_simple
  public :: calcp_piso


contains

 ! I will use simple include statements to include specific files, you can look into any of these for more details
 ! I use C style include with #include and not fortran style 'include' to avoid warnings related to preprocessor
 ! directives in the two files below.
#include "calcp_simple.f90" 

#include "calcp_piso.f90"


end module