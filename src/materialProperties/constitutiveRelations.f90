module constitutiveRelations

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing built-in linear solvers library.
!
!  Discussion:
!
!    Say something...
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

  private 

  public :: ..., 

contains

subroutine CarreauYassudaModel()
!
!  Purpose: 
!
!    This subroutine ...
!
!  Discussion:
!
!    Say something...
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    dd Month yyy
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    A. Ref, Title, Publication, pp.1-2, year.
!
!  Parameters:
!
!    Input, integer, name( size ), description.
!    Input/Output, type, name( size ), description.
!    Output, type, name( size ), description.
!
  use amodule, only: ...

  implicit none
!
! Parameters
!
  A list ...
!
! Local variables
!
  real(dp) :: mu0 = 22e-3
  real(dp) :: muinf = 2.2e-3
  real(dp) :: nexp = 0.392_dp
  real(dp) :: aexp = 0.644_dp
  real(dp) :: lambdcy = 0.110_dp

end subroutine

end module XXX