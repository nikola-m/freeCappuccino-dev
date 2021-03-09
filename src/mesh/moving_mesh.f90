module moving_mesh

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing data and functions to support dynamic mesh flow solver.
!
!  Discussion:
!
!    What ever is the way you move and modify the mesh, fluid part needs to handle it. This module contains data and functions
!    to support it. 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    NOTE: STILL BEING DEVELOPED!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types

  implicit none

  logical :: DynamicMesh = .False.              ! Turn on/off dynamic mesh flow solution
  real(dp), dimension(:), allocatable :: flgrid ! Moving grid fluxes
  real(dp), dimension(:), allocatable :: volo   ! Last time-step cell volumes - maybe here maybe in geometry module?

  private 

  public :: 

contains

subroutine ZZZ()
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
  A list...


end subroutine

end module