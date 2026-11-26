subroutine modify_viscosity_inlet
!
! Purpose:
!   Eddy viscosity at inlet boundaries, for the model selected in 'turbulence_registry'.
!
  use turbulence_registry, only: turb_update_inlet

  implicit none

  if ( associated(turb_update_inlet) ) call turb_update_inlet ! otherwise laminar flow or algebraic SGS model

end subroutine
