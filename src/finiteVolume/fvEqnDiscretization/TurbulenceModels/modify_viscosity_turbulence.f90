subroutine modify_viscosity_turbulence
!
! Purpose:
!   Called once per outer iteration (time step). Prepares the quantities every turbulence model needs
!   (velocity gradients, Reynolds-stress-like terms for post-processing, strain and vorticity)
!   and then calls the model selected in 'turbulence_registry'.
!
  use variables, only: u,v,w,dUdxi,dVdxi,dWdxi
  use gradients, only: grad_gauss
  use turbulence_registry, only: turb_update

  implicit none

  ! Velocity gradients: 
  call grad_gauss(U,dUdxi)
  call grad_gauss(V,dVdxi)
  call grad_gauss(W,dWdxi)

  call calcstress

  call calc_strain_and_vorticity

  if ( associated(turb_update) ) call turb_update ! otherwise laminar flow

end subroutine
