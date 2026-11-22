module allspeeds
!
!  Purpose:
!
!    Thermodynamic state of the all-speeds (compressible) block-coupled p-U solver. Thin layer between the coupled solver
!    and the energy equation module (energy_equation, ScalarEqns/energy.f90), which owns the energy equation, the
!    thermodynamic properties and the time levels (density, pressure, energy variable, BDF coefficients).
!
!    * ideal gas:   rho = p / (R T),   R = Rair (parameters),   cp = cp(T) (energy_equation: cpModel)
!    * time levels: shifted once per time step (energy_time_levels), BDF1/BDF2 coefficients tc0, tc1, tc2
!    * the energy equation is solved after the coupled p-U step by calc_energy (energy_equation), see there for
!      the formulations (total enthalpy by default), convection schemes and boundary conditions
!    * optional Sutherland viscosity law (namelist: cpl_sutherland, cpl_Tref, cpl_Ssuth)
!
  use types
  use parameters
  use variables
  use energy_equation, only: energy_init, energy_time_levels, update_thermodynamic_properties, speed_of_sound_sq, &
                             den_o, den_oo, p_o, p_oo, tc0, tc1, tc2, energy_ready, &
                             sutherland, mu_Tref, mu_S

  implicit none

  private

  public :: as_init, as_time_levels, update_density, update_viscosity, speed_of_sound_sq
  public :: sutherland, mu_Tref, mu_S
  public :: den_o, den_oo, p_o, p_oo, tc0, tc1, tc2, as_ready

  logical, save :: as_ready = .false.

contains


subroutine update_density
!
!  Ideal gas density from the current pressure and temperature, cells and boundary faces.
!
  implicit none
  den = p/( Rair*T + small )
end subroutine


subroutine update_viscosity
!
!  Sutherland's law  mu = mu_ref (T/T_ref)^(3/2) (T_ref + S)/(T + S), cells and boundary faces. Laminar flow: vis = mu.
!
  implicit none
  if ( .not. sutherland ) return
  vis = viscos*( T/mu_Tref )**1.5_dp*( mu_Tref + mu_S )/( T + mu_S )
end subroutine


subroutine as_init
!
!  Called once at the first iteration: make the state consistent (density, viscosity, properties, old levels).
!
  implicit none

  call update_density
  call update_viscosity
  call energy_init

  as_ready = .true.

end subroutine


subroutine as_time_levels
!
!  Once per time step (detected by the change of itime). See energy_time_levels.
!
  implicit none

  call energy_time_levels

end subroutine

end module allspeeds
