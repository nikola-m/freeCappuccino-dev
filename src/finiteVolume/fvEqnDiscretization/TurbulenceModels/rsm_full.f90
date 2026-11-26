module rsm_full
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    High-Reynolds-number Reynolds stress transport model: SSG pressure-strain, epsilon equation, wall functions.
!
!  Discussion:
!
!    All the physics is in 'rsm_common' (stress equations, SSG, epsilon equation, momentum coupling through 'aij').
!    Here: alpha = 1 (no wall blending), wall-adjacent cells hold the log-layer stresses (from the wall shear stress),
!    epsilon in the wall cell is Cmu^3/4 k^3/2/(kappa y), wall shear stress from the standard wall function.
!    Mesh: first cell centre in the log layer, 30 < y+ < 300.
!    Initial/inlet stresses are isotropic, 2/3 k delta_ij, with k and epsilon read like for the k-epsilon models.
!
!  Reference:
!    [1] C.G. Speziale, S. Sarkar, T.B. Gatski, J. Fluid Mech. 227, 1991.   [2] Launder, Reece, Rodi, J. Fluid Mech. 68, 1975.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use rsm_common

  implicit none

  private

  public :: modify_viscosity_rsm_full
  public :: modify_viscosity_inlet_rsm_full

contains

subroutine modify_viscosity_rsm_full
  implicit none

  if ( .not.allocated(rs) ) call rsm_allocate( .false. )

  call rsm_shift_time
  call rsm_sources( .false. )
  call rsm_solve_stresses( .false. )
  call rsm_update_k
  call rsm_solve_epsilon( .false. )
  call rsm_effective_viscosity( .false. )

end subroutine

subroutine modify_viscosity_inlet_rsm_full
  implicit none
  call rsm_inlet
end subroutine

end module
