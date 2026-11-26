module turbulence_registry
!
! Purpose:
!   Single place where a turbulence model is selected by name.
!
!   Before: the model name was compared (select case) in three places - the flag setter in
!   TurbModelData, 'modify_viscosity_turbulence' (every outer iteration) and 'modify_viscosity_inlet'.
!   Adding a model meant editing three files and was easy to get half-done.
!
!   Now: 'select_turbulence_model' is called once after the input file is read. It
!     - points the procedure pointers 'turb_update' and 'turb_update_inlet' to the model routines,
!     - sets the solution flags solveTKE, solveEpsilon, solveOmega, solveNuTilda and the coupling flag
!       'lasm' (excess anisotropy 'aij' of an algebraic or transported Reynolds stress model is passed to the
!       momentum equations as an explicit source, see asm_stress_terms). Unknown names stop the run.
!   A new model is added by one 'case' block below.
!
use types, only: dp
use parameters, only: lasm
use TurbModelData
use k_epsilon_std
use k_epsilon_std_2lewt
use k_epsilon_rng
use k_epsilon_rlzb
use k_epsilon_rlzb_2lewt
use k_epsilon_zeta_f
use k_omega_SST
use transition_sst
use easm_rumsey_gatski_komega
use earsm_gatski_speziale
use rsm_full
use ebrsm
use sgs_one_equation_dynamic
use DDES_k_omega_SST
use IDDES_k_omega_SST
use k_eqn_eddy
use spalart_allmaras
use wale_sgs
use vortexID_sgs
use sigmaSGS
use vremanSGS

implicit none

abstract interface
  subroutine turb_proc_interface
  end subroutine
end interface

! Model entry points, set by select_turbulence_model.
procedure(turb_proc_interface), pointer :: turb_update       => null() ! Solve model eqs. and update eddy viscosity
procedure(turb_proc_interface), pointer :: turb_update_inlet => null() ! Eddy viscosity at inlet boundaries

private
public :: turb_update, turb_update_inlet, select_turbulence_model

contains

subroutine select_turbulence_model
  implicit none

  select case ( TurbModel%name )

    case ('none')
      turb_update       => null()
      turb_update_inlet => null()

    case ('k_epsilon_std')
      turb_update       => modify_viscosity_k_epsilon_std
      turb_update_inlet => modify_viscosity_inlet_k_epsilon_std
      solveTKE = .true.;  solveEpsilon = .true.

    case ('k_epsilon_rng')
      turb_update       => modify_viscosity_k_epsilon_rng
      turb_update_inlet => modify_viscosity_inlet_k_epsilon_rng
      solveTKE = .true.;  solveEpsilon = .true.

    case ('k_epsilon_std_2lewt')
      turb_update       => modify_viscosity_k_epsilon_std_2lewt
      turb_update_inlet => modify_viscosity_inlet_k_epsilon_std_2lewt
      solveTKE = .true.;  solveEpsilon = .true.

    case ('k_epsilon_rlzb')
      turb_update       => modify_viscosity_k_epsilon_rlzb
      turb_update_inlet => modify_viscosity_inlet_k_epsilon_rlzb
      solveTKE = .true.;  solveEpsilon = .true.

    case ('k_epsilon_rlzb_2lewt')
      turb_update       => modify_viscosity_k_epsilon_rlzb_2lewt
      turb_update_inlet => modify_viscosity_inlet_k_epsilon_rlzb_2lewt
      solveTKE = .true.;  solveEpsilon = .true.

    case ('k_epsilon_zeta_f')
      turb_update       => modify_viscosity_k_epsilon_zeta_f
      turb_update_inlet => modify_viscosity_inlet_k_epsilon_zeta_f
      solveTKE = .true.;  solveEpsilon = .true.

    case ('k_omega_sst')
      turb_update       => modify_viscosity_k_omega_sst
      turb_update_inlet => modify_viscosity_inlet_k_omega_sst
      solveTKE = .true.;  solveOmega = .true.

    case ('transition_sst')
      turb_update       => modify_viscosity_transition_sst
      turb_update_inlet => modify_viscosity_inlet_transition_sst
      solveTKE = .true.;  solveOmega = .true.

    case ('easm_rumsey_gatski_komega')
      turb_update       => modify_viscosity_easm_rumsey_gatski_komega
      turb_update_inlet => modify_viscosity_inlet_easm_rumsey_gatski_komega
      solveTKE = .true.;  solveOmega = .true.
      lasm = .true.

    case ('earsm_gatski_speziale')
      turb_update       => modify_viscosity_earsm_gatski_speziale
      turb_update_inlet => modify_viscosity_inlet_earsm_gatski_speziale
      solveTKE = .true.;  solveEpsilon = .true.
      lasm = .true.

    case ('rsm_full')
      turb_update       => modify_viscosity_rsm_full
      turb_update_inlet => modify_viscosity_inlet_rsm_full
      solveTKE = .true.;  solveEpsilon = .true.
      lasm = .true.

    case ('ebrsm')
      turb_update       => modify_viscosity_ebrsm
      turb_update_inlet => modify_viscosity_inlet_ebrsm
      solveTKE = .true.;  solveEpsilon = .true.
      lasm = .true.

    case ('sgs_one_equation_dynamic')
      turb_update       => modify_viscosity_sgs_one_equation_dynamic
      turb_update_inlet => modify_viscosity_inlet_sgs_one_equation_dynamic
      solveTKE = .true.

    case ('DDES_k_omega_sst')
      turb_update       => modify_viscosity_DDES_k_omega_sst
      turb_update_inlet => modify_viscosity_inlet_DDES_k_omega_sst
      solveTKE = .true.;  solveOmega = .true.

    case ('IDDES_k_omega_sst')
      turb_update       => modify_viscosity_IDDES_k_omega_sst
      turb_update_inlet => modify_viscosity_inlet_IDDES_k_omega_sst
      solveTKE = .true.;  solveOmega = .true.

    case ('Spalart_Allmaras')
      turb_update       => modify_viscosity_spalart_allmaras
      turb_update_inlet => modify_viscosity_inlet_spalart_allmaras
      solveNuTilda = .true.

    case ('k_eqn_eddy')
      turb_update       => modify_viscosity_k_eqn_eddy
      turb_update_inlet => modify_viscosity_inlet_k_eqn_eddy
      solveTKE = .true.

    ! Algebraic SGS models - no transport equations, no inlet treatment.
    case ('WALE')
      turb_update => modify_viscosity_wale_sgs
    case ('vortexIDsgs')
      turb_update => modify_viscosity_vortexID_sgs
    case ('sigmaSGS')
      turb_update => modify_viscosity_sigma_sgs
    case ('vremanSGS')
      turb_update => modify_viscosity_vreman_sgs

    case default
      write(*,'(3a)') ' *** Unknown turbulence model: "', trim(TurbModel%name), '". See turbulence_registry.f90 for the list.'
      stop

  end select

end subroutine

end module
