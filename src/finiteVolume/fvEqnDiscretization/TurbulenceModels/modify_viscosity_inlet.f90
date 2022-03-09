subroutine modify_viscosity_inlet

  use TurbModelData, only : TurbModel
  use k_epsilon_std
  use k_epsilon_std_2lewt  
  use k_epsilon_rng
  use k_epsilon_rlzb  
  use k_epsilon_rlzb_2lewt 
  use k_omega_SST
  use DDES_k_omega_SST
  use IDDES_k_omega_SST
  use k_eqn_eddy
  use spalart_allmaras
  ! use smagorinsky
  ! use wale_sgs
  
  implicit none

!+----------------------------------------------------------------------------+!

  select case ( TurbModel%name )

    case ('k_epsilon_std')
      call modify_viscosity_inlet_k_epsilon_std

    case ('k_epsilon_rng')
      call modify_viscosity_inlet_k_epsilon_rng

    case ('k_omega_sst') 
      call modify_viscosity_inlet_k_omega_sst

    case ('Spalart_Allmaras')
      call modify_viscosity_inlet_spalart_allmaras

    case ('k_eqn_eddy')
      call modify_viscosity_inlet_k_eqn_eddy

    case ('k_epsilon_std_2lewt')
      call modify_viscosity_inlet_k_epsilon_std_2lewt

    case ('k_epsilon_rlzb')
      call modify_viscosity_inlet_k_epsilon_rlzb

    case ('k_epsilon_rlzb_2lewt')
      call modify_viscosity_inlet_k_epsilon_rlzb_2lewt

    case ('DDES_k_omega_sst')
      call modify_viscosity_inlet_DDES_k_omega_sst

    case ('IDDES_k_omega_sst')
      call modify_viscosity_inlet_IDDES_k_omega_sst

    ! case ('Smagorinsky')
    !   call modify_viscosity_inlet_smagorinsky

    ! case ('WALE')
    !   call modify_viscosity_inlet_wale_sgs

    case default ! e.g. laminar flow
      return
      
  end select


end subroutine