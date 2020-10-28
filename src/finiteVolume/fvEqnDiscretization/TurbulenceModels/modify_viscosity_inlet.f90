subroutine modify_viscosity_inlet()
  use parameters, only: TurbModel
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

  implicit none

  select case (TurbModel)
    case (1)
      call modify_viscosity_inlet_k_epsilon_std()
    case (2)
      call modify_viscosity_inlet_k_epsilon_rng()
    case (3)
      call modify_viscosity_inlet_k_omega_sst()
    case (4)
      call modify_viscosity_inlet_k_omega_sst()
    case (5)
      call modify_viscosity_inlet_spalart_allmaras()
    case (6)
      call modify_viscosity_inlet_k_eqn_eddy()
    case (7)
      call modify_viscosity_inlet_k_epsilon_std_2lewt()    
    case (8)
      call modify_viscosity_inlet_k_epsilon_rlzb() 
    case (9)
      call modify_viscosity_inlet_k_epsilon_rlzb_2lewt()   
    case (10)
      call modify_viscosity_inlet_DDES_k_omega_sst
    case (11)
      call modify_viscosity_inlet_IDDES_k_omega_sst          
    case default
  end select

end subroutine