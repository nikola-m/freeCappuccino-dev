subroutine modify_viscosity_turbulence

  use TurbModelData, only: TurbModel
  use variables, only: u,v,w,dUdxi,dVdxi,dWdxi
  use gradients, only: grad_gauss
  use k_epsilon_std
  use k_epsilon_std_2lewt  
  use k_epsilon_rng
  use k_epsilon_rlzb  
  use k_epsilon_rlzb_2lewt 
  use k_epsilon_zeta_f
  use k_omega_SST
  use DDES_k_omega_SST
  use IDDES_k_omega_SST
  use k_eqn_eddy
  use spalart_allmaras
  use wale_sgs
  use vortexID_sgs
  use sigmaSGS
  use vremanSGS

  implicit none

!+----------------------------------------------------------------------------+!


  ! Velocity gradients: 
  call grad_gauss(U,dUdxi)
  call grad_gauss(V,dVdxi)
  call grad_gauss(W,dWdxi)

  call calcstress

  call calc_strain_and_vorticity

  select case ( TurbModel%name ) 

    case ('k_epsilon_std')
      call modify_viscosity_k_epsilon_std

    case ('k_epsilon_rng')
      call modify_viscosity_k_epsilon_rng

    case ('k_omega_sst') 
      call modify_viscosity_k_omega_sst

    case ('Spalart_Allmaras')
      call modify_viscosity_spalart_allmaras

    case ('k_eqn_eddy')
      call modify_viscosity_k_eqn_eddy

    case ('k_epsilon_std_2lewt')
      call modify_viscosity_k_epsilon_std_2lewt

    case ('k_epsilon_rlzb')
      call modify_viscosity_k_epsilon_rlzb

    case ('k_epsilon_rlzb_2lewt')
      call modify_viscosity_k_epsilon_rlzb_2lewt

    case ('k_epsilon_zeta_f')
      call modify_viscosity_k_epsilon_zeta_f

    case ('DDES_k_omega_sst')
      call modify_viscosity_DDES_k_omega_sst

    case ('IDDES_k_omega_sst')
      call modify_viscosity_IDDES_k_omega_sst

    case ('WALE')
      call modify_viscosity_wale_sgs

    case('vortexIDsgs')
      call modify_viscosity_vortexID_sgs

    case('sigmaSGS')
      call modify_viscosity_sigma_sgs

    case('vremanSGS')
      call modify_viscosity_vreman_sgs

    case default ! e.g. laminar flow
      return
      
  end select


end subroutine