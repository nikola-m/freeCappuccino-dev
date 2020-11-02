subroutine modify_viscosity
  use parameters, only: TurbModel
  use variables, only: u,v,w,dUdxi,dVdxi,dWdxi
  use gradients, only: grad
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
  use smagorinsky

  implicit none

!+----------------------------------------------------------------------------+!


  ! Velocity gradients: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  call calcstress

  call calc_strain_and_vorticity

  select case (TurbModel)

    case(0)
      return

    case (1)
      call modify_viscosity_k_epsilon_std()

    case (2)
      call modify_viscosity_k_epsilon_rng()

    case (3) 
      LowRe = .false.
      call modify_viscosity_k_omega_sst()

    case (4)
      LowRe = .true.
      call modify_viscosity_k_omega_sst()

    case (5)
      call modify_viscosity_spalart_allmaras()

    case (6)
      call modify_viscosity_k_eqn_eddy()

    case (7)
      call modify_viscosity_k_epsilon_std_2lewt()

    case (8)
      call modify_viscosity_k_epsilon_rlzb()

    case (9)
      call modify_viscosity_k_epsilon_rlzb_2lewt()

    case (10)
      call modify_viscosity_DDES_k_omega_sst

    case (11)
      call modify_viscosity_IDDES_k_omega_sst

    case (12)
      call modify_viscosity_smagorinsky

    case default
      
  end select



end subroutine

