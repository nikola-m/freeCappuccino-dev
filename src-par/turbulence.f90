module turbulence
!
! Module that contains driver subroutines for turbulence models, 
! where each turbulence model is in separate module and is included here.
!

  use parameters, only: TurbModel
  use variables, only: u,v,w,dudxi,dvdxi,dwdxi
  use gradients
  use k_epsilon_std
  use k_epsilon_rng
  use k_omega_sst
  use k_eqn_eddy
  use spalart_allmaras


  implicit none

  public


contains


subroutine correct_turbulence()

  implicit none

  ! Velocity gradients: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  call calcstress
  
  call calc_strain_and_vorticity

  select case (TurbModel)

    case (1)
      call correct_turbulence_k_epsilon_std()
    case (2)
      call correct_turbulence_k_epsilon_rng()
    case (3)
      LowRe = .false.
      call correct_turbulence_k_omega_sst()
    case (4)
      LowRe = .true.
      call correct_turbulence_k_omega_sst()
    case (5)
      call correct_turbulence_spalart_allmaras()
    case (6)
      call correct_turbulence_k_eqn_eddy()

    case default
      
  end select

end subroutine


subroutine correct_turbulence_inlet()

  implicit none

  select case (TurbModel)
    case (1)
      call correct_turbulence_inlet_k_epsilon_std()
    case (2)
      call correct_turbulence_inlet_k_epsilon_rng()
    case (3)
      call correct_turbulence_inlet_k_omega_sst()
    case (4)
      call correct_turbulence_inlet_k_omega_sst()
    case (5)
      call correct_turbulence_inlet_spalart_allmaras()
    case (6)
      call correct_turbulence_inlet_k_eqn_eddy()
    case default
  end select

end subroutine


end module