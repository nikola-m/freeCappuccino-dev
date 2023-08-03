module TurbModelData
!
! Purpose:
!   This module contains the solution parameters that will be passed to turbulence models.
!
!   If we have e.g. two-equation turbulence model like k-epsilon, it uses
!   parameters such as under-relaxation factor 'urf', deferred correction parameter 'gam', etc.,
!   for two scalars k, and epsilon taking arrays values of these parameters at ite=1 and ied=2.
!
!   Spallart-Allmaras is one equation model, in that module 
!
!   Convection scheme 'cScheme', diffusion scheme 'dScheme', face-gradient non-orthogonal
!   correction relaxation parameter 'nrelax', and linear solver 'lSolver' have same setting for 
!   all turbulence scalar fields that are solved simultaneously, to make things less complicated.
!
use types

implicit none

! Napravi TurbModel koji ima i niz Skalara TurbModel%Scalar(1)%lSolver


type TurbScalar
  character (len=20) :: name = ' '                ! Scalar field identifyer
  real(dp) :: urf  = 0.7_dp                       ! Under-relaxation factors.
  real(dp) :: gds = 1.0_dp                        ! Deferred correction factor.
  character( len=30 ) :: cScheme = 'linearUpwind' ! Convection scheme - default is second order upwind.
  character( len=12 ) :: dScheme = 'skewness'     ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
  integer :: nrelax = 0                           ! Relaxation parameter non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed, 1/0/-1.
  character( len=12 ) :: lSolver = 'bicgstab'     ! Linear algebraic solver.
  integer :: maxiter = 10                         ! Max number of iterations in linear solver.
  real(dp) :: tolAbs = 1e-10                      ! Absolute residual level.
  real(dp) :: tolRel = 0.01_dp                    ! Relative drop in residual to exit linear solver.
  real(dp) :: resor                               ! Residual norm - e.g. L1 norm of initial residual
end type

!
! Example: type(Scalar) :: epsilon = Scalar('epsilon', 0.7, 1.0, 'linearUpwind', 'skewness', 0, 'bicgstab', 10, 1e-10, 0.01 )
!

type TurbulenceModel
  character (len=20) :: name = 'none' ! Turbulence model name - unique identifier.
  type(TurbScalar), dimension(4) :: Scalar ! Data for solution of discretized PDE for specific scalar
  real(dp) :: urfVis = 1.0_dp ! Under-relaxation facotr for eddy-viscosity (if the model is eddy-viscosity model)
end type


type( TurbulenceModel ) :: TurbModel

logical :: SolveTKE = .False.
logical :: solveEpsilon = .False.
logical :: solveOmega = .False.
logical :: solveNuTilda = .False.


public :: TurbModel ! params
public :: solveTKE, solveEpsilon, solveOmega, solveNuTilda

contains


subroutine set_turb_scalar_solution_flags
!
! Parameters used in postprocessing phase
!

implicit none


  select case ( TurbModel%name )

    case ('k_epsilon_std')
      solveTKE = .true.
      solveEpsilon = .true.

    case ('k_epsilon_rng')
        solveTKE = .true.
        solveEpsilon = .true.

    case ('k_omega_sst') 
        solveTKE = .true.
        solveOmega = .true.

    case ('Spalart_Allmaras')
        solveNuTilda = .true.

    case ('k_eqn_eddy')
        solveTKE = .true.

    case ('k_epsilon_std_2lewt')
      solveTKE = .true.
      solveEpsilon = .true.

    case ('k_epsilon_rlzb')
      solveTKE = .true.
      solveEpsilon = .true.

    case ('k_epsilon_rlzb_2lewt')
      solveTKE = .true.
      solveEpsilon = .true.

    case ('k_epsilon_zeta_f')
        solveTKE = .true.
        solveEpsilon = .true.

    case ('k_omega_EARSM_WJ') 
        solveTKE = .true.
        solveOmega = .true.

    case ('DDES_k_omega_sst')
        solveTKE = .true.
        solveOmega = .true.

    case ('IDDES_k_omega_sst')
        solveTKE = .true.
        solveOmega = .true.

    case default
      return
      
  end select

end subroutine



end module

