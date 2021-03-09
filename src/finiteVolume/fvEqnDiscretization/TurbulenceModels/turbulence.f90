module turbulence
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

! Number of turbulence scalar equations - determines the size of solution parameter arrays.
integer, parameter :: nsc = 8 

! To access solution parameters, fields need their number, 
! eddy viscosity is allways the last
integer, parameter :: ite=1
integer, parameter :: ied=2
integer, parameter :: inutild = 1
integer, parameter :: ivis=8

!
! Model and solution parameters - modified trough input.nml file
!
character (len=20) :: TurbModel = 'none'

logical :: SolveTKE = .False.
logical :: solveEpsilon = .False.
logical :: solveOmega = .False.
logical :: solveNuTilda = .False.

!
! We set default values for solution parmeters for all turbulence sclars, change them trough input.nml dictionary
!
real(dp), dimension(nsc) :: urf  = 0.5 ! Under-relaxation factors.
real(dp), dimension(nsc) :: gds = 1.0 ! Deferred correction factor.
character( len=30 ) :: cScheme = 'linearUpwind' ! Convection scheme - default is second order upwind.
character( len=12 ) :: dScheme = 'skewness' ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
integer :: nrelax = 0 ! Relaxation parameter non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed, 1/0/-1.
character( len=12 ) :: lSolver = 'bicgstab'  ! Linear algebraic solver.
integer :: maxiter = 5       ! Max number of iterations in linear solver.
real(dp) :: tolAbs = 1e-13   ! Absolute residual level.
real(dp) :: tolRel = 0.01    ! Relative drop in residual to exit linear solver.


public :: TurbModel, urf, gds, cScheme, dScheme, nrelax, lSolver, maxiter, tolAbs, tolRel ! params

contains


subroutine set_turb_scalar_solution_flags
!
! Parameters used in postprocessing phase
!

implicit none


  select case (TurbModel)

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

