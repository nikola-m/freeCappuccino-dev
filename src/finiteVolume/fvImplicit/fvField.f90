module fvField
!
! Definition of a datatype - a finite volume field, which is a composite type
! formed of discrete tensor field and a finite volume equation, along with
! various parameters pertinent to equation discretization and solution.
! fvField is usually a tensor field for which we need to write and solve 
! finite volume discretization of governing equation to obtain its solution.
!
use types
use geometry
use sparse_matrix
use tensorFields

implicit none

public

contains

!
! The fvEquation derived data type (all objects from fvImplicit belong to this type.)
!
type, extends(volScalarField) :: fvField
  real(dp), dimension(:,:), allocatable :: Grad ! Gradient field, eg. phi%Grad(1,inp)
  real(dp), dimension(:), allocatable :: o      ! n-1 (previous) timestep value, eg. phi%o(inp)
  real(dp), dimension(:), allocatable :: oo     ! n-2 timestep value, eg. phi%oo(inp)
  real(dp), dimension(:), allocatable :: aver   ! Average value, eg. phi%aver(inp)
  
  type(fvEquation) :: Eqn    ! eg. phi%Eqn; phi%Eqn%res(inp)
  real(dp) :: urf            ! Underrelaxation factor eg. phi%urf
  real(dp) :: gds            ! Gamma blending factor [0,1] for deffered correction
  character( len = 20 ) :: gradient_scheme
  character( len = 20 ) :: gradient_limiter
  character( len = 20 ) :: convection_scheme
  character( len = 20 ) :: diffusion_scheme

! Funkcije:
! -funkcija za alokaciju ili to moze u mainu...
! -func za init gde se inicijalnizuje polje i gde se podesava tip BC i inicijalne vrednosti boundary fejsova
! -func sa diskretizovanom jedn. koja se zove calc pa nesto, recimo calcScalar
! -func za flukseve koju poziva calcScalar
! -fun za eksplicitnu korekciju granicnih uslova posle solve, ali to mozda moze globalno?

end type fvField

type, extends(volScalarField) :: fvVectorField
  real(dp), dimension(:,:), allocatable :: xGrad ! Gradient field, eg. U%xGrad(1,inp)
  real(dp), dimension(:,:), allocatable :: yGrad ! Gradient field, eg. U%yGrad(1,inp)
  real(dp), dimension(:,:), allocatable :: zGrad ! Gradient field, eg. U%zGrad(1,inp)    
  real(dp), dimension(:), allocatable :: xo      ! n-1 (previous) timestep value, eg.U%xo(icell) 
  real(dp), dimension(:), allocatable :: yo      ! n-1 (previous) timestep value
  real(dp), dimension(:), allocatable :: zo      ! n-1 (previous) timestep value    
  real(dp), dimension(:), allocatable :: xoo     ! n-2 timestep value
  real(dp), dimension(:), allocatable :: yoo     ! n-2 timestep value  
  real(dp), dimension(:), allocatable :: zoo     ! n-2 timestep value  
  real(dp), dimension(:), allocatable :: xaver   ! Average value, eg. U%xAver(inp) 
  real(dp), dimension(:), allocatable :: yaver   ! Average value, eg. U%yAver(inp)
  real(dp), dimension(:), allocatable :: zaver   ! Average value, eg. U%zAver(inp)

  ! Finite volume equation and realted parameters
  type(fvEquation) :: Eqn                         ! eg. phi%Eqn%res(inp)
  real(dp) :: urf                                 ! Underrelaxation factor
  real(dp) :: gds                                 ! Gamma blending factor [0,1] for deffered correction 
  ! FVM Discretisation parameters
  character( len = 20 ) :: gradient_scheme
  character( len = 20 ) :: gradient_limiter
  character( len = 20 ) :: convection_scheme
  character( len = 20 ) :: diffusion_scheme


  ! integer :: n_sample ! broj semplovanih vrednosti za statitiku, eg. U%n_sample

! Funkcije:
! -funkcija za alokaciju ili to moze u mainu...
! -func za init gde se inicijalnizuje polje i gde se podesava tip BC i inicijalne vrednosti boundary fejsova
! -func sa diskretizovanom jedn. koja se zove calc pa nesto, recimo calcScalar
! -func za flukseve koju poziva calcScalar
! -fun za eksplicitnu korekciju granicnih uslova posle solve.

end type fvField

end module