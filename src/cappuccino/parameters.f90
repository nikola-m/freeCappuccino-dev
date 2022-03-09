module parameters  
  use types

  real(dp), parameter :: one = 1.0_dp
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: small = 1e-20
  real(dp), parameter :: great = 1e+20
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: twothirds = 2./3._dp
  real(dp), parameter :: onethird = 1./3._dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: slarge = 1e+30 

  ! Law of the wall parameters
  real(dp), parameter :: CAPPA = 0.41_dp 
  real(dp), parameter :: Rair = 287.058  ! J/(kg K)
  real(dp), parameter :: ELOG = 8.432_dp
  real(dp), parameter :: ctrans = 11.63


  character( len=10 ) :: mesh_format   ! Mesh format type
  integer :: pRefCell  ! Pressure reference cell
  integer :: mpoints   ! No. of monitoring points
  real(dp) :: flomas   ! mass flow at inlet
  real(dp) :: densit   ! Fluid density
  real(dp) :: viscos   ! Molecular dynamic viscosity
  real(dp) :: tolerance ! Residual toleance for SIMPLE
  real(dp) :: facnap   ! Underelaxation factor for Reynolds stresses
  real(dp) :: facflx   ! Underelaxation factor for turbulent heat fluxes
  real(dp) :: pranl     ! (= 0.7_dp for air, 7.0_dp for water, read it from input file.)
  logical :: lbuoy      ! Bouyancy effects - are they included in momentum and turbulence eqns. If yes we calculate heat fluxes.
  logical :: boussinesq ! Is Boussinesq hypothesis evoked yes=1/no=0, read from input file.
  real(dp) :: tref      ! Reference temperature, read from input file
  real(dp) :: beta      ! Thermal expansion coefficient, read from input file
  real(dp) :: gravx,gravy,gravz ! Components of gravity acceleration vector, read from input file
  real(dp) :: facvis            ! Under-relaxation for viscosity
  integer :: numstep            ! Total number of timesteps
  integer :: itime,itimes, itimee  ! Current timestep value, first and last timestep value in the current simulation
  integer :: nzapis    ! After how many timesteps we write backup and output files
  integer :: maxit     ! Maximum number of iterations in timestep, also max. number of iterations for SIMPLE iteration
  real(dp) :: timestep ! Timestep size
  real(dp) :: time     ! Total elapsed time
  logical :: CoNumFix         ! Is Courant no. fixed during time-stepping
  real(dp) :: CoNumFixValue   ! Fixed value for Courant number - set in modinp for now - may be read in input
  real(dp) :: CoNum ! Courant number (max value). 
  real(dp) :: meanCoNum  ! Courant number (mean value). 
  ! character(len=9) :: timechar! A char string to write current timestep

    
  ! Logicals, mostly read from simulation-input file:
  logical :: lturb,lread,lwrite,ltest   ! turbulent simulation, read restart file, write restart file, print residual of the linear solver,.,..      
  logical :: ltransient                 ! LTRANSIENT is TRUE for transient (non-stationary) simulations              
  logical :: lasm,lles,lsgdh,lggdh,lafm ! eddy-viscosity, algebraic stress model or LES, simple gradient or generalized gradient hypothesis, algerbaic flux model
  logical :: bdf,bdf2,bdf3,cn           ! control for the time-stepping algorithm
  logical :: simple,piso                ! control for the velocity-pressure coupling algorithm
  logical :: compressible               ! If true together with SIMPLE, then pressure based algorithm for all speeds is activated.
  logical :: const_mflux                ! control for constant flow rate 


  integer :: ncorr                   ! PISO control parameter: no. of Piso corrections.
  integer :: icorr                   ! PISO iteration no.: icorr=1..ncorr
  integer :: npcor                   ! No. of pressure-corrections; non-orthogonality correctors
  integer :: ipcorr                  ! Iteration no.: ipcorr=1..npcor
  integer :: nigrad = 1              ! No. of iters. for iterative cell-centered gradient calculation
  integer, parameter :: nipgrad = 2  ! No. of stages for 'pressure extrapolation at boundary + pressure gradient update' loop

  ! logical :: roughWall          ! Is aerodynamically rough wall assumed
  ! real(dp) :: erough            ! E parameter for rough wall
  ! real(dp) :: zzero             ! z0 - wall roughness [m] for aerodynamically rough walls

  real(dp) :: magUbar, gradPcmf ! Magnitude of the bulk velocity,and pressure grad that will drive the constant mass-flux flow (cmf)
  real(dp) :: sumLocalContErr, globalContErr, cumulativeContErr, res5Mass   ! Continuity errors

  real(dp), dimension(12) :: resor ! Normalized residuals for monitoring convergence of SIMPLE iterations.

  character(len=70) :: title         ! Case title-short description for monitor file.

  character(len=100) :: input_file,grid_file,monitor_file,restart_file

  
end module parameters