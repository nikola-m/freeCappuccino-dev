subroutine read_input_file
!
! Purpose:
!   Open, Read and Process Input File
!
! Reference:
!  Reading namelist file for input parameters is based on a code of Hiroaki Nishikawa NIA
!
  use types
  use parameters
  use gradients, only: lstsq, lstsq_qr, lstsq_dm, gauss, limiter
  use velocity, only: calcU, urfU, gdsU, cSchemeU, dSchemeU, nrelaxU, lSolverU, maxiterU, tolAbsU, tolRelU
  use pressure, only: calcP, urfP, lSolverP, maxiterP, tolAbsP, tolRelP
  use nablap, only: pscheme
  use temperature
  use energy
  use concentration
  use mhd
  use TurbModelData, only: TurbModel,set_turb_scalar_solution_flags
  use rheology, only: npow, Consst, shearmin, Tau_0, megp,  &
                      muplastic, muzero, muinfty, lamtime, acy,&    
                      calcVis, urfVis, non_newtonian_model
  use monitors


  implicit none

  integer :: os

  ! 
  ! > Parameters namelist for simulation control
  !
  namelist / input_parameters /   &
        title, &              ! Descriptive name of the case which will be written in monitor file.                             
        mesh_format, &        ! Mesh format 'nativeMesh' for native polyMesh format or 'foamMesh' for OpenFOAM polyMesh.
        lread, &              ! Read restart file - continue simulation from saved state?  True/False.
        ltest, &              ! Verbosity for linear solver convergence (for troubleshooting) True/False. 
        !
        ! Approximation details for velocity field.
        !
        calcU, &              ! Activate Velocity field caculation? True/False.
        urfU, &               ! Under-relaxation factor for momentum eq
        gdsU, &               ! gamma-deferred correction parameter for velocity
        cSchemeU, &           ! Convection scheme for momentun eq
        dSchemeU, &           ! Difussion scheme for momentum eq
        nrelaxU, &            ! Integer parameter for difussion scheme - advanced
        lSolverU, &           ! Linear algebraic solver for momentum eq
        maxiterU, &           ! Max no of iterations for linear U-V-W equations
        tolAbsU, &            ! Absolute tolerance level for residual for linear U-V-W equations
        tolRelU, &            ! Relative tolerance level for residual for linear U-V-W equations
        pscheme,&             ! pressure interpolation scheme 
        !
        ! Approximation details for pressure/pressure correction.
        !
        calcP, &              ! Activate Pressure field caculation? True/False.
        urfP, &               ! Under-relaxation factor for pressure
        lSolverP, &           ! Linear algebraic solver for pressure/pressure correction
        maxiterP, &           ! Max no of iterations for pressure/pressure correction
        tolAbsP, &            ! Absolute tolerance level for residual for pressure/pressure correction
        tolRelP, &            ! Relative tolerance level for residual for pressure/pressure correction  
        !
        ! Physical model: Viscous flows - Laminar/Turbulent, Activation and model details.
        !
        TurbModel, &          ! Turbulence model. Now a number, maybe char string in the future.    
        !
        ! Physical model: Rheology - Non-Newtonian fluids - model and solution details
        ! 
        calcVis, &            ! Activate dynamic viscosity recalculation (Temperature dependance, Non-Newtonian flows)
        non_newtonian_model, &! Self-explanatory
        npow, &               ! Exponent for power-law fluids
        Consst, &             ! Consistency index for power-law fluids
        shearmin, &           ! Lower limit for the shear rate magnitude for power-law fluids
        Tau_0, &              ! Yield stress
        megp,  &              ! Exponential growth parameter
        muplastic, &          ! Plastic viscosity
        muzero, &             ! Zero shear rate viscosity
        muinfty, &            ! Infinity shear rate viscosity
        lamtime, &            ! Natural time - time parameter
        acy, &                ! Exponent a in Carreau-Yasuda model
        urfVis, &             ! Under-relaxation parameter
        !
        ! Physical model: Heat transfer, Temperature as energy equation, Activation and model details.
        !
        calcT, &              ! Activate Temperature equation caculation? True/False.
        urfT, &              ! Under-relaxation factors.
        gdsT, &              ! Deferred correction factor.
        cSchemeT, &          ! Convection scheme - default is second order upwind.
        dSchemeT, &          ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
        nrelaxT, &           ! Relaxation parameter non-orthogonal correction for face gradient.
        lSolverT, &          ! Linear algebraic solver.
        maxiterT, &          ! Max number of iterations in linear solver.
        tolAbsT, &           ! Absolute residual level.
        tolRelT, &           ! Relative drop in residual to exit linear solver.
        sigt, &             ! sigma_t
        !
        ! Physical model: Energy equation, Activation and model details.
        !
        calcEn, &             ! Activate Energy equation caculation? True/False.
        urfEn, &              ! Under-relaxation factors.
        gdsEn, &              ! Deferred correction factor.
        cSchemeEn, &          ! Convection scheme - default is second order upwind.
        dSchemeEn, &          ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
        nrelaxEn, &           ! Relaxation parameter non-orthogonal correction for face gradient.
        lSolverEn, &          ! Linear algebraic solver.
        maxiterEn, &          ! Max number of iterations in linear solver.
        tolAbsEn, &           ! Absolute residual level.
        tolRelEn, &           ! Relative drop in residual to exit linear solver.
        sigtEn, &             ! sigma_t
        solveTotalEnergy, &   !@ What we solve here - default is total energy, change it 
        solveInternalEnergy, &!@ in input.nml file.
        solveEnthalpy, &      !@
        addViscDiss, &        ! Add viscous dissipation term? T/F
        !
        ! Physical model: Buoyancy 
        !
        lbuoy, &              ! Buoyancy activated True/False
        boussinesq, &         ! Bousinesq approximation for buoyancy
        tref, &               ! Reference temperature for buoyant flows
        gravx, &              !g Three components of gravity vector
        gravy, &              !g -
        gravz, &              !g -
        !
        ! Physical model: Dispersion of a passive scalar: Activation and solution details
        !
        calcCon, &            ! Activate passive scalar concentration field caculation? True/False.
        urfCon, &             ! Under-relaxation factors.
        gdsCon, &             ! Deferred correction factor.
        cSchemeCon, &         ! Convection scheme - default is second order upwind.
        dSchemeCon, &         ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
        nrelaxCon, &          ! Type of non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed. 1/0/-1.
        lSolverCon, &         ! Linear algebraic solver.
        maxiterCon, &         ! Max number of iterations in linear solver.
        tolAbsCon, &          ! Absolute residual level.
        tolRelCon, &          ! Relative drop in residual to exit linear solver.
        sigCon, &             ! Prandtl-Schmidt
        !
        ! Physical model: Magnetohydrodynamics using Electric potential: Activation and solution details
        !
        calcEpot, &           ! Activate Electric potential field caculation? True/False.
        urfEpot, &            ! Under-relaxation factor.
        gdsEpot, &            ! Deferred correction factor.
        lSolverEpot, &        ! Linear algebraic solver.
        maxiterEpot, &        ! Max number of iterations in linear solver.
        tolAbsEpot, &         ! Absolute residual level.
        tolRelEpot, &         ! Relative drop in residual to exit linear solver.
        !
        ! Definition of physical properties of fluid.
        !
        densit, &             ! Fluid density [kg/m3]
        viscos, &             ! Molecular dynamic viscosity [Pa.s]       
        pranl, &              ! Prandtl coefficient for specific fluid
        beta, &               ! Thermal expansion coefficient
        !
        ! Turbulent Heat flux model and underrelaxation for turb heat fluxes and RST.
        !
        lsgdh, &              !h Simple gradient hypothesis for heat-fluxes, or...
        lggdh, &              !h Generalized gradient hypothesis for heat fluxes, or...
        lafm, &               !h Algebraic flux modelling for heat fluxes.
        facflx, &             !h Underrelaxation factor for heat fluxes
        facnap, &             ! Underrelaxation factor for Reynolds stress tensor calculation.
        !
        ! Unsteady simulation
        !
        ltransient, &         !% |Unsteady simulation True/False and chose ONE algorithm  below
        bdf, &                !% |Backward-Euler; First-Order Implicit, or...
        bdf2, &               !% |Second-Order Backward Euler; Second-Order Implicit, or...
        bdf3, &               !% |Third-order backard, or...
        CN, &                 !% |Crank-Nicolson.
        !
        ! Numerical approximation: Gradient approximation and limiting
        !
        lstsq, &              !@ |Gradient approximation, chose ONE: Unweighted Least-Square gradient approximation, or...
        lstsq_qr, &           !@ |Least square gradients based on thin QR decomposition, or...
        lstsq_dm, &           !@ |Distance-square weighted version of Least-Square Gradient, or...
        gauss, &              !@ |Cell based Gauss gradient with simple linear interpolation.
        nigrad, &             ! Number of iterations for Gauss gradient approximation.   
        limiter, &            ! Gradient limiter - used for all fields
        !
        ! Solution method: pressure-velocity coupling
        !             
        SIMPLE, &             !# |Pressure-velocity coupling method - SIMPLE, or...
        PISO, &               !# |Pressure-velocity coupling method - PISO.
        compressible, &       !# SIMPLE for compressible flows. Activate Energy eqn.
        ncorr, &              ! Number of PISO corrections - only relevant for PISO.
        npcor, &              ! Number of iterations for pressure/pressure correction equation, i.e. Number of Nonorthogonal corrections.
        pRefCell, &           ! Reference cell for setting pressure level (since we have pure Neumann problem)
        tolerance, &          ! Lower tolerance bound for outer (SIMPLE/PISO) iterations
        !
        ! Simulation run details
        !
        numstep, &            ! Total number of timesteps or (SIMPLE steps if steady case).                   
        timestep, &           ! Timestep size  - also known as dt. 
        nzapis, &             ! Program writes output files every NZAPIS timesteps.
        maxit, &              ! Number of SIMPLE/PISO iterations in every timestep (SIMPLE steps if steady case and numstep=1).
        CoNumFix, &           !Co  |Adaptive timestep size based on target Courant number, True/False.
        CoNumFixValue,&       !Co  |If CoNumFix=True then set target maximum Courant number here.
        !
        ! Constant-mass flux: Activation and prescription of target bulk velocity.
        !
        const_mflux, &        !mdot | Do we have constant mass flow in the domain True/False.
        magUbar               !mdot | Target bulk velocity for constant mass flow situation.


!
! > Read the input parameters, defined in the file named as eg. 'input.nml'.
!

  write(*,*) "**************************************************************"
  write(*,*) " List of simulation control parameters and their values"
  write(*,*)

  open(unit=10,file=input_file,form='formatted',status='old',iostat=os)
  read(unit=10,nml=input_parameters)

  write(*,nml=input_parameters) ! Print the namelist variables.
  close(10)
  write(*,*)


  ! Turbulent flow computation condition:
  lturb =  .False.
  if( TurbModel%name /= 'none') then
    lturb = .True. 
    call  set_turb_scalar_solution_flags
  endif


  ! > Open files for data at monitoring points // this functionality is now based on function from 'monitors' module.
  call define_monitors


end subroutine