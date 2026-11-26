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
  use gradients, only: lstsq, lstsq_qr, lstsq_dm, lstsq_fw, gauss, node_gauss, limiter
  use velocity, only: calcU, urfU, gdsU, cSchemeU, dSchemeU, nrelaxU, lSolverU, maxiterU, tolAbsU, tolRelU
  use pressure, only: calcP, urfP, lSolverP, maxiterP, tolAbsP, tolRelP
  use nablap, only: pscheme
  use implicit_bc_compressible, only: cplbc_name, cplbc_kind, cplbc_p0, cplbc_T0, cplbc_mdot, cplbc_dirx, cplbc_diry, &
                                      cplbc_dirz, cplbc_urf
  use simple_allspeeds, only: asimple, asimple_urfP, asimple_beta, asimple_lSolverP, asimple_maxitP, &
                              asimple_tolRelP, asimple_tolAbsP, asimple_pseudoCFL, asimple_cflRamp, asimple_wallFunc, &
                              pamg_theta, pamg_coarse_size, &
                              pamg_nu1, pamg_nu2, pamg_max_levels, pamg_restart, pamg_verbose
  use coupled_pU, only: coupled, cpl_maxit, cpl_restart, cpl_tolRel, cpl_tolAbs, cpl_urfP, cpl_rebuild, cpl_contScale, &
                        cpl_rhoBlend, cpl_pseudoCFL, cpl_cflRamp, &
                        cpl_sutherland, cpl_Tref, cpl_Ssuth, &
                        amg_verbose, amg_smoother, amg_cycle, amg_theta, amg_coarse_size, amg_nu1, amg_nu2, amg_ccf
  use temperature
  use energy_equation, only: calcEn, urfEn, gdsEn, cSchemeEn, limiterEn, dSchemeEn, nrelaxEn, lSolverEn, maxiterEn, &
                             tolAbsEn, tolRelEn, sigtEn, solveEnthalpy, solveInternalEnergy, solveTotalEnergy, &
                             solveTotalEnthalpy, addViscDiss, conservativeEn, cpModel, cp_const, cp_Tref, cp_poly, &
                             nasa_lo, nasa_hi, nasa_R, shomate, shomate_M, sutherland, mu_Tref, mu_S
  use concentration
  use mhd
  use TurbModelData, only: TurbModel
  use turbulence_registry, only: select_turbulence_model
  use rheology, only: npow, Consst, shearmin, Tau_0, megp,  &
                      muplastic, muzero, muinfty, lamtime, acy,&    
                      calcVis, urfVis, non_newtonian_model
  use monitors
  use clsvof, only: twoPhase, cl_rho1, cl_rho2, cl_mu1, cl_mu2, cl_sigma, cl_eps, cl_band, cl_nReinit, cl_cfl, cl_fmin, &
                    cl_ab2, cl_harmonic, cl_ic, cl_inside, cl_center, cl_radius, cl_flow, cl_vel, cl_omega, cl_vcenter, &
                    cl_T, cl_inflowF
  use timescale, only: autotime, lenscale_option
  use moving_mesh, only: DynamicMesh, mm_motion, mm_patch, mm_amp, mm_freq, mm_vel, mm_center, mm_axis, mm_omega, &
                         mm_rbf, mm_radius, mm_epsilon


  implicit none

  integer :: os

  ! 
  ! > Parameters namelist for simulation control
  !
  namelist / input_parameters /   &
        title, &              ! Descriptive name of the case which will be written in monitor file.                             
        mesh_format, &        ! Mesh format 'nativeMesh' for native polyMesh format or 'foamMesh' for OpenFOAM polyMesh.
        lread, &              ! Read restart file - continue simulation from saved state?  True/False.
        lwrite,&              ! Control writing of the certain outut files. True/False.
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
        TurbModel, &          ! Turbulence model.   
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
        solveTotalEnthalpy, & !@ What we solve here - default is total enthalpy, change it
        solveTotalEnergy, &   !@ in input.nml file (priority: enthalpy, internal energy, total energy, total enthalpy).
        solveInternalEnergy, &!@
        solveEnthalpy, &      !@
        addViscDiss, &        ! Add viscous terms (dissipation / work of viscous stresses)? T/F
        limiterEn, &          ! Limiter of the MUSCL scheme: none, vanLeer, minmod, superbee, vanAlbada.
        conservativeEn, &     ! Conservative form of the unsteady and convective terms? T/F
        cpModel, &            ! Specific heat model: constant, polynomial, nasa, shomate.
        cp_const, &           ! Specific heat at constant pressure [J/kgK] (cpModel = constant).
        cp_Tref, &            ! Reference temperature of the enthalpy [K].
        cp_poly, &            ! Polynomial cp = c0 + c1 T + c2 T^2 + c3 T^3.
        nasa_lo, nasa_hi, &   ! NASA 7-coefficient polynomials (T < 1000 K, T >= 1000 K).
        nasa_R, &             ! Specific gas constant of the NASA polynomials.
        shomate, shomate_M, & ! Shomate equation coefficients A..F, molar mass [kg/mol].
        sutherland, &         ! Sutherland viscosity law (also cpl_sutherland of the coupled solver).
        mu_Tref, mu_S, &      ! Sutherland reference temperature and temperature S [K].
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
        autotime, &           !# Automatic timestep for pseudotransient runs,
        lenscale_option, &    !# length scale option for pseudotransient runs.
        !
        ! Numerical approximation: Gradient approximation and limiting
        !
        lstsq, &              !@ |Gradient approximation, chose ONE: Unweighted Least-Square gradient approximation, or...
        lstsq_qr, &           !@ |Least square gradients based on thin QR decomposition, or...
        lstsq_dm, &           !@ |Distance-square weighted version of Least-Square Gradient, or...
        lstsq_fw, &           !@ |Distance-square weighted LSQ with precomputed face weights (any cell shape), or...
        node_gauss, &         !@ |Node based Gauss gradient using pseudolaplacian interpolation, or...
        gauss, &              !@ |Cell based Gauss gradient with simple linear interpolation.
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
        !
        ! Alternative all-speeds SIMPLE (module simple_allspeeds) and its AMG for the pressure correction
        !
        asimple, &            !# |All-speeds SIMPLE (Przulj-Basara type) instead of SIMPLE/PISO, True/False.
        asimple_urfP, &       ! Under-relaxation of the pressure update.
        asimple_beta, &       ! Weight of the central part of the face density (0 = upwind density).
        asimple_lSolverP, &   ! Solver for p': 'amg' (own AMG+FGMRES) or any name of csrsolve.
        asimple_maxitP, &     ! Max. no. of iterations of the p' solver.
        asimple_tolRelP, &    ! Relative residual drop of the p' solver.
        asimple_tolAbsP, &    ! Absolute residual level of the p' solver.
        asimple_pseudoCFL, &  ! Steady compressible: pseudo-time CFL of the density term in the p' equation (0 = off).
        asimple_cflRamp, &    ! Iterations of the geometric ramp of the pseudo CFL from 1 to asimple_pseudoCFL.
        asimple_wallFunc, &    ! Wall function: none, standard, reichardt, spl, compressible.
        pamg_theta, &          ! AMG: strength-of-connection threshold.
        pamg_coarse_size, &    ! AMG: direct solution if the number of unknowns is <= this.
        pamg_nu1, &            ! AMG: pre-smoothing sweeps.
        pamg_nu2, &            ! AMG: post-smoothing sweeps.
        pamg_max_levels, &     ! AMG: max. no. of levels.
        pamg_restart, &        ! AMG: FGMRES restart length.
        pamg_verbose, &        ! AMG: print hierarchy.
        !
        ! Block-coupled pressure-velocity solver and block AMG (used if coupled=T)
        !
        coupled, &            !# |Block-coupled p-U solver instead of SIMPLE/PISO, True/False.
        cpl_maxit, &          ! Max. no. of FGMRES iterations per outer iteration.
        cpl_restart, &        ! FGMRES restart length.
        cpl_tolRel, &         ! Relative residual reduction of the linear solver per outer iteration.
        cpl_tolAbs, &         ! Absolute residual level of the linear solver.
        cpl_urfP, &           ! Under-relaxation of the pressure update (1 = none).
        cpl_rebuild, &        ! Rebuild AMG aggregates every cpl_rebuild outer iterations.
        cpl_contScale, &      ! Scaling of the continuity rows.
        cpl_rhoBlend, &       ! All speeds: weight of the deferred central density correction (0 = upwind density).
        cpl_pseudoCFL, &      ! All speeds: pseudo-time CFL in the continuity equation (0 = off).
        cpl_cflRamp, &        ! All speeds: iterations of the geometric ramp of the pseudo CFL from 1 to cpl_pseudoCFL.
        cpl_sutherland, &     ! All speeds: Sutherland viscosity law (viscos = value at cpl_Tref).
        cpl_Tref, &           ! All speeds: Sutherland reference temperature.
        cpl_Ssuth, &          ! All speeds: Sutherland temperature S [K].
        cplbc_name, &         ! All speeds: patch names with a compressible boundary condition (implicit_bc_compressible).
        cplbc_kind, &         ! Kinds: compinlet, massflow, totalpressure, stagnation, riemann, compoutlet, open.
        cplbc_p0, &           ! Total pressure (static pressure for compoutlet/open if > 0) [Pa].
        cplbc_T0, &           ! Total temperature [K].
        cplbc_mdot, &         ! Mass flow rate through the patch [kg/s] (massflow).
        cplbc_dirx, cplbc_diry, cplbc_dirz, & ! Inflow direction (default: opposite to the patch normal).
        cplbc_urf, &          ! Under-relaxation of the boundary pressure update.
        amg_smoother, &       ! AMG smoother: 'ilu0' or 'gs'.
        amg_cycle, &          ! AMG cycle: 1=V, 2=W.
        amg_theta, &          ! AMG strength-of-connection threshold.
        amg_coarse_size, &    ! AMG: size (cells) of the coarsest level solved by dense LU.
        amg_nu1, &            ! AMG: pre-smoothing sweeps.
        amg_nu2, &            ! AMG: post-smoothing sweeps.
        amg_ccf, &            ! AMG: coarse-grid correction factor.
        amg_verbose, &        ! AMG: print hierarchy information (0/1).
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
        magUbar, &            !mdot | Target bulk velocity for constant mass flow situation.
        !
        ! Two-phase flow: CLSVOF (coupled level set and volume of fluid) method
        !
        twoPhase, &           ! Activate the two-phase solver (clsvof module) True/False.
        cl_rho1, cl_rho2, &   ! Densities of phase 1 (vof=1, lsf>0) and phase 2 [kg/m3].
        cl_mu1, cl_mu2, &     ! Dynamic viscosities of phase 1 and phase 2 [Pa s].
        cl_sigma, &           ! Surface tension coefficient [N/m].
        cl_eps, &             ! Half-width of the smoothed Heaviside function in cell sizes.
        cl_band, &            ! Half-width of the narrow band with the distance function [cells].
        cl_nReinit, &         ! Sweeps of PDE reinitialization of the level set outside of the narrow band.
        cl_cfl, &             ! Max Courant number of one substep of interface advection.
        cl_fmin, &            ! Volume fractions below cl_fmin (above 1-cl_fmin) are treated as 0 (1).
        cl_ab2, &             ! Adams-Bashforth extrapolation of the volume flux used for the advection.
        cl_harmonic, &        ! Harmonic mixing of viscosities.
        cl_ic, &              ! Initial interface: sphere, cylinder, plane.
        cl_inside, &          ! Phase inside of the sphere/cylinder: 1 (drop) or 2 (bubble).
        cl_center, &          ! Center of the initial sphere/cylinder (x,y,z).
        cl_radius, &          ! Radius of the initial sphere/cylinder.
        cl_flow, &            ! solve, or prescribed flow for advection tests: translate, rotate, vortex.
        cl_vel, &             ! Velocity of the translation test.
        cl_omega, &           ! Angular velocity of the rotation test.
        cl_vcenter, &         ! Center of the rotation test (x,y).
        cl_T, &               ! Period of the reversing vortex test.
        cl_inflowF, &         ! Volume fraction entering by open boundaries (advection tests only).
        !
        ! Moving (dynamic) mesh: deformation by radial basis functions (RBF), prescribed motion of the wall patches
        !
        DynamicMesh, &        ! Activate the moving mesh True/False.
        mm_patch, &           ! Names of the wall patches that move.
        mm_motion, &          ! Motion of the patches: oscillate, translate, rotate, user (edit mm_user_displacement).
        mm_amp, mm_freq, &    ! oscillate: amplitude vector and frequency, d = amp*sin(2 pi freq t).
        mm_vel, &             ! translate: velocity vector.
        mm_center, mm_axis, & ! rotate: point on the axis and direction of the axis.
        mm_omega, &           ! rotate: angular velocity [rad/s].
        mm_rbf, &             ! RBF: wendland, gaussian, multiquadric, inverse_multiquadric, thin_plate.
        mm_radius, &          ! wendland: support radius (<= 0: diagonal of the bounding box).
        mm_epsilon, &         ! Shape parameter of gaussian, multiquadric, inverse_multiquadric.
        !
        ! Reference values for reporting forces, fluxes, and coefficients.
        !
        rhoref, &  ! Reference density
        Uref, &    ! Reference velocity
        Aref, &    ! Reference area 
        pref, &    ! Reference pressure 
        lref       ! Reference length


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
    call select_turbulence_model
  endif


  ! > Open files for data at monitoring points // this functionality is now based on function from 'monitors' module.
  call define_monitors


end subroutine