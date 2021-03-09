!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
program cappuccino
!
!   __               _____                                  _             
!  / _|             /  __ \                                (_)            
! | |_ _ __ ___  ___| /  \/ __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___  
! |  _| '__/ _ \/ _ \ |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \ 
! | | | | |  __/  __/ \__/\ (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |
! |_| |_|  \___|\___|\____/\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/ 
!                               | |   | |                                 
!                               |_|   |_|                                 
! 
!
! Description:
!  A 3D unstructured finite volume solver for Computational Fluid Dynamics.   
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  use types
  use parameters
  use geometry
  use variables
  use fieldManipulation
  use sparse_matrix
  use velocity
  use pressure
  use temperature
  use energy
  use concentration
  use mhd
  use rheology
  use statistics
  use monitors
  use utils, only: show_logo

  implicit none

  integer :: iter
  integer :: itimes, itimee
  real(dp):: source
  real :: start, finish
  character( len = 9) :: timechar
!                                                                       
!******************************************************************************
!

!  Check command line arguments
  ! narg=command_argument_count()
  call get_command_argument(1,input_file)
  call get_command_argument(2,monitor_file)
  call get_command_argument(3,restart_file)

  ! Open simulation log file
  open(unit=6,file=monitor_file)
  rewind 6

  ! Print Cappuccino logo in the header of monitor file
  call show_logo


!-----------------------------------------------------------
!  Initialization, mesh definition, sparse matrix allocation
!-----------------------------------------------------------
  call read_input_file

  call read_mesh

  call create_CSR_matrix

  call create_fields

  call init

  !
  !===============================================
  ! Time loop: 
  !===============================================

  write(6,'(a)') ' '
  write(6,'(a)') '  Start iteration!'
  write(6,'(a)') ' '

  itimes = itime+1
  itimee = itime+numstep

  time_loop: do itime=itimes,itimee


    ! Update variables - shift in time: 
    call time_shift
   
    ! Set inlet boundary conditions
    if(itime.eq.itimes) call bcin

    ! Courant number report:
    call CourantNo


    iteration_loop: do iter=1,maxit 

      write(*,'(2(2x,a,i0))') 'Timestep: ',itime,',Iteration: ',iter

      call cpu_time(start)

      !---------------------------------------------------------------

      ! Calculate velocities.
      if(calcU) call calcuvw  

      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(calcP) then

        if(SIMPLE)   call calcp_simple
        if(PISO)     call calcp_piso

      endif 

      ! Non-Newtonian fluid - modification of modify_viscosity
      if( calcVis ) call modifyViscosityNonNewtonianFluid

      ! Turbulence
      if( lturb )    call modify_viscosity_turbulence

      ! Temperature or Total energy/Internal energy/Enhtalpy as energy variables
      if( calcT )   call calculate_temperature
      if( calcEn )  call calc_energy
      
      ! Concentration of the passive scalar
      if( calcCon )  call calculate_concentration_field

      ! Electric potential field for mhd flows.
      if( calcEpot )  call calculate_electric_potential

      !---------------------------------------------------------------

      !  Observe change in values during simulation
      call log_monitored_values

      call cpu_time(finish)
      write(timechar,'(f9.5)') finish-start
      write(6,'(3a)') '  ExecutionTime = ',adjustl(timechar),' s'
      write(6,*)

      !---------------------------------------------------------------
      !  Simulation management - residuals, loop control, output, etc.
      !---------------------------------------------------------------

      ! Check residual and stop program if residuals diverge
      source = max(resor(iu),resor(iv),resor(iw)) 

      if( source.gt.slarge ) then
        write(6,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
        stop
      endif

      ! If residuals fall to level below tolerance level - simulation is finished.
      if( .not.ltransient .and. source.lt.tolerance ) then

        call write_restart_files
        call writefiles
        write(6,"(//,10x,a)") "*** Successful end -  iterations converged ***" 
        exit time_loop

      end if

      if(ltransient) then 

        ! Has converged within timestep or has reached maximum no. of SIMPLE iterations per timetstep:
        if( source.lt.tolerance .or. iter.ge.maxit ) then 

          ! Correct driving force for a constant mass flow rate simulation:
          if(const_mflux) call constant_mass_flow_forcing

          ! Write values at monitoring points // comment out for now.
          ! call writehistory

          ! Recalculate time-average values for statistics:
          call calc_statistics 

          ! Write field values after nzapis iterations or at the end of time-dependent simulation:
          if( mod(itime,nzapis).eq.0  .or. itime.eq.numstep ) then

            call write_restart_files
            call writefiles

          endif

          cycle time_loop
       
        endif

      end if 

    end do iteration_loop

    ! Write field values after nzapis iterations or at the end of false-time-stepping simulation:
    if(.not.ltransient .and. ( mod(itime,nzapis).eq.0 .or. (itime-itimes+1).eq.numstep ) ) then

      call write_restart_files
      call writefiles
      
    endif

    if(ltransient) call flush(6)

  end do time_loop

end program