!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
program cappuccino
!                               
! 
!   __                _____                                 _             
!  / _|              / ____|                               (_)            
! | |_ _ __ ___  ___| |     __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___  
! |  _| '__/ _ \/ _ \ |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \.
! | | | | |  __/  __/ |___| (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |
! |_| |_|  \___|\___|\_____\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/ 
!                               | |   | |                                 
!                               |_|   |_|                                 
! '
!             MM          MMM
!            NMN.           MM
!           MMM..   OMM.    MMM.
!          .MMO.    MMM      MMM                         ...
!          MMM:     MMM      MMM                       .. .?DMD7..
!          MMM      MMM      MMM                      ..MMMMMMMMMMM~
!          MMM:     .M       MMM.                    .NMMMMMMMMMMMMMM
!          .MMM     .M.     .MMM.    ... .?MMMMM .. .MMMMMMMMMMMMMMMMM
!           MMM.    .M .    MMM. ....MMMMMMMMMMMM$..MMMMMM?     .7MMMMM.
!              MM.  M  M   MM  ..MMMMMMM~.. .MMMMM7MMMMMM.        7MMMMN
!                            =MMMMMM7..  ..MMM.MMMMMMMMM+         .MMMMM.
!                         DMMMMMD..  . . MMMM. .MMMMMMMM.        ..MMMMM.
!                    ..=MMMMMZ.     ...MMMM      MMMMMMMI.        .MMMM$.
!                    MMMMM8.   ..  .NMMMM..      :MMMMMMM         :MMMM.
!                 :MMMMM.......  ~MMMMN           MMMMMMMMO.      MMMM+
!             . ?MMMM:. ....  ,MMMMM              .MMMMMMMMM.    :MMMM.
!           ..IMMMM  .  .  =MMMMMI                ..MMMMMMMM.    MMMM=
!           .MMMM. .... DMMMMM?                     8MMMMM=     NMMMM
!          +MMM.   .~MMMMMM~                        .MMMM.     ,MMMM.
!         ~MM?.~MMMMMMM$                              MMMD    .MMMMM
!        .DMMMMMMM$                                   MMMM.  .MMMMM
!         .MMMM.                                      .MMMN..MMMMM,
!         ..MMMM.                                     .MMMM.MMMMMM
!           =MMMZ..   .=. ..         ..      =. ,:     .MMMMMMMMM
!           .MMMM~..  7  .MM.M   M   MM.    M.  ..  M  .MMMMMMMM+
!             MMMM....Z. .MM.M...M...MM..O:.  M . ~.M...ZMMMMMMM
!              MMMM,. ., :  ,:   :  :  :    .,  ,,  ,  . MMMMMMM
!               MMMM7. ..................................MMMMMM
!                MMMM8 ..................................MMMMMM
!                 NMMMM  ................................MMMMD
!                  7MMMM ............................... MMMM
!                    MMMMO............................. :MMM
!                     MMMMM..............................MMM
!                       MMMMM..........................MMMMM
!                        ZMMMMD.......................MMMMM
!                          MMMMMM.................. MMMMM
!                            MMMMMM. ............OMMMMMZ
!                              NMMMMM8.......=MMMMMMMI
!                                =MMMMMMMMMMMMMMMMM
!                                   NMMMMMMMMMM
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

    if (ltransient) then

      ! Update variables - shift in time: 
      call time_shift

      ! Courant number report:
      call CourantNo
    
    endif 


    iteration_loop: do iter=1,maxit 

      write(*,'(2x,2(a,i0))') 'Timestep: ',itime,', Iteration: ',iter

      call cpu_time(start)

      !--- Sequential equation solution ------------------------------------------------------------

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

      !--- END: Sequential equation solution -------------------------------------------------------

      !  Observe change in values during simulation
      call log_monitored_values

      call cpu_time(finish)
      write(timechar,'(f9.5)') finish-start
      write(6,'(3a)') '  ExecutionTime = ',adjustl(timechar),' s'
      write(6,*)

      !---------------------------------------------------------------
      !  Simulation management - residuals, loop control, output, etc.
      !---------------------------------------------------------------

      ! Check normalized residual and stop program if residuals diverge
      ! We check here only residual norms of u-mom, v-mom, w-mom 
      source = max( resor(1),resor(2),resor(3) ) 
      
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

      if (ltransient) then 

        ! Has converged within timestep or has reached maximum no. of SIMPLE iterations per timestep:
        if( source.lt.tolerance .or. iter.ge.maxit ) then 

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