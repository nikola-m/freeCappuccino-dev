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
  use title_mod
  use fieldManipulation
  use sparse_matrix
  use temperature
  use concentration
  use mhd
  use utils, only: show_logo

  implicit none

  integer :: k
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

      ! Calculate velocities.
      if(lcal(iu)) call calcuvw  

      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(lcal(ip)) then

        if(SIMPLE)   call calcp_simple
        if(PISO)     call calcp_piso

      endif 

      ! Turbulence
      if(lturb)    call modify_viscosity

      !Scalars: Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calculate_temperature_field

      ! if(lcal(ivart)) call calculate_temperature_variance_field
      
      if(lcal(icon))  call calculate_concentration_field

      if(lcal(iep))  call calculate_electric_potential

     
      ! Log scaled residuals
      write(6,'(2x,a)') 'Scaled residuals:'
      write(6,'(2x,11(a,4x))')     (chvarSolver(k), k=1,nphi)
      write(6,'(a,11(1PE9.3,2x))') '> ',(resor(k), k=1,nphi)

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
      if( .not.ltransient .and. source.lt.sormax ) then
        call write_restart_files
        call writefiles
        exit time_loop
      end if

      if(ltransient) then 

        ! Has converged within timestep or has reached maximum no. of SIMPLE iterations per timetstep:
        if( source.lt.sormax .or. iter.ge.maxit ) then 

          ! Correct driving force for a constant mass flow rate simulation:
          if(const_mflux) then
            call constant_mass_flow_forcing
          endif

          ! Write values at monitoring points and recalculate time-average values for statistics:
          call writehistory
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