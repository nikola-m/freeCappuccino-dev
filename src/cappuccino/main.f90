!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
program cappuccino
!
!******************************************************************************
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

  implicit none

  integer :: iter
  integer :: narg
  integer :: itimes, itimee
  real(dp):: source
  real :: start, finish
!                                                                       
!******************************************************************************
!

!  Check command line arguments
  narg=command_argument_count()
  call get_command_argument(1,input_file)
  call get_command_argument(2,monitor_file)
  call get_command_argument(3,restart_file)
  call get_command_argument(4,out_folder_path)

!-----------------------------------------------------------
!  Initialization, mesh definition, sparse matrix allocation
!-----------------------------------------------------------
  call openfiles

  call read_input_file

  call read_mesh

  call create_CSR_matrix

  call allocate_arrays

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

  time_loop: do itime=itimes,itimee ! time_loop 


    ! Update variables - shift in time: 
    call time_shift
   
    ! Set inlet boundary conditions
    if(itime.eq.itimes) call bcin

    ! Courant number report:
    call CourantNo

    ! 
    !===============================================
    ! SIMPLE iteration loop:
    !===============================================
    !
    iteration_loop: do iter=1,maxit 

      write(*,'(2x,a,i0)') 'Iter. ',iter

      call cpu_time(start)

      ! Calculate velocities.
      call calcuvw  

      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(SIMPLE)   call CALCP
      if(PISO)     call PISO_multiple_correction

      ! Turbulence
      if(lturb)    call correct_turbulence()

      !Scalars: Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calculate_temperature_field()

      ! if(lcal(ivart)) call calculate_temperature_variance_field()
      
      if(lcal(icon))  call calculate_concentration_field()


      call cpu_time(finish)
      write(timechar,'(f9.3)') finish-start
      write(6,'(3a)') '  ExecutionTime = ',trim(adjustl(timechar)),' s'
      write(6,*)

      !---------------------------------------------------------------
      !  Simulation management - residuals, loop control, output, etc.
      !---------------------------------------------------------------

      ! Check residual and stop program if residuals diverge
      source=max(resor(iu),resor(iv),resor(iw),resor(ip)) 
      if(source.gt.slarge) then
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
          if(source.lt.sormax.or.iter.ge.maxit) then 

            ! Correct driving force for a constant mass flow rate simulation:
            if(const_mflux) call constant_mass_flow_forcing

            ! Write field values after nzapis iterations or at the end of time-dependent simulation:
            if(mod(itime,nzapis).eq.0  .or. itime.eq.numstep) then
              call write_restart_files
              call writefiles
            endif

            ! Write values at monitoring points and recalculate time-average values for statistics:
            call writehistory
            call calc_statistics 

            cycle time_loop
         
          endif

      end if 

    end do iteration_loop

    ! Write field values after nzapis iterations or at the end of false-time-stepping simulation:
    if(.not.ltransient .and. ( mod(itime,nzapis).eq.0 .or. itime.eq.numstep ) ) then
        call write_restart_files
        call writefiles
    endif

    if(ltransient) call flush(6)

  end do time_loop

end program