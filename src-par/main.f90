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
  use sparse_matrix
  use variables
  use title_mod
  use temperature
  use concentration
  use fieldManipulation
  use periodicity
  use utils
  use mpi

  implicit none

  integer :: iter
  integer :: itimes, itimee
  real(dp):: source
  real :: start, finish
!                                                                       
!******************************************************************************
!

!----------------------------------------------------------------------
! MPI start up

  call MPI_INIT(ierr)                   
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )  
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )  
  
  this = myid + 1

  if(nproc .eq. 1) then
    nproc = 0
    this = 0
  endif

!----------------------------------------------------------------------


!  Check command line arguments
  call get_command_argument(1,input_file)
  call get_command_argument(2,monitor_file)
  call get_command_argument(3,restart_file)

!-----------------------------------------------
!  Initialization, grid definition 
!-----------------------------------------------

  if (myid .eq. 0) then

    ! Simulation log file, monitor file for residuals
    open(unit=6,file=monitor_file)
    rewind 6
    
    ! Print cappuccino logo to log file.
    call show_logo

    write(*,'(a)') ' '
    write(*,'(a,i2)') '  Parallel run. Number of processes, np = ', nproc
    write(*,'(a)') ' '

  endif

  ! Read input file
  call read_input_file

  ! Read and process grid files
  call read_mesh

  ! Create sparse matrix data structure (CSR format)
  call create_CSR_matrix
 
  ! Allocate working arrays
  call allocate_arrays

  ! Initialisation of fields
  call init
 
  ! call add_random_noise_to_field(U,40)
  
  call face_mapping
!
!===============================================
!     T i m e   l o o p : 
!===============================================

  if (myid .eq. 0) then
    write(6,'(a)') ' '
    write(6,'(a)') '  Start iteration!'
    write(6,'(a)') ' '
  endif

  itimes = itime+1
  itimee = itime+numstep

  time_loop: do itime=itimes,itimee ! time_loop 

    ! Update variables - shift in time: 
    call time_shift

    ! Set inlet boundary conditions at every timestep
    if(itime.eq.itimes) call bcin

    ! Courant number report:
    call CourantNo

    ! if ( mod(itime,45).eq.0 ) call add_random_noise_to_field( U, 10 )

! 
!===============================================
!.....ITERATION loop
!===============================================
!
    iteration_loop: do iter=1,maxit 

      if (myid .eq. 0) then
        call cpu_time(start)
      endif

      ! Calculate velocities.
      call calcuvw 

      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(SIMPLE)   call calcp_simple
      if(PISO)     call calcp_piso

      ! Turbulence
      if(lturb)    call correct_turbulence()

      !Scalars: Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calculate_temperature_field

      ! if(lcal(ivart)) call calculate_temperature_variance_field()

      if(lcal(icon))  call calculate_concentration_field()


      if (myid .eq. 0) then
        call cpu_time(finish)
        write(timechar,'(f9.3)') finish-start
        write(6,'(3a)') '  ExecutionTime = ',trim(timechar),' s'
        write(6,*)
      endif

      ! Check residual
      source = max( resor(iu),resor(iv),resor(iw),resor(ip) ) 

      ! Now find global maximum
      call global_max( source )

      ! If iteration diverges
      if(source.gt.slarge) then
        if ( myid .eq. 0 ) write(6,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
        call abort_mission
      endif

      ! If residuals fall to level below tolerance level - simulation is finished.
      if( .not.ltransient  .and. source.lt.sormax ) then
        call write_restart_files
        call writefiles        
        exit time_loop
      end if

      if(ltransient) then 

        ! Has converged within timestep or has reached maximum no. of SIMPLE iterations per timetstep:
        if(source.lt.sormax.or.iter.ge.maxit) then 

          ! Correct driving force for a constant mass flow rate simulation:
          if(const_mflux) then
            call constant_mass_flow_forcing
            call recirculate_flow
          endif

          ! Write values at monitoring points and recalculate time-average values for statistics:
          call writehistory 
          call calc_statistics 

          ! Write field values after nzapis iterations or at the end of time-dependent simulation:
          if( mod(itime,nzapis).eq.0 .or. itime.eq.numstep ) then
            call write_restart_files
            call writefiles
          endif

          cycle time_loop
       
        endif

      end if 

    end do iteration_loop

    ! Write field values after nzapis iterations or at the end of false-time-stepping simulation:
    if(.not.ltransient .and. ( mod(itime,nzapis).eq.0 .or. itime.eq.numstep ) ) then
      call write_restart_files
      call writefiles
    endif


    if( ltransient .and. myid.eq.0 ) call flush(6)

  end do time_loop

  ! MPI final call
  call MPI_Finalize(ierr)
      
end program