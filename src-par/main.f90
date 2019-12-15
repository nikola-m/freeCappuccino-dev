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
  use mpi

  implicit none

  integer :: iter, i, ijp, ijn, inp, ib, iface
  integer :: narg
  integer :: itimes, itimee
  real(dp):: magUbarStar, rUAw, gragPplus, flowDirection
  real(dp):: source
  real(dp):: suma,dt
  real :: start, finish
!                                                                       
!******************************************************************************
!

!----------------------------------------------------------------------
! MPI start up

  call MPI_INIT(ierr)                   

  call MPI_COMM_SIZE(MPI_COMM_WORLD, & 
                               nproc,&  
                               ierr  )  

  call MPI_COMM_RANK(MPI_COMM_WORLD, & 
                               myid, &  
                               ierr  )  
  
  this = myid + 1

  if(nproc .eq. 1) then
    nproc = 0
    this = 0
  endif

  ! write(*,'(2(a,i2))') ' np = ', nproc, ' myid = ', myid
!----------------------------------------------------------------------


!  Check command line arguments
  narg=command_argument_count()
  if (narg==0.or.narg<4) write(*,'(a)') 'Usage: '&
  &'./cappuccino <input_file> <monitor_file> <restart_file> <out_folder_path>'
  call get_command_argument(1,input_file)
  call get_command_argument(2,monitor_file)
  call get_command_argument(3,restart_file)
  call get_command_argument(4,out_folder_path)

!-----------------------------------------------
!  Initialization, grid definition 
!-----------------------------------------------

  if (myid .eq. 0) then

    ! Open files
    call openfiles

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
    include 'CourantNo.h'

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
  !  ! MPI final call
  ! call MPI_Finalize(ierr)
  ! stop  
      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(SIMPLE)   call CALCP
      if(PISO)     call PISO_multiple_correction

      ! Turbulence
      if(lturb)    call correct_turbulence()

      !Scalars: Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calculate_temperature_field()

      ! if(lcal(ivart)) call calculate_temperature_variance_field()

      if(lcal(icon))  call calculate_concentration_field()


      if (myid .eq. 0) then
        call cpu_time(finish)
        write(timechar,'(f9.3)') finish-start
        write(6,'(3a)') '  ExecutionTime = ',trim(adjustl(timechar)),' s'
        write(6,*)
      endif


      ! ! Residual normalization, convergence check  
      ! do i=1,nphi
      !   resor(i)=resor(i)*rnor(i)
      ! end do


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
      if(.not.ltransient  .and. source.lt.sormax ) then
          call write_restart_files
          call writefiles        
          exit time_loop
      end if

      if(ltransient) then 

          ! Has converged within timestep or has reached maximum no. of SIMPLE iterations per timetstep:
          if(source.lt.sormax.or.iter.ge.maxit) then 

            if(const_mflux) then
              ! Correct driving force for a constant mass flow rate.
              include 'constant_mass_flow_forcing.f90'
            endif

            ! Write field values after nzapis iterations or at the end of time-dependent simulation:
            if( mod(itime,nzapis).eq.0 .and. itime.ne.numstep ) then
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


    if( ltransient .and. myid.eq.0 ) call flush(6)

  end do time_loop

  ! MPI final call
  call MPI_Finalize(ierr)
      
end program