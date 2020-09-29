module my_mpi_module
!
! Module that contains data and functions to perform inter-process exchange
! of data in buffer cells, between processes, in parallel execution of the program.
! 
  use types

  implicit none

    ! MPI related 
    integer :: lenbuf ! Buffer size, total no of faces that divide this and other domains
    integer :: numConnections ! broj konektovanih domena na ovaj trenutni
    integer, dimension(:), allocatable :: neighbProcNo       ! size[1,num_connections]; Value: self descriptive.
    integer, dimension(:), allocatable :: neighbProcOffset   ! size[1,num_connections+1]; Value: Where, in arrays of size lenbuf, do the faces of certain conection start.
    integer, dimension(:), allocatable :: bufind             ! size[1,len_buffer]; Represents indexes of cell at process boudary.
    real(dp), dimension(:), allocatable :: buffer            ! size[1,len_buffer]; Stores buffer values of type real.

  public

contains

!
!***********************************************************************
!
  subroutine global_isum(i) 
!
!***********************************************************************
!
!   Estimates global sum among all processors of an integer variable.                    
!
!***********************************************************************
!
  use types
  use mpi

  implicit none

  integer :: i

  integer :: isum
  integer :: ierr

  call mpi_allreduce      &               
   (i,                    & ! send buffer
    isum,                 & ! recv buffer 
    1,                    & ! length     
    mpi_integer,          & ! datatype  
    mpi_sum,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  i = isum

  end subroutine global_isum

!
!***********************************************************************
!
  subroutine global_sum(phi) 
!
!***********************************************************************
!
!   Estimates global sum among all processors.
!   Used e.g. in dot product. Every process creates it own sum
!   end stores it in phi. Than it all gets gathered and summed and sent
!   back to each process.                       
!
!***********************************************************************
!
  use types
  use mpi
  
  implicit none

  real(dp), intent(inout) :: phi

  integer :: ierr
  real(dp) :: phisum

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  CALL MPI_ALLREDUCE      &               
   (phi,                  & ! send buffer
    phisum,               & ! recv buffer 
    1,                    & ! length     
    MPI_DOUBLE_PRECISION, & ! datatype  mpi_double_precision
    MPI_SUM,              & ! operation 
    MPI_COMM_WORLD,       & ! communicator            
    ierr) 

  phi = phisum

  end subroutine global_sum

!
!***********************************************************************
!
  subroutine global_max(phi) 
!
!***********************************************************************
!
!   Estimates global maximum value of phi among all processes.
!   Every process creates it own maximum end stores it phi.
!   Than it all gets gathered and summed and sent back to each process.                       
!
!***********************************************************************
!
  use types
  use mpi
  
  implicit none


  real(dp) :: phi

  real(dp) :: phimax
  integer :: ierr

  call mpi_allreduce      &               
   (phi,                  & ! send buffer
    phimax,               & ! recv buffer 
    1,                    & ! length     
    mpi_double_precision, & ! datatype  
    mpi_max,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  phi = phimax

  end subroutine global_max

!
!***********************************************************************
!
  subroutine global_min(phi) 
!
!***********************************************************************
!
!   Estimates global minimum value of phi among all processes.
!   Every process creates it own minimum end stores it phi.
!   Than it all gets gathered and summed and sent back to each process.                       
!
!***********************************************************************
!
  use types
  use mpi 

  implicit none

  real(dp) :: phi

  real(dp) :: phimin
  integer :: ierr

  call mpi_allreduce      &               
   (phi,                  & ! send buffer
    phimin,               & ! recv buffer 
    1,                    & ! length     
    mpi_double_precision, & ! datatype  
    mpi_min,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  phi = phimin

  end subroutine global_min

!
!***********************************************************************
!
  subroutine synchronize 
!
!***********************************************************************
!
!  Calls mpi_barrier and syncronizes all processes at an end of a task.                      
!
!***********************************************************************
!
  use mpi

  implicit none

  integer :: ierr

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  end subroutine

!
!***********************************************************************
!
  subroutine abort_mission
!
!***********************************************************************
!
!   Finalize parallel communication and  program execution.                    
!
!***********************************************************************
!
  use parameters
  use mpi
  
  implicit none

  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_finalize(ierr)
  stop

  end subroutine


end module