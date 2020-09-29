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