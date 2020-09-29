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