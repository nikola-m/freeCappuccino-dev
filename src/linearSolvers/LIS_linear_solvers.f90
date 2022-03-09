module LIS_linear_solvers
!
! Description:
!   Interfaces LIS linear solver library.
!
!  Author:
!    Nikola Mirkov, nmirkov@vinca.rs
!
!  Date:
!    25/11/2015
!

#include "lisf.h"

  implicit none


  public


contains



subroutine lis_spsolve( lis_solver_options, fi, rhs, chvar )
!
! A subroutine which handles linear system solution using LIS library.
!
  use types
  use parameters, only: zero
  use geometry, only: numCells
  use sparse_matrix, only: nnz, ioffset, ja, a

  implicit none

  character( len=* ), intent(in) :: lis_solver_options
  real(dp), dimension(numCells), intent(inout) :: fi      ! On input-current field; on output-the solution vector
  real(dp), dimension(numCells), intent(in) :: rhs        ! The right hand side of the linear system
  character( len=1), intent(in) :: chvar                  ! Character string containing name of the solved field, printed on stdout

!
! > Local variables
!

  integer, parameter :: numthrd = 1 ! hard coded number of threads

  LIS_MATRIX :: Alis
  LIS_VECTOR :: b,x
  LIS_SOLVER :: solver
  LIS_INTEGER :: ierr
  LIS_INTEGER :: i,k,row
  LIS_SCALAR :: xval
  LIS_INTEGER :: iter,iter_double,iter_quad
  LIS_REAL :: resid,res0


! > 1. Initialization

  call lis_initialize(ierr)

  call omp_set_num_threads(numthrd)


! > 2. Matrix creation

  call lis_matrix_create(0,Alis,ierr)
  call lis_matrix_set_size(Alis,0,numCells,ierr)

  row = 0
  k = 1
  do i=1,nnz
    if( i == ioffset(k) ) then
      row  = row+1 
      k = k + 1
    endif
    call lis_matrix_set_value(LIS_INS_VALUE,row,ja(i),a(i),Alis,ierr)
  enddo
  call lis_matrix_set_type(Alis,LIS_MATRIX_CSR,ierr)
  call lis_matrix_assemble(Alis,ierr)



! > 3. Vector creation

! RHS vector
  call lis_vector_create(0,b,ierr)
  call lis_vector_set_size( b,0,numCells,ierr )
  do i=1,numCells
   call lis_vector_set_value( LIS_INS_VALUE,i,rhs(i),b,ierr )
  enddo

! Solution vector - initial value
  call lis_vector_create( 0,x,ierr )
  call lis_vector_set_size( x,0,numCells,ierr )
  do i=1,numCells
   call lis_vector_set_value( LIS_INS_VALUE,i,fi(i),x,ierr )
  enddo

! Calculate initial residual norm
  res0 = 0.0_dp
  do i=1,numCells
    resid = rhs(i) 
    do k = ioffset(i), ioffset(i+1)-1
      resid = resid -  a(k) * fi( ja(k) ) 
    enddo
    res0 = res0 + abs(resid) 
  enddo


! > 4. Solver creation and setting solver options
  call lis_solver_create( solver,ierr )

  call lis_solver_set_option( trim( lis_solver_options ), solver, ierr )

  call lis_solver_set_optionC(solver,ierr)


! > 5. Solver execution

  call lis_solve(Alis,b,x,solver,ierr)

  ! > Write it in our solution vector
  do i=1,numCells
   call lis_vector_get_value(x, i, xval, ierr)
   fi(i) = xval
  enddo


! > Output iterative solution summary
  call lis_solver_get_iterex(solver,iter,iter_double,iter_quad,ierr)
  call lis_solver_get_residualnorm(solver,resid,ierr)


  ! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  LIS['//trim(lis_solver_options)//']:  Solving for ', &
  trim( chvar ),', Initial residual = ',res0,', Final residual = ',resid,', No Iterations ',iter

! > 6. Finalization

 call lis_solver_destroy(solver,ierr)
 call lis_matrix_destroy(Alis,ierr)
 call lis_vector_destroy(x,ierr)
 call lis_vector_destroy(b,ierr)

 call lis_finalize(ierr)

end subroutine


end module