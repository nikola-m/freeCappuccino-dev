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


! subroutine lis_spsolve( ioffset, ja, a, fi, rhs, itr_max, tol_abs, tol_rel, chvar)
subroutine lis_spsolve( lis_solver_options, fi, rhs, chvar)
!
! A subroutine which handles linear system solution using LIS library.
!
  use types
  use geometry, only: numCells
  use sparse_matrix, only: nnz, ioffset, ja, a

  implicit none

  ! integer, dimension(numCells+1), intent(in) :: ioffset ! The offsets of each row in coef array
  ! integer, dimension(nnz), intent(in) :: ja             ! Columns array
  ! real(dp), dimension(nnz), intent(in) :: a             ! Coefficient array
  ! integer, dimension(n), intent(in) :: diag             ! Position of diagonal elements in coeff. array
  character( len=* ), intent(in) :: lis_solver_options
  real(dp), dimension(numCells), intent(inout) :: fi           ! On input-current field; on output-the solution vector
  real(dp), dimension(numCells), intent(in) :: rhs             ! The right hand side of the linear system
  ! integer, intent(in) :: itr_max                        ! Maximum number of iterations
  ! real(dp), intent(in) :: tol_abs                       ! Absolute tolerance level for residual
  ! real(dp), intent(in) :: tol_rel                       ! Relative tolerance level for residual
  character( len=1), intent(in) :: chvar               ! Character string containing name of the solved field, printed on stdout


  LIS_MATRIX :: Alis
  LIS_VECTOR :: b,x
  LIS_SOLVER :: solver
  LIS_INTEGER :: ierr
  LIS_INTEGER :: i,nz_num,cell_num,row,icell
  LIS_SCALAR :: xval
  LIS_INTEGER :: iter,iter_double,iter_quad
  LIS_REAL :: resid,res0
  ! LIS_INTEGER, allocatable :: ptr(:),index(:)
  ! LIS_SCALAR, dimension(:), allocatable :: value
  ! LIS_INTEGER :: nsol
  ! real(dp) :: time,itime,ptime,p_c_time,p_i_time
  ! character(len=256) :: resname!,solname ! filenames for solution vector and residual vector
  ! character(len=10) :: solver_type
  ! character(len=10) :: preconditioner
  real(dp), parameter :: zero = 0.0d0
  integer, parameter :: numthrd = 1 ! hard coded number of threads

 
  ! Sizes
  nz_num = nnz
  cell_num = numCells


! > 1. Initialization

  call lis_initialize(ierr)

  ! call omp_set_num_threads(numthrd)


! > 2. Matrix creation

 call lis_matrix_create(0,Alis,ierr)
 call lis_matrix_set_size(Alis,0,cell_num,ierr)

 row = 0
 icell = 1
 do i=1,nz_num
   if( i == ioffset(icell) ) then
   row  = row+1 
   icell = icell + 1
   endif
   call lis_matrix_set_value(LIS_INS_VALUE,row,ja(i),a(i),Alis,ierr)
 enddo
 call lis_matrix_set_type(Alis,LIS_MATRIX_CSR,ierr)
 call lis_matrix_assemble(Alis,ierr)



! > 3. Vector creation

! rhs vector
 call lis_vector_create(0,b,ierr)
 call lis_vector_set_size( b,0,cell_num,ierr )
 do i=1,cell_num
   call lis_vector_set_value( LIS_INS_VALUE,i,rhs(i),b,ierr )
 enddo

! solution vector
 call lis_vector_create( 0,x,ierr )
 call lis_vector_set_size( x,0,cell_num,ierr )
 do i=1,cell_num
   call lis_vector_set_value( LIS_INS_VALUE,i,fi(i),x,ierr )
 enddo



! > 4. Solver creation and setting solver options
 call lis_solver_create( solver,ierr )

 call lis_solver_set_option( trim( lis_solver_options ),solver,ierr )

 call lis_solver_set_optionC(solver,ierr)


! > 5. Solver execution

  call lis_solver_get_residualnorm(solver,res0,ierr)

  call lis_solve(Alis,b,x,solver,ierr)

  ! write solution 
  ! solname = "rjesenje.txt"
  ! call lis_output_vector(x,LIS_FMT_MM,solname,ierr);

  ! write residual history
  ! resname = "residual.txt"
  ! call lis_solver_output_rhistory(solver, resname, ierr)

  ! > Write it in our solution vector
  do i=1,cell_num
   call lis_vector_get_value(x, i, xval, ierr)
   fi(i) = xval
  enddo


! > Output iterative solution summary
      call lis_solver_get_iterex(solver,iter,iter_double,iter_quad,ierr)
      call lis_solver_get_residualnorm(solver,resid,ierr)
      ! call lis_solver_get_timeex(solver,time,itime,ptime,p_c_time,p_i_time,ierr)
      ! call lis_solver_get_solver(solver,nsol,ierr)
      ! call lis_solver_get_solvername(nsol,solvername,ierr)

        ! write(6,'(a,i0)') '  LIS['//trim(lis_solver_options)//']: number of iterations = ',iter
        ! write(6,'(a,es11.4)') '  LIS['//trim(lis_solver_options)//']:   linear solver      = ',itime
        ! write(6,'(a,es11.4)') '  LIS['//trim(lis_solver_options)//']: residual             = ',resid


  ! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  LIS['//trim(lis_solver_options)//']:  Solving for ', &
  trim( chvar ),', Initial residual = ',res0,', Final residual = ',resid,', No Iterations ',iter

! > 6. Finalization

 call lis_solver_destroy(solver,ierr)
 call lis_matrix_unset(Alis,ierr)
 call lis_matrix_destroy(Alis,ierr)
 call lis_vector_destroy(x,ierr)
 call lis_vector_destroy(b,ierr)

 call lis_finalize(ierr)

end subroutine lis_spsolve

end module LIS_linear_solvers
