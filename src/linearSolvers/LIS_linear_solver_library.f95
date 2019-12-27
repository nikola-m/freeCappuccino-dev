module LIS_linear_solver_library
!
! Description:
!   Interfaces LIS linear solver library.
!
! Usage:
!
! > 1. Initialization
!  call lis_initialize(ierr)
!
! > 2. Matrix creation
!
! LIS_INTEGER i,n
! LIS_MATRIX A
! n = XXXX
! call lis_matrix_create(0,A,ierr)
! call lis_matrix_set_size(A,0,n,ierr)
! do i=1,n
!   if( i>1 ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,1.0d0,A,ierr)
!   if( i<n ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,1.0d0,A,ierr)
!   call lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0d0,A,ierr)
! enddo
! call lis_matrix_set_type(A,LIS_MATRIX_CSR,ierr)
! call lis_matrix_assemble(A,ierr)
!
! Use following LIS integer TAGS:
!   LIS_INS_VALUE or LIS_INS_VALUE
! Destroy:
!   call lis_matrix_destroy(LIS_MATRIX A, LIS_INTEGER ierr)
! Already known format:
!   call lis_matrix_set_csr(LIS_INTEGER nnz, LIS_INTEGER ptr(),
!   LIS_INTEGER index(), LIS_SCALAR value(), LIS_MATRIX A, LIS_INTEGER ierr)
! eg.
!   n = 4; nnz = 10; k = 0;
!   call lis_matrix_set_csr(n,nnz,ptr,index,value,ierr)
!   call lis_matrix_create(0,A,ierr)
!   call lis_matrix_set_size(A,0,n,ierr)
!   ...populate ptr,index,value
!   call lis_matrix_set_csr(nnz,ptr,index,value,A, ierr)
!   call lis_matrix_assemble(A,ierr)
! where ptr <=> ioffset (samo sto pocinje od 0) array; index <=> ja array (samo sto pocinje od 0); value <=> amat (samo sto pocinje od 0)!
!
! > 3. Vector creation
!
! LIS_INTEGER i,n
! LIS_VECTOR v
! n = 4
! call lis_vector_create(0,v,ierr)
! call lis_vector_set_size(v,0,n,ierr)
! do i=1,n
!   call lis_vector_set_value(LIS_INS_VALUE,i,DBLE(i),v,ierr)
! enddo
!
! Use following LIS integer TAGS:
!   LIS_INS_VALUE or LIS_INS_VALUE
! Also: to create vectro with same information (as opposed to lis_vector_copy)
!   call lis_vector_duplicate(LIS_VECTOR vin, LIS_VECTOR vout, LIS_INTEGER ierr)
! Destroy:
!   call lis_vector_destroy(LIS_VECTOR v,LIS_INTEGER ierr)
!
! > 4. Solver creation
! LIS_MATRIX A
! LIS_VECTOR b,x
! LIS_SOLVER solver
! /* Create matrix and vector */
! call lis_solver_create(solver,ierr)
!
! > 5. Value assignment for matrices and vectors
! ...
!
! > 6. Solver assignment
! call lis_solver_set_option(’-i bicg -p none’,solver,ierr)
! call lis_solver_set_option(’-tol 1.0e-12’,solver,ierr)
!      call lis_solver_set_option("-print mem",solver,ierr);
!      call lis_solver_set_optionC(solver,ierr)
!
! > 7. Solver execution
! call lis_solve(A,b,x,solver,ierr)
!
! > 8. Finalization
!
!  Author:
!    Nikola Mirkov, nmirkov@vinca.rs
!
!  Date:
!    25/11/2015
!


  use types
  use geometry, only: numCells,numTotal
  use sparse_matrix, only: nnz, ioffset, ja, a, diag

  #include "lisf.h"

  implicit none

  character(len=70) :: lis_solver_options

  public


contains


! subroutine solve_csr(numCells,nnz,ioffset,ja,aval,phi,su)
subroutine lis_spsolve( numCells, nnz, ioffset, ja, a, fi, rhs, itr_max, tol_abs, tol_rel, chvar)
!
! A subroutine which handles linear system solution using LIS library.
!
  use types
  implicit none

  ! integer, intent(in) :: numCells
  ! integer, intent(in) :: nnz
  ! integer, dimension(numCells+1), intent(in) :: ioffset
  ! integer, dimension(nnz), intent(in) :: ja
  ! real(dp), dimension(nnz), intent(in) :: aval
  ! real(dp), dimension(numCells), intent(in) :: su
  ! real(dp), dimension(numCells), intent(inout) :: phi
  integer, intent(in) :: n                              ! Number of unknowns, length of a solution vector
  integer, intent(in) :: nnz                            ! Number of non-zero elements in sparse matrix
  integer, dimension(n+1), intent(in) :: ioffset        ! The offsets of each row in coef array
  integer, dimension(nnz), intent(in) :: ja             ! Columns array
  real(dp), dimension(nnz), intent(in) :: a             ! Coefficient array
  integer, dimension(n), intent(in) :: diag             ! Position of diagonal elements in coeff. array
  real(dp), dimension(n), intent(inout) :: fi           ! On input-current field; on output-the solution vector
  real(dp), dimension(n), intent(in) :: rhs             ! The right hand side of the linear system
  integer, intent(in) :: itr_max                        ! Maximum number of iterations
  real(dp), intent(in) :: tol_abs                       ! Absolute tolerance level for residual
  real(dp), intent(in) :: tol_rel                       ! Relative tolerance level for residual
  character( len=* ), intent(in) :: chvar               ! Character string containing name of the solved field, printed on stdout


  LIS_MATRIX :: A
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
  character(len=10) :: solver_type
  character(len=10) :: preconditioner
  real(dp), parameter :: zero = 0.0d0
  integer, parameter :: numthrd = 2 ! hard coded number of threads

 

! > 1. Initialization

  call lis_initialize(ierr)

  call omp_set_num_threads(numthrd)

  ! Sizes
  nz_num = nnz
  cell_num = n

! > 2. Matrix creation

 call lis_matrix_create(0,A,ierr)
 call lis_matrix_set_size(A,0,cell_num,ierr)

 row = 0
 icell = 1
 do i=1,nz_num
   if( i.eq.ioffset(icell) ) then
   row  = row+1 
   icell = icell + 1
   endif
   call lis_matrix_set_value(LIS_INS_VALUE,row,ja(i),a(i),A,ierr)
 enddo
 call lis_matrix_set_type(A,LIS_MATRIX_CSR,ierr)
 call lis_matrix_assemble(A,ierr)

!  Matrix in CSR format:
!  allocate(ptr(cell_num+1))
!  allocate(index(nz_num))
!  allocate(value(nz_num))
!...or new version of Fortran will allocate uppon the assignement
!  ptr = ioffset
!  index = ja
!  value = a
!
! call lis_matrix_set_csr(nz_num,ptr,index,value,A, ierr)
! call lis_matrix_set_type(A,LIS_MATRIX_CSR,ierr)
! call lis_matrix_assemble(A,ierr)

! > 3. Vector creation

! rhs vector
 call lis_vector_create(0,b,ierr)
 call lis_vector_set_size( b,0,cell_num,ierr )
 do i=1,n
   call lis_vector_set_value( LIS_INS_VALUE,i,rhs(i),b,ierr )
 enddo

! solution vector
 call lis_vector_create( 0,x,ierr )
 call lis_vector_set_size( x,0,cell_num,ierr )
 ! Initialize solutio vector
 do i=1,n
   call lis_vector_set_value( LIS_INS_VALUE,i,fi(i),x,ierr )
 enddo



! > 6. Solver creation and setting solver options
 call lis_solver_create( solver,ierr )

 solver_type = 'cg' ! 'gs',...
 preconditioner = 'ssor'

 ! call lis_solver_set_option("-print mem",solver,ierr)
 call lis_solver_set_option( trim( lis_solver_options ),solver,ierr )
 call lis_solver_set_option( "-i cg -p ssor",solver,ierr)
 ! call lis_solver_set_option("-i gs",solver,ierr)
 ! call lis_solver_set_option("-maxiter 1000",solver,ierr)

 call lis_solver_set_optionC(solver,ierr)

! > 7. Solver execution
  ! write(6,*) ' '

  call lis_solver_get_residualnorm(solver,res0,ierr)

  res_init = res0

  call lis_solve(A,b,x,solver,ierr)

  ! write solution 
  ! solname = "rjesenje.txt"
  ! call lis_output_vector(x,LIS_FMT_MM,solname,ierr);

  ! write residual history
  ! resname = "residual.txt"
  ! call lis_solver_output_rhistory(solver, resname, ierr)

  ! > Write it in our solution vector
  do i=1,n
   call lis_vector_get_value(x, i, xval, ierr)
   fi(i) = xval
  enddo


! > Output iterative solution summary
      call lis_solver_get_iterex(solver,iter,iter_double,iter_quad,ierr)
      ! call lis_solver_get_timeex(solver,time,itime,ptime,p_c_time,p_i_time,ierr)
      call lis_solver_get_residualnorm(solver,resid,ierr)
      ! call lis_solver_get_solver(solver,nsol,ierr)
      ! call lis_solver_get_solvername(nsol,solvername,ierr)

      ! if( debug ) then
        ! write(6,'(2a,i4)')     trim( solvername ),': number of iterations = ',iter
        ! write(6,'(2a,es11.4)') trim( solvername ),': elapsed time         = ',time
        ! write(6,'(2a,es11.4)') trim( solvername ),':   preconditioner     = ',ptime
        ! write(6,'(2a,es11.4)') trim( solvername ),':     matrix creation  = ',p_c_time
        ! write(6,'(2a,es11.4)') trim( solvername ),':   linear solver      = ',itime
        ! write(6,'(2a,es11.4)') trim( solvername ),': residual             = ',resid
      ! endif

  ! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  LIS['//trim(solver_type)//'('//trim(preconditioner)//')]:  Solving for ', &
  trim( chvar ),', Initial residual = ',res0,', Final residual = ',resid,', No Iterations ',iter

! > 8. Finalization

 call lis_solver_destroy(solver,ierr)
 call lis_matrix_unset(A,ierr)
 call lis_matrix_destroy(A,ierr)
 call lis_vector_destroy(x,ierr)
 call lis_vector_destroy(b,ierr)

! deallocate(ptr)
! deallocate(index)
! deallocate(value)

 call lis_finalize(ierr)

end subroutine lis_spsolve

end module LIS_linear_solver_library
