program test_sparse
!
! Purpose:
!   Test iterative solvers from linear_solvers module.
! Description:
!   The program computes the solution to the system of linear
!   equations with a square matrix A and
!   right-hand side B, where A is the coefficient matrix.
!
! Test 1 - Real nonsymmetric 5x5 system :
!
!   6.80  -6.05  -0.45   8.32  -9.67
!  -2.11  -3.30   2.58   2.71  -5.14
!   5.66   5.36  -2.70   4.35  -7.26
!   5.97  -4.44   0.27  -7.17   6.08
!   8.23   1.08   9.04   2.14  -6.87
!
! and B is the right-hand side:
!
!   4.02
!   6.19 
!  -8.22
!  -7.57
!  -3.03
!
! Solution
!  -0.80
!  -0.70
!   0.59
!   1.32 
!   0.57
!
!
! Test 2 - symmetric positive definite 5x5 system!
!
!    3.14   0.17  -0.90   1.65  -0.72
!    0.17   0.79   0.83  -0.65   0.28
!   -0.90   0.83   4.53  -3.70   1.60
!    1.65  -0.65  -3.70   5.32  -1.37
!   -0.72   0.28   1.60  -1.37   1.98
!
!  and B is the right-hand side:
!
!   -7.29
!    9.25
!    5.99
!   -1.94
!   -8.30 
!
! Solution
!  -6.02
!  15.62
!   3.02
!   3.25
!  -8.78
!
  use types, only: dp
  implicit none

  integer :: i,ioffset(6),diag(5)
  integer :: ja(25)
  integer :: itr_max
  real(dp) :: a(25), b(5), x(5), xa(5)
  real(dp) :: res0, tol_abs, tol_rel

  itr_max = 5
  tol_rel = 1e-7
  tol_abs = 1e-10

!
! Test 1 - Real nonsymmetric 5x5 system
!
  write(*,'(a)') ' '
  write(*,'(a)') ' Test 1 - Real nonsymmetric 5x5 system.'
  write(*,'(a)') ' '

  a = (/ &
  6.80, -6.05,  -0.45,   8.32,  -9.67,  &
  -2.11,  -3.30,   2.58,   2.71, -5.14, &
  5.66,  5.36,  -2.70,   4.35,  -7.26,  &
  5.97,  -4.44,   0.27,  -7.17,   6.08, &
  8.23,   1.08,   9.04,   2.14,  -6.87 /)

  b = (/ 4.02, 6.19,-8.22,-7.57,-3.03 /)

  ja = (/ &
      1,2,3,4,5,&
      1,2,3,4,5,&
      1,2,3,4,5,&
      1,2,3,4,5,&
      1,2,3,4,5 &
      /)

  ioffset = (/ 1, 6, 11, 16, 21, 26 /)

  diag = (/ 1, 7, 13, 19, 25 /)

  xa = (/ -0.80, -0.70, 0.59, 1.32, 0.57 /)

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The ILU(0) preconditioned BiCGStab solver:'

  call spsolve('bicgstab_ilu', x, b, res0, itr_max, tol_abs, tol_rel, 'x')

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The ILU(0) preconditioned Restarted GMRES solver:'

  call spsolve('pmgmres_ilu', x, b, res0, itr_max, tol_abs, tol_rel, 'x')

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

!
! Test 2 - symmetric positive definite 5x5 system
!
  write(*,'(a)') ' '
  write(*,'(a)') ' Test 2 - symmetric positive definite 5x5 system.'
  write(*,'(a)') ' '

  a = (/ &
    3.14,   0.17,  -0.90,   1.65,  -0.72, &
    0.17,   0.79,   0.83,  -0.65,   0.28, &
   -0.90,   0.83,   4.53,  -3.70,   1.60, &
    1.65,  -0.65,  -3.70,   5.32,  -1.37, &
   -0.72,   0.28,   1.60,  -1.37,   1.98 /)

  b = (/ -7.29, 9.25, 5.99, -1.94, -8.30 /)

  xa = (/  -6.02,15.62,3.02,3.25,-8.78 /)

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The Incomplete Cholesky preconditioned CG solver:'

  call spsolve('iccg', x, b, res0, itr_max, tol_abs, tol_rel, 'x')

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The Diagonally preconditioned CG solver:'

  call spsolve('dpcg', x, b, res0, itr_max, tol_abs, tol_rel, 'x')

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The Gauss-Seidel solver:'

  call spsolve('gauss-seidel', x, b, res0, 100, tol_abs, tol_rel, 'x')

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

end program

