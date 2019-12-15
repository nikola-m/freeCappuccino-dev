program steadyHeatEq
!***********************************************************************
!
! A solver for stationary heat equation.
! -div k grad (T) = Q
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use sparse_matrix, only: create_CSR_matrix,su
  use field_initialization, only: initialize_scalar_field
  use utils, only: show_logo, get_unit
  use output
  ! use LIS_linear_solver_library

  implicit none
!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: it
  integer :: output_file
  integer :: monitor_file
  real(dp), dimension(:), allocatable:: k 
  real(dp), dimension(:), allocatable:: T
  real(dp), dimension(:,:), allocatable:: dTdxi

!
!***********************************************************************
!

!
! > Open log file
! 
  monitor_file = 6
  open(unit=monitor_file,file='log')
  rewind monitor_file
! 
! > Print code logo and timestamp in log file
!
  call show_logo

!
! >  Open & Read mesh file, calculate mesh geometrical quantities, allocate arrays
!
  call read_mesh

  allocate(k(numCells))
  allocate(T(numTotal))
  allocate(dTdxi(3,numTotal))

!
! >  Index arrays of matrix elements stored in CSR format
!
  call create_CSR_matrix


! 
! > Temperature field initialization.
! 
  call initialize_scalar_field(T,dTdxi,'T')

  ! call get_unit( output_file )
  ! open(unit=output_file,file='temperature.vtu')

  ! call vtu_write_scalar_field ( output_file, 'Temp', T )

  ! close( output_file)

!
! > Discretisation and solution of the problem
!

  write(*,'(a)') ' '
  write(*,'(a)') '  Assembling linear system!'  

  ! Source term
  su(1:numCells) = 0.0_dp

    su(8232:numCells) = 1.0 !Vol(8232:numCells)
  
  ! Add volumetric source term
  !su(1:numCells) = 8*pi**2*sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells))*Vol(1:numCells)

  !  Coefficient array for Laplacian
  k = -1.0_dp       

  ! Laplacian operator and BCs         
  call laplacian(k,T) 

  it = 7
  sor(it) = 1e-16
  nsw(it) = 5000

  ! Solve system

  write(*,'(a)') ' '
  write(*,'(a)') '  Solving linear system!'  
  write(*,'(a)') ' '  

  ! 1)  Incomplete Cholesky Conjugate Gradient
  call iccg(T,it) 
  ! 2) Gauss-Seidel
  ! call gaussSeidel(T,it)
  ! 3) Iterative solver from LIS library
  ! call solve_csr(numCells,nnz,ioffset,ja,a,su,p) 


  !
  ! > Output
  !
  call get_unit( output_file )
  open(unit=output_file,file='temperature.vtu')

  call vtu_write_scalar_field ( output_file, 'Temp', T )

  close( output_file)

  write(*,'(a)') ' '
  write(*,'(a)') '  Output file written - simulation ended correctly.'  

end program