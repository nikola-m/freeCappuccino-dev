program poisson
!***********************************************************************
!
! Test discretisation of the Lapalcian operator on example problem:
! Laplacian(p) = -8 pi^2 sin(2 pi x) sin(2 pi y),
! with exact solution:
! p_exact = sin(2 pi x) sin(2 pi y)
! Test of meshes with different resolution and check the rate of grid 
! convergence.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use sparse_matrix, only: create_CSR_matrix,su
  use utils, only: show_logo
  use my_mpi_module
  use mpi

  implicit none


!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: ijp, ijn
  real(dp) :: lh
  real(dp), dimension(:), allocatable:: p

!
!***********************************************************************
!

! MPI start up

  call MPI_INIT(ierr)                   

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr  )  

  call MPI_COMM_RANK(MPI_COMM_WORLD, myid,  ierr  )  

  this = myid + 1


  ! > Create log file, print code logo and timestamp in log file
  if (myid .eq. 0) then

    open(unit=6,file='log_poisson')
    rewind 6

    call show_logo

  endif


!
! >  Open & Read mesh file, calculate mesh geometrical quantities, allocate arrays
!
  call read_mesh


!
! >  Index arrays of matrix elements stored in CSR format
!
  call create_CSR_matrix

!
! > Discretisation and solution of the problem
!

  allocate( p(numTotal) )

  ! Source term
  su(1:numCells) = 8*pi**2*sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells))*Vol(1:numCells)

  ! Initialize solution
  p = 0.0_dp     

  ! Laplacian operator and BCs         
  call laplacian(p) 

  sor(ip) = 1e-10
  nsw(ip) = 500

  ! Solve system
  ! 1)
  call dpcg(p,ip) 
  ! 2)
  ! call gaussSeidel(p,ip)
  ! 3)
  ! call solve_csr(numCells,nnz,ioffset,ja,a,su,p) 

  ! do i=1,numCells
  !   write(6,'(es11.4)') p(i)
  ! enddo 

  ! Cell size
  ijp = owner(1)
  ijn  = neighbour(1)
  lh = abs(xc(ijp)-xc(ijn))
  
  ! L_inf error norm: 
  if( myid .eq. 0 ) then
    write(*,'(a)') ' '
    write(*,'(2es11.4)') lh , maxval( abs( p(1:numCells)-sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells)) ) )
  endif

  ! MPI final call
  call MPI_Finalize(ierr)

end program