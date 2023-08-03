module sparse_matrix
!
! Module defining CSR (Compressed Sparse Row) sparse matrix data type.
!
  use types
  use geometry, only: numCells,numInnerFaces,owner,neighbour,numPeriodic,numBoundaries,bctype,nfaces,startFace,startFaceTwin
  use utils, only: csr_to_k, find_index_position, find_main_diag_element_positions, i4vec2_sort_a, i4vec_print, i4vec_print2

  ! Matrix in sparse format CSR(ia,ja,a) and COO(iacoo,ja,a)
  integer :: nnz                                ! No. of nonzeros in sparse system matrix.
  integer, dimension(:), allocatable :: ia      ! Positions where is the start of the next row in matrix 'a'; size[1:ncell+1] CSR format.
  integer, dimension(:), allocatable :: iacoo   ! Rows - coordinate format [1:nnz].
  integer, dimension(:), allocatable :: ja      ! Columns [1:nnz].
  integer, dimension(:), allocatable :: diag    ! Position of diagonal elements [1:ncell].
  real(dp), dimension(:), allocatable :: a      ! Coefficient matrix [1:nnz].

  integer, dimension(:), allocatable :: icell_jcell_csr_index !(i,j) matrix element transfered to a position in an array of length (1:nnz)
  integer, dimension(:), allocatable :: jcell_icell_csr_index !(j,i) matrix element transfered to a position in an array of length (1:nnz)

  ! Coefficients resulting form fvm discretization:
  ! real(dp), dimension(:), allocatable :: res           ! Residual vector for linear solvers
  real(dp), dimension(:), allocatable :: spu,spv,sp    ! Source terms for the left hand side
  real(dp), dimension(:), allocatable :: su, sv, sw    ! Source terms for the left hand side of equation
  real(dp), dimension(:), allocatable :: apu, apv, apw ! Reciprocal values of diagonal coefficients

  real(dp), dimension(:), allocatable :: h,rU,rV,rW ! Arrays used in piso algorithm.
  
  !
  ! The CSR matrix derived type
  !
  type csrMatrix
    integer, dimension(:), allocatable :: ia
    integer, dimension(:), allocatable :: ja
    integer, dimension(:), allocatable :: diag   
    real(dp), dimension(:), allocatable :: a
  end type


public

contains

!
! > Create new CSR matrix object for given nnz, numCells, ia and ja data arrays.
!
function new_csrMatrix( ) result(csr)
  implicit none
  integer :: i
  type(csrMatrix) :: csr
  

  ! Option 1.

  allocate(csr%ia ( numCells+1 ))
  allocate(csr%ja ( nnz ))
  allocate(csr%diag ( numCells ))
  allocate(csr%a ( nnz ))

  do i=1,numCells+1
    csr % ia(i) = ia(i)
  enddo

  do i=1,numCells
    csr % diag(i) = diag(i)
  enddo
  
  do i=1,nnz
    csr % ja(i) = ja(i)
  enddo

  !
  ! Option 2. Instead of all this above - allocate on asignement
  !
  ! csr % ia = ia
  ! csr % diag = diag
  ! csr % ja = ja
  ! csr % a  = a
  

end function

!
! > Create new CSR matrix object for given mesh data.
!

subroutine create_CSR_matrix
!
! Define sparsity pattern according to given mesh connectivity data,
! and allocate arays representing system matrix in CSR format.
!

  implicit none
  
  integer :: i
  integer :: icell,ijp,ijn
  integer :: istart, iend
  integer :: k,ib,iPer


  ! Number of non-zero elements in sparse matrix: nnz
  ! We brake it down like this:
  ! If we have cell P and cell N as its neighbour, we will have:
  ! In Pth row, matrix elements: PP, PN 
  ! In Nth row, matrix elements: NP, NN
  ! So finding one inner face (or periodic boundary face, since we treat them as inner faces)
  ! we will have two off-diagonal elements, therefore 2*numInnerFaces (2*numPeriodic).
  ! That was about off diagonal elements, we will also need elements on the main diagonal
  ! of the matrix, and we have numCells of them, since there is numCells linear equations
  ! in a system on numCells unknowns. Nikola
  nnz = 2*(numInnerFaces+numPeriodic) + numCells

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of nonzero coefficients in sparse matrix, nnz = ', nnz

  allocate ( iacoo(nnz) ) 
  allocate ( ja(nnz) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Creating CSR matrix sparsity pattern based on mesh data.'

!
!  > Populate sparsity arrays for CSR format
!

  ! This will be postions of diagonal elements, but later we'll sort ia,ja.
  do icell = 1,numCells
    iacoo(icell) = icell
       ja(icell) = icell
  enddo

  do i = 1,numInnerFaces
    iacoo(numCells+i) = owner(i)
       ja(numCells+i) = neighbour(i) 
  enddo

  do i = 1,numInnerFaces
    iacoo(numCells+numInnerFaces+i) = neighbour(i) 
       ja(numCells+numInnerFaces+i) = owner(i)
  enddo

  ! Matrix elements resulting from connecting cells of the periodic boundaries
  
  iPer = 0 ! counts every instance of finding periodic boundary
  k = numCells+2*numInnerFaces ! start counting from this number, see why few lines above, wehere ia,ja are populated.

  do ib=1,numBoundaries

    if (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1 ! found one pair of periodic boundaries, increase counter by one

      ! Loop trough faces
      do i=1,nfaces(ib)

        ijp = owner( startFace(ib) + i )
        ijn = owner( startFaceTwin(iPer) + i )

        k = k + 1 ! next
        iacoo(k) = ijp
           ja(k) = ijn

        k = k + 1
        iacoo(k) = ijn
           ja(k) = ijp


      enddo

    endif

  enddo
!
!  > Lexically sort the ia, ja values.
!
  call i4vec2_sort_a ( nnz, iacoo, ja )

  ! call i4vec_print2 ( 10, iacoo, ja, '  First 10 lines of Sorted IA and JA arrays:' )


!
! > Find positions of diagonal elements in COO matrix format
!
 allocate ( diag(numCells) )

 call find_main_diag_element_positions ( iacoo, ja, nnz, diag, numCells )

 ! call i4vec_print ( 10, diag, '  First 20 lines of Diagonal adjacency vector:' )


!
! > Find positions of row starting matrix elements in COO format and put these positions in 'ia' array.
!
 allocate ( ia(numCells+1) )

 istart = 1
 iend = 1
 do icell = 1,numCells
   call find_index_position( icell, istart, iend, iacoo, nnz, ia(icell) )
   istart = ia(icell)+1
   iend = nnz
 enddo
 ia(numCells+1) = nnz+1 ! poslednji element

 ! call i4vec_print ( 10, ia, '  First 10 lines of ia vector:' )


!
! > We do not need row indices - information on rows is in ia (CSR format).
!
  deallocate ( iacoo )

!
! > Allocate array for sparse matrix element values (matrix is in CSR format).
!
  allocate( a(nnz) )

!
! > Source vectors
!
  allocate( su(numCells) ) 
  allocate( sv(numCells) )  
  allocate( sw(numCells)) 

! > Sources moved to LHS, i.e. to main diagonal.
!
  allocate( spu(numCells) )  
  allocate( spv(numCells) ) 
  allocate( sp(numCells) ) 

!
! > Residual vector
! 
  ! allocate(res(numCells) ) 

!
! > The inverses of main diagonal elements for u,v and w momentum equations.
!
  allocate(apu(numCells) )
  allocate(apv(numCells) ) 
  allocate(apw(numCells) ) 

!
! > Allocate array storing positions in the array for sparse matrix elements where specific cell pair coefs are stored.
!
  allocate( icell_jcell_csr_index( numInnerFaces+numPeriodic ) )
  allocate( jcell_icell_csr_index( numInnerFaces+numPeriodic ) )
  
!
! >Position arrays of matrix elements stored in CSR format
!
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    ! Index of the (icell,jcell) matrix element:
    icell_jcell_csr_index(i) = csr_to_k( ijp, ijn, ia, ja ) 

    ! Index of the (jcell,icell) matrix element:
    jcell_icell_csr_index(i) = csr_to_k( ijn, ijp, ia, ja ) 
  enddo

  ! Matrix elements resulting from connecting cells of the periodic boundaries

  iPer = 0 ! counts every instance of finding periodic boundary

  k = numInnerFaces ! counter for arrays filled below

  do ib=1,numBoundaries

    if (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        ijp = owner( startFace(ib) + i )
        ijn = owner( startFaceTwin(iPer) + i )

        k = k + 1 ! next

        ! Index of the (icell,jcell) matrix element:
        icell_jcell_csr_index(k) = csr_to_k( ijp, ijn, ia, ja ) 

        ! Index of the (jcell,icell) matrix element:
        jcell_icell_csr_index(k) = csr_to_k( ijn, ijp, ia, ja ) 


      enddo

    endif

  enddo


end subroutine


end module sparse_matrix