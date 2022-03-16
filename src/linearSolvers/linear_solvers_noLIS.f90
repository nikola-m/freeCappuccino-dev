module linear_solvers

!********************************************************************************************************************************132
!
!  Purpose: 
!
!    A module containing built-in linear solvers library for sparse linear systems.
!
!  Discussion:
!
!    Here we provide a general subroutine 'solve' which can receive as na argument, either a derived data type fvEquation,
!    where all the solution parameters are provided within the data type, or it can accept a CSR matrix with other
!    required parameters listed in a subroutine call. 
!    The list of available linear solvers inlcude: 
!      - Gauss-Seidel
!      - Incomplete Cholesky Preconditioned Conjugate Gradient - iccg
!      - Diagonally Preconditioned CG - dpcg
!      - Restarted GMRES with ILU preconditioner  - pmgmres_ilu,
!      - BiConjugate Gradient Stabilised - BiCGStab with ILU preconditioner.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use parameters
  use geometry, only: numCells,numTotal
  use sparse_matrix, only: nnz, ioffset, ja, a, diag

  implicit none

  public

contains


subroutine csrsolve(solver, fi, rhs, res0, itr_max, tol_abs, tol_rel, chvar)
!
!  Purpose:
!   Main routine for solution of sparse linear systems, adjusted to use global variables of freeCappuccino.
!
!  Discussion:
!    Call specific algorithm for iterative solution of linear systems based on given input.
!
  implicit none
!
! Parameters
!
  character( len = *) :: solver                       ! Char string for linear solver name.
  real(dp), dimension(numTotal), intent(inout) :: fi  ! On input-current field; on output-the solution vector
  real(dp), dimension(numCells), intent(in) :: rhs    ! The right hand side of the linear system
  real(dp), intent(out) :: res0                       ! Output - initial residual of the linear system in L1 norm 
  integer, intent(in) :: itr_max                      ! Maximum number of iterations
  real(dp), intent(in) :: tol_abs                     ! Absolute tolerance level for residual
  real(dp), intent(in) :: tol_rel                     ! Relative tolerance level for residual
  character( len=* ), intent(in) :: chvar             ! Character string containing name of the solved field, printed on stdout



  if( solver .eq. 'gauss-seidel') then

    call GaussSeidel( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), rhs, res0, itr_max, tol_abs, tol_rel, chvar, ltest )

  elseif( solver .eq. 'dpcg' ) then

    call dpcg( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), rhs, res0, itr_max, tol_abs, tol_rel, chvar, ltest )

  elseif( solver .eq. 'iccg' ) then 

    call iccg( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), rhs, res0, itr_max, tol_abs, tol_rel, chvar, ltest )

  elseif( solver .eq. 'bicgstab' ) then  

    call bicgstab( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), rhs, res0, itr_max, tol_abs, tol_rel, chvar, ltest ) 

  elseif( solver .eq. 'pmgmres' ) then 

    ! ** NOTE for GMRES(m) **
    ! Here we have hardcoded the restart parameter - m in restarted GMRES algorithm - GMRES(m), to m=4. 
    ! Play with this number if you want, the greater like m=20, better the convergence, but much more memory is required.

    call pmgmres_ilu( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), rhs, res0, itr_max, 4, tol_abs, tol_rel, chvar, ltest )

  else

    write(*,*) "  Error in linear solvers module: non-existing linear solver!"
    stop 


  endif


end subroutine


!***********************************************************************
!
subroutine GaussSeidel( n, nnz, ioffset, ja, a, diag, fi, rhs, resor, itr_max, tol_abs, tol_rel, chvar, verbose )
!
!***********************************************************************
!
!    This routine incorporates the  
!    Gauss-Seidel solver for sparse matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vinca.rs
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: n                              ! Number of unknowns, length of a solution vector
  integer, intent(in) :: nnz                            ! Number of non-zero elements in sparse matrix
  integer, dimension(n+1), intent(in) :: ioffset        ! The offsets of each row in coef array
  integer, dimension(nnz), intent(in) :: ja             ! Columns array
  real(dp), dimension(nnz), intent(in) :: a             ! Coefficient array
  integer, dimension(n), intent(in) :: diag             ! Position of diagonal elements in coeff. array
  real(dp), dimension(n), intent(inout) :: fi           ! On input-current field; on output-the solution vector
  real(dp), dimension(n), intent(in) :: rhs             ! The right hand side of the linear system
  real(dp), intent(out) :: resor                        ! Output - initial residual of the linear system in L1 norm
  integer, intent(in) :: itr_max                        ! Maximum number of iterations
  real(dp), intent(in) :: tol_abs                       ! Absolute tolerance level for residual
  real(dp), intent(in) :: tol_rel                       ! Relative tolerance level for residual
  character( len=* ), intent(in) :: chvar               ! Character string containing name of the solved field, printed on stdout
  logical, intent(in) :: verbose                        ! Boolean saying do we print residual of each iteration to stdout.

!
! Local variables
!
  integer :: i, k, l, itr_used
  real(dp) :: res0, rsm, resl, factor
  real(dp), dimension(:), allocatable :: res                         ! Residual vector

!
! Start iterations
!

  allocate( res(n))

  itr_used = 0
  factor = 0.0

  do l=1,itr_max

!
! Calculate initial residual vector and update the solution right away
!
  do i=1,n
    res(i) = rhs(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi( ja(k) ) 
    enddo
    fi(i) = fi(i) + res(i)/(a( diag(i) )+small)   
  enddo


! L1-norm of residual
  if(l.eq.1)  then
    res0=sum( abs(res) )

    if( res0.lt.tol_abs ) then
      write(*,'(3a,1PE10.3,a,1PE10.3,a)') '  Gauss-Seidel:  Solving for ',trim( chvar ), &
      ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations 1'
      return
    endif  

  endif

  if( verbose .and. l.eq.1 ) write(*,'(20x,a,1pe10.3)') 'res0 = ',res0  

  ! L1-norm of residual
  resl = sum( abs(res) )

  itr_used = itr_used + 1
  
!
! Check convergence
!
  if (l.eq.1) then
    factor = sum( abs( a(diag(1:n)) * fi(1:n) )) + small
    ! Initial normalized residual - for convergence report
    resor = res0/factor
  endif

  rsm = resl/(res0+small)

  if( verbose ) write(*,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',trim( chvar ),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if( rsm.lt.tol_rel .or. resl.lt.tol_abs ) exit ! Criteria for exiting iterations.


!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  Gauss-Seidel:  Solving for ',trim( chvar ), &
  ', Initial residual = ',resor,', Final residual = ',resl/factor,', No Iterations ',itr_used 

  deallocate( res )

end subroutine


!***********************************************************************
!
subroutine dpcg( n, nnz, ioffset, ja, a, diag, fi, rhs, resor, itr_max, tol_abs, tol_rel, chvar, verbose )
!
!***********************************************************************
!
!    This routine incorporates the diagonally preconditioned 
!    Conjugate Gradient solver for symmetric matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vin.bg.ac.rs
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: n                              ! Number of unknowns, length of a solution vector
  integer, intent(in) :: nnz                            ! Number of non-zero elements in sparse matrix
  integer, dimension(n+1), intent(in) :: ioffset        ! The offsets of each row in coef array
  integer, dimension(nnz), intent(in) :: ja             ! Columns array
  real(dp), dimension(nnz), intent(in) :: a             ! Coefficient array
  integer, dimension(n), intent(in) :: diag             ! Position of diagonal elements in coeff. array
  real(dp), dimension(n), intent(inout) :: fi           ! On input-current field; on output-the solution vector
  real(dp), dimension(n), intent(in) :: rhs             ! The right hand side of the linear system
  real(dp), intent(out) :: resor                        ! Output - initial residual of the linear system in L1 norm
  integer, intent(in) :: itr_max                        ! Maximum number of iterations
  real(dp), intent(in) :: tol_abs                       ! Absolute tolerance level for residual
  real(dp), intent(in) :: tol_rel                       ! Relative tolerance level for residual
  character( len=* ), intent(in) :: chvar               ! Character string containing name of the solved field, printed on stdout
  logical, intent(in) :: verbose                        ! Boolean saying do we print residual of each iteration to stdout.
!
! Local variables
!
  integer :: i, k, l, itr_used
  real(dp) :: res0, rsm, resl, factor
  real(dp) :: s0, sk, alf, bet, pkapk
  real(dp), dimension(:), allocatable :: res,pk,zk

!
! Initalize working arrays
!

  allocate( res(n), pk(n), zk(n) )

  factor = 0.0_dp
  pk = 0.0_dp
  zk = 0.0_dp
  res = 0.0_dp

!
! Calculate initial residual vector and the norm
!
  do i=1,n
    res(i) = rhs(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi( ja(k) ) 
    enddo
  enddo

  ! L1-norm of the residual
  res0=sum( abs(res) )
  
    if(res0.lt.tol_abs) then
      write(*,'(3a,1PE10.3,a,1PE10.3,a)') '  PCG(Jacobi):  Solving for ',trim( chvar ), &
      ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations 0'
      return
    endif
!
! If verbose==true, print the norm 
!
  if( verbose ) write(*,'(20x,a,1pe10.3)') 'res0 = ',res0

  s0=1.e20

  itr_used = 0

!
! Start iterations
!

  do l=1,itr_max
!
! Solve for zk(ijk) -- diagonal preconditioner
!
  do i=1,n
    zk(i) = res(i) / a( diag(i) )
  enddo
  
  ! Inner product
  sk = sum( res*zk ) !..or  dot_product(res,zk)

!
! Calculate beta
!
  bet=sk/s0

!
! Calculate new search vector pk
!
  pk = zk + bet*pk

!
! Calculate scalar product (pk.a pk) and alpha (overwrite zk)
!
  do i=1,n
    zk(i) = 0.0_dp 
    do k = ioffset(i),ioffset(i+1)-1
      zk(i) = zk(i) + a(k) * pk( ja(k) ) 
    enddo
  enddo

  ! Inner product
  pkapk=sum( pk*zk )

  alf=sk/pkapk

  ! Update solution vector
  fi(1:n) = fi(1:n) + alf*pk

  ! Update residual vector
  res = res - alf*zk

  ! L1-norm of residual - current value
  resl = sum( abs( res ) )

  s0=sk

  itr_used = itr_used + 1
  
!
! Check convergence
!
  if (l.eq.1) then
    factor = sum( abs( a(diag(1:n)) * fi(1:n) )) + small
    ! Initial normalized residual - for convergence report
    resor = res0/factor
  endif

  rsm = resl/(res0+small) ! Relative residual - current value

  if( verbose ) write(*,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',trim( chvar) ,' sweep = ',l,' resl = ',resl,' rsm = ',rsm

  if( rsm.lt.tol_rel .or. resl.lt.tol_abs ) exit ! Criteria for exiting iterations.

!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(Jacobi):  Solving for ',trim(chvar),', Initial residual = ',resor, &
  ', Final residual = ',resl/factor,', No Iterations ',itr_used

  deallocate( res, pk, zk )

end subroutine


!***********************************************************************
!
subroutine iccg( n, nnz, ioffset, ja, a, diag, fi, rhs, resor, itr_max, tol_abs, tol_rel, chvar, verbose )
!
!***********************************************************************
!
!    This routine incorporates the incomplete Cholesky preconditioned 
!    Conjugate Gradient solver for symmetric matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vinca.rs
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: n                              ! Number of unknowns, length of a solution vector
  integer, intent(in) :: nnz                            ! Number of non-zero elements in sparse matrix
  integer, dimension(n+1), intent(in) :: ioffset        ! The offsets of each row in coef array
  integer, dimension(nnz), intent(in) :: ja             ! Columns array
  real(dp), dimension(nnz), intent(in) :: a             ! Coefficient array
  integer, dimension(n), intent(in) :: diag             ! Position of diagonal elements in coeff. array
  real(dp), dimension(n), intent(inout) :: fi           ! On input-current field; on output-the solution vector
  real(dp), dimension(n), intent(in) :: rhs             ! The right hand side of the linear system
  real(dp), intent(out) :: resor                        ! Output - initial residual of the linear system in L1 norm
  integer, intent(in) :: itr_max                        ! Maximum number of iterations
  real(dp), intent(in) :: tol_abs                       ! Absolute tolerance level for residual
  real(dp), intent(in) :: tol_rel                       ! Relative tolerance level for residual
  character( len=* ), intent(in) :: chvar               ! Character string containing name of the solved field, printed on stdout
  logical, intent(in) :: verbose                        ! Boolean saying do we print residual of each iteration to stdout.

!
! Local variables
!
  integer :: i, k, l, itr_used
  real(dp) :: res0, rsm, resl, factor
  real(dp) :: s0, sk, alf, bet, pkapk
  real(dp), dimension(:), allocatable :: pk,zk,d
  real(dp), dimension(:), allocatable :: res                         ! Residual vector

!
! Initalize working arrays
!

  allocate( res(n), pk(n), zk(n), d(n) )

  factor = 0.0_dp
  pk = 0.0_dp
  zk = 0.0_dp
  d = 0.0_dp
  res = 0.0_dp

!
! Calculate initial residual vector and the norm
!

  do i=1,n
    res(i) = rhs(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi( ja(k) ) 
    enddo
  enddo

  ! L^1-norm of residual
  res0=sum( abs(res) )
  ! res0 = sqrt( sum( (res**2)/dble(n)) )

    ! Initially converged solution
    if( res0.lt.tol_abs ) then
      write(*,'(3a,1PE10.3,a,1PE10.3,a)') '  PCG(IC0):  Solving for ',trim( chvar ), &
      ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations 0'
      return
    endif  

  if( verbose ) write(*,'(20x,a,1pe10.3)') 'res0 = ',res0

!
! Calculate elements of diagonal preconditioning matrix
!
  do i=1,n
    d(i) = a( diag(i) )
    do k = ioffset(i), diag(i)-1
      d(i) = d(i) - a( k )**2 * d( ja( k )) 
    end do
    d(i) =  1.0_dp / d(i)
  enddo

  s0=1.e20

  itr_used = 0
  
!
! Start iterations
!
  do l=1,itr_max
!
! Solve for zk(ijk) -- forward substitution
!
  do i=1,n
    zk(i) = res(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo

  zk = zk/(d+small)     
!
! Backward substitution
!
  do i=n,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo
  
  ! Inner product
  sk = sum(res*zk) !..or  dot_product(res,zk)

!
! Calculate beta
!
  bet=sk/s0

!
! Calculate new search vector pk
!
  pk = zk + bet*pk

!
! Calculate scalar product (pk.a pk) and alpha (overwrite zk)
!
  do i=1,n
    zk(i) = 0.0_dp 
    do k = ioffset(i),ioffset(i+1)-1
      zk(i) = zk(i) + a(k) * pk( ja(k) ) 
    enddo
  enddo

  ! Inner product
  pkapk=sum(pk*zk)

  alf=sk/pkapk

  ! Update solution vector
  fi(1:n) = fi(1:n) + alf*pk

  ! Update residual vector
  res = res - alf*zk

  ! L1-norm of residual
  resl = sum( abs(res) )

  s0=sk

  itr_used = itr_used + 1
  
!
! Check convergence
!

  if (l.eq.1) then
    factor = sum( abs( a(diag(1:n)) * fi(1:n) )) + small
    ! Initial normalized residual - for convergence report
    resor = res0/factor
  endif

  rsm = resl/(res0+small)

  if( verbose ) write(*,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',trim( chvar),' sweep = ',l,' resl = ',resl,' rsm = ',rsm

  if( rsm.lt.tol_rel .or. resl.lt.tol_abs ) exit ! Criteria for exiting iterations.

!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(IC0):  Solving for ',trim( chvar ), &
  ', Initial residual = ',resor,', Final residual = ',resl/factor,', No Iterations ',itr_used

  deallocate( res, pk, zk, d )

end subroutine


subroutine bicgstab( n, nnz, ioffset, ja, a, diag, fi, rhs, resor, itr_max, tol_abs, tol_rel, chvar, verbose )
!
!***********************************************************************
!
! This routine incorporates the 
! BiCGStab - Bi-Conjugate Gradient Stabilised (Van der Vorst 1990)
! for sparse linear systems with matrices in CSR format.   
! The precodnitioning is done by incomplete LU factorization.
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: n                              ! Number of unknowns, length of a solution vector
  integer, intent(in) :: nnz                            ! Number of non-zero elements in sparse matrix
  integer, dimension(n+1), intent(in) :: ioffset        ! The offsets of each row in coef array
  integer, dimension(nnz), intent(in) :: ja             ! Columns array
  real(dp), dimension(nnz), intent(in) :: a             ! Coefficient array
  integer, dimension(n), intent(in) :: diag             ! Position of diagonal elements in coeff. array
  real(dp), dimension(n), intent(inout) :: fi           ! On input-current field; on output-the solution vector
  real(dp), dimension(n), intent(in) :: rhs             ! The right hand side of the linear system
  real(dp), intent(out) :: resor                        ! Output - initial residual of the linear system in L1 norm
  integer, intent(in) :: itr_max                        ! Maximum number of iterations
  real(dp), intent(in) :: tol_abs                       ! Absolute tolerance level for residual
  real(dp), intent(in) :: tol_rel                       ! Relative tolerance level for residual
  character( len=* ), intent(in) :: chvar               ! Character string containing name of the solved field, printed on stdout
  logical, intent(in) :: verbose                        ! Boolean saying do we print residual of each iteration to stdout.

!
!     Local variables
!
  integer :: i, k, l, itr_used
  real(dp) :: res0, rsm, resl, factor
  real(dp) :: alf, beto, gam, bet, om, ukreso
  real(dp), dimension(:), allocatable :: reso,pk,uk,zk,vk,d
  real(dp), dimension(:), allocatable :: res                 ! Residual vector


  allocate( res(n), reso(n), pk(n), uk(n), zk(n), vk(n), d(n) )

!
! Calculate initial residual vector and the norm
!

  do i=1,n
    res(i) = rhs(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi( ja(k) ) 
    enddo
  enddo

  ! L1-norm of residual
  res0=sum(abs(res))
  ! res0 = sqrt( sum( (res**2)/dble(n)) )

    if(res0.lt.tol_abs) then
      write(*,'(3a,1PE10.3,a,1PE10.3,a)') '  BiCGStab(ILU(0)):  Solving for ',trim(chvar), &
      ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations 0'
      return
    endif

  if( verbose ) write(*,'(20x,a,1pe10.3)') 'res0 = ',res0

!
! Calculate elements of diagonal preconditioning matrix
!
  do i=1,n
    d(i) = a( diag(i) )
    do k = ioffset(i), diag(i)-1
      do l = diag( ja(k) ), ioffset( ja( k )+1 )-1
        ! kolona u kojoj je trenutni dijagonalni element je i
        ! kada je pronadje izlazi sa poslednjom vrednosti l indeksa:
        if ( ja( l ) == i ) exit 
      end do
      d(i) = d(i) - a( k ) * d( ja( k )) * a( l )
    end do
    d(i) = 1.0_dp / d(i)
  enddo 

!
! Initialize working arrays and constants
!
  reso  = res    
  pk = 0.0_dp
  uk = 0.0_dp
  zk = 0.0_dp
  vk = 0.0_dp
  factor = 0.0_dp

  alf = 1.0_dp
  beto = 1.0_dp
  gam = 1.0_dp

  itr_used = 0

!
! Start iterations
!

  do l=1,itr_max

!
! Calculate beta and omega
!
  bet = sum(res*reso)
  om = bet*gam/(alf*beto+small)
  beto = bet

!
! Calculate pk
!
  pk = res + om*(pk -alf*uk)


!
! Solve (M ZK = PK) - forward substitution
!
  do i=1,n
    zk(i) = pk(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo

  zk = zk/(d+small)


!
! Backward substitution
!
  do i=n,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo


!
! Matvec 1: Uk = A*pk
!
  do i=1,n
    uk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      uk(i) = uk(i) +  a(k) * zk( ja(k) ) 
    enddo
  enddo


!
! Calculate scalar product uk*reso, and gamma
!
  ukreso = sum(uk*reso)
  gam = bet/ukreso


!
! Update 'fi' and calculate 'w' (overwrite 'res; - it is res-update)
!
  fi(1:n) = fi(1:n) + gam*zk
  res     = res     - gam*uk   ! <- W


!
! Solve (M Y = W); Y overwrites zk; forward substitution
!
  do i=1,n
    zk(i) = res(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k ) * zk( ja( k ) )
    end do
    zk(i) = zk(i)*d(i)
  enddo

  zk = zk/(d+small)

!
! Backward substitution
!
  do i=n,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ) )
    end do
    zk(i) = zk(i)*d(i)
  enddo

!
! Matvec 2: v = A*y (vk = A*zk); vk = csrMatVec(a,zk)
!
  do i=1,n
    vk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      vk(i) = vk(i) +  a(k) * zk( ja(k) ) 
    enddo
  enddo  

!
! Calculate alpha (alf)
!
  alf = sum(vk*res) / (sum(vk*vk)+small)

  ! Update solution vector
  fi(1:n) = fi(1:n) + alf*zk

  ! Update residual vector
  res = res - alf*vk

  ! L1-norm of residual
  resl = sum(abs(res))

  itr_used = itr_used + 1
  
!
! Check convergence
!
  if (l.eq.1) then
    factor = sum( abs( a(diag(1:n)) * fi(1:n) )) + small
    ! Initial normalized residual - for convergence report
    resor = res0/factor
  endif

  rsm = resl/(res0+small)

  if( verbose ) write(*,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',trim(chvar),' sweep = ',l,' resl = ',resl,' rsm = ',rsm

  if( rsm.lt.tol_rel .or. resl.lt.tol_abs ) exit ! Criteria for exiting iterations.

!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(*,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  BiCGStab(ILU(0)):  Solving for ',trim(chvar), &
  ', Initial residual = ',resor,', Final residual = ',resl/factor,', No Iterations ',itr_used

  deallocate( res, reso, pk, uk, zk, vk, d )

end subroutine


subroutine pmgmres_ilu ( n, nz_num, ia, ja, a, ua, x, rhs, resor, itr_max, mr, &
  tol_abs, tol_rel, chvar, verbose )

!*****************************************************************************80
!
!! PMGMRES_ILU applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2017
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!    Modification for usage in Cappuccino CFD code by Nikola Mirkov
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) AU(N), integer pointer to the matrix values at
!    diagonal.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Output, real ( kind = 8 ) RESOR, initial normalized residual (L1) norm.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
!    Input, character chvarSolver, a string with a field name which is solved
!    this is written in monitor file along number of iterations and relative
!    residual drop.
!
!    Input, logical VERBOSE, a boolen deciding wether residual history and  
!    iteration details will be printed.
  
  implicit none

  character( len = * ) :: chvar

  integer ( kind = 4 ), intent(in) :: mr
  integer ( kind = 4 ), intent(in) :: n
  integer ( kind = 4 ), intent(in) :: nz_num

  real ( kind = 8 ), intent(in) :: a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(mr+1)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(mr+1)
  real ( kind = 8 ) h(mr+1,mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ), intent(in) :: ia(n+1)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ), intent(in) :: itr_max
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ), intent(in) :: ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) l(ia(n+1)+1)
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rho
  real ( kind = 8 ), intent(out) :: resor
  real ( kind = 8 ) factor
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ), intent(in) :: rhs(n)
  real ( kind = 8 ) s(mr+1)
  real ( kind = 8 ), intent(in) :: tol_abs
  real ( kind = 8 ), intent(in) :: tol_rel
  integer ( kind = 4 ), intent(in) :: ua(n)
  real ( kind = 8 ) v(n,mr+1)
  logical verbose
  real ( kind = 8 ), intent(inout) :: x(n)
  real ( kind = 8 ) y(mr+1)


  itr_used = 0
  rho_tol = 0. ! To eliminate 'uninitialized' warning.
  k_copy = 0   ! To eliminate 'uninitialized' warning.
  factor = 0.

  ! NOTE: In our case the elements are already aranged
  ! call rearrange_cr ( n, nz_num, ia, ja, a )

  ! NOTE: We provide diagonal as input argument
  ! call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )


  iter_loop: do itr = 1, itr_max

    ! do i=1,n
    !   r(i) = rhs(i) 
    !   do k = ia(i),ia(i+1)-1
    !     r(i) = r(i) -  a(k) * x( ja(k) ) 
    !   enddo
    ! enddo

    do i=1,n
      r(i) = rhs(i) - sum( a( ia(i) : (ia(i+1)-1) ) * x( ja( ia(i) : (ia(i+1)-1) ) ) )
    end do

    ! call ax_cr ( n, nz_num, ia, ja, a, x, r )
    ! r(1:n) = rhs(1:n) - r(1:n)

    ! ! rho = sqrt ( dot_product ( r, r ) )
    rho = sum ( abs ( r ) )

    if ( itr == 1 ) then

      ! Normalization factor for scaled residuals
      factor = sum( abs( a( ua(1:n) ) * x(1:n) )) + 1e-20

      resor = rho / factor

      rho_tol = rho * tol_rel

    endif

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    mr_loop: do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .or. rho <= tol_abs ) then
        exit mr_loop
      end if

    end do mr_loop

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .or. rho <= tol_abs ) then
      exit iter_loop
    end if

  end do iter_loop

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  ! Write linear solver report:
  write(*,'(a,i0,3a,1PE10.3,a,1PE10.3,a,I0)') '  PMGMRES_ILU(',mr,'):  Solving for ',trim(chvar), &
  ', Initial residual = ',resor,', Final residual = ',rho/factor,', No Iterations ',itr_used


  return
end

subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

  return
end

subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end

subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

!*****************************************************************************80
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) tl
  integer ( kind = 4 ) ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

!*****************************************************************************80
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end

subroutine mult_givens ( c, s, k, g )

!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer ( kind = 4 ) K, indicates the location of the first
!    vector entry.
!
!    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = 8 ) c
  real ( kind = 8 ) g(1:k+1)
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end

subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

!*****************************************************************************80
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end
subroutine rearrange_cr ( n, nz_num, ia, ja, a )

!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
!
!    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) i4temp
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end

end module linear_solvers