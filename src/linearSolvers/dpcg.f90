!***********************************************************************
!
subroutine dpcg(fi,ifi)
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
  use types 
  use parameters
  use geometry, only: numCells, numTotal
  use sparse_matrix
  use title_mod

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(dp), dimension(numTotal), intent(inout) :: fi 

!
! Local variables
!
  integer :: i, k, ns, l, itr_used
  real(dp), dimension(numCells) :: pk,zk
  real(dp) :: rsm, resmax, res0, resl
  real(dp) :: s0, sk, alf, bet, pkapk, tol
  real(dp) :: factor ! normalization factor for scaled residuals

! residual tolerance
  resmax = sor(ifi)
  tol=1e-13

  itr_used = 0

!
! Initalize working arrays
!
  pk = 0.0_dp
  zk = 0.0_dp
  res = 0.0_dp

!
! Calculate initial residual vector and the norm
!
  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi(ja(k)) 
    enddo
  enddo

  ! Normalization factor for scaled residuals
  factor = sum( abs( a(diag(1:numCells)) * fi(1:numCells) ))

  ! L1-norm of the residual
  res0=sum(abs(res))

  ! Initial normalized residual - for convergence report
  resor(ifi) = res0/(factor+small)
  
  if(res0.lt.tol) then
    write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(Jacobi):  Solving for ',trim(chvarSolver(ifi)), &
    ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations ',0
    return
  endif
!
! If ltest=true, print the norm 
!
  if(ltest) write(66,'(20x,a,1pe10.3)') 'res0 = ',res0

  s0=1.e20
!
! Start iterations
!
  ns=nsw(ifi)

  do l=1,ns
!
! Solve for zk(ijk) -- diagonal preconditioner
!
  do i=1,numCells
    zk(i) = res(i) / a( diag(i) )
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
  do i=1,numCells
    zk(i) = 0.0_dp 
    do k = ioffset(i),ioffset(i+1)-1
      zk(i) = zk(i) + a(k) * pk( ja(k) ) 
    enddo
  enddo

  ! Inner product
  pkapk=sum(pk*zk)

  alf=sk/pkapk

  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*pk

  ! Update residual vector
  res = res - alf*zk

  ! L^1-norm of residual
  resl = sum(abs(res))

  s0=sk

  itr_used = itr_used + 1
  
!
! Check convergence
!

  rsm = resl/(res0+small)
  if(ltest) write(6,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit

!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
  '  PCG(Jacobi):  Solving for ',trim(chvarSolver(ifi)), &
  ', Initial residual = ',resor(ifi),', Final residual = ',resl/(factor+small),', No Iterations ',itr_used

end subroutine
