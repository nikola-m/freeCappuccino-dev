!***********************************************************************
!
subroutine iccg(fi,ifi)
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
  use types 
  use parameters
  use geometry
  use my_mpi_module
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
  integer :: i, k, ns, l, itr_used, ib, iface, ijn, ipro
  real(dp), dimension(numCells) :: zk
  real(dp), dimension(numTotal) :: pk,d
  real(dp) :: rsm, resmax, res0, resl, tol
  real(dp) :: s0, sk, alf, bet, pkapk

! residual tolerance
  resmax = sor(ifi)
  tol = 1e-13

  itr_used = 0

!
! Initalize working arrays
!
  pk = 0.0_dp
  zk = 0.0_dp
  d = 0.0_dp
  res = 0.0_dp
!
! Calculate initial residual vector and the norm
!

  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi( ja(k) ) 
    enddo
  enddo

  ! Residual contribution from cells at processor boundary.
  ipro = 0
  do ib=1,numBoundaries
    if ( bctype(ib) == 'process' ) then
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        k = owner(iface)
        ijn = iBndValueStart(ib) + i
        ipro = ipro + 1

        res( k ) = res( k ) - apr( ipro )*fi( ijn )

      enddo
    endif
  enddo

  ! L1-norm of residual
  res0=sum(abs(res))

  call global_sum(res0)

  if(res0.lt.tol) then
    if (myid .eq. 0) write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(IC0):  Solving for ',trim(chvarSolver(ifi)), &
    ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations ',0
    return
  endif

  if(ltest) write(6,'(20x,a,1pe10.3)') 'res0 = ',res0

!
! Calculate elements of preconditioning matrix diagonal
!
  do i=1,numCells
    d(i) = a( diag(i) )
    do k = ioffset(i), diag(i)-1
      d(i) = d(i) - a( k ) * d( ja( k )) * a( k )
    end do
  enddo



  ! ! Contribution from cells at processor boundary.
  ! ipro = 0
  ! call exchange(d)
  
  ! do ib=1,numBoundaries

  !   if ( bctype(ib) == 'process' ) then

  !     do i=1,nfaces(ib)
        
  !       iface = startFace(ib) + i
  !       k = owner(iface)
  !       ijn = iBndValueStart(ib) + i
  !       ipro = ipro + 1

  !       d( k ) = d( k ) - apr( ipro )*d( ijn )*apr( ipro )

  !     enddo

  !   endif

  ! enddo

  do i=1,numCells
    d(i) =  1.0_dp / (d(i) + small)
  enddo


  s0=1.e20
!
! Start iterations
!
  ns=nsw(ifi)

  do l=1,ns
!
! Solve for zk(ijk) -- forward substitution
!
  do i=1,numCells
    zk(i) = res(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo


  do i=1,numCells
    zk(i) = zk(i) / (d(i)+small) 
  enddo 
!
! Backward substitution
!
  do i=numCells,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo

  
  ! Inner product
  sk = sum(res*zk) !..or  dot_product(res,zk)

  call global_sum(sk)

!
! Calculate beta
!
  bet=sk/s0

!
! Calculate new search vector pk
!
  pk(1:numCells) = zk + bet*pk(1:numCells)



!
! Calculate scalar product (pk.a pk) and alpha (overwrite zk)
!

  ! Matvec:
  do i=1,numCells
    zk(i) = 0.0_dp 
    do k = ioffset(i),ioffset(i+1)-1
      zk(i) = zk(i) + a(k) * pk( ja(k) ) 
    enddo
  enddo

  ! ! Processor boundaries
  ! ipro = 0

  !  call exchange( pk )

  ! do ib=1,numBoundaries

  !   if ( bctype(ib) == 'process' ) then

  !     do i=1,nfaces(ib)
  !       iface = startFace(ib) + i
  !       k = owner(iface)
  !       ijn = iBndValueStart(ib) + i
  !       ipro = ipro + 1

  !       zk( k ) = zk( k ) + apr( ipro ) * pk( ijn )

  !     enddo

  !   endif

  ! enddo

  ! Inner product
  pkapk=sum(pk*zk)

  call global_sum(pkapk)

  alf=sk/pkapk

  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*pk(1:numCells)

  ! Update residual vector
  res = res - alf*zk


  ! L1-norm of residual
  resl = sum(abs(res))
  
  call global_sum(resl)
  
  s0=sk

  itr_used = itr_used + 1

!
! Check convergence 
!
  if(l.eq.1) resor(ifi) = res0
  rsm = resl/(resor(ifi)+small)
  if(ltest) write(6,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit

!
! End of iteration loop
!
  end do

  ! MPI exchange
  call exchange( fi )

! Write linear solver report:
  if ( myid.eq.0 ) write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(IC0):  Solving for ',trim(chvarSolver(ifi)), &
  ', Initial residual = ',res0,', Final residual = ',resl,', No Iterations ',itr_used

end subroutine
