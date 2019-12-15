module linear_solvers

!***********************************************************************
!
subroutine linear_solver(solver, fi, itr_max, tol_abs, tol_rel, chvar)
!
!***********************************************************************
!
!
!
  use types
  use parameters
  use geometry, only: numCells
  use sparse_matrix, only: nnz, ioffset, ia, a, diag

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(dp), dimension(numTotal), intent(inout) :: fi 
!
! Local variables
!

  if( solver .eq. xx ) then

    call dpcg( numCells, nnz, ioffset, ja, a, diag, fi, su, itr_max, tol_abs, tol_rel, chvar )

  elseif( solver .eq. xx ) then 

    call iccg( numCells, nnz, ioffset, ja, a, diag, fi, su, itr_max, tol_abs, tol_rel, chvar )

  elseif( solver .eq. xx ) then  

    call bicgstab( numCells, nnz, ioffset, ja, a, diag, fi, su, itr_max, tol_abs, tol_rel, chvar ) 

  elseif( solver .eq. xx ) then 

    call solve_csr( numCells, nnz, ioffset, ja, a, su, fi)

  elseif( solver .eq. xx ) then 

    call pmgmres_ilu ( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), ip, su, 100, 4, 1e-8, sor(ip) )
    
  endif


end subroutine

end module linear_solvers