module asimple_amg
!
!  Purpose:
!
!    Scalar algebraic multigrid (AMG) solver for the pressure-correction equation of the all-speeds SIMPLE algorithm
!    (module simple_allspeeds). Works on the CSR matrix of freeCappuccino (ia, ja, a), which is NOT symmetric for
!    compressible flows because of the convective (density-pressure) part of the equation.
!
!  Discussion:
!
!    * Coarsening : aggregation, strength of connection |a_ij| >= theta max_k |a_ik| (k /= i), greedy aggregates made of
!                   a seed and its unaggregated strong neighbours, leftovers join the strongest neighbouring aggregate.
!                   Piecewise constant prolongation P (P_iJ = 1 if cell i belongs to aggregate J), restriction R = P^T.
!    * Coarse op. : Galerkin  A_c = P^T A P, assembled by summing the fine entries into (agg(i), agg(j)).
!    * Smoother   : Gauss-Seidel, forward before and backward after the coarse grid correction (symmetric sweep).
!    * Cycle      : V(nu1,nu2); coarsest level (n <= coarse_size) solved directly by dense LU with partial pivoting.
!    * Outer      : right-preconditioned flexible GMRES(m), one V-cycle as preconditioner. Plain aggregation V-cycle is a
!                   poor stand-alone solver (piecewise constant P), as a preconditioner of a Krylov method it is robust.
!
!    The hierarchy is rebuilt at every call: the pressure-correction matrix changes in each SIMPLE iteration and the setup
!    cost is small compared to momentum assembly. Matrix of the finest level is copied, the arrays of the code are not touched.
!
!    This module is independent of block_amg (which is left untouched).
!
  use types

  implicit none

  private

  ! Parameters - modified through input.nml (asimple_amg_*)
  real(dp) :: pamg_theta       = 0.25_dp   ! Strength threshold
  integer  :: pamg_coarse_size = 40        ! Direct solve if number of unknowns <= this
  integer  :: pamg_nu1         = 1         ! Pre-smoothing sweeps
  integer  :: pamg_nu2         = 1         ! Post-smoothing sweeps
  integer  :: pamg_max_levels  = 25
  integer  :: pamg_restart     = 20        ! GMRES restart length
  logical  :: pamg_verbose     = .false.

  public :: pamg_theta, pamg_coarse_size, pamg_nu1, pamg_nu2, pamg_max_levels, pamg_restart, pamg_verbose
  public :: asimple_amg_solve

  type level_t
    integer :: n = 0, nnz = 0, nc = 0
    integer,  allocatable :: ia(:), ja(:), dg(:)  ! CSR, dg(i) = position of the diagonal element of row i
    real(dp), allocatable :: a(:)
    integer,  allocatable :: agg(:)               ! aggregate (coarse unknown) of every unknown of this level
    real(dp), allocatable :: x(:), b(:), r(:)
  end type level_t

  type(level_t), allocatable, save :: lev(:)
  integer, save :: nlev = 0
  real(dp), allocatable, save :: lu(:,:)          ! dense LU of the coarsest matrix
  integer,  allocatable, save :: piv(:)

contains


subroutine asimple_amg_solve( n, ia, ja, a, rhs, x, tolRel, tolAbs, maxit, res0, res, its )
!
!  Solve A x = rhs. On entry x is the initial guess. res0 and res are L1 norms of the initial and final residual.
!
  implicit none
  integer,  intent(in) :: n, maxit
  integer,  intent(in) :: ia(n+1), ja(*)
  real(dp), intent(in) :: a(*), rhs(n), tolRel, tolAbs
  real(dp), intent(inout) :: x(n)
  real(dp), intent(out) :: res0, res
  integer,  intent(out) :: its

  call amg_setup( n, ia, ja, a )
  call fgmres( rhs, x, tolRel, tolAbs, maxit, res0, res, its )

end subroutine


!=======================================================================================================================
!  SETUP
!=======================================================================================================================

subroutine amg_setup( n, ia, ja, a )
  implicit none
  integer,  intent(in) :: n, ia(n+1), ja(*)
  real(dp), intent(in) :: a(*)

  integer :: l, i, k, nnz

  if ( allocated(lev) ) deallocate( lev )
  allocate( lev(pamg_max_levels) )

  ! Finest level: copy of the code's matrix
  nnz = ia(n+1)-1
  call alloc_level( lev(1), n, nnz )
  lev(1)%ia = ia(1:n+1)
  lev(1)%ja = ja(1:nnz)
  lev(1)%a  = a(1:nnz)
  call find_diag( lev(1) )

  nlev = 1
  do l = 1, pamg_max_levels-1
    if ( lev(l)%n <= pamg_coarse_size ) exit
    call aggregate( lev(l) )
    ! Stop coarsening if it stalls
    if ( lev(l)%nc >= 0.9_dp*lev(l)%n .or. lev(l)%nc < 1 ) exit
    call galerkin( lev(l), lev(l+1) )
    nlev = l+1
  enddo

  ! Dense LU of the coarsest level
  i = lev(nlev)%n
  if ( allocated(lu) ) deallocate( lu, piv )
  allocate( lu(i,i), piv(i) )
  lu = 0.0_dp
  do i = 1, lev(nlev)%n
    do k = lev(nlev)%ia(i), lev(nlev)%ia(i+1)-1
      lu( i, lev(nlev)%ja(k) ) = lev(nlev)%a(k)
    enddo
  enddo
  call dense_lu( lev(nlev)%n, lu, piv )

  if ( pamg_verbose ) then
    write(*,'(a,i0,a)') '  asimple AMG: ', nlev, ' levels, unknowns per level:'
    write(*,'(2x,12i9)') ( lev(l)%n, l=1,nlev )
  endif

end subroutine


subroutine alloc_level( lv, n, nnz )
  implicit none
  type(level_t), intent(inout) :: lv
  integer, intent(in) :: n, nnz
  lv%n = n;  lv%nnz = nnz
  allocate( lv%ia(n+1), lv%ja(nnz), lv%dg(n), lv%a(nnz) )
  allocate( lv%x(n), lv%b(n), lv%r(n) )
  lv%x = 0.0_dp;  lv%b = 0.0_dp;  lv%r = 0.0_dp
end subroutine


subroutine find_diag( lv )
  implicit none
  type(level_t), intent(inout) :: lv
  integer :: i, k
  lv%dg = 0
  do i = 1, lv%n
    do k = lv%ia(i), lv%ia(i+1)-1
      if ( lv%ja(k) == i ) then
        lv%dg(i) = k
        exit
      endif
    enddo
  enddo
end subroutine


subroutine aggregate( lv )
!
!  Greedy aggregation based on strength of connection. Result in lv%agg and lv%nc.
!
  implicit none
  type(level_t), intent(inout) :: lv

  integer :: i, j, k, nc, jbest
  real(dp) :: amax, abest
  logical, allocatable :: strong(:)
  integer, allocatable :: nstrong(:)

  if ( allocated(lv%agg) ) deallocate( lv%agg )
  allocate( lv%agg(lv%n) )
  allocate( strong(lv%nnz) )
  lv%agg = 0

  ! Strong connections
  do i = 1, lv%n
    amax = 0.0_dp
    do k = lv%ia(i), lv%ia(i+1)-1
      if ( lv%ja(k) /= i ) amax = max( amax, abs(lv%a(k)) )
    enddo
    do k = lv%ia(i), lv%ia(i+1)-1
      strong(k) = ( lv%ja(k) /= i .and. abs(lv%a(k)) >= pamg_theta*amax .and. amax > 0.0_dp )
    enddo
  enddo

  ! Pass 1: seed + strong unaggregated neighbours (seeds whose strong neighbours are all free)
  nc = 0
  do i = 1, lv%n
    if ( lv%agg(i) /= 0 ) cycle
    ! seed only if none of its strong neighbours is already aggregated: keeps aggregates compact
    k = 0
    do j = lv%ia(i), lv%ia(i+1)-1
      if ( strong(j) .and. lv%agg( lv%ja(j) ) /= 0 ) k = k + 1
    enddo
    if ( k > 0 ) cycle
    nc = nc + 1
    lv%agg(i) = nc
    do j = lv%ia(i), lv%ia(i+1)-1
      if ( strong(j) ) lv%agg( lv%ja(j) ) = nc
    enddo
  enddo

  ! Pass 2: leftovers join the strongest neighbouring aggregate, isolated ones make their own
  do i = 1, lv%n
    if ( lv%agg(i) /= 0 ) cycle
    jbest = 0;  abest = -1.0_dp
    do j = lv%ia(i), lv%ia(i+1)-1
      k = lv%ja(j)
      if ( k /= i .and. lv%agg(k) /= 0 .and. abs(lv%a(j)) > abest ) then
        abest = abs(lv%a(j));  jbest = lv%agg(k)
      endif
    enddo
    if ( jbest == 0 ) then
      nc = nc + 1;  lv%agg(i) = nc
    else
      lv%agg(i) = jbest
    endif
  enddo

  lv%nc = nc
  deallocate( strong )
  if ( allocated(nstrong) ) deallocate( nstrong )

end subroutine


subroutine galerkin( fine, coarse )
!
!  Coarse matrix A_c = P^T A P for piecewise constant P: entries of fine rows are summed into (agg(i), agg(j)).
!
  implicit none
  type(level_t), intent(inout) :: fine, coarse

  integer :: i, k, I_, J_, nc, cnt, pos, m
  integer, allocatable :: first(:), members(:), marker(:), cols(:)
  real(dp), allocatable :: vals(:)
  integer, allocatable :: rowcount(:)

  nc = fine%nc

  ! Fine rows grouped by aggregate (counting sort)
  allocate( first(nc+1), members(fine%n), rowcount(nc) )
  rowcount = 0
  do i = 1, fine%n
    rowcount( fine%agg(i) ) = rowcount( fine%agg(i) ) + 1
  enddo
  first(1) = 1
  do I_ = 1, nc
    first(I_+1) = first(I_) + rowcount(I_)
  enddo
  rowcount = 0
  do i = 1, fine%n
    I_ = fine%agg(i)
    members( first(I_) + rowcount(I_) ) = i
    rowcount(I_) = rowcount(I_) + 1
  enddo

  ! Pass 1: count unique columns per coarse row
  allocate( marker(nc) )
  marker = 0
  cnt = 0
  do I_ = 1, nc
    do m = first(I_), first(I_+1)-1
      i = members(m)
      do k = fine%ia(i), fine%ia(i+1)-1
        J_ = fine%agg( fine%ja(k) )
        if ( marker(J_) /= I_ ) then
          marker(J_) = I_
          cnt = cnt + 1
        endif
      enddo
    enddo
  enddo

  call alloc_level( coarse, nc, cnt )

  ! Pass 2: fill. marker(J) holds the position of column J in the current row.
  allocate( cols(nc), vals(nc) )
  marker = 0
  coarse%ia(1) = 1
  pos = 0
  do I_ = 1, nc
    cnt = 0
    do m = first(I_), first(I_+1)-1
      i = members(m)
      do k = fine%ia(i), fine%ia(i+1)-1
        J_ = fine%agg( fine%ja(k) )
        if ( marker(J_) == 0 ) then
          cnt = cnt + 1
          marker(J_) = cnt
          cols(cnt) = J_
          vals(cnt) = fine%a(k)
        else
          vals( marker(J_) ) = vals( marker(J_) ) + fine%a(k)
        endif
      enddo
    enddo
    ! sort the columns of this row (insertion sort, rows are short) and store
    call sort_row( cnt, cols, vals )
    do m = 1, cnt
      pos = pos + 1
      coarse%ja(pos) = cols(m)
      coarse%a(pos)  = vals(m)
      marker( cols(m) ) = 0
    enddo
    coarse%ia(I_+1) = pos + 1
  enddo

  call find_diag( coarse )

  deallocate( first, members, rowcount, marker, cols, vals )

end subroutine


subroutine sort_row( n, cols, vals )
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: cols(n)
  real(dp), intent(inout) :: vals(n)
  integer :: i, j, c
  real(dp) :: v
  do i = 2, n
    c = cols(i);  v = vals(i)
    j = i - 1
    do while ( j >= 1 )
      if ( cols(j) <= c ) exit
      cols(j+1) = cols(j);  vals(j+1) = vals(j)
      j = j - 1
    enddo
    cols(j+1) = c;  vals(j+1) = v
  enddo
end subroutine


!=======================================================================================================================
!  CYCLE
!=======================================================================================================================

subroutine gs_forward( lv, nsweep )
  implicit none
  type(level_t), intent(inout) :: lv
  integer, intent(in) :: nsweep
  integer :: i, k, s
  real(dp) :: sum
  do s = 1, nsweep
    do i = 1, lv%n
      sum = lv%b(i)
      do k = lv%ia(i), lv%ia(i+1)-1
        if ( lv%ja(k) /= i ) sum = sum - lv%a(k)*lv%x( lv%ja(k) )
      enddo
      lv%x(i) = sum / lv%a( lv%dg(i) )
    enddo
  enddo
end subroutine


subroutine gs_backward( lv, nsweep )
  implicit none
  type(level_t), intent(inout) :: lv
  integer, intent(in) :: nsweep
  integer :: i, k, s
  real(dp) :: sum
  do s = 1, nsweep
    do i = lv%n, 1, -1
      sum = lv%b(i)
      do k = lv%ia(i), lv%ia(i+1)-1
        if ( lv%ja(k) /= i ) sum = sum - lv%a(k)*lv%x( lv%ja(k) )
      enddo
      lv%x(i) = sum / lv%a( lv%dg(i) )
    enddo
  enddo
end subroutine


recursive subroutine vcycle( l )
  implicit none
  integer, intent(in) :: l
  integer :: i, k
  real(dp) :: s

  if ( l == nlev ) then
    ! Coarsest level: direct solution
    lev(l)%x = lev(l)%b
    call dense_lu_solve( lev(l)%n, lu, piv, lev(l)%x )
    return
  endif

  ! Pre-smoothing from zero initial guess
  lev(l)%x = 0.0_dp
  call gs_forward( lev(l), pamg_nu1 )

  ! Residual restricted to the coarse level: b_c = P^T ( b - A x )
  lev(l+1)%b = 0.0_dp
  do i = 1, lev(l)%n
    s = lev(l)%b(i)
    do k = lev(l)%ia(i), lev(l)%ia(i+1)-1
      s = s - lev(l)%a(k)*lev(l)%x( lev(l)%ja(k) )
    enddo
    lev(l+1)%b( lev(l)%agg(i) ) = lev(l+1)%b( lev(l)%agg(i) ) + s
  enddo

  call vcycle( l+1 )

  ! Prolongation and correction
  do i = 1, lev(l)%n
    lev(l)%x(i) = lev(l)%x(i) + lev(l+1)%x( lev(l)%agg(i) )
  enddo

  ! Post-smoothing
  call gs_backward( lev(l), pamg_nu2 )

end subroutine


!=======================================================================================================================
!  FLEXIBLE GMRES
!=======================================================================================================================

subroutine matvec( x, y )
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: y(:)
  integer :: i, k
  real(dp) :: s
  do i = 1, lev(1)%n
    s = 0.0_dp
    do k = lev(1)%ia(i), lev(1)%ia(i+1)-1
      s = s + lev(1)%a(k)*x( lev(1)%ja(k) )
    enddo
    y(i) = s
  enddo
end subroutine


subroutine precond( r, z )
  implicit none
  real(dp), intent(in)  :: r(:)
  real(dp), intent(out) :: z(:)
  lev(1)%b = r
  call vcycle( 1 )
  z = lev(1)%x
end subroutine


subroutine fgmres( rhs, x, tolRel, tolAbs, maxit, res0, res, its )
  implicit none
  real(dp), intent(in) :: rhs(:), tolRel, tolAbs
  real(dp), intent(inout) :: x(:)
  integer,  intent(in) :: maxit
  real(dp), intent(out) :: res0, res
  integer,  intent(out) :: its

  integer :: n, m, i, j, k, iter
  real(dp), allocatable :: V(:,:), Z(:,:), H(:,:), cs(:), sn(:), g(:), y(:), r(:), w(:)
  real(dp) :: beta, tmp, target

  n = size(x)
  m = pamg_restart
  allocate( V(n,m+1), Z(n,m), H(m+1,m), cs(m), sn(m), g(m+1), y(m), r(n), w(n) )

  call matvec( x, w )
  r = rhs - w
  res0 = sum( abs(r) )
  res  = res0
  its  = 0
  target = max( tolAbs, tolRel*res0 )
  if ( res0 <= tolAbs ) then
    deallocate( V, Z, H, cs, sn, g, y, r, w )
    return
  endif

  iter = 0
  outer: do
    beta = sqrt( sum(r*r) )
    if ( beta <= 0.0_dp ) exit
    V(:,1) = r/beta
    g = 0.0_dp;  g(1) = beta
    H = 0.0_dp

    k = 0
    do j = 1, m
      iter = iter + 1
      k = j
      call precond( V(:,j), Z(:,j) )
      call matvec( Z(:,j), w )
      ! Modified Gram-Schmidt
      do i = 1, j
        H(i,j) = dot_product( w, V(:,i) )
        w = w - H(i,j)*V(:,i)
      enddo
      H(j+1,j) = sqrt( sum(w*w) )
      if ( H(j+1,j) > 1.0e-300_dp ) V(:,j+1) = w/H(j+1,j)
      ! Previous Givens rotations
      do i = 1, j-1
        tmp      =  cs(i)*H(i,j) + sn(i)*H(i+1,j)
        H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j)
        H(i,j)   = tmp
      enddo
      ! New rotation
      tmp = sqrt( H(j,j)**2 + H(j+1,j)**2 )
      if ( tmp > 0.0_dp ) then
        cs(j) = H(j,j)/tmp;  sn(j) = H(j+1,j)/tmp
      else
        cs(j) = 1.0_dp;  sn(j) = 0.0_dp
      endif
      H(j,j) = tmp;  H(j+1,j) = 0.0_dp
      g(j+1) = -sn(j)*g(j)
      g(j)   =  cs(j)*g(j)
      ! L2 estimate is used inside the loop; the L1 norm is evaluated at the end of the restart cycle
      if ( abs(g(j+1)) <= target/sqrt(real(n,dp)) .or. iter >= maxit ) exit
    enddo

    ! Solve upper triangular system and update
    do i = k, 1, -1
      y(i) = g(i)
      do j = i+1, k
        y(i) = y(i) - H(i,j)*y(j)
      enddo
      y(i) = y(i)/H(i,i)
    enddo
    do i = 1, k
      x = x + y(i)*Z(:,i)
    enddo

    call matvec( x, w )
    r = rhs - w
    res = sum( abs(r) )
    its = iter
    if ( res <= target .or. iter >= maxit ) exit outer
  enddo outer

  deallocate( V, Z, H, cs, sn, g, y, r, w )

end subroutine


!=======================================================================================================================
!  DENSE LU (coarsest level)
!=======================================================================================================================

subroutine dense_lu( n, A, piv )
  implicit none
  integer, intent(in) :: n
  real(dp), intent(inout) :: A(n,n)
  integer, intent(out) :: piv(n)
  integer :: i, j, k, p
  real(dp) :: amax, tmp, scale

  ! Pure Neumann problems (closed domain, incompressible) give a singular coarse matrix with a constant null space:
  ! a vanishing pivot is replaced by a huge one, which pins the corresponding unknown to zero.
  scale = 0.0_dp
  do k = 1, n
    scale = max( scale, abs(A(k,k)) )
  enddo

  do k = 1, n
    p = k;  amax = abs(A(k,k))
    do i = k+1, n
      if ( abs(A(i,k)) > amax ) then
        amax = abs(A(i,k));  p = i
      endif
    enddo
    piv(k) = p
    if ( p /= k ) then
      do j = 1, n
        tmp = A(k,j);  A(k,j) = A(p,j);  A(p,j) = tmp
      enddo
    endif
    if ( abs(A(k,k)) < 1.0e-10_dp*scale ) A(k,k) = 1.0e30_dp*scale
    do i = k+1, n
      A(i,k) = A(i,k)/A(k,k)
      do j = k+1, n
        A(i,j) = A(i,j) - A(i,k)*A(k,j)
      enddo
    enddo
  enddo
end subroutine


subroutine dense_lu_solve( n, A, piv, b )
  implicit none
  integer, intent(in) :: n, piv(n)
  real(dp), intent(in) :: A(n,n)
  real(dp), intent(inout) :: b(n)
  integer :: i, j
  real(dp) :: tmp
  do i = 1, n
    if ( piv(i) /= i ) then
      tmp = b(i);  b(i) = b(piv(i));  b(piv(i)) = tmp
    endif
    do j = 1, i-1
      b(i) = b(i) - A(i,j)*b(j)
    enddo
  enddo
  do i = n, 1, -1
    do j = i+1, n
      b(i) = b(i) - A(i,j)*b(j)
    enddo
    b(i) = b(i)/A(i,i)
  enddo
end subroutine

end module asimple_amg
