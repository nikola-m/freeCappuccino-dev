module block_amg
!
! Aggregation-based algebraic multigrid for 4x4 block-structured sparse systems
! (u,v,w,p unknowns per cell), used as preconditioner of a flexible GMRES
! iteration for the block-coupled pressure-velocity solver.
!
! * Matrix storage       : block CSR, 4x4 dense blocks, pattern taken from the scalar CSR
!                          of the mesh (sparse_matrix module) - columns are sorted per row.
! * Coarsening           : greedy aggregation on strength of connection measured with the
!                          Frobenius norm of D_i^{-1} A_ij  (scale invariant, so momentum and
!                          continuity rows can be scaled differently).
! * Transfer operators   : piecewise-constant prolongation (each aggregate -> 1 coarse block-row),
!                          restriction is its transpose, coarse operator by Galerkin sums
!                          A_c = R A P  - computed block-wise.
! * Smoothers            : symmetric block Gauss-Seidel, or block ILU(0).
! * Coarsest level      : dense LU with partial pivoting.
! * Cycle                : V (cycle_type=1) or W (cycle_type=2); an optional over-correction
!                          factor on the coarse-grid correction compensates for the poor
!                          approximation of plain aggregation.
! * Outer iteration      : right-preconditioned flexible GMRES(m).
!
  implicit none

  integer, parameter, private :: dp = kind(1.0d0)
  integer, parameter :: bs = 4          ! block size

  type bcsr
    integer :: n = 0
    integer :: nnz = 0
    integer, allocatable :: ia(:)        ! row pointers (n+1)
    integer, allocatable :: ja(:)        ! column (block) indices, sorted in each row
    integer, allocatable :: dg(:)        ! position of diagonal block in each row
    real(dp), allocatable :: val(:,:,:)  ! (bs,bs,nnz)
  end type

  type amg_level
    type(bcsr) :: A
    real(dp), allocatable :: dinv(:,:,:)    ! inverse of diagonal blocks (GS) or of U diagonal (ILU0)
    real(dp), allocatable :: lu(:,:,:)      ! ILU(0) factors, same pattern as A
    integer, allocatable :: agg(:)          ! fine row -> coarse row  (all but coarsest level)
    integer, allocatable :: emap(:)         ! fine nnz -> coarse nnz  (all but coarsest level)
    integer :: nc = 0
    real(dp), allocatable :: x(:,:), b(:,:), r(:,:), y(:,:)   ! work arrays (bs,n)
    real(dp), allocatable :: dense(:,:)     ! coarsest level LU
    integer, allocatable :: piv(:)
  end type

  type amg_hierarchy
    logical :: built = .false.
    integer :: nlev = 0
    integer, allocatable :: perm(:)
    type(amg_level), allocatable :: lev(:)
  end type

  ! Parameters (can be set from outside, defaults are reasonable)
  integer  :: amg_max_levels  = 30
  integer  :: amg_coarse_size = 16       ! switch to dense LU when n <= this
  real(dp) :: amg_theta       = 0.25_dp  ! strength threshold
  integer  :: amg_nu1         = 1        ! pre-smoothing sweeps
  integer  :: amg_nu2         = 1        ! post-smoothing sweeps
  integer  :: amg_cycle       = 2        ! 1=V, 2=W
  real(dp) :: amg_ccf         = 1.0_dp   ! coarse-grid correction factor
  character(len=8) :: amg_smoother = 'ilu0'   ! 'gs' or 'ilu0'
  integer  :: amg_verbose     = 0

  private :: lu4_inv, strength_aggregate, build_coarse_pattern, smooth_gs, smooth_ilu0, &
             amg_cycle_rec, bmv, residual, factor_ilu0, factor_diag, coarse_factor, coarse_solve, &
             bcsr_from_pattern, galerkin, alloc_work, mv4, mm4

contains

!=============================================================================
! 4x4 helpers
!=============================================================================

!=============================================================================
! 4x4 helpers (explicitly unrolled - far faster than the intrinsic matmul on tiny arrays)
!=============================================================================
  pure function mv4(A, x) result(y)
    real(dp), intent(in) :: A(bs,bs), x(bs)
    real(dp) :: y(bs)
    y(1) = A(1,1)*x(1) + A(1,2)*x(2) + A(1,3)*x(3) + A(1,4)*x(4)
    y(2) = A(2,1)*x(1) + A(2,2)*x(2) + A(2,3)*x(3) + A(2,4)*x(4)
    y(3) = A(3,1)*x(1) + A(3,2)*x(2) + A(3,3)*x(3) + A(3,4)*x(4)
    y(4) = A(4,1)*x(1) + A(4,2)*x(2) + A(4,3)*x(3) + A(4,4)*x(4)
  end function

  pure function mm4(A, B) result(C)
    real(dp), intent(in) :: A(bs,bs), B(bs,bs)
    real(dp) :: C(bs,bs)
    integer :: j
    do j=1,bs
      C(1,j) = A(1,1)*B(1,j) + A(1,2)*B(2,j) + A(1,3)*B(3,j) + A(1,4)*B(4,j)
      C(2,j) = A(2,1)*B(1,j) + A(2,2)*B(2,j) + A(2,3)*B(3,j) + A(2,4)*B(4,j)
      C(3,j) = A(3,1)*B(1,j) + A(3,2)*B(2,j) + A(3,3)*B(3,j) + A(3,4)*B(4,j)
      C(4,j) = A(4,1)*B(1,j) + A(4,2)*B(2,j) + A(4,3)*B(3,j) + A(4,4)*B(4,j)
    end do
  end function

  subroutine lu4_inv(a, ainv)
    ! Gauss-Jordan inverse with partial pivoting. Singular -> zero pivot guarded.
    real(dp), intent(in)  :: a(bs,bs)
    real(dp), intent(out) :: ainv(bs,bs)
    real(dp) :: m(bs,2*bs), t, f
    integer :: i, j, k, p
    m = 0.0_dp
    m(:,1:bs) = a
    do i=1,bs
      m(i,bs+i) = 1.0_dp
    end do
    do k=1,bs
      p = k
      do i=k+1,bs
        if (abs(m(i,k)) > abs(m(p,k))) p = i
      end do
      if (p /= k) then
        do j=1,2*bs
          t = m(k,j); m(k,j) = m(p,j); m(p,j) = t
        end do
      end if
      if (abs(m(k,k)) < 1.0e-300_dp) m(k,k) = 1.0e-300_dp
      t = 1.0_dp/m(k,k)
      m(k,:) = m(k,:)*t
      do i=1,bs
        if (i /= k) then
          f = m(i,k)
          if (f /= 0.0_dp) m(i,:) = m(i,:) - f*m(k,:)
        end if
      end do
    end do
    ainv = m(:,bs+1:2*bs)
  end subroutine

!=============================================================================
! Build block CSR from scalar CSR pattern; values supplied by caller later.
! ia,ja are 1-based; ja need not be sorted. Result has sorted columns.
!=============================================================================
  subroutine bcsr_from_pattern(A, n, ia, ja, perm)
    type(bcsr), intent(out) :: A
    integer, intent(in) :: n, ia(n+1), ja(:)
    integer, intent(out), optional :: perm(:)   ! perm(old position) = new (sorted) position
    integer, allocatable :: orig(:)
    integer :: i, k, l, t, to
    A%n = n
    A%nnz = ia(n+1) - 1
    allocate(A%ia(n+1), A%ja(A%nnz), A%dg(n), A%val(bs,bs,A%nnz))
    A%ia = ia
    A%ja = ja(1:A%nnz)
    allocate(orig(A%nnz))
    do k=1,A%nnz
      orig(k) = k
    end do
    do i=1,n
      ! insertion sort of columns in row i
      do k=ia(i)+1, ia(i+1)-1
        t = A%ja(k); to = orig(k); l = k-1
        do while (l >= ia(i))
          if (A%ja(l) <= t) exit
          A%ja(l+1) = A%ja(l); orig(l+1) = orig(l); l = l-1
        end do
        A%ja(l+1) = t; orig(l+1) = to
      end do
      A%dg(i) = 0
      do k=ia(i), ia(i+1)-1
        if (A%ja(k) == i) A%dg(i) = k
      end do
    end do
    A%val = 0.0_dp
    if (present(perm)) then
      do k=1,A%nnz
        perm(orig(k)) = k
      end do
    end if
  end subroutine

!=============================================================================
! y = A x
!=============================================================================
  subroutine bmv(A, x, y)
    type(bcsr), intent(in) :: A
    real(dp), intent(in)  :: x(bs,A%n)
    real(dp), intent(out) :: y(bs,A%n)
    integer :: i, k, j
    do i=1,A%n
      y(:,i) = 0.0_dp
      do k=A%ia(i), A%ia(i+1)-1
        j = A%ja(k)
        y(:,i) = y(:,i) + mv4(A%val(:,:,k), x(:,j))
      end do
    end do
  end subroutine

  subroutine residual(A, x, b, r)
    type(bcsr), intent(in) :: A
    real(dp), intent(in)  :: x(bs,A%n), b(bs,A%n)
    real(dp), intent(out) :: r(bs,A%n)
    integer :: i, k, j
    do i=1,A%n
      r(:,i) = b(:,i)
      do k=A%ia(i), A%ia(i+1)-1
        j = A%ja(k)
        r(:,i) = r(:,i) - mv4(A%val(:,:,k), x(:,j))
      end do
    end do
  end subroutine

!=============================================================================
! Aggregation
!=============================================================================
  subroutine strength_aggregate(A, theta, agg, nc)
    type(bcsr), intent(in) :: A
    real(dp), intent(in) :: theta
    integer, intent(out) :: agg(A%n), nc
    real(dp), allocatable :: s(:)          ! strength per nnz
    logical, allocatable :: strong(:)
    real(dp) :: dinv(bs,bs), smax, blk(bs,bs)
    integer :: i, k, j, best
    logical :: ok

    allocate(s(A%nnz), strong(A%nnz))
    s = 0.0_dp; strong = .false.
    do i=1,A%n
      call lu4_inv(A%val(:,:,A%dg(i)), dinv)
      smax = 0.0_dp
      do k=A%ia(i), A%ia(i+1)-1
        if (A%ja(k) == i) cycle
        blk = mm4(dinv, A%val(:,:,k))
        s(k) = sqrt(sum(blk*blk))
        smax = max(smax, s(k))
      end do
      do k=A%ia(i), A%ia(i+1)-1
        if (A%ja(k) == i) cycle
        if (s(k) >= theta*smax .and. s(k) > 0.0_dp) strong(k) = .true.
      end do
    end do

    agg = 0; nc = 0
    ! pass 1: seeds with all strong neighbours free
    do i=1,A%n
      if (agg(i) /= 0) cycle
      ok = .true.
      do k=A%ia(i), A%ia(i+1)-1
        if (strong(k) .and. agg(A%ja(k)) /= 0) then
          ok = .false.; exit
        end if
      end do
      if (.not. ok) cycle
      nc = nc + 1
      agg(i) = nc
      do k=A%ia(i), A%ia(i+1)-1
        if (strong(k)) agg(A%ja(k)) = nc
      end do
    end do
    ! pass 2: leftovers joined to strongest aggregated neighbour
    do i=1,A%n
      if (agg(i) /= 0) cycle
      best = 0; smax = -1.0_dp
      do k=A%ia(i), A%ia(i+1)-1
        j = A%ja(k)
        if (j == i .or. agg(j) == 0) cycle
        if (s(k) > smax) then
          smax = s(k); best = agg(j)
        end if
      end do
      if (best /= 0) agg(i) = best
    end do
    ! pass 3: isolated points become singleton aggregates
    do i=1,A%n
      if (agg(i) == 0) then
        nc = nc + 1; agg(i) = nc
      end if
    end do
    deallocate(s, strong)
  end subroutine

!=============================================================================
! Coarse pattern of Galerkin product (piecewise constant P) + fine->coarse map.
!=============================================================================
  subroutine build_coarse_pattern(Af, agg, nc, Ac, emap)
    type(bcsr), intent(in) :: Af
    integer, intent(in) :: agg(Af%n), nc
    type(bcsr), intent(out) :: Ac
    integer, intent(out) :: emap(Af%nnz)
    integer, allocatable :: cnt(:), mark(:), cols(:), tmpia(:)
    integer :: i, k, I_, J_, ncnnz, p, q
    integer, allocatable :: first(:), members(:)

    ! members of each aggregate (counting sort)
    allocate(first(nc+1), members(Af%n))
    first = 0
    do i=1,Af%n
      first(agg(i)+1) = first(agg(i)+1) + 1
    end do
    first(1) = 1
    do I_=1,nc
      first(I_+1) = first(I_) + first(I_+1)
    end do
    allocate(cnt(nc))
    cnt = 0
    do i=1,Af%n
      members(first(agg(i)) + cnt(agg(i))) = i
      cnt(agg(i)) = cnt(agg(i)) + 1
    end do

    ! pass 1: count distinct coarse columns per coarse row
    allocate(mark(nc), tmpia(nc+1))
    mark = 0
    tmpia(1) = 1
    do I_=1,nc
      ncnnz = 0
      do p=first(I_), first(I_+1)-1
        i = members(p)
        do k=Af%ia(i), Af%ia(i+1)-1
          J_ = agg(Af%ja(k))
          if (mark(J_) /= I_) then
            mark(J_) = I_; ncnnz = ncnnz + 1
          end if
        end do
      end do
      tmpia(I_+1) = tmpia(I_) + ncnnz
    end do

    allocate(cols(tmpia(nc+1)-1))
    mark = 0
    do I_=1,nc
      q = tmpia(I_)
      do p=first(I_), first(I_+1)-1
        i = members(p)
        do k=Af%ia(i), Af%ia(i+1)-1
          J_ = agg(Af%ja(k))
          if (mark(J_) /= I_) then
            mark(J_) = I_; cols(q) = J_; q = q + 1
          end if
        end do
      end do
    end do

    call bcsr_from_pattern(Ac, nc, tmpia, cols)

    ! emap: position of (I,J) in coarse row I
    mark = 0
    do I_=1,nc
      do q=Ac%ia(I_), Ac%ia(I_+1)-1
        mark(Ac%ja(q)) = q
      end do
      do p=first(I_), first(I_+1)-1
        i = members(p)
        do k=Af%ia(i), Af%ia(i+1)-1
          emap(k) = mark(agg(Af%ja(k)))
        end do
      end do
      do q=Ac%ia(I_), Ac%ia(I_+1)-1
        mark(Ac%ja(q)) = 0
      end do
    end do
    deallocate(first, members, cnt, mark, tmpia, cols)
  end subroutine

!=============================================================================
! Smoothers (operate on L%x with right-hand side L%b)
!=============================================================================
  subroutine factor_diag(L)
    type(amg_level), intent(inout) :: L
    integer :: i
    if (.not. allocated(L%dinv)) allocate(L%dinv(bs,bs,L%A%n))
    do i=1,L%A%n
      call lu4_inv(L%A%val(:,:,L%A%dg(i)), L%dinv(:,:,i))
    end do
  end subroutine

  subroutine factor_ilu0(L)
    ! Block ILU(0): L_ik = A_ik*Dinv_k (unit-lower in the scaled sense), U_ij updated A_ij,
    ! Dinv_i = inverse of the updated diagonal block.  Needs sorted columns.
    type(amg_level), intent(inout) :: L
    integer, allocatable :: pos(:)
    integer :: i, k, kk, j, jj, n
    real(dp) :: lik(bs,bs)
    n = L%A%n
    if (.not. allocated(L%lu)) allocate(L%lu(bs,bs,L%A%nnz))
    if (.not. allocated(L%dinv)) allocate(L%dinv(bs,bs,n))
    L%lu = L%A%val
    allocate(pos(n)); pos = 0
    do i=1,n
      do kk=L%A%ia(i), L%A%ia(i+1)-1
        pos(L%A%ja(kk)) = kk
      end do
      do kk=L%A%ia(i), L%A%ia(i+1)-1
        k = L%A%ja(kk)
        if (k >= i) exit
        lik = mm4(L%lu(:,:,kk), L%dinv(:,:,k))
        L%lu(:,:,kk) = lik
        do jj=L%A%dg(k)+1, L%A%ia(k+1)-1
          j = L%A%ja(jj)
          if (pos(j) /= 0) L%lu(:,:,pos(j)) = L%lu(:,:,pos(j)) - mm4(lik, L%lu(:,:,jj))
        end do
      end do
      call lu4_inv(L%lu(:,:,L%A%dg(i)), L%dinv(:,:,i))
      do kk=L%A%ia(i), L%A%ia(i+1)-1
        pos(L%A%ja(kk)) = 0
      end do
    end do
    deallocate(pos)
  end subroutine

  subroutine smooth_gs(L, nsweep)
    type(amg_level), intent(inout) :: L
    integer, intent(in) :: nsweep
    integer :: s, i, k, j
    real(dp) :: t(bs)
    do s=1,nsweep
      do i=1,L%A%n
        t = L%b(:,i)
        do k=L%A%ia(i), L%A%ia(i+1)-1
          j = L%A%ja(k)
          if (j /= i) t = t - mv4(L%A%val(:,:,k), L%x(:,j))
        end do
        L%x(:,i) = mv4(L%dinv(:,:,i), t)
      end do
      do i=L%A%n,1,-1
        t = L%b(:,i)
        do k=L%A%ia(i), L%A%ia(i+1)-1
          j = L%A%ja(k)
          if (j /= i) t = t - mv4(L%A%val(:,:,k), L%x(:,j))
        end do
        L%x(:,i) = mv4(L%dinv(:,:,i), t)
      end do
    end do
  end subroutine

  subroutine smooth_ilu0(L, nsweep)
    type(amg_level), intent(inout) :: L
    integer, intent(in) :: nsweep
    integer :: s, i, k, j
    real(dp) :: t(bs)
    do s=1,nsweep
      call residual(L%A, L%x, L%b, L%r)
      do i=1,L%A%n
        t = L%r(:,i)
        do k=L%A%ia(i), L%A%dg(i)-1
          t = t - mv4(L%lu(:,:,k), L%y(:,L%A%ja(k)))
        end do
        L%y(:,i) = t
      end do
      do i=L%A%n,1,-1
        t = L%y(:,i)
        do k=L%A%dg(i)+1, L%A%ia(i+1)-1
          j = L%A%ja(k)
          t = t - mv4(L%lu(:,:,k), L%r(:,j))
        end do
        L%r(:,i) = mv4(L%dinv(:,:,i), t)
      end do
      L%x = L%x + L%r
    end do
  end subroutine

!=============================================================================
! Dense coarsest-level LU
!=============================================================================
  subroutine coarse_factor(L)
    type(amg_level), intent(inout) :: L
    integer :: n, i, k, j, p, nn
    real(dp) :: t, f
    n = L%A%n; nn = bs*n
    if (allocated(L%dense)) deallocate(L%dense)
    if (.not. allocated(L%piv)) allocate(L%piv(nn))
    allocate(L%dense(nn,nn)); L%dense = 0.0_dp
    do i=1,n
      do k=L%A%ia(i), L%A%ia(i+1)-1
        j = L%A%ja(k)
        L%dense(bs*(i-1)+1:bs*i, bs*(j-1)+1:bs*j) = L%A%val(:,:,k)
      end do
    end do
    do k=1,nn
      p = k
      do i=k+1,nn
        if (abs(L%dense(i,k)) > abs(L%dense(p,k))) p = i
      end do
      L%piv(k) = p
      if (p /= k) then
        do j=1,nn
          t = L%dense(k,j); L%dense(k,j) = L%dense(p,j); L%dense(p,j) = t
        end do
      end if
      if (abs(L%dense(k,k)) < 1.0e-300_dp) L%dense(k,k) = 1.0e-300_dp
      do i=k+1,nn
        f = L%dense(i,k)/L%dense(k,k)
        L%dense(i,k) = f
        if (f /= 0.0_dp) L%dense(i,k+1:nn) = L%dense(i,k+1:nn) - f*L%dense(k,k+1:nn)
      end do
    end do
  end subroutine

  subroutine coarse_solve(L)
    type(amg_level), intent(inout) :: L
    real(dp) :: v(bs*L%A%n), t
    integer :: nn, k, i
    nn = bs*L%A%n
    v = reshape(L%b, [nn])
    do k=1,nn
      if (L%piv(k) /= k) then
        t = v(k); v(k) = v(L%piv(k)); v(L%piv(k)) = t
      end if
      do i=k+1,nn
        v(i) = v(i) - L%dense(i,k)*v(k)
      end do
    end do
    do k=nn,1,-1
      v(k) = (v(k) - dot_product(L%dense(k,k+1:nn), v(k+1:nn)))/L%dense(k,k)
    end do
    L%x = reshape(v, [bs, L%A%n])
  end subroutine

!=============================================================================
! Hierarchy construction
!=============================================================================
  subroutine alloc_work(L)
    type(amg_level), intent(inout) :: L
    allocate(L%x(bs,L%A%n), L%b(bs,L%A%n), L%r(bs,L%A%n), L%y(bs,L%A%n))
    L%x = 0.0_dp; L%b = 0.0_dp; L%r = 0.0_dp; L%y = 0.0_dp
  end subroutine

  subroutine galerkin(Lf, Lc)
    ! coarse values = sums of fine blocks via emap
    type(amg_level), intent(in) :: Lf
    type(amg_level), intent(inout) :: Lc
    integer :: k
    Lc%A%val = 0.0_dp
    do k=1,Lf%A%nnz
      Lc%A%val(:,:,Lf%emap(k)) = Lc%A%val(:,:,Lf%emap(k)) + Lf%A%val(:,:,k)
    end do
  end subroutine

  subroutine amg_setup(h, n, ia, ja, val, rebuild)
    ! Build (rebuild=.true. or first call) or update (values only) the hierarchy.
    ! n,ia,ja : scalar CSR pattern, 1-based, ja may be unsorted;  val(bs,bs,nnz) in ja order.
    type(amg_hierarchy), intent(inout) :: h
    integer, intent(in) :: n, ia(n+1), ja(:)
    real(dp), intent(in) :: val(bs,bs,*)
    logical, intent(in) :: rebuild
    integer :: l, nc, k, nnz
    integer, allocatable :: agg(:), emap(:)

    nnz = ia(n+1) - 1
    if (rebuild .or. .not. h%built) then
      if (allocated(h%lev)) deallocate(h%lev)
      allocate(h%lev(amg_max_levels))
      if (allocated(h%perm)) deallocate(h%perm)
      allocate(h%perm(nnz))
      call bcsr_from_pattern(h%lev(1)%A, n, ia, ja, h%perm)
      h%built = .false.
    end if
    do k=1,nnz
      h%lev(1)%A%val(:,:,h%perm(k)) = val(:,:,k)
    end do

    if (.not. h%built) then
      l = 1
      do
        call alloc_work(h%lev(l))
        if (h%lev(l)%A%n <= amg_coarse_size .or. l == amg_max_levels) exit
        allocate(agg(h%lev(l)%A%n))
        call strength_aggregate(h%lev(l)%A, amg_theta, agg, nc)
        if (nc > 0.85_dp*h%lev(l)%A%n .or. nc < 1) then      ! coarsening stalled
          deallocate(agg); exit
        end if
        allocate(emap(h%lev(l)%A%nnz))
        call build_coarse_pattern(h%lev(l)%A, agg, nc, h%lev(l+1)%A, emap)
        call move_alloc(agg, h%lev(l)%agg)
        call move_alloc(emap, h%lev(l)%emap)
        h%lev(l)%nc = nc
        call galerkin(h%lev(l), h%lev(l+1))
        l = l + 1
      end do
      h%nlev = l
      h%built = .true.
    else
      do l=1,h%nlev-1
        call galerkin(h%lev(l), h%lev(l+1))
      end do
    end if

    do l=1,h%nlev-1
      if (amg_smoother == 'ilu0') then
        call factor_ilu0(h%lev(l))
      else
        call factor_diag(h%lev(l))
      end if
    end do
    call coarse_factor(h%lev(h%nlev))
    if (amg_verbose > 0 .and. rebuild) call amg_report(h)
  end subroutine

  subroutine amg_report(h)
    type(amg_hierarchy), intent(in) :: h
    integer :: l
    real(dp) :: oc
    oc = 0.0_dp
    do l=1,h%nlev
      oc = oc + real(h%lev(l)%A%nnz, dp)
    end do
    write(*,'(a,i0,a,f6.3)') '  block AMG: levels = ', h%nlev, ',  operator complexity = ', oc/real(h%lev(1)%A%nnz,dp)
    do l=1,h%nlev
      write(*,'(a,i2,a,i9,a,i10)') '    level', l, '  rows', h%lev(l)%A%n, '  blocks', h%lev(l)%A%nnz
    end do
  end subroutine

!=============================================================================
! Cycle:  x <- x + coarse correction + smoothing, on level l (x,b in lev(l)%x,%b)
!=============================================================================
  recursive subroutine amg_cycle_rec(h, l)
    type(amg_hierarchy), intent(inout) :: h
    integer, intent(in) :: l
    integer :: i, g, nrep

    if (l == h%nlev) then
      call coarse_solve(h%lev(l))
      return
    end if
    if (amg_smoother == 'ilu0') then
      call smooth_ilu0(h%lev(l), amg_nu1)
    else
      call smooth_gs(h%lev(l), amg_nu1)
    end if
    call residual(h%lev(l)%A, h%lev(l)%x, h%lev(l)%b, h%lev(l)%r)
    h%lev(l+1)%b = 0.0_dp
    do i=1,h%lev(l)%A%n
      h%lev(l+1)%b(:,h%lev(l)%agg(i)) = h%lev(l+1)%b(:,h%lev(l)%agg(i)) + h%lev(l)%r(:,i)
    end do
    h%lev(l+1)%x = 0.0_dp
    nrep = 1
    if (amg_cycle >= 2 .and. l+1 < h%nlev) nrep = 2
    do g=1,nrep
      call amg_cycle_rec(h, l+1)
    end do
    do i=1,h%lev(l)%A%n
      h%lev(l)%x(:,i) = h%lev(l)%x(:,i) + amg_ccf*h%lev(l+1)%x(:,h%lev(l)%agg(i))
    end do
    if (amg_smoother == 'ilu0') then
      call smooth_ilu0(h%lev(l), amg_nu2)
    else
      call smooth_gs(h%lev(l), amg_nu2)
    end if
  end subroutine

  subroutine amg_precond(h, v, z)
    ! z = M^{-1} v  (one cycle from zero initial guess)
    type(amg_hierarchy), intent(inout) :: h
    real(dp), intent(in)  :: v(bs,h%lev(1)%A%n)
    real(dp), intent(out) :: z(bs,h%lev(1)%A%n)
    h%lev(1)%x = 0.0_dp
    h%lev(1)%b = v
    call amg_cycle_rec(h, 1)
    z = h%lev(1)%x
  end subroutine

!=============================================================================
! Right-preconditioned flexible GMRES(m) with AMG cycle as preconditioner
!   x  : initial guess on entry, solution on exit
!   returns initial/final residual 2-norms and iteration count
!=============================================================================
  subroutine amg_fgmres(h, x, b, tol_rel, tol_abs, maxit, restart, res0, resf, iters)
    type(amg_hierarchy), intent(inout) :: h
    real(dp), intent(inout) :: x(:,:)
    real(dp), intent(in)  :: b(:,:)
    real(dp), intent(in)  :: tol_rel, tol_abs
    integer, intent(in)   :: maxit, restart
    real(dp), intent(out) :: res0, resf
    integer, intent(out)  :: iters
    real(dp), allocatable :: V(:,:,:), Z(:,:,:), r(:,:), w(:,:)
    real(dp), allocatable :: H_(:,:), cs(:), sn(:), g(:), yv(:)
    real(dp) :: beta, tmp, target_res
    integer :: n, m, j, i, kk, itot

    n = h%lev(1)%A%n
    m = max(1, restart)
    allocate(V(bs,n,m+1), Z(bs,n,m), r(bs,n), w(bs,n))
    allocate(H_(m+1,m), cs(m), sn(m), g(m+1), yv(m))

    call residual(h%lev(1)%A, x, b, r)
    beta = sqrt(sum(r*r))
    res0 = beta
    resf = beta
    iters = 0
    target_res = max(tol_abs, tol_rel*res0)
    if (beta <= target_res) return

    itot = 0
    outer: do
      V(:,:,1) = r/beta
      g = 0.0_dp; g(1) = beta
      H_ = 0.0_dp
      kk = 0
      do j=1,m
        call amg_precond(h, V(:,:,j), Z(:,:,j))
        call bmv(h%lev(1)%A, Z(:,:,j), w)
        do i=1,j                                     ! modified Gram-Schmidt
          H_(i,j) = sum(w*V(:,:,i))
          w = w - H_(i,j)*V(:,:,i)
        end do
        H_(j+1,j) = sqrt(sum(w*w))
        if (H_(j+1,j) > 1.0e-300_dp) V(:,:,j+1) = w/H_(j+1,j)
        do i=1,j-1                                   ! apply previous Givens rotations
          tmp      =  cs(i)*H_(i,j) + sn(i)*H_(i+1,j)
          H_(i+1,j)= -sn(i)*H_(i,j) + cs(i)*H_(i+1,j)
          H_(i,j)  = tmp
        end do
        tmp = sqrt(H_(j,j)**2 + H_(j+1,j)**2)
        if (tmp < 1.0e-300_dp) then
          cs(j) = 1.0_dp; sn(j) = 0.0_dp
        else
          cs(j) = H_(j,j)/tmp; sn(j) = H_(j+1,j)/tmp
        end if
        H_(j,j)   = cs(j)*H_(j,j) + sn(j)*H_(j+1,j)
        H_(j+1,j) = 0.0_dp
        g(j+1) = -sn(j)*g(j)
        g(j)   =  cs(j)*g(j)
        itot = itot + 1
        kk = j
        resf = abs(g(j+1))
        if (resf <= target_res .or. itot >= maxit) exit
      end do
      ! solve upper-triangular system and update x
      do i=kk,1,-1
        yv(i) = g(i)
        do j=i+1,kk
          yv(i) = yv(i) - H_(i,j)*yv(j)
        end do
        yv(i) = yv(i)/H_(i,i)
      end do
      do j=1,kk
        x = x + yv(j)*Z(:,:,j)
      end do
      call residual(h%lev(1)%A, x, b, r)
      beta = sqrt(sum(r*r))
      resf = beta
      if (beta <= target_res .or. itot >= maxit) exit outer
    end do outer
    iters = itot
    deallocate(V, Z, r, w, H_, cs, sn, g, yv)
  end subroutine

  subroutine amg_matvec(h, x, y)
    type(amg_hierarchy), intent(in) :: h
    real(dp), intent(in)  :: x(bs,h%lev(1)%A%n)
    real(dp), intent(out) :: y(bs,h%lev(1)%A%n)
    call bmv(h%lev(1)%A, x, y)
  end subroutine

end module block_amg
