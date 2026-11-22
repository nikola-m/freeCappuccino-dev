program test_amg
  use block_amg
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer :: nx, ny, n, nnz, i, j, c, k, it
  integer, allocatable :: ia(:), ja(:)
  real(dp), allocatable :: val(:,:,:), x(:,:), b(:,:), xe(:,:)
  type(amg_hierarchy) :: h
  real(dp) :: r0, rf, cc
  integer :: dx(5), dy(5), cnt
  character(len=8) :: sm
  call get_command_argument(1, sm)
  if (len_trim(sm)>0) amg_smoother = sm
  amg_verbose = 1
  do nx = 32, 128, 96
    ny = nx; n = nx*ny
    allocate(ia(n+1), ja(5*n), val(4,4,5*n), x(4,n), b(4,n), xe(4,n))
    dx = [0,-1,1,0,0]; dy = [0,0,0,-1,1]
    cc = 1.0_dp
    nnz = 0; val = 0.0_dp
    do j=1,ny; do i=1,nx
      c = i + (j-1)*nx; ia(c) = nnz+1
      do k=1,5
        if (i+dx(k)<1 .or. i+dx(k)>nx .or. j+dy(k)<1 .or. j+dy(k)>ny) cycle
        nnz = nnz+1; ja(nnz) = c + dx(k) + dy(k)*nx
        if (k==1) then
          val(1,1,nnz)=4.1; val(2,2,nnz)=4.1; val(3,3,nnz)=4.1; val(4,4,nnz)=4.0*cc*0.25+1e-3
        else
          val(1,1,nnz)=-1; val(2,2,nnz)=-1; val(3,3,nnz)=-1; val(4,4,nnz)=-cc*0.25
          ! pressure gradient in momentum, divergence in continuity (D=-G^T)
          if (dx(k)/=0) then
            val(1,4,nnz)=0.5*dx(k); val(4,1,nnz)=0.5*dx(k)
          else
            val(2,4,nnz)=0.5*dy(k); val(4,2,nnz)=0.5*dy(k)
          end if
        end if
      end do
    end do; end do
    ia(n+1) = nnz+1
    call random_number(xe); b = 0
    ! b = A xe
    do c=1,n
      do k=ia(c),ia(c+1)-1
        b(:,c) = b(:,c) + matmul(val(:,:,k), xe(:,ja(k)))
      end do
    end do
    call amg_setup(h, n, ia, ja, val, .true.)
    x = 0
    call amg_fgmres(h, x, b, 1.0d-10, 1.0d-14, 200, 30, r0, rf, it)
    print '(a,i0,a,i0,a,es9.2,a,es9.2,a,es9.2)','n=',n,' fgmres its=',it,' res0=',r0,' resf=',rf,' err=',maxval(abs(x-xe))
    deallocate(ia,ja,val,x,b,xe)
    h%built = .false.
  end do
end program
