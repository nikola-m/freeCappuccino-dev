module gradients
!
! Module for cell centered gradients, gradient limiters and surface normal gradients.
!

use types
use geometry
use sparse_matrix, only: ia,ja,diag

implicit none

logical :: lstsq, lstsq_qr, lstsq_dm, gauss, node_gauss  ! Gradient discretization approach
logical :: lstsq_fw = .false.   ! Weighted LSQ with precomputed face weights (any cell shape)
character(len=20) :: limiter  ! Gradient limiter. Options: none, Barth-Jespersen, Venkatakrishnan, R3, multidimensional

real(dp),dimension(:,:), allocatable ::  Dmat      !  Dmat(6,numCells): unique entries (11,12,13,22,23,33) of the inverse of the symmetric LSQ matrix
real(dp),dimension(:,:), allocatable ::  lsq_cp, lsq_cn, lsq_cb  ! lstsq_fw: face weights for owner, neighbour, boundary face
real(dp),dimension(:,:,:), allocatable ::  D       !  when using qr version of the subroutine size(3,6,nxyz)!

real(dp), parameter :: small = 1e-20
real(dp), parameter :: zero = 0.0_dp



interface grad
  module procedure grad_scalar_field
  module procedure grad_vector_field
  module procedure grad_scalar_field_w_option  
end interface

interface sngrad
  module procedure sngrad
  module procedure sngrad2
  module procedure sngrad3
end interface

private

public :: lstsq, lstsq_qr, lstsq_dm, lstsq_fw, gauss, node_gauss, limiter
public :: grad, grad_gauss, grad_gauss_corrected, sngrad, create_lsq_grad_matrix



contains


!***********************************************************************
!
subroutine allocate_lsq_grad_matrix
!
!***********************************************************************
!
implicit none
  
  integer :: ierr

  if( lstsq .or. lstsq_dm ) then

    allocate( Dmat(6,numCells), stat=ierr ) 
      if(ierr /= 0)write(*,*)"allocation error: Dmat"

  elseif( lstsq_qr ) then

    allocate( D(3,6,numCells), stat=ierr ) 
      if(ierr /= 0)write(*,*)"allocation error: D"
      
  endif

end subroutine


!***********************************************************************
!
subroutine create_lsq_grad_matrix
!
!***********************************************************************
!
!  Discussion:
!    Prepare System Matrix For Least-Squares Gradient Calculation.
!    It is done only once at the beggiing of the calculation.
!
!***********************************************************************
!

  implicit none

  call allocate_lsq_grad_matrix

  if (lstsq) then

    call create_matrix_lsq

  elseif (lstsq_qr) then 

    call create_matrix_lsq_qr

  elseif (lstsq_dm) then 

    call create_matrix_lsq_dm

  elseif (lstsq_fw) then 

    call create_lsq_face_weights

  endif 

end subroutine


!***********************************************************************
!
subroutine grad_scalar_field(phi,dPhidxi)
!
!***********************************************************************
!

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(3,numTotal), intent(inout) :: dPhidxi

  dPhidxi = 0.0_dp
  
  if (lstsq) then 

    call grad_lsq(phi,dPhidxi)

  elseif (lstsq_qr) then 

    call grad_lsq_qr(phi,dPhidxi)

  elseif (lstsq_dm) then 

    call grad_lsq_dm(phi,dPhidxi)

  elseif (lstsq_fw) then 

    call grad_lsq_fw(phi,dPhidxi)

  elseif ( node_gauss ) then

    call grad_node_based_gauss(phi,dPhidxi)

  else ! if(gauss) then

    call grad_gauss(phi,dPhidxi)

  endif 

  
  ! > Gradient limiter:

  select case ( limiter ) 

    case( 'Barth-Jespersen' )
      call slope_limiter_Barth_Jespersen(phi, dPhidxi)

    case( 'Venkatakrishnan' )
      call slope_limiter_Venkatakrishnan(phi, dPhidxi)

    case( 'R3' )
      call slope_limiter_r345(phi, dPhidxi)

    case( 'multidimensional' )
      call slope_limiter_multidimensional(phi, dPhidxi)

    case default 
      ! 'no-limit' or 'none'
      
  end select


end subroutine


!***********************************************************************
!
subroutine grad_vector_field(U,V,W,dUdxi,dVdxi,dWdxi)
!
!***********************************************************************
!

  implicit none

  real(dp), dimension(numTotal), intent(in) :: U,V,W
  real(dp), dimension(3,numTotal), intent(inout) :: dUdxi,dVdxi,dWdxi

  dUdxi=0.0_dp
  dVdxi=0.0_dp
  dWdxi=0.0_dp

  if (lstsq) then

    call grad_lsq(U,dUdxi)
    call grad_lsq(V,dVdxi)
    call grad_lsq(W,dWdxi)

  elseif (lstsq_qr) then

    call grad_lsq_qr(U,dUdxi)
    call grad_lsq_qr(V,dVdxi)
    call grad_lsq_qr(W,dWdxi)

  elseif (lstsq_dm) then

    call grad_lsq_dm(U,dUdxi)
    call grad_lsq_dm(V,dVdxi)
    call grad_lsq_dm(W,dWdxi)

  elseif (lstsq_fw) then 

    call grad_lsq_fw(U,dUdxi)
    call grad_lsq_fw(V,dVdxi)
    call grad_lsq_fw(W,dWdxi)

  elseif ( node_gauss ) then

    call grad_node_based_gauss(u,dUdxi)
    call grad_node_based_gauss(v,dVdxi)
    call grad_node_based_gauss(w,dWdxi)

  else ! gauss gradient

    call grad_gauss(U,dUdxi)
    call grad_gauss(V,dVdxi)
    call grad_gauss(W,dWdxi)

  endif

end subroutine


!***********************************************************************
!
subroutine grad_scalar_field_w_option(phi,dPhidxi,option,option_limiter)
!
!***********************************************************************
!
! The main reason why we write this subroutine is to correct velocities
! in SIMPLE algorithm with conservative gradients, which is possible
! with Gauss rule. 
! We noticed it is better for calculation precision.
! But calling gradients with option may be nice anyhow.
!
!***********************************************************************
!

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(3,numTotal), intent(inout) :: dPhidxi
  character( len=* ), intent(in) :: option
  character( len=* ), intent(in), optional :: option_limiter  ! default: no limiting

  ! 'gauss_corrected' starts from the gradient that is passed in
  if ( option /= 'gauss_corrected' ) dPhidxi = 0.0_dp
  
  if ( option == 'lsq' ) then 

    call grad_lsq(phi,dPhidxi)

  elseif ( option == 'lsq_qr' ) then 

    call grad_lsq_qr(phi,dPhidxi)

  elseif ( option == 'wlsq' ) then 

    call grad_lsq_dm(phi,dPhidxi)

  elseif ( option == 'lsq_fw' ) then

    call grad_lsq_fw(phi,dPhidxi)

  elseif ( option == 'gauss' ) then

    call grad_gauss(phi,dPhidxi)

  elseif ( option == 'gauss_corrected' ) then

    call grad_gauss_corrected(phi,dPhidxi)

  endif 

  !
  ! Gradient limiter:
  !
  if ( .not. present(option_limiter) ) return

  select case ( option_limiter ) 

    case( 'Barth-Jespersen' )
      call slope_limiter_Barth_Jespersen(phi, dPhidxi)

    case( 'Venkatakrishnan' )
      call slope_limiter_Venkatakrishnan(phi, dPhidxi)

    case( 'R3' )
      call slope_limiter_r345(phi, dPhidxi)

    case( 'multidimensional' )
      call slope_limiter_multidimensional(phi, dPhidxi)

    case default 
      ! 'no-limit' or 'none'
      
  end select

end subroutine




!***********************************************************************
!
subroutine slope_limiter_Barth_Jespersen(phi, dPhidxi)
!
!***********************************************************************
!
!     Calculates slope limiter and appiies to scalar gradient:
!     Barth and Jespersen slope limiter:
!
!     AIAA-89-0366, The design and application of upwind schemes
!     on unstructured meshes, T.J.Barth, D.C.Jespersen, 1989.
!
!***********************************************************************
!

  implicit none

  ! Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numTotal) :: dPhidxi


  ! Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  ! real(dp) :: fimax,fimin
  real(dp) :: deltamax,deltamin


  !fimin = minval(phi(1:numCells))
  !fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ia(inp) ))
    phi_min = phi(ja( ia(inp) ))

    do k=ia(inp)+1, ia(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_min, phi(ja(k)) )      
    enddo

    deltamax = phi_max - phi(inp)
    deltamin = phi_min - phi(inp)

    slopelimit = 1.0_dp

    do k=ia(inp), ia(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      delta_face=dPhidxi(1,ijp)*(xc(ijn)-xc(ijp))+dPhidxi(2,ijp)*(yc(ijn)-yc(ijp))+dPhidxi(3,ijp)*(zc(ijn)-zc(ijp)) 


     if( abs(delta_face) < 1.e-6 )then
       r = 1.0_dp
     else if( delta_face > 0.0 )then
       r = deltamax/delta_face
     else
       r = deltamin/delta_face
     endif

     slopelimit = min( slopelimit , r )

    enddo

    dPhidxi(:,inp) = slopelimit*dPhidxi(:,inp)

  enddo

end subroutine




!***********************************************************************
!
subroutine slope_limiter_Venkatakrishnan(phi, dPhidxi)
!
!***********************************************************************
!
!     Calculates slope limiter and applies it to the scalar gradient:
!     Venkatakrishnan slope limiter:
!
!    AIAA-93-0880, On the accuracy of limiters and convergence
!    to steady state solutions, V.Venkatakrishnan, 1993
!
!     Differences to the previous version:
!     - bounds are the LOCAL extrema over the cell and its neighbours
!       (the old code used global field extrema, so it practically never
!       limited, and computed the local ones without using them),
!     - the smooth form with eps^2 is used, so the limiter is differentiable
!       near extrema (better iterative convergence),
!     - eps^2 = (K*h/L)^3 * (phi_max-phi_min)^2 with h = vol^(1/3) and
!       L^3 = total volume: invariant to units and to the size of the
!       domain. The old absolute threshold |delta|<1e-6 is gone,
!     - face value is estimated at the midpoint of the line P-N.
!
!     K is a tuning constant; larger K limits less (more accurate, less
!     monotone). K=0.01..0.2 are recommended values by Wang (https://doi.org/10.2514/6.1996-2091)
!
!***********************************************************************
!
  implicit none

  real(dp), dimension(numTotal) :: phi
  real(dp), dimension(3,numTotal) :: dPhidxi

  real(dp), parameter :: K = 0.01_dp ! values in the range [0.01 .. 0.2_dp]

  integer  :: inp, ijn, k_nb
  real(dp) :: phi_p, phi_max, phi_min
  real(dp) :: fimax, fimin, range, voltot
  real(dp) :: d1, d2, eps2, philim, slopelimit

  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))
  range = fimax - fimin

  ! Uniform field: nothing to limit
  if ( range <= 0.0_dp ) return

  voltot = sum( vol(1:numCells) )

  do inp = 1, numCells

    phi_p   = phi(inp)
    phi_max = phi_p
    phi_min = phi_p

    do k_nb = ia(inp), ia(inp+1)-1
      phi_max = max( phi_max, phi(ja(k_nb)) )
      phi_min = min( phi_min, phi(ja(k_nb)) )
    enddo

    eps2 = (K**3) * ( vol(inp)/voltot ) * range**2

    slopelimit = 1.0_dp

    do k_nb = ia(inp), ia(inp+1)-1

      if ( k_nb == diag(inp) ) cycle
      ijn = ja(k_nb)

      ! Unlimited change of phi from the centre to the face (midpoint)
      d2 = 0.5_dp*( dPhidxi(1,inp)*(xc(ijn)-xc(inp)) &
                  + dPhidxi(2,inp)*(yc(ijn)-yc(inp)) &
                  + dPhidxi(3,inp)*(zc(ijn)-zc(inp)) )

      ! Allowed change towards the local maximum / minimum
      if ( d2 > 0.0_dp ) then
        d1 = phi_max - phi_p
      else
        d1 = phi_min - phi_p
      endif

      ! (d1^2+eps2)*d2 + 2*d2^2*d1 over d2*(d1^2+2*d2^2+d1*d2+eps2); d2 cancels
      philim = ( d1*d1 + eps2 + 2.0_dp*d1*d2 ) / ( d1*d1 + 2.0_dp*d2*d2 + d1*d2 + eps2 )

      slopelimit = min( slopelimit, philim )

    enddo

    dPhidxi(:,inp) = slopelimit*dPhidxi(:,inp)

  enddo

end subroutine


!***********************************************************************
!
subroutine slope_limiter_r345(phi, dPhidxi)
!
!***********************************************************************
!
!     Calculates slope limiter and applies to scalar gradient:
!     R3-R5 slope limiters:
!
!     Based on paper: AIAA2022-1374 by Hiroaki Nishikawa
!
!***********************************************************************
!

  implicit none

  ! Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numTotal) :: dPhidxi


  ! Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  ! real(dp) :: fimax,fimin
  real(dp) :: deltamax,deltamin


  ! fimin = minval(phi(1:numCells))
  ! fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ia(inp) ))
    phi_min = phi(ja( ia(inp) ))

    do k=ia(inp)+1, ia(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_min, phi(ja(k)) )      
    enddo


    deltamax = phi_max - phi(inp)
    deltamin = phi_min - phi(inp)

    slopelimit = 1.0_dp

    do k=ia(inp), ia(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      delta_face=dPhidxi(1,ijp)*(xc(ijn)-xc(ijp))+dPhidxi(2,ijp)*(yc(ijn)-yc(ijp))+dPhidxi(3,ijp)*(zc(ijn)-zc(ijp)) 


     if( abs(delta_face) < 1.e-6 )then
       r = 1.0_dp
     else if( delta_face > 0.0 )then
       r = deltamax/delta_face
     else
       r = deltamin/delta_face
     endif

     ! slopelimit = min( slopelimit , (r**3+4*r)/(r**3+r**2+r+4.0) ) ! R3
        slopelimit = min( slopelimit , (r**4+2*r**3-4*r**2+8*r)/(r**4+r**3+2*r**2-4*r+8.0) ) ! R4
          ! slopelimit = min( slopelimit , (r**5+8*r**3-16*r**2+16*r)/(r**5+r**4+8*r**2-16*r+16) ) ! R5    

    enddo

    dPhidxi(:,inp) = slopelimit*dPhidxi(:,inp)

  enddo

end subroutine



!***********************************************************************
!
subroutine slope_limiter_multidimensional(phi, dPhidxi)
!
!***********************************************************************
!
!  Calculates slope limiter and applies to scalar gradient:
!  Multidimensional slope limiter
!  Ref.: SE Kim, B Makarov, D Caraeni - A Multi-Dimensional Linear 
!  Reconstruction Scheme for Arbitrary Unstructured Meshes, AIAA 2003-3990.
!  The same slope limiter is used in Fluent.
!
!***********************************************************************
!

  implicit none

  ! Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numTotal) :: dPhidxi


  ! Locals
  integer :: inp,ijp,k,iface
  real(dp) :: dPhi,dPhimax,dPhimin
  real(dp) :: gtx,gty,gtz
  real(dp) :: gx,gy,gz,gn
  real(dp) :: xpn,ypn,zpn,dpn
  real(dp) :: nx,ny,nz
  real(dp), dimension(:), allocatable :: phimax,phimin

  allocate(phimax(numCells),phimin(numCells))

  ! Find phi_max and phi_min in  neighbours of Co, including itself.
  do inp = 1,numCells

    ! max and min values over current cell and neighbors
    phimax(inp) = phi(ja( ia(inp) ))
    phimin(inp) = phi(ja( ia(inp) ))

    do k=ia(inp)+1, ia(inp+1)-1
      phimax(inp) = max( phimax(inp), phi(ja(k)) )
      phimin(inp) = min( phimin(inp), phi(ja(k)) )      
    enddo

  enddo


  ! Loop over inner faces:
  do iface=1,numInnerFaces

    do k=1,2 ! Do the same for both owner and neighbour cell

      if (k==1) then
        ijp = owner(iface)
      else
        ijp = neighbour(iface)
      endif

      ! Initialize gradient vector with current unlimited value
      gx = dPhidxi(1,ijp)
      gy = dPhidxi(2,ijp)
      gz = dPhidxi(3,ijp)

      ! Distance vector between cell center and face center
      xpn=xf(iface)-xc(ijp)
      ypn=yf(iface)-yc(ijp)
      zpn=zf(iface)-zc(ijp)

      dpn=sqrt(xpn**2+ypn**2+zpn**2)

      nx = xpn/dpn
      ny = ypn/dpn
      nz = zpn/dpn

      gn = gx*nx+gy*ny+gz*nz

      gtx = gx - gn*nx
      gty = gy - gn*ny
      gtz = gz - gn*nz

      ! Increment from cell center to face centroid f
      dPhi=gx*xpn+gy*ypn+gz*zpn 

      dPhimax = phimax(ijp)-phi(ijp)
      dPhimin = phimin(ijp)-phi(ijp)

      ! Check for overshoots and undershoots and correct accordingly.

      ! Overshoot:
      if ( phimax(ijp) > phi(ijp) .and. dPhi > dPhimax ) then

        gx = gtx + nx*dPhimax
        gy = gty + ny*dPhimax
        gz = gtz + nz*dPhimax

      endif

      ! Undershoot:
      if ( phimin(ijp) < phi(ijp) .and. dPhi < dPhimin ) then

        gx = gtx + nx*dPhimin
        gy = gty + ny*dPhimin
        gz = gtz + nz*dPhimin

      endif

      dPhidxi(1,ijp) = gx
      dPhidxi(2,ijp) = gy
      dPhidxi(3,ijp) = gz

    enddo

  enddo

  deallocate(phimax,phimin)

end subroutine



!***********************************************************************
!
subroutine invert_sym3( a )
!
!***********************************************************************
!
!  Purpose:
!    In-place (pseudo-)inverse of a symmetric positive semi-definite 3x3
!    matrix stored as a(1:6) = (a11,a12,a13,a22,a23,a33).
!
!  Method:
!    1. Diagonal (Jacobi) equilibration M = S*A*S, S = diag(1/sqrt(a_ii)).
!       M has unit diagonal, so the result no longer depends on cell size,
!       units or aspect ratio, and no absolute 'small' number is needed.
!    2. Eigen-decomposition of M by cyclic Jacobi rotations. Accurate even
!       for tiny eigenvalues and needs only a few sweeps for a 3x3 matrix.
!    3. Eigenvalues below rank_tol*max(eigenvalue) are discarded (rank-
!       revealing generalised inverse). A cell whose stencil does not span
!       all three directions (2D meshes, collapsed cells) therefore gets a
!       finite answer instead of 1e+20: the gradient component along an
!       undetermined coordinate axis is exactly zero, and the components
!       that the stencil does determine are still reproduced exactly.
!       (If the degenerate direction is not a coordinate axis the
!       undetermined component is bounded but arbitrary, since this is not
!       the minimum-norm solution.)
!    4. A^+ = S * V * diag(1/lambda) * V^T * S
!
!  It is called once per cell at start-up, so its cost is irrelevant.
!
!***********************************************************************
!
  implicit none

  real(dp), intent(inout) :: a(6)

  real(dp), parameter :: rank_tol = 1.0e-8_dp
  integer,  parameter :: ip(3) = (/ 1, 1, 2 /)
  integer,  parameter :: iq(3) = (/ 2, 3, 3 /)

  real(dp) :: s(3), m(3,3), v(3,3), jr(3,3), lam(3), linv(3)
  real(dp) :: theta, t, c, sn, off, lmax
  integer  :: sweep, k, p, q

  ! 1) Equilibrate. A zero diagonal means no stencil extent in that direction.
  s = 0.0_dp
  if ( a(1) > 0.0_dp ) s(1) = 1.0_dp/sqrt(a(1))
  if ( a(4) > 0.0_dp ) s(2) = 1.0_dp/sqrt(a(4))
  if ( a(6) > 0.0_dp ) s(3) = 1.0_dp/sqrt(a(6))

  m(1,1) = a(1)*s(1)*s(1)
  m(1,2) = a(2)*s(1)*s(2)
  m(1,3) = a(3)*s(1)*s(3)
  m(2,2) = a(4)*s(2)*s(2)
  m(2,3) = a(5)*s(2)*s(3)
  m(3,3) = a(6)*s(3)*s(3)
  m(2,1) = m(1,2)
  m(3,1) = m(1,3)
  m(3,2) = m(2,3)

  ! 2) Cyclic Jacobi eigenvalue iteration
  v = 0.0_dp
  v(1,1) = 1.0_dp
  v(2,2) = 1.0_dp
  v(3,3) = 1.0_dp

  do sweep = 1, 16

    off = abs(m(1,2)) + abs(m(1,3)) + abs(m(2,3))
    if ( off < 1.0e-15_dp ) exit

    do k = 1, 3
      p = ip(k)
      q = iq(k)
      if ( abs(m(p,q)) < 1.0e-300_dp ) cycle

      theta = ( m(q,q) - m(p,p) ) / ( 2.0_dp*m(p,q) )
      if ( abs(theta) > 1.0e150_dp ) then
        t = 0.5_dp/theta
      else
        t = sign(1.0_dp,theta) / ( abs(theta) + sqrt(theta*theta + 1.0_dp) )
      endif
      c  = 1.0_dp/sqrt(t*t + 1.0_dp)
      sn = t*c

      jr = 0.0_dp
      jr(1,1) = 1.0_dp
      jr(2,2) = 1.0_dp
      jr(3,3) = 1.0_dp
      jr(p,p) = c
      jr(q,q) = c
      jr(p,q) = sn
      jr(q,p) = -sn

      m = matmul( transpose(jr), matmul(m,jr) )
      v = matmul( v, jr )
    enddo

  enddo

  ! 3) Pseudo-inverse of the eigenvalues
  do k = 1, 3
    lam(k) = m(k,k)
  enddo
  lmax = maxval(lam)

  linv = 0.0_dp
  do k = 1, 3
    if ( lmax > 0.0_dp .and. lam(k) > rank_tol*lmax ) linv(k) = 1.0_dp/lam(k)
  enddo

  ! 4) Undo equilibration and store the 6 unique entries
  do q = 1, 3
    do p = 1, 3
      m(p,q) = sum( v(p,:)*linv(:)*v(q,:) )
    enddo
  enddo

  a(1) = s(1)*s(1)*m(1,1)
  a(2) = s(1)*s(2)*m(1,2)
  a(3) = s(1)*s(3)*m(1,3)
  a(4) = s(2)*s(2)*m(2,2)
  a(5) = s(2)*s(3)*m(2,3)
  a(6) = s(3)*s(3)*m(3,3)

end subroutine



!***********************************************************************
!
subroutine grad_gauss_corrected(u,dudxi)
!
!***********************************************************************
!
!   Purpose:
!     Gauss gradient with non-orthogonality (skewness) correction.
!
!   Description:
!     The face value is the linear interpolate between cell centres, taken
!     at the point where the line P-N crosses the face plane, corrected to
!     the face centroid with a gradient that is interpolated to that point
!     from the INCOMING gradient:
!
!       phi_f = phi_lin + grad_i . (x_f - x_i)
!
!     On entry dudxi must hold an initial gradient (e.g. from a least-squares
!     routine); on exit it holds the corrected, conservative Gauss gradient.
!     Port of the routine of the same name in fvxGradient.f90.
!
!***********************************************************************
!
  implicit none

  real(dp), dimension(numTotal), intent(in)      :: u
  real(dp), dimension(3,numTotal), intent(inout) :: dudxi

  integer  :: i, ijp, ijn, ijb, iface
  real(dp) :: fxn, fxp, xi, yi, zi, gx, gy, gz, fie, volr
  real(dp), dimension(:,:), allocatable :: g0

  allocate( g0(3,numCells) )
  g0 = dudxi(:,1:numCells)

  dudxi = 0.0_dp

  do i = 1, numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    fxn = facint(i)
    fxp = 1.0_dp - fxn

    ! Point on the line P-N where the linear interpolation is exact
    xi = xc(ijp)*fxp + xc(ijn)*fxn
    yi = yc(ijp)*fxp + yc(ijn)*fxn
    zi = zc(ijp)*fxp + zc(ijn)*fxn

    gx = g0(1,ijp)*fxp + g0(1,ijn)*fxn
    gy = g0(2,ijp)*fxp + g0(2,ijn)*fxn
    gz = g0(3,ijp)*fxp + g0(3,ijn)*fxn

    fie = u(ijp)*fxp + u(ijn)*fxn + gx*(xf(i)-xi) + gy*(yf(i)-yi) + gz*(zf(i)-zi)

    dudxi(1,ijp) = dudxi(1,ijp) + fie*arx(i)
    dudxi(2,ijp) = dudxi(2,ijp) + fie*ary(i)
    dudxi(3,ijp) = dudxi(3,ijp) + fie*arz(i)

    dudxi(1,ijn) = dudxi(1,ijn) - fie*arx(i)
    dudxi(2,ijn) = dudxi(2,ijn) - fie*ary(i)
    dudxi(3,ijn) = dudxi(3,ijn) - fie*arz(i)

  enddo

  do i = 1, numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijb = numCells + i
    call gradbc( arx(iface), ary(iface), arz(iface), u(ijb), dudxi(1,ijp), dudxi(2,ijp), dudxi(3,ijp) )
  enddo

  do ijp = 1, numCells
    volr = 1.0_dp/vol(ijp)
    dudxi(:,ijp) = dudxi(:,ijp)*volr
  enddo

  deallocate( g0 )

end subroutine



!***********************************************************************
!
subroutine create_lsq_face_weights
!
!***********************************************************************
!
!  Purpose:
!    Precompute, for every face, the vectors that turn a difference of the
!    field across the face into a contribution to the weighted least-squares
!    gradient of the adjacent cells:
!
!        grad(phi)_P = sum_f  c_f * ( phi_N(f) - phi_P )
!        c_f = w_f * Ainv_P * d_f,   w_f = 1/|d_f|^2,   d_f = x_N - x_P
!
!    Ainv_P is the inverse LSQ matrix of cell P (see invert_sym3).
!    Boundary faces use the face centre as the 'neighbour' position.
!
!  Why:
!    At run time the gradient is only a sweep over faces: no per-cell
!    solve and no dependence on the number of faces of a cell, so it is
!    valid for any polyhedral cell (unlike the QR variant, which is
!    dimensioned for 6 neighbours). Results equal those of lstsq_dm.
!
!***********************************************************************
!
  implicit none

  integer  :: i, ijp, ijn, iface, ierr
  real(dp) :: w, Dx, Dy, Dz
  logical  :: had_dmat

  allocate( lsq_cp(3,numInnerFaces), lsq_cn(3,numInnerFaces), &
            lsq_cb(3,numBoundaryFaces), stat=ierr )
  if (ierr /= 0) write(*,*) "allocation error: lsq face weights"

  ! Reuse the weighted LSQ matrix builder; it leaves the inverse in Dmat.
  had_dmat = allocated(Dmat)
  if (.not. had_dmat) allocate( Dmat(6,numCells), stat=ierr )
  call create_matrix_lsq_dm

  do i = 1, numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    Dx = xc(ijn)-xc(ijp)
    Dy = yc(ijn)-yc(ijp)
    Dz = zc(ijn)-zc(ijp)
    w  = 1.0_dp/(Dx*Dx + Dy*Dy + Dz*Dz)

    lsq_cp(1,i) = w*( Dmat(1,ijp)*Dx + Dmat(2,ijp)*Dy + Dmat(3,ijp)*Dz )
    lsq_cp(2,i) = w*( Dmat(2,ijp)*Dx + Dmat(4,ijp)*Dy + Dmat(5,ijp)*Dz )
    lsq_cp(3,i) = w*( Dmat(3,ijp)*Dx + Dmat(5,ijp)*Dy + Dmat(6,ijp)*Dz )

    ! Neighbour sees offset -d and difference (phi_P - phi_N): same product
    ! d*(phi_N-phi_P), but with its own inverse matrix.
    lsq_cn(1,i) = w*( Dmat(1,ijn)*Dx + Dmat(2,ijn)*Dy + Dmat(3,ijn)*Dz )
    lsq_cn(2,i) = w*( Dmat(2,ijn)*Dx + Dmat(4,ijn)*Dy + Dmat(5,ijn)*Dz )
    lsq_cn(3,i) = w*( Dmat(3,ijn)*Dx + Dmat(5,ijn)*Dy + Dmat(6,ijn)*Dz )

  enddo

  do i = 1, numBoundaryFaces

    iface = numInnerFaces + i
    ijp = owner(iface)

    Dx = xf(iface)-xc(ijp)
    Dy = yf(iface)-yc(ijp)
    Dz = zf(iface)-zc(ijp)
    w  = 1.0_dp/(Dx*Dx + Dy*Dy + Dz*Dz)

    lsq_cb(1,i) = w*( Dmat(1,ijp)*Dx + Dmat(2,ijp)*Dy + Dmat(3,ijp)*Dz )
    lsq_cb(2,i) = w*( Dmat(2,ijp)*Dx + Dmat(4,ijp)*Dy + Dmat(5,ijp)*Dz )
    lsq_cb(3,i) = w*( Dmat(3,ijp)*Dx + Dmat(5,ijp)*Dy + Dmat(6,ijp)*Dz )

  enddo

  if (.not. had_dmat) deallocate( Dmat )

end subroutine


!***********************************************************************
!
subroutine grad_lsq_fw(Phi,dPhidxi)
!
!***********************************************************************
!
!  Weighted least-squares gradient using precomputed face weights
!  (see create_lsq_face_weights).
!
!***********************************************************************
!
  implicit none

  real(dp), dimension(numTotal), intent(in)      :: Phi
  real(dp), dimension(3,numTotal), intent(inout) :: dPhidxi

  integer  :: i, ijp, ijn, iface
  real(dp) :: delta

  dPhidxi = 0.0_dp

  do i = 1, numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    delta = Phi(ijn)-Phi(ijp)
    dPhidxi(:,ijp) = dPhidxi(:,ijp) + lsq_cp(:,i)*delta
    dPhidxi(:,ijn) = dPhidxi(:,ijn) + lsq_cn(:,i)*delta
  enddo

  do i = 1, numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells + i
    dPhidxi(:,ijp) = dPhidxi(:,ijp) + lsq_cb(:,i)*( Phi(ijn)-Phi(ijp) )
  enddo

end subroutine



subroutine create_matrix_lsq
!
!***********************************************************************
!
!  Purpose:
!  Create coefficient matrix for the function which
!  calculates cell-centered gradients using UNWEIGHTED Least-Squares approach.
!
!  Description:
!  Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.,
!  also in Muzaferija and Gossman JCP paper from 1995.
!
!
!***********************************************************************
!

  implicit none

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: Dx,Dy,Dz
!
!***********************************************************************
!
  ! Coefficient matrix - should be calculated only once 

  ! Initialize Dmat matrix:
  Dmat = 0.0_dp

  ! Inner faces:                                             
  do i=1,numInnerFaces   

    ijp = owner(i)
    ijn = neighbour(i)

    Dx = xc(ijn)-xc(ijp)
    Dy = yc(ijn)-yc(ijp)
    Dz = zc(ijn)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + Dx*Dx
    Dmat(1,ijn) = Dmat(1,ijn) + Dx*Dx 

    Dmat(4,ijp) = Dmat(4,ijp) + Dy*Dy
    Dmat(4,ijn) = Dmat(4,ijn) + Dy*Dy

    Dmat(6,ijp) = Dmat(6,ijp) + Dz*Dz
    Dmat(6,ijn) = Dmat(6,ijn) + Dz*Dz  

    Dmat(2,ijp) = Dmat(2,ijp) + Dx*Dy
    Dmat(2,ijn) = Dmat(2,ijn) + Dx*Dy

    Dmat(3,ijp) = Dmat(3,ijp) + Dx*Dz
    Dmat(3,ijn) = Dmat(3,ijn) + Dx*Dz

    Dmat(5,ijp) = Dmat(5,ijp) + Dy*Dz
    Dmat(5,ijn) = Dmat(5,ijn) + Dy*Dz

  enddo     
 


  ! Boundary faces:
  do i=1,numBoundaryFaces

    iface = numInnerFaces + i
    ijp = owner(iface)

    Dx = xf(iface)-xc(ijp)
    Dy = yf(iface)-yc(ijp)
    Dz = zf(iface)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + Dx*Dx
    Dmat(4,ijp) = Dmat(4,ijp) + Dy*Dy
    Dmat(6,ijp) = Dmat(6,ijp) + Dz*Dz
    Dmat(2,ijp) = Dmat(2,ijp) + Dx*Dy 
    Dmat(3,ijp) = Dmat(3,ijp) + Dx*Dz 
    Dmat(5,ijp) = Dmat(5,ijp) + Dy*Dz 

  end do


  ! Invert the symmetric 3x3 matrix of every cell once and store its 6 unique
  ! entries. See invert_sym3 for how scaling and rank-deficiency are treated.
  do inp=1,numCells
    call invert_sym3( Dmat(1:6,inp) )
  enddo

end subroutine


subroutine grad_lsq(Phi,dPhidxi)
!
!***********************************************************************
!
!  Purpose:
!  Calculates cell-centered gradients using UNWEIGHTED Least-Squares approach.
!
!  Description:
!  Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.,
!  also in Muzaferija and Gossman JCP paper from 1995.
!
!  Arguments:
!
!  Phi - field variable which gradient we look for.
!  DPhiDXi - cell centered gradient - a three component gradient vector.
!
!  Example call:
!  CALL GRADFI_LSQ(U,DUDXI)
!
!***********************************************************************
!

  implicit none

  real(dp),dimension(numTotal), intent(in)   :: Phi
  real(dp),dimension(3,numTotal), intent(inout) :: dPhidxi

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: b1,b2,b3 
  real(dp) :: delta,Dx,Dy,Dz
!
!***********************************************************************
!

  !
  ! *** COMMENT: *** 
  ! We want to save space, so we will reuse dPhidx,dPhidy,dPhidz 
  ! to store rhs vectors. These were, in earlier version, denoted
  ! with 'b1','b2', 'b3' respectively.
  ! This way it is maybe hardrer to read code but introduces great savings,
  ! because each of 'b1','b2', 'b3' should be numCells long!
  !

  ! Initialize rhs vector
  dPhidxi = 0.0_dp

  ! Inner faces:

  do i=1,numInnerFaces    

    ijp = owner(i)
    ijn = neighbour(i)

    delta = Phi(ijn)-Phi(ijp)  

    Dx = ( xc(ijn)-xc(ijp) ) * delta
    Dy = ( yc(ijn)-yc(ijp) ) * delta
    Dz = ( zc(ijn)-zc(ijp) ) * delta

    dPhidxi(1,ijp) = dPhidxi(1,ijp) + Dx
    dPhidxi(2,ijp) = dPhidxi(2,ijp) + Dy
    dPhidxi(3,ijp) = dPhidxi(3,ijp) + Dz 

    dPhidxi(1,ijn) = dPhidxi(1,ijn) + Dx
    dPhidxi(2,ijn) = dPhidxi(2,ijn) + Dy 
    dPhidxi(3,ijn) = dPhidxi(3,ijn) + Dz   

  enddo     

  ! Boundary faces:

  do i=1,numBoundaryFaces

    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells + i

    delta = Phi(ijn)-Phi(ijp)

    Dx = ( xf(iface)-xc(ijp) ) * delta 
    Dy = ( yf(iface)-yc(ijp) ) * delta
    Dz = ( zf(iface)-zc(ijp) ) * delta 

    dPhidxi(1,ijp) = dPhidxi(1,ijp) + Dx
    dPhidxi(2,ijp) = dPhidxi(2,ijp) + Dy
    dPhidxi(3,ijp) = dPhidxi(3,ijp) + Dz

  enddo


  !
  ! Solve the system A*X = B.
  ! 

  do inp=1,numCells 

    b1 = dPhidxi(1,inp)
    b2 = dPhidxi(2,inp)
    b3 = dPhidxi(3,inp)
            
    ! Dmat holds the symmetric inverse as (11,12,13,22,23,33):
    dPhidxi(1,inp) = Dmat(1,inp)*b1 + Dmat(2,inp)*b2 + Dmat(3,inp)*b3
    dPhidxi(2,inp) = Dmat(2,inp)*b1 + Dmat(4,inp)*b2 + Dmat(5,inp)*b3
    dPhidxi(3,inp) = Dmat(3,inp)*b1 + Dmat(5,inp)*b2 + Dmat(6,inp)*b3

  enddo 

  
end subroutine



! Least square gradients via QR decomposition
!***********************************************************************
!
subroutine create_matrix_lsq_qr
!
!***********************************************************************
!
!      Purpose:
!      Forms a system matrix D (that is R^(-1)*Q^T from it's QR factorisation) 
!      which is used for cell-centered gradients using Least-Squares approach.
!
!      Description:
!      Uses QR decomposition of system matrix via Householder or via
!      Gramm-Schmidt.
!      QR decomposition is precomputed and R^(-1)*Q^T is stored in 
!      D array for every cell.
!
!***********************************************************************
!

  use matrix_module

  implicit none

  !
  ! Locals
  !
  integer, parameter :: n=3, m=6  ! m is the number of neighbours, e.g. for hex cell it is 6.
  integer :: i,l,ijp,ijn,inp,iface

  integer, dimension(:), allocatable :: neighbour_index  

  real(dp), dimension(m,n) :: Dtmp
  real(dp), dimension(n,m) :: Dtmpt

#ifndef LAPACK

  ! REAL(dp), DIMENSION(m,n) :: R
  REAL(dp), DIMENSION(m,n) :: Q1
  REAL(dp), DIMENSION(n,n) :: R1
  REAL(dp), DIMENSION(n,m) :: Q1t

#else

  integer :: k
  INTEGER :: INFO
  REAL(dp), DIMENSION(n) :: TAU
  INTEGER, DIMENSION(n) :: WORK
  REAL(dp), DIMENSION(m) :: v1,v2,v3
  REAL(dp), DIMENSION(m,m) :: H1,H2,H3,Ieye
  REAL(dp), DIMENSION(n,n) :: R
  REAL(dp), DIMENSION(m,m) :: Q

#endif

 

! Coefficient matrix - should be calculated only once 
!**************************************************************************************************
    Dtmp = 0.0_dp

    allocate(neighbour_index(numCells))
    neighbour_index = 0

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

      neighbour_index(ijp) = neighbour_index(ijp) + 1 
      l = neighbour_index(ijp)
      D(1,l,ijp) = xc(ijn)-xc(ijp)
      D(2,l,ijp) = yc(ijn)-yc(ijp)
      D(3,l,ijp) = zc(ijn)-zc(ijp)
  
      neighbour_index(ijn) = neighbour_index(ijn) + 1 
      l = neighbour_index(ijn)   
      D(1,l,ijn) = xc(ijp)-xc(ijn)
      D(2,l,ijn) = yc(ijp)-yc(ijn)
      D(3,l,ijn) = zc(ijp)-zc(ijn)
                                                            
  enddo     

  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ! ijb = numCells + i
      neighbour_index(ijp) = neighbour_index(ijp) + 1
      l = neighbour_index(ijp)
      D(1,l,ijp) = xf(iface)-xc(ijp)
      D(2,l,ijp) = yf(iface)-yc(ijp)
      D(3,l,ijp) = zf(iface)-zc(ijp)
  end do

  ! Form system matrix using QR decomposition:

  do inp=1,numCells

  l = neighbour_index(inp)

  Dtmpt = D(:,:,inp)
  Dtmp = transpose(Dtmpt)


#ifndef LAPACK

  ! 1) Decompose A=QR using Householder
  ! call householder_qr(Dtmp, m, n, Q, R)
  ! 2) Decompose A=QR using Gram-Schmidt
  call mgs_qr(Dtmp, m, n, Q1, R1)

  Q1t = transpose(Q1)
  ! Q1t = Q(1:n,1:m)      ! NOTE: A=Q1R1 is so-called 'thin QR factorization' - see Golub & Van Loan
                        ! Here Q1 is actually Q1^T a transpose of Q1(thin Q - Q with m-n column stripped off)
  ! R1 = R(1:n,1:n)       ! our Q1 is thin transpose of original Q.
  R1 = inv(R1)          ! inv is a function in matrix_module, now works only for 3x3 matrices.
  Q1t  = matmul(R1,Q1t) ! this is actually R^(-1)*Q^T - a matrix of size n x m.
  D(:,:,INP) = Q1t      ! Store it for later.


#else

  ! LAPACK routine DGEQRF
  CALL DGEQRF( l, N, Dtmp, M, TAU, WORK, N, INFO )  

  ! Upper triangular matrix R
  R(1:n,1:n)=Dtmp(1:n,1:n)

  ! Create reflectors
  !H(i) = I - TAU * v * v'
  Ieye=eye(l)
  !v(1:i-1) = 0. and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i)
  v1(1) = 1.; v1(2:l)=Dtmp(2:l,1)
  H1 = rank_one_update(Ieye,l,l,v1,v1,-TAU(1))
  v2(1) = 0.; v2(2) = 1.; v2(3:l)=Dtmp(3:l,2)
  H2 = rank_one_update(Ieye,l,l,v2,v2,-TAU(2))
  v3(1:2) = 0.; v3(3) = 1.; v3(4:l)=Dtmp(4:l,3)
  H3 = rank_one_update(Ieye,l,l,v3,v3,-TAU(3))
  ! The matrix Q is represented as a product of elementary reflectors H1, H2, ..., Hn
  Q=matmul(H1,H2)
  Q=matmul(Q,H3)

  ! Form R_1^(-1)*Q_1^T explicitly:
  do k=1,neighbour_index(inp)
    d(1,k,inp) = q(k,1)/r(1,1) - (r(1,2)*q(k,2))/(r(1,1)*r(2,2)) + (q(k,3)*(r(1,2)*r(2,3) - r(1,3)*r(2,2)))/(r(1,1)*r(2,2)*r(3,3))
    d(2,k,inp) = q(k,2)/r(2,2) - (r(2,3)*q(k,3))/(r(2,2)*r(3,3))
    d(3,k,inp) = q(k,3)/r(3,3)
  enddo


#endif  


  enddo

  deallocate(neighbour_index)

end subroutine



! Least square gradients via QR decomposition
!***********************************************************************
!
subroutine grad_lsq_qr(Phi,dPhidxi)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using Least-Squares approach.
!
!      Description:
!      Uses QR decomposition of system matrix via Householder or via
!      Gramm-Schmidt.
!      QR decomposition is precomputed and R^(-1)*Q^T is stored in 
!      D array for every cell.
!
!      Arguments:
!
!      Phi - field variable
!      dPhidxi - reconstructed cell centered gradient vector.
!
!      Example call:
!      CALL GRAD_LSTSQ_QR(U,dUdxi)
!
!***********************************************************************
!

  use matrix_module

  implicit none

  integer, parameter :: n=3, m=6  ! m is the number of neighbours, e.g. for structured 3D mesh it's 6

  real(dp), dimension(numTotal), intent(in)   :: Phi
  real(dp), dimension(3,numCells), intent(inout) :: dPhidxi

  !
  !    Locals
  !
  integer ::  i,l,ijp,ijn,inp,iface

  integer, dimension(:), allocatable :: neighbour_index  
  real(dp), dimension(:,:), allocatable  :: b ! rhs vector


  allocate( neighbour_index(numCells), b(m,numCells) )
  neighbour_index = 0
  b = 0.0_dp

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

        neighbour_index(ijp) = neighbour_index(ijp) + 1
        l = neighbour_index(ijp)
        b(l,ijp) = Phi(ijn)-Phi(ijp)

        neighbour_index(ijn) = neighbour_index(ijn) + 1        
        l = neighbour_index(ijn)   
        b(l,ijn) = Phi(ijp)-Phi(ijn)
                                                                                                                   
  enddo     

  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells + i

        neighbour_index(ijp) = neighbour_index(ijp) + 1
        l = neighbour_index(ijp)
        b(l,ijp) = Phi(ijn)-Phi(ijp)

  enddo

! Solve overdetermined system in least-square sense

  ! Cell loop
  do inp=1,numCells

    l = neighbour_index(inp)     

    !  ...using precomputed QR factorization and storing R^(-1)*Q^T in D
    dPhidxi(1,INP) = sum(D(1,1:l,inp)*b(1:l,inp))
    dPhidxi(2,INP) = sum(D(2,1:l,inp)*b(1:l,inp))
    dPhidxi(3,INP) = sum(D(3,1:l,inp)*b(1:l,inp))

  enddo


end subroutine



! Weighted least square gradients
subroutine create_matrix_lsq_dm
!
!***********************************************************************
!
!      Purpose:
!      Creates and stores system matrix based on geometry data for 
!      cell-centered gradient reconstruction using WEIGHTED Least-Squares approach.
!
!      Description:
!      Approach taken from a paper:
!      Dimitry Mavriplis, "Revisiting the Least Squares procedure for Gradient 
!      Reconstruction on Unstructured Meshes." NASA/CR-2003-212683.
!      
!      Weights based on inverse cell-centers distance is added to improve conditioning 
!      of the system matrix for skewed meshes.
!      Reduces storage requirements compared to QR subroutine.
!      System matrix is symmetric and can be solved efficiently using 
!      Cholesky decomposition or by matrix inversion. 
!
!***********************************************************************
!

  implicit none

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface
  real(dp) :: w
  real(dp) :: Dx,Dy,Dz
!
!***********************************************************************
!
  ! Coefficient matrix - should be calculated only once 

  ! Initialize dmat matrix:
  Dmat = 0.0d0

  ! Inner faces:     

  do i=1,numInnerFaces                                                       
    ijp = owner(i)

    !***------------------------------------
    ijn = neighbour(i)

    w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

    Dx = xc(ijn)-xc(ijp)
    Dy = yc(ijn)-yc(ijp)
    Dz = zc(ijn)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
    Dmat(1,ijn) = Dmat(1,ijn) + w*Dx*Dx 

    Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
    Dmat(4,ijn) = Dmat(4,ijn) + w*Dy*Dy

    Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
    Dmat(6,ijn) = Dmat(6,ijn) + w*Dz*Dz  

    Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy
    Dmat(2,ijn) = Dmat(2,ijn) + w*Dx*Dy

    Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz
    Dmat(3,ijn) = Dmat(3,ijn) + w*Dx*Dz

    Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz
    Dmat(5,ijn) = Dmat(5,ijn) + w*Dy*Dz     
    !***------------------------------------

  ! !
  ! ! **Extended interpolation molecule: neighbours of neighbours**
  ! !
  !   nb_loop: do k = ia( neighbour(i) ), ia( neighbour(i)+1 )-1

  !     ijn = ja(k)

  !     if (ijn == ijp) cycle nb_loop ! nemoj ownera dirati

  !     ! Paznja kada je ja(k) = diag( neighbour(i) ) ijn ce uzeti samog cell_neighboura, 
  !     ! tako da ovaj loop moze uraditi isto sto i gornji, pa je u toj verziji
  !     ! racunanja gradijenata ovaj loop gore nepotreban. Eto simplifikacije koda...

  !     w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

  !     Dx = xc(ijn)-xc(ijp)
  !     Dy = yc(ijn)-yc(ijp)
  !     Dz = zc(ijn)-zc(ijp)

  !     Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
  !     Dmat(1,ijn) = Dmat(1,ijn) + w*Dx*Dx 

  !     Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
  !     Dmat(4,ijn) = Dmat(4,ijn) + w*Dy*Dy

  !     Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
  !     Dmat(6,ijn) = Dmat(6,ijn) + w*Dz*Dz  

  !     Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy
  !     Dmat(2,ijn) = Dmat(2,ijn) + w*Dx*Dy

  !     Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz
  !     Dmat(3,ijn) = Dmat(3,ijn) + w*Dx*Dz
 
  !     Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz
  !     Dmat(5,ijn) = Dmat(5,ijn) + w*Dy*Dz 

  !   enddo nb_loop
                                                           
  enddo     

 
  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)

    w = 1.0_dp/((xf(iface)-xc(ijp))**2+(yf(iface)-yc(ijp))**2+(zf(iface)-zc(ijp))**2)

    Dx = xf(iface)-xc(ijp)
    Dy = yf(iface)-yc(ijp)
    Dz = zf(iface)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
    Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
    Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
    Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy 
    Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz 
    Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz 
  end do


  ! Invert the symmetric 3x3 matrix of every cell once and store its 6 unique
  ! entries. See invert_sym3 for how scaling and rank-deficiency are treated.
  do inp=1,numCells
    call invert_sym3( Dmat(1:6,inp) )
  enddo

  
end subroutine




! Weighted least square gradients
subroutine grad_lsq_dm(Phi,dPhidxi)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using WEIGHTED Least-Squares approach.
!
!      Description:
!      Approach taken from a paper:
!      Dimitry Mavriplis, "Revisiting the Least Squares procedure for Gradient 
!      Reconstruction on Unstructured Meshes." NASA/CR-2003-212683.
!      
!      Weights based on inverse cell-centers distance is added to improve 
!      conditioning of the system matrix for skewed meshes.
!      Reduces storage requirements compared to QR subroutine.
!      System matrix is symmetric and can be solved efficiently using 
!      Cholesky decomposition or by matrix inversion. 
!
!      Arguments:
!
!      Phi - field variable which gradient we look for.
!      dPhidxi- cell centered gradient - a three component gradient vector.
!
!      Uses module variable Dmat - LSQ matrix with geometry data.
!
!      Example call:
!      CALL GRADFI_LSQ_DM(U,DUDXI)
!
!***********************************************************************
!

  implicit none

  real(dp),dimension(numTotal), intent(in)   :: Phi
  real(dp),dimension(3,numTotal), intent(inout) :: dPhidxi

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: w
  real(dp) :: Dx,Dy,Dz
  real(dp) :: b1,b2,b3

!
!***********************************************************************
!

  !
  ! *** COMMENT: *** 
  ! We want to save space, so we will reuse dPhidx,dPhidy,dPhidz 
  ! to store rhs vectors. These were, in earlier version, denoted
  ! with 'b1','b2', 'b3' respectively.
  ! This way it is maybe hardrer to read code but introduces great savings,
  ! because each of 'b1','b2', 'b3' should be numCells long!
  !

  ! Initialize rhs vector
  dPhidxi = 0.0_dp

  ! Inner faces:

  do i=1,numInnerFaces                                                       
    ijp = owner(i)

    !***------------------------------------
    ijn = neighbour(i)

    w = ( Phi(ijn)-Phi(ijp) ) &
       /( (xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2 )

    Dx = w * ( xc(ijn)-xc(ijp) ) 
    Dy = w * ( yc(ijn)-yc(ijp) )
    Dz = w * ( zc(ijn)-zc(ijp) )

    dPhidxi(1,ijp) = dPhidxi(1,ijp) + Dx 
    dPhidxi(2,ijp) = dPhidxi(2,ijp) + Dy
    dPhidxi(3,ijp) = dPhidxi(3,ijp) + Dz 

    dPhidxi(1,ijn) = dPhidxi(1,ijn) + Dx
    dPhidxi(2,ijn) = dPhidxi(2,ijn) + Dy 
    dPhidxi(3,ijn) = dPhidxi(3,ijn) + Dz      
    !***------------------------------------
    
  ! !
  ! ! **Extended interpolation molecule: neighbours of neighbours**
  ! !
  !   nb_loop2: do k = ia( neighbour(i) ), ia( neighbour(i)+1 )-1

  !     ijn = ja(k)

  !     if (ijn == ijp) cycle nb_loop2 ! nemoj ownera dirati

  !     ! Paznja kada je ja(k) = diag( neighbour(i) ) ijn ce uzeti samo cell_neighboura, 
  !     ! tako da ovaj loop moze uraditi isto sto i gornji, pa je u toj verziji
  !     ! racunanja gradijenata ovaj loop gore nepotreban. Eto simplifikacije koda...

  !     w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

  !     Dx = w * ( xc(ijn)-xc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
  !     Dy = w * ( yc(ijn)-yc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
  !     Dz = w * ( zc(ijn)-zc(ijp) ) * ( Phi(ijn)-Phi(ijp) )

    ! dPhidxi(1,ijp) = dPhidxi(1,ijp) + Dx 
    ! dPhidxi(2,ijp) = dPhidxi(2,ijp) + Dy
    ! dPhidxi(3,ijp) = dPhidxi(3,ijp) + Dz 

    ! dPhidxi(1,ijn) = dPhidxi(1,ijn) + Dx
    ! dPhidxi(2,ijn) = dPhidxi(2,ijn) + Dy 
    ! dPhidxi(3,ijn) = dPhidxi(3,ijn) + Dz

  !   enddo nb_loop2
                                                            
  enddo     


  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells+i

    w = ( Phi(ijn)-Phi(ijp) ) &
       /( (xf(iface)-xc(ijp))**2+(yf(iface)-yc(ijp))**2+(zf(iface)-zc(ijp))**2 )

    Dx = w*(xf(iface)-xc(ijp)) 
    Dy = w*(yf(iface)-yc(ijp))
    Dz = w*(zf(iface)-zc(ijp)) 

    dPhidxi(1,ijp) = dPhidxi(1,ijp) + Dx
    dPhidxi(2,ijp) = dPhidxi(2,ijp) + Dy
    dPhidxi(3,ijp) = dPhidxi(3,ijp) + Dz

  end do

  ! Calculate gradient

  do inp=1,numCells 

    b1 = dPhidxi(1,inp)
    b2 = dPhidxi(2,inp)
    b3 = dPhidxi(3,inp)
            
    ! Dmat holds the symmetric inverse as (11,12,13,22,23,33):
    dPhidxi(1,inp) = Dmat(1,inp)*b1 + Dmat(2,inp)*b2 + Dmat(3,inp)*b3
    dPhidxi(2,inp) = Dmat(2,inp)*b1 + Dmat(4,inp)*b2 + Dmat(5,inp)*b3
    dPhidxi(3,inp) = Dmat(3,inp)*b1 + Dmat(5,inp)*b2 + Dmat(6,inp)*b3

  enddo 

  
end subroutine



! Node based Gauss gradients
subroutine grad_node_based_gauss(u,dudxi)
!
!  Purpose: 
!
!    This subroutine implements so called Green-Gauss Node Based Gradients.
!
!  Discussion:
!
!    First interpolated values at nodes are found using pseudolaplacian 
!    technique (from node_interpolation module). Then face centroid values
!    are found using node (vertex) averaging. After that we proceed with 
!    conventional Gauss gradient approach.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    19th September 2022
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    S.-E.Kim,B.Makarov, and D.Caraeni, A Multi-Dimensional Linear Reconstruction Scheme for Arbitrary Unstructured Grids, AIAA 2003-3990, 2003.
!
!  Parameters:
!
!    Input, real, u( numTotal ), a scalar field for which we calculate gradients.
!    Input/Output, real, dudxi( 3,numTotal ), gradient vector.
!
  use node_interpolation, only: interpolate_to_nodes

  implicit none
!
! Parameters
!
  real(dp), dimension(numTotal), intent(in) :: u
  real(dp), dimension(3,numTotal), intent(inout) :: dudxi
!
! Local variables
!
  integer :: iface,k,ino,ijp,ijn
  real(dp) :: fie, dfxe, dfye, dfze, volr
  real(dp), dimension(:), allocatable :: uN


  allocate( uN(numNodes) )

  ! Interpolate to nodes using interpolation weights generated at init stage.
  call interpolate_to_nodes(u, uN)
   
  ! Initialize new gradient vector field
  dudxi = 0.0_dp

  ! Loop over inner faces
  do iface=1,numFaces

    ! Initialize face centroid interpolated value
    fie = 0.0_dp

    ! Loop over face nodes
    do k=1,nnodes(iface)

      ! Take global index of the node
      ino = node(k,iface)

      ! Accumulate to the sum the field value at the specific node
      fie = fie + uN(ino)

    enddo

    ! Make arithmetic average based on number of nodes at face
    fie = fie / dble( nnodes(iface) )

    ! (interpolated mid-face value)x(area)
    dfxe = fie*arx(iface)
    dfye = fie*ary(iface)
    dfze = fie*arz(iface)

    ! Get indices of the owner and neighbor cells, and
    ! accumulate contribution at owner and neighbour cell center.

    ijp = owner(iface)
    dudxi(1,ijp) = dudxi(1,ijp) + dfxe
    dudxi(2,ijp) = dudxi(2,ijp) + dfye
    dudxi(3,ijp) = dudxi(3,ijp) + dfze
  
    if (iface .le. numInnerFaces ) then

    ijn = neighbour(iface)   
    dudxi(1,ijn) = dudxi(1,ijn) - dfxe
    dudxi(2,ijn) = dudxi(2,ijn) - dfye
    dudxi(3,ijn) = dudxi(3,ijn) - dfze

    endif

  enddo


  ! Calculate gradient components at cv-centers
  do ijp=1,numCells
    volr = 1./vol(ijp)
    dudxi(:,ijp) = dudxi(:,ijp)*volr
  enddo

  deallocate( uN )

end subroutine



! Gauss gradients
subroutine grad_gauss(u,dudxi)
!
!***********************************************************************
!
!   Purpose:
!     Calculates cell centered gradient using Gauss theorem.
!   Parameters:
!     u - the field for which we compute the gradient
!     dudxi - gradient vector field
!
!     gauss gradient rule:
!     ------->                                 ->
!     grad(u) = 1/vol * sum_{i=1}^{i=nf} (u)_f*sf
!     where:
!     grad(u) - cell centered gradient vector
!     (u)_f   - face interpolated value of scalar u
!     vol     - cell volume
!     sf      - cell face area vector
!     nf      - number of faces in a cell
!
!***********************************************************************
!

  implicit none

  ! Arguments
  real(dp), dimension(numTotal), intent(in) :: u
  real(dp), dimension(3,numTotal), intent(inout) :: dudxi

  ! Local
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: volr,dfxe,dfye,dfze,fie

  
  ! Initialize new gradient
  dudxi = 0.0_dp

  ! Calculate terms integrated over surfaces

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)

    ! Linear interpolation to face
    fie = u(ijp) + (u(ijn)-u(ijp))*facint(i)

    ! (interpolated mid-face value)x(area)
    dfxe = fie*arx(i)
    dfye = fie*ary(i)
    dfze = fie*arz(i)


    ! Accumulate contribution at cell center and neighbour
    dudxi(1,ijp) = dudxi(1,ijp) + dfxe
    dudxi(2,ijp) = dudxi(2,ijp) + dfye
    dudxi(3,ijp) = dudxi(3,ijp) + dfze
     
    dudxi(1,ijn) = dudxi(1,ijn) - dfxe
    dudxi(2,ijn) = dudxi(2,ijn) - dfye
    dudxi(3,ijn) = dudxi(3,ijn) - dfze

  enddo

  ! Contribution from boundaries

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijb = numCells + i
    
    call gradbc( arx(iface), ary(iface), arz(iface), u(ijb), dudxi(1,ijp), dudxi(2,ijp), dudxi(3,ijp) )

  enddo


  ! Calculate gradient components at cv-centers
  do ijp=1,numCells

    volr = 1.0_dp/vol(ijp)

    dudxi(:,ijp) = dudxi(:,ijp)*volr

  enddo


end subroutine


subroutine gradbc(sx,sy,sz,fi,dfx,dfy,dfz)
!***********************************************************************
!
! This routine calculates the contribution of a 
! boundary cell face to the gradient at CV-center.
!
!***********************************************************************
  use types

  implicit none

  real(dp), intent(in) :: sx,sy,sz
  real(dp), intent(in) :: fi
  real(dp), intent(inout)  :: dfx,dfy,dfz

  dfx = dfx + fi*sx
  dfy = dfy + fi*sy
  dfz = dfz + fi*sz

end subroutine


!******************************************************************************
!
subroutine sngrad(i, ijp, ijn, arx, ary, arz, lambda, Phi, dPhidxi, &
                  dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii)
!
!******************************************************************************
!
!  Surface normal gradient with non-orthogonal correction done in two
!  possible ways - either by skewness correction of intersection point
!  offset.
!
!  Check out reference paper: 
!    Mirkov, Rasuo, Kenjeres, JCP, Vol. 287, 2015.
!
!******************************************************************************
!
  implicit none

  integer, intent(in) :: i, ijp, ijn
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numTotal), intent(in) :: Phi
  real(dp), dimension(3,numTotal), intent(in) :: dPhidxi
  real(dp), intent(out) :: dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii
!
! Locals
!
  real(dp) :: fxp,fxn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: vole

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! dpn * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Interpolate gradients defined at CV centers to faces
  dfixi = dPhidxi(1,ijp)*fxp+dPhidxi(1,ijn)*fxn
  dfiyi = dPhidxi(2,ijp)*fxp+dPhidxi(2,ijn)*fxn
  dfizi = dPhidxi(3,ijp)*fxp+dPhidxi(3,ijn)*fxn

  ! du/dx_i interpolated at cell face:
  dfixii = dfixi + arx/vole*( Phi(ijn)-Phi(ijp)-dfixi*xpn-dfiyi*ypn-dfizi*zpn ) 
  dfiyii = dfiyi + ary/vole*( Phi(ijn)-Phi(ijp)-dfixi*xpn-dfiyi*ypn-dfizi*zpn ) 
  dfizii = dfizi + arz/vole*( Phi(ijn)-Phi(ijp)-dfixi*xpn-dfiyi*ypn-dfizi*zpn )

  ! Create non-orthogonal correction to gradient only
  dfixi = dfixi*(arx-Df(i)*xpn)
  dfiyi = dfiyi*(ary-Df(i)*ypn)
  dfizi = dfizi*(arz-Df(i)*zpn)

end subroutine


!******************************************************************************
!
subroutine sngrad2(ijp, ijn, arx, ary, arz, lambda, Phi, dPhidxi, &
                  dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii)
!
!******************************************************************************
!
!  Surface normal gradient with non-orthogonal correction done in two
!  possible ways - either by skewness correction of intersection point
!  offset.
!
!  Check out reference paper: 
!    Mirkov, Rasuo, Kenjeres, JCP, Vol. 287, 2015.
!
!******************************************************************************
!
  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numTotal), intent(in) :: Phi
  real(dp), dimension(3,numTotal), intent(in) :: dPhidxi
  real(dp), intent(out) :: dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii
!
! Locals
!
  real(dp) :: fxp,fxn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: dpn,are,costheta
  real(dp) :: vole,costn
  real(dp) :: d1x,d1y,d1z
  real(dp) :: d2x,d2y,d2z


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are = sqrt(arx**2+ary**2+arz**2)

  ! Angle between unit vectors of the face normal and i_xi along dpn direction - we need cosine
  costheta = (arx*xpn+ary*ypn+arz*zpn)/(are*dpn)

  ! Relaxation factor for higher-order cell face gradient
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! if(nrelax == 1) then
  ! ! Minimal correction: nrelax = +1 :
     costn = costheta
  ! elseif(nrelax == 0) then
  ! ! Orthogonal correction: nrelax =  0 : 
    ! costn = 1.0_dp
  ! elseif(nrelax == -1) then
  ! ! Over-relaxed approach: nrelax = -1 :
    ! costn = 1.0_dp/costheta  
  ! endif

  ! dpn * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Interpolate gradients defined at CV centers to faces
  dfixi = dPhidxi(1,ijp)*fxp+dPhidxi(1,ijn)*fxn
  dfiyi = dPhidxi(2,ijp)*fxp+dPhidxi(2,ijn)*fxn
  dfizi = dPhidxi(3,ijp)*fxp+dPhidxi(3,ijn)*fxn


  ! Overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn

  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn

  !.....du/dx_i interpolated at cell face:
  dfixii = dfixi*d1x + arx/vole*( Phi(ijn)-Phi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfiyii = dfiyi*d1y + ary/vole*( Phi(ijn)-Phi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfizii = dfizi*d1z + arz/vole*( Phi(ijn)-Phi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z )


end subroutine



!******************************************************************************
!
subroutine sngrad3(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, Phi, dPhidxi, &
                  dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii)
!
!******************************************************************************
!
!  Surface normal gradient with non-orthogonal correction done in two
!  possible ways - either by skewness correction of intersection point
!  offset.
!
!  Check out reference paper: 
!    Mirkov, Rasuo, Kenjeres, JCP, Vol. 287, 2015.
!
!******************************************************************************
!
  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf, yf, zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numTotal), intent(in) :: Phi
  real(dp), dimension(3,numTotal), intent(in) :: dPhidxi
  real(dp), intent(out) :: dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii
!
! Locals
!
  real(dp) :: fxp,fxn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz
  ! real(dp) :: dpn
  real(dp) :: are
  ! real(dp) :: costheta
  real(dp) :: costn
  real(dp) :: vole
  real(dp) :: d1x,d1y,d1z
  real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  real(dp) :: gfxdnn1,gfxdpp1


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! ! Distance from P to neighbor N
  ! dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! ! Angle between unit vectors of the face normal and i_xi along dpn direction - we need cosine
  ! costheta=arx*xpn+ary*ypn+arz*zpn/(are*dpn)

  ! Relaxation factor for higher-order cell face gradient
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! if(nrelax == 1) then
  !   ! Minimal correction: nrelax = +1 :
  !   costn = costheta
  ! elseif(nrelax == 0) then
    ! Orthogonal correction: nrelax =  0 : 
    costn = 1.0_dp
  ! elseif(nrelax == -1) then
  !   ! Over-relaxed approach: nrelax = -1 :
  !   costn = 1.0_dp/costheta  
  ! endif

  ! dpn * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Interpolate gradients defined at CV centers to faces
  dfixi = dPhidxi(1,ijp)*fxp+dPhidxi(1,ijn)*fxn
  dfiyi = dPhidxi(2,ijp)*fxp+dPhidxi(2,ijn)*fxn
  dfizi = dPhidxi(3,ijp)*fxp+dPhidxi(3,ijn)*fxn



  ! |-- Intersection point offset and skewness correction -->

    nxx = arx/are 
    nyy = ary/are 
    nzz = arz/are

    ! Find points P' and Pj'
    xpp=xf-(xf-xc(ijp))*nxx
    ypp=yf-(yf-yc(ijp))*nyy 
    zpp=zf-(zf-zc(ijp))*nzz

    xep=xf-(xf-xc(ijn))*nxx 
    yep=yf-(yf-yc(ijn))*nyy 
    zep=zf-(zf-zc(ijn))*nzz     

    xpnp = xep-xpp 
    ypnp = yep-ypp 
    zpnp = zep-zpp

    volep = arx*xpnp+ary*ypnp+arz*zpnp

   ! Overrelaxed correction vector d2, where S=dpn+d2
    d1x = costn
    d1y = costn
    d1z = costn
    
    xpnp = xpnp*costn
    ypnp = ypnp*costn
    zpnp = zpnp*costn

    ! The cell face interpolated gradient (d phi / dx_i)_j:
    ! Nonorthogonal corrections:      ___
    ! gfxdnn1 =>> dot_product(dPhidxi,dNN')
    ! And:                            ___
    ! gfxdnn1 =>> dot_product(dPhidxi,dPP')
    gfxdnn1 = dPhidxi(1,ijn)*(xep-xc(ijn))+dPhidxi(2,ijn)*(yep-yc(ijn))+dPhidxi(3,ijn)*(zep-zc(ijn))
    gfxdpp1 = dPhidxi(1,ijp)*(xpp-xc(ijp))+dPhidxi(2,ijp)*(ypp-yc(ijp))+dPhidxi(3,ijp)*(zpp-zc(ijp))

    dfixii = dfixi*d1x + arx/volep*( Phi(ijn)+gfxdnn1-Phi(ijp)-gfxdpp1-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
    dfiyii = dfiyi*d1y + ary/volep*( Phi(ijn)+gfxdnn1-Phi(ijp)-gfxdpp1-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
    dfizii = dfizi*d1z + arz/volep*( Phi(ijn)+gfxdnn1-Phi(ijp)-gfxdpp1-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 

end subroutine


end module gradients