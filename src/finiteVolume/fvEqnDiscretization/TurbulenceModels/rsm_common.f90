module rsm_common
!********************************************************************************************************************************132
!
!  Purpose:
!
!    Common core of the Reynolds-stress transport models 'rsm_full' (high-Re, wall functions) and 'ebrsm'
!    (elliptic blending, low-Re), on top of 'turb_transport'.
!
!  Discussion:
!
!    Six transport equations for the stresses R_ij = <u_i u_j> (stored as rho-free kinematic values in rs(:,1:6),
!    order xx,xy,xz,yy,yz,zz - the same as 'aij'), plus the dissipation rate epsilon in the global array 'ed':
!
!      DR_ij/Dt = P_ij + Phi_ij - eps_ij + div[ (nu + Cs (2/3) k T) grad R_ij ]
!      Deps/Dt  = ( Ceps1 P - Ceps2 eps ) / T + div[ (nu + Ceps (2/3) k T) grad eps ]
!
!      P_ij   = - R_ik dU_j/dx_k - R_jk dU_i/dx_k,        P = P_kk/2
!      Phi_ij = (1-alpha^2) Phi^w_ij + alpha^2 Phi^h_ij   (alpha = 1 everywhere in 'rsm_full')
!      eps_ij = (1-alpha^2) (R_ij/k) eps + alpha^2 (2/3) eps delta_ij
!
!    Phi^h is the SSG pressure-strain model [1] written with a_ij = R_ij/k - 2/3 delta_ij (constants of Manceau [2], these are
!    the SSG constants C1=3.4, C1*=1.8, C2=4.2, C4=1.25, C5=0.4 divided for the change b_ij = a_ij/2).
!    Phi^w = -5 (eps/k) ( R_ik n_j n_k + R_jk n_i n_k - 1/2 R_kl n_k n_l (n_i n_j + delta_ij) ), n = grad alpha/|grad alpha|.
!
!    Simplification with respect to the literature: the Daly-Harlow (generalised gradient) diffusion tensor Cs T R_kl is
!    contracted with the direction of the gradient of the diffused quantity (function 'ggdh'), giving the scalar
!    diffusivity the engine needs. Exact for 1D gradients (channel, boundary layer), approximate otherwise. Isotropic
!    'nu + Cs (2/3) k T' over-diffuses in the wall-normal direction (R_nn < 2/3 k) and gave Cf 25% low in a channel.
!
!    Coupling to the momentum equations - the whole difficulty of segregated RSM solvers: an explicit divergence of the
!    Reynolds stress makes the coupling weak and leads to checkerboarding/divergence. We use the standard remedy:
!      T_ij = -rho R_ij = 2 mu_a S_ij - 2/3 rho k delta_ij - rho k a^ex_ij,
!      mu_a = rho C_a k T (artificial viscosity, in 'vis' - implicit),   a^ex_ij = R_ij/k - 2/3 delta_ij + 2 C_a T S_ij
!    where a^ex is the excess anisotropy used by 'asm_stress_terms' (explicit, activated by 'lasm').
!
!  Reference:
!    [1] C.G. Speziale, S. Sarkar, T.B. Gatski, J. Fluid Mech. 227, 1991, pp. 245-272.
!    [2] R. Manceau, Int. J. Heat Fluid Flow 51, 2015, pp. 395-400.  R. Manceau, K. Hanjalic, Phys. Fluids 14, 2002.
!    [3] B.E. Launder, G.J. Reece, W. Rodi, J. Fluid Mech. 68, 1975, pp. 537-566.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
  use types
  use parameters
  use geometry
  use variables
  use gradients
  use TurbModelData, only: TurbModel
  use turb_transport

  implicit none

  ! SSG pressure-strain constants for a_ij = R_ij/k - 2/3 delta_ij
  real(dp), parameter :: C1 = 1.7_dp, C1S = 0.9_dp, C2 = 1.05_dp, C3 = 0.8_dp, C3S = 0.65_dp, C4 = 0.625_dp, C5 = 0.2_dp
  real(dp), parameter :: CEPS1 = 1.44_dp, CEPS2 = 1.83_dp, A1EB = 0.03_dp
  real(dp), parameter :: CS = 0.22_dp, CEPS = 0.18_dp, CA = 0.1_dp         ! diffusion, artificial viscosity coefficients
  real(dp), parameter :: MUA_MAX = 1.0e3_dp                                   ! bound of mu_a/mu
  real(dp), parameter :: CMU = 0.09_dp, CMU25 = sqrt(sqrt(CMU)), CMU75 = CMU25**3

  integer, parameter :: ixx=1, ixy=2, ixz=3, iyy=4, iyz=5, izz=6
  integer, parameter, dimension(6) :: ii = [1,1,1,2,2,3], jj = [1,2,3,2,3,3]
  integer, parameter :: iEps = 7                        ! slot of epsilon in TurbModel%Scalar (stresses use 1..6)

  real(dp), dimension(:,:), allocatable :: rs, rso, rsoo       ! stresses (numTotal,6) and old time levels
  real(dp), dimension(:,:,:), allocatable :: drs              ! their gradients (3,numTotal,6)
  real(dp), dimension(:,:), allocatable :: dAlphadxi ! gradient of the blending function (3,numTotal)
  real(dp), dimension(:,:), allocatable :: fixvalR            ! wall-function values of the stresses in wall cells
  real(dp), dimension(:,:), allocatable :: sx                 ! explicit sources (numCells,6), kinematic
  real(dp), dimension(:), allocatable :: cx                   ! implicit coefficient (kinematic, eps/k type)
  real(dp), dimension(:), allocatable :: Tt, alpha, pkin      ! time scale, blending function, production P
  real(dp), dimension(:), allocatable :: dif, suv, spv
  logical, dimension(:), allocatable :: fixed
  real(dp), dimension(:), allocatable :: fixval
  integer, dimension(6) :: stamps = -1

  private

  public :: rs, alpha, dAlphadxi, Tt, ixx, ixy, ixz, iyy, iyz, izz, iEps, CMU, CMU25, CMU75, CEPS1, CEPS2, A1EB
  public :: rsm_allocate, rsm_shift_time, rsm_sources, rsm_solve_stresses, rsm_solve_epsilon, rsm_update_k
  public :: rsm_effective_viscosity, rsm_inlet

contains


subroutine rsm_allocate( lowRe )
!
! Allocate the stress fields and start from isotropic turbulence R_ij = 2/3 k delta_ij (walls: zero if lowRe).
!
  implicit none
  logical, intent(in) :: lowRe
  integer :: ib, i, ijb, m

  allocate( rs(numTotal,6), rso(numTotal,6), rsoo(numTotal,6), drs(3,numTotal,6), sx(numCells,6) )
  allocate( cx(numCells), Tt(numCells), alpha(numTotal), pkin(numCells), dAlphadxi(3,numTotal), fixvalR(numCells,6) )
  dAlphadxi = 0.0_dp;  fixvalR = 0.0_dp
  allocate( dif(numTotal), suv(numCells), spv(numCells), fixed(numCells), fixval(numCells) )

  rs = 0.0_dp
  rs(:,ixx) = twothirds*te
  rs(:,iyy) = twothirds*te
  rs(:,izz) = twothirds*te
  alpha = 1.0_dp
  if ( lowRe ) then
    do ib=1,numBoundaries
      if ( bctype(ib) /= 'wall' ) cycle
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        rs(ijb,:) = 0.0_dp
        alpha(ijb) = 0.0_dp
      enddo
    enddo
  endif
  rso = rs;  rsoo = rs
  drs = 0.0_dp
  Tt = 1.0_dp
  do m=1,6
    call grad( rs(:,m), drs(:,:,m) )
  enddo

  write(*,'(a)') ' '
  write(*,'(a)') '  **Allocated and initialised Reynolds stress fields.'
  write(*,'(a)') ' '

end subroutine


subroutine rsm_shift_time
  implicit none
  integer :: m
  do m=1,6
    call shift_time_levels( rs(:,m), rso(:,m), rsoo(:,m), stamps(m) )
  enddo
end subroutine


subroutine rsm_update_k
!
! k = R_kk/2 with realizability limits on the diagonal and the Schwarz inequality on the shear stresses.
!
  implicit none
  integer :: inp
  real(dp) :: lim

  do inp=1,numCells
    rs(inp,ixx) = max( rs(inp,ixx), small )
    rs(inp,iyy) = max( rs(inp,iyy), small )
    rs(inp,izz) = max( rs(inp,izz), small )
    lim = sqrt( rs(inp,ixx)*rs(inp,iyy) );  rs(inp,ixy) = max( -lim, min( lim, rs(inp,ixy) ) )
    lim = sqrt( rs(inp,ixx)*rs(inp,izz) );  rs(inp,ixz) = max( -lim, min( lim, rs(inp,ixz) ) )
    lim = sqrt( rs(inp,iyy)*rs(inp,izz) );  rs(inp,iyz) = max( -lim, min( lim, rs(inp,iyz) ) )
    te(inp) = 0.5_dp*( rs(inp,ixx) + rs(inp,iyy) + rs(inp,izz) )
  enddo
  call updateBoundary( te )

end subroutine


subroutine rsm_sources( lowRe )
!
! Explicit sources sx = P_ij + Phi_ij - eps_ij (evaluated with the stresses of the previous iteration), implicit
! coefficient cx (a positive multiple of eps/k that is added to both sides: it makes the matrix diagonally dominant
! without changing the converged solution), production P, time scale T, and the blending function derivatives.
!
  implicit none
  logical, intent(in) :: lowRe
  integer :: inp, m, i, j, kk, l
  real(dp), dimension(3,3) :: R, g, S, W, a, P, Phih, Phiw, epsij, Phi
  real(dp), dimension(3) :: n
  real(dp) :: k, eps, al2, a2, pk, trasS, rnn, nu, gn, dd

  nu = viscos/densit

  do inp=1,numCells

    k = max( te(inp), small )
    eps = max( ed(inp), small )

    do m=1,6
      R(ii(m),jj(m)) = rs(inp,m);  R(jj(m),ii(m)) = rs(inp,m)
    enddo
    g(1,:) = dUdxi(:,inp);  g(2,:) = dVdxi(:,inp);  g(3,:) = dWdxi(:,inp)
    S = 0.5_dp*( g + transpose(g) )
    W = 0.5_dp*( g - transpose(g) )

    ! Production
    P = 0.0_dp
    do i=1,3
      do j=1,3
        do kk=1,3
          P(i,j) = P(i,j) - R(i,kk)*g(j,kk) - R(j,kk)*g(i,kk)
        enddo
      enddo
    enddo
    pk = 0.5_dp*( P(1,1) + P(2,2) + P(3,3) )
    pkin(inp) = pk
    gen(inp) = den(inp)*pk

    ! Anisotropy
    a = R/k
    do i=1,3
      a(i,i) = a(i,i) - twothirds
    enddo
    a2 = sum( a*a )
    trasS = sum( a*S )

    ! Homogeneous (SSG) pressure-strain
    Phih = -( C1 + C1S*pk/eps )*eps*a + ( C3 - C3S*sqrt(a2) )*k*S
    do i=1,3
      do j=1,3
        do kk=1,3
          Phih(i,j) = Phih(i,j) + C2*eps*a(i,kk)*a(kk,j) &
                    + C4*k*( a(i,kk)*S(j,kk) + a(j,kk)*S(i,kk) ) &
                    + C5*k*( a(i,kk)*W(j,kk) + a(j,kk)*W(i,kk) )
        enddo
      enddo
      Phih(i,i) = Phih(i,i) - C2*eps*a2/3.0_dp - twothirds*C4*k*trasS
    enddo

    ! Wall (elliptic blending) part and blended dissipation
    al2 = alpha(inp)**2
    epsij = al2*twothirds*eps*reshape( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    Phiw = 0.0_dp
    if ( lowRe ) then
      gn = sqrt( sum( dAlphadxi(:,inp)**2 ) )
      if ( gn > small ) then
        n = dAlphadxi(:,inp)/gn
        rnn = 0.0_dp
        do i=1,3
          do j=1,3
            rnn = rnn + R(i,j)*n(i)*n(j)
          enddo
        enddo
        do i=1,3
          do j=1,3
            dd = 0.0_dp;  if (i==j) dd = 1.0_dp
            Phiw(i,j) = -5.0_dp*eps/k*( dot_product(R(i,:),n)*n(j) + dot_product(R(j,:),n)*n(i) - 0.5_dp*rnn*( n(i)*n(j) + dd ) )
          enddo
        enddo
      endif
      epsij = epsij + (1.0_dp-al2)*eps/k*R
    endif
    Phi = (1.0_dp-al2)*Phiw + al2*Phih

    do m=1,6
      i = ii(m);  j = jj(m)
      sx(inp,m) = P(i,j) + Phi(i,j) - epsij(i,j)
    enddo

    ! Implicit coefficient (positive): dissipation, slow return-to-isotropy, wall relaxation
    cx(inp) = eps/k*( (1.0_dp-al2)*6.0_dp + al2*( C1 + C1S*max(pk,0.0_dp)/eps ) )   ! 6 = 1 (dissipation) + 5 (Phi^w)

    ! Time scale of the model
    if ( lowRe ) then
      Tt(inp) = max( k/eps, 6.0_dp*sqrt(nu/eps) )
    else
      Tt(inp) = k/eps
    endif

  enddo

end subroutine


function ggdh( inp, gr, c )
!
! Gradient-aligned generalised gradient diffusion: Cs T R_kl g_k g_l/|g|^2 (the Daly-Harlow diffusivity tensor contracted
! with the direction of the gradient g of the diffused quantity, bounded by the isotropic value). Scalar form for the engine.
!
  implicit none
  integer, intent(in) :: inp
  real(dp), dimension(3), intent(in) :: gr
  real(dp), intent(in) :: c
  real(dp) :: ggdh, g2, rgg

  g2 = sum(gr*gr)
  if ( g2 > small ) then
    rgg = rs(inp,ixx)*gr(1)**2 + rs(inp,iyy)*gr(2)**2 + rs(inp,izz)*gr(3)**2 &
        + 2.0_dp*( rs(inp,ixy)*gr(1)*gr(2) + rs(inp,ixz)*gr(1)*gr(3) + rs(inp,iyz)*gr(2)*gr(3) )
    rgg = rgg/g2
  else
    rgg = twothirds*te(inp)
  endif
  ggdh = c*Tt(inp)*max( rgg, 0.1_dp*te(inp) )

end function


subroutine rsm_solve_stresses( lowRe )
!
! The six stress equations: sources from 'rsm_sources', with cx R_old added to both sides. Low-Re: R_w = 0 (Dirichlet).
! High-Re: wall cells hold the log-layer values (wall_stresses).
!
  implicit none
  logical, intent(in) :: lowRe
  integer :: inp, m

  fixed = .false.
  if ( .not.lowRe ) call wall_stresses

  do m=1,6

    do inp=1,numCells
      dif(inp) = viscos + den(inp)*ggdh( inp, drs(:,inp,m), CS )
    enddo
    call copy_owner_to_boundary( dif )

    do inp=1,numCells
      suv(inp) = den(inp)*( sx(inp,m) + cx(inp)*rs(inp,m) )*vol(inp)
      spv(inp) = den(inp)*cx(inp)*vol(inp)
    enddo
    fixval = fixvalR(:,m)

    call transport_scalar( rs(:,m), drs(:,:,m), rso(:,m), rsoo(:,m), m, 'R_ij', suv, spv, dif, &
                           wallDirichlet=lowRe, fixed=fixed, fixval=fixval, positive=.false. )

  enddo

end subroutine


subroutine wall_stresses
!
! Log-layer values of the stresses in wall-adjacent cells, from the wall shear stress of the previous iteration:
!   k_w = u_tau^2/sqrt(Cmu),  R_ij = k_w ( 1.214 t_i t_j + 0.238 n_i n_j + 0.548 s_i s_j + sqrt(Cmu) (t_i n_j + n_i t_j) )
! (ratios 5.1 : 1.0 : 2.3 : 1.0 of the log layer, normalised to 4.2 = tr R/2 / u_tau^2; n points out of the fluid).
!
  implicit none
  integer :: ib, i, iface, ijp, iWall, m
  real(dp) :: utau2, kw, nxf, nyf, nzf, vtp, xtp, ytp, ztp
  real(dp), dimension(3) :: nn, tt3, ss3
  real(dp), dimension(3,3) :: R

  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      call wall_normal_tangent( iface, ijp, nxf, nyf, nzf, vtp, xtp, ytp, ztp )
      if ( vtp <= small ) cycle                       ! no wall-parallel flow: leave the cell free
      nn = [nxf,nyf,nzf];  tt3 = [xtp,ytp,ztp]
      ss3 = [ nn(2)*tt3(3)-nn(3)*tt3(2), nn(3)*tt3(1)-nn(1)*tt3(3), nn(1)*tt3(2)-nn(2)*tt3(1) ]
      utau2 = max( abs(tau(iWall))/den(ijp), small )
      kw = utau2/sqrt(CMU)
      R = 0.0_dp
      do m=1,3
        R(m,:) = kw*( 1.214_dp*tt3(m)*tt3 + 0.238_dp*nn(m)*nn + 0.548_dp*ss3(m)*ss3 + 0.30_dp*( tt3(m)*nn + nn(m)*tt3 ) )
      enddo
      fixed(ijp) = .true.
      do m=1,6
        fixvalR(ijp,m) = R(ii(m),jj(m))
      enddo
    enddo
  enddo

end subroutine


subroutine rsm_solve_epsilon( lowRe )
!
! Dissipation rate equation. High-Re: wall-cell value eps = Cmu^3/4 k^3/2/(kappa y). Low-Re: eps = 2 nu k/y^2.
!
  implicit none
  logical, intent(in) :: lowRe
  integer :: inp, ib, i, iface, ijp, iWall
  real(dp) :: ce1, kk

  do inp=1,numCells
    kk = max(te(inp),small)
    ce1 = CEPS1
    if ( lowRe ) ce1 = CEPS1*( 1.0_dp + A1EB*( 1.0_dp - alpha(inp)**3 )*max(pkin(inp),0.0_dp)/max(ed(inp),small) )
    suv(inp) = den(inp)*ce1*max(pkin(inp),0.0_dp)*ed(inp)/kk*vol(inp)
    if ( lowRe ) suv(inp) = den(inp)*ce1*max(pkin(inp),0.0_dp)/Tt(inp)*vol(inp)
    spv(inp) = den(inp)*CEPS2*ed(inp)/kk*vol(inp)
    if ( lowRe ) spv(inp) = den(inp)*CEPS2/Tt(inp)*vol(inp)
    dif(inp) = viscos + den(inp)*ggdh( inp, dEDdxi(:,inp), CEPS )
  enddo
  call copy_owner_to_boundary( dif )

  fixed = .false.
  iWall = 0
  do ib=1,numBoundaries
    if ( bctype(ib) /= 'wall' ) cycle
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)
      iWall = iWall + 1
      fixed(ijp) = .true.
      if ( lowRe ) then
        fixval(ijp) = 2.0_dp*viscos/den(ijp)*te(ijp)/dnw(iWall)**2
      else
        fixval(ijp) = CMU75*te(ijp)**1.5_dp/(cappa*dnw(iWall))
      endif
    enddo
  enddo

  call transport_scalar( ed, dEDdxi, edo, edoo, iEps, 'epsilon', suv, spv, dif, fixed=fixed, fixval=fixval )

  ! Floor: nu_t/nu <= 1e4. An epsilon overshoot to ~0 gives T -> infinity and destroys the momentum solution.
  do inp=1,numCells
    ed(inp) = max( ed(inp), CMU*te(inp)**2/(1.0e4_dp*viscos/densit) )
  enddo
  call updateBoundary( ed )

end subroutine


subroutine rsm_effective_viscosity( lowRe )
!
! mu_a = rho CA k T goes to 'vis' (implicit), the excess anisotropy to 'aij' (explicit source in the momentum equations).
! Walls: wall functions (high-Re) or molecular viscosity (low-Re).
!
  implicit none
  logical, intent(in) :: lowRe
  integer :: inp, ib, i, iface, ijp, ijb, iWall
  real(dp) :: k, mua, urf, visold, nxf, nyf, nzf, vtp, xtp, ytp, ztp, ut2, utau
  real(dp), dimension(3,3) :: g, S

  urf = TurbModel%urfVis

  do inp=1,numCells
    k = max( te(inp), small )
    visold = vis(inp)
    mua = min( CA*k*Tt(inp), MUA_MAX*viscos/densit )   ! kinematic artificial viscosity, bounded
    vis(inp) = viscos + den(inp)*mua
    vis(inp) = urf*vis(inp) + (1.0_dp-urf)*visold

    g(1,:) = dUdxi(:,inp);  g(2,:) = dVdxi(:,inp);  g(3,:) = dWdxi(:,inp)
    S = 0.5_dp*( g + transpose(g) )
    if ( .not.lasm ) cycle
    aij(1,inp) = rs(inp,ixx)/k - twothirds + 2.0_dp*mua/k*S(1,1)
    aij(2,inp) = rs(inp,ixy)/k             + 2.0_dp*mua/k*S(1,2)
    aij(3,inp) = rs(inp,ixz)/k             + 2.0_dp*mua/k*S(1,3)
    aij(4,inp) = rs(inp,iyy)/k - twothirds + 2.0_dp*mua/k*S(2,2)
    aij(5,inp) = rs(inp,iyz)/k             + 2.0_dp*mua/k*S(2,3)
    aij(6,inp) = rs(inp,izz)/k - twothirds + 2.0_dp*mua/k*S(3,3)
  enddo

  call updateBoundary( vis )

  if ( lowRe ) then
    iWall = 0
    do ib=1,numBoundaries
      if ( bctype(ib) /= 'wall' ) cycle
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1
        call wall_normal_tangent( iface, ijp, nxf, nyf, nzf, vtp, xtp, ytp, ztp )
        ut2 = abs( (u(ijb)-u(ijp))*xtp + (v(ijb)-v(ijp))*ytp + (w(ijb)-w(ijp))*ztp )
        tau(iWall) = viscos*ut2/dnw(iWall)
        utau = sqrt( tau(iWall)/den(ijb) )
        ypl(iWall) = den(ijb)*utau*dnw(iWall)/viscos
        visw(iWall) = viscos
        vis(ijb) = viscos
      enddo
    enddo
  else
    call wall_viscosity_loglaw( CMU25 )
  endif

end subroutine


subroutine rsm_inlet
!
! Inlet: isotropic stresses from k, viscosity from k and epsilon.
!
  implicit none
  integer :: i, ib, ijb

  if ( .not.allocated(rs) ) return

  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        rs(ijb,:) = 0.0_dp
        rs(ijb,ixx) = twothirds*te(ijb);  rs(ijb,iyy) = twothirds*te(ijb);  rs(ijb,izz) = twothirds*te(ijb)
        vis(ijb) = viscos + den(ijb)*CMU*te(ijb)**2/(ed(ijb)+small)
      end do
    endif
  enddo

end subroutine

end module
