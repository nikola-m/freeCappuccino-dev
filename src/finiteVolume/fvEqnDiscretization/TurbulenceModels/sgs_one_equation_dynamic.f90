module sgs_one_equation_dynamic
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Localized dynamic one-equation subgrid-scale model for LES (Kim & Menon 1995, Menon & Kim 1996), on unstructured meshes.
!
!  Discussion:
!
!    Transport equation for the SGS kinetic energy, on top of 'turb_transport' (scalar slot 1, k = te):
!
!      Dk/Dt = P - Ceps k^3/2 / Delta + div[ (nu + nut) grad k ],     P = nut 2 S_ij S_ij,     nut = Ck Delta sqrt(k)
!
!    Delta = V^(1/3). Ck and Ceps are computed locally and dynamically with a test filter (hat):
!
!      L_ij = (u_i u_j)^ - u^_i u^_j,           k_test = L_kk/2
!      M_ij = -2 [ Delta^ sqrt(k_test) S^_ij - Delta ( sqrt(k) S_ij )^ ]               (deviatoric part)
!      Ck   = < L^d_ij M_ij > / < M_ij M_ij >
!      Ceps = nu [ ( |grad u|^2 )^ - |grad u^|^2 ] / [ k_test^3/2 / Delta^ - ( k^3/2 / Delta )^ ]       (Menon-Kim, simplified)
!
!    The test filter is FACE BASED: phi^_P = ( phi_P + sum_nb phi_nb ) / ( 1 + n_nb ), no directions, no structured mesh
!    required. The ratio Delta^/Delta is evaluated per cell from the second moment of the filter kernel,
!    Delta^2 = Delta^2 + 12 sigma^2,  sigma^2 = (1/3) sum_nb |x_nb - x_P|^2 /(1+n_nb)   (2.1 for a uniform hex mesh: 7-point stencil).
!    The brackets <.> in Ck, Ceps (needed because the local formula is ill-conditioned) are the same face-based filter, followed by
!    clipping to [0, CK_MAX], [CEPS_MIN, CEPS_MAX] and relaxation in time. Negative Ck (backscatter) is clipped: no
!    backscatter model. No wall damping: use on wall-resolved meshes (nut ~ sqrt(k) -> 0 at the wall as k is zero there),
!    wall-function treatment is not provided. Not tested against DNS databases - see the technical documentation.
!
!  Reference:
!    [1] W.-W. Kim, S. Menon, A new dynamic one-equation subgrid-scale model for large eddy simulations, AIAA 95-0356, 1995.
!    [2] S. Menon, W.-W. Kim, High Reynolds number flow simulations using the localized dynamic subgrid-scale model, AIAA 96-0425.
!    [3] M. Germano, U. Piomelli, P. Moin, W. Cabot, Phys. Fluids A 3(7), 1991.
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

  real(dp), parameter :: CK_MAX = 0.2_dp
  real(dp), parameter :: CEPS_MIN = 0.3_dp, CEPS_MAX = 2.5_dp
  real(dp), parameter :: CK0 = 0.07_dp, CEPS0 = 1.0_dp          ! initial values
  real(dp), parameter :: URF_C = 0.5_dp                         ! relaxation of the dynamic coefficients

  integer, parameter, dimension(6) :: ii6 = [1,1,1,2,2,3], jj6 = [1,2,3,2,3,3]
  real(dp), dimension(:), allocatable :: delta, ratio, nbcount, ck, ceps, dif, suv, spv
  logical :: first = .true.

  private

  public :: modify_viscosity_sgs_one_equation_dynamic
  public :: modify_viscosity_inlet_sgs_one_equation_dynamic

contains


subroutine modify_viscosity_sgs_one_equation_dynamic
  implicit none

  if ( first ) call allocate_and_initialise

  call dynamic_coefficients
  call calcsc_ksgs
  call modify_mu_eff

end subroutine


subroutine allocate_and_initialise
!
! Filter width Delta, number of neighbours and per-cell ratio Delta^/Delta from the second moment of the face-based filter.
!
  implicit none
  integer :: i, ijp, ijn, inp
  real(dp), dimension(:), allocatable :: m2
  real(dp) :: d2

  allocate( delta(numCells), ratio(numCells), nbcount(numCells), ck(numCells), ceps(numCells) )
  allocate( dif(numTotal), suv(numCells), spv(numCells), m2(numCells) )

  nbcount = 1.0_dp
  m2 = 0.0_dp
  do i=1,numInnerFaces
    ijp = owner(i);  ijn = neighbour(i)
    d2 = (xc(ijn)-xc(ijp))**2 + (yc(ijn)-yc(ijp))**2 + (zc(ijn)-zc(ijp))**2
    nbcount(ijp) = nbcount(ijp) + 1.0_dp;  nbcount(ijn) = nbcount(ijn) + 1.0_dp
    m2(ijp) = m2(ijp) + d2;  m2(ijn) = m2(ijn) + d2
  enddo

  do inp=1,numCells
    delta(inp) = vol(inp)**(1.0_dp/3.0_dp)
    ratio(inp) = sqrt( 1.0_dp + 12.0_dp*( m2(inp)/(3.0_dp*nbcount(inp)) )/delta(inp)**2 )
  enddo

  ck = CK0
  ceps = CEPS0
  first = .false.

  write(*,'(a,2(1x,f7.3))') '  **SGS dynamic model: test filter width ratio, min and max over cells:', minval(ratio), maxval(ratio)

end subroutine


subroutine test_filter( phi, phih )
!
! Face-based test filter of a cell field: average of the cell and its face neighbours.
!
  implicit none
  real(dp), dimension(numCells), intent(in) :: phi
  real(dp), dimension(numCells), intent(out) :: phih
  integer :: i, ijp, ijn

  phih = phi
  do i=1,numInnerFaces
    ijp = owner(i);  ijn = neighbour(i)
    phih(ijp) = phih(ijp) + phi(ijn)
    phih(ijn) = phih(ijn) + phi(ijp)
  enddo
  phih = phih/nbcount

end subroutine


subroutine dynamic_coefficients
!
! Local Ck and Ceps from the Germano identity with the test filter.
!
  implicit none
  integer :: inp, m
  real(dp), dimension(:), allocatable :: uh, vh, wh, s2, s2h, g2, g2h, rk, rkh, wk
  real(dp), dimension(:,:), allocatable :: Lij, Sh, SqS, SqSh, Mij, LM, MM
  real(dp), dimension(3,3) :: g, S, Ld, Mt
  real(dp) :: ktest, trL, num, den_, nu, dtest, ceps_new, ck_new

  allocate( uh(numCells), vh(numCells), wh(numCells), s2(numCells), s2h(numCells), g2(numCells), g2h(numCells) )
  allocate( rk(numCells), rkh(numCells), wk(numCells) )
  allocate( Lij(numCells,6), Sh(numCells,6), SqS(numCells,6), SqSh(numCells,6), Mij(numCells,6), LM(numCells,1), MM(numCells,1) )

  nu = viscos/densit

  ! Resolved velocity products, strain components and gradient magnitude, test filtered.
  call test_filter( u(1:numCells), uh );  call test_filter( v(1:numCells), vh );  call test_filter( w(1:numCells), wh )

  call test_filter( u(1:numCells)*u(1:numCells), Lij(:,1) );  Lij(:,1) = Lij(:,1) - uh*uh
  call test_filter( u(1:numCells)*v(1:numCells), Lij(:,2) );  Lij(:,2) = Lij(:,2) - uh*vh
  call test_filter( u(1:numCells)*w(1:numCells), Lij(:,3) );  Lij(:,3) = Lij(:,3) - uh*wh
  call test_filter( v(1:numCells)*v(1:numCells), Lij(:,4) );  Lij(:,4) = Lij(:,4) - vh*vh
  call test_filter( v(1:numCells)*w(1:numCells), Lij(:,5) );  Lij(:,5) = Lij(:,5) - vh*wh
  call test_filter( w(1:numCells)*w(1:numCells), Lij(:,6) );  Lij(:,6) = Lij(:,6) - wh*wh

  rk = sqrt( max(te(1:numCells),0.0_dp) )

  do inp=1,numCells
    g(1,:) = dUdxi(:,inp);  g(2,:) = dVdxi(:,inp);  g(3,:) = dWdxi(:,inp)
    S = 0.5_dp*(g + transpose(g))
    g2(inp) = sum(g*g)
    do m=1,6
      SqS(inp,m) = S(ii6(m),jj6(m))
    enddo
  enddo
  do m=1,6
    call test_filter( SqS(:,m), Sh(:,m) )                        ! S^_ij
    call test_filter( rk*SqS(:,m), SqSh(:,m) )                   ! ( sqrt(k) S_ij )^
  enddo
  call test_filter( g2, g2h )                                    ! ( |grad u|^2 )^
  call test_filter( rk, rkh )

  ! Local Ck: numerator <L^d M>, denominator <M M>
  do inp=1,numCells
    trL = ( Lij(inp,1) + Lij(inp,4) + Lij(inp,6) )
    ktest = max( 0.5_dp*trL, 0.0_dp )
    dtest = delta(inp)*ratio(inp)

    do m=1,6
      Mij(inp,m) = -2.0_dp*( dtest*sqrt(ktest)*Sh(inp,m) - delta(inp)*SqSh(inp,m) )
    enddo
    ! deviatoric parts
    trL = trL/3.0_dp
    Ld = 0.0_dp
    do m=1,6
      Ld(ii6(m),jj6(m)) = Lij(inp,m);  Ld(jj6(m),ii6(m)) = Lij(inp,m)
    enddo
    Ld(1,1) = Ld(1,1) - trL;  Ld(2,2) = Ld(2,2) - trL;  Ld(3,3) = Ld(3,3) - trL
    Mt = 0.0_dp
    do m=1,6
      Mt(ii6(m),jj6(m)) = Mij(inp,m);  Mt(jj6(m),ii6(m)) = Mij(inp,m)
    enddo
    trL = (Mt(1,1)+Mt(2,2)+Mt(3,3))/3.0_dp
    Mt(1,1) = Mt(1,1) - trL;  Mt(2,2) = Mt(2,2) - trL;  Mt(3,3) = Mt(3,3) - trL

    LM(inp,1) = sum(Ld*Mt)
    MM(inp,1) = sum(Mt*Mt)
  enddo

  call test_filter( LM(:,1), wk );  LM(:,1) = wk              ! brackets <.>
  call test_filter( MM(:,1), wk );  MM(:,1) = wk

  ! Local Ceps: viscous dissipation of the band between the filters over the model dissipation difference
  do inp=1,numCells
    ktest = max( 0.5_dp*( Lij(inp,1) + Lij(inp,4) + Lij(inp,6) ), 0.0_dp )
    dtest = delta(inp)*ratio(inp)
    rkh(inp) = ktest**1.5_dp/dtest - ( rk(inp)**3/delta(inp) )     ! unfiltered k^3/2/Delta, filtered below
  enddo
  ! numerator and denominator are filtered (brackets) before the division
  ! nu [ (|grad u|^2)^ - |grad u^|^2 ]: |grad u^|^2 is approximated by the filtered-field gradient energy of the cell field
  s2 = nu*max( g2h - g2, 0.0_dp )
  call test_filter( s2, s2h )
  call test_filter( rkh, wk )

  do inp=1,numCells
    if ( MM(inp,1) > small ) then
      ck_new = LM(inp,1)/MM(inp,1)
    else
      ck_new = CK0
    endif
    ck_new = max( 0.0_dp, min( CK_MAX, ck_new ) )

    if ( abs(wk(inp)) > small ) then
      ceps_new = s2h(inp)/wk(inp)
    else
      ceps_new = CEPS0
    endif
    ceps_new = max( CEPS_MIN, min( CEPS_MAX, ceps_new ) )

    ck(inp)   = URF_C*ck_new   + (1.0_dp-URF_C)*ck(inp)
    ceps(inp) = URF_C*ceps_new + (1.0_dp-URF_C)*ceps(inp)
  enddo

  write(*,'(2x,a,f7.4,a,f7.4,a,f7.3,a,f7.3)') 'Ck: mean ', sum(ck)/numCells, ' max ', maxval(ck), '   Ceps: mean ', &
        sum(ceps)/numCells, ' max ', maxval(ceps)

  deallocate( uh, vh, wh, s2, s2h, g2, g2h, rk, rkh, wk, Lij, Sh, SqS, SqSh, Mij, LM, MM )

end subroutine


subroutine calcsc_ksgs
  implicit none
  integer :: inp

  do inp=1,numCells
    gen(inp) = (vis(inp)-viscos)*magStrain(inp)**2
    suv(inp) = max(gen(inp),zero)*vol(inp)
    spv(inp) = den(inp)*ceps(inp)*sqrt(max(te(inp),zero))/delta(inp)*vol(inp) - min(gen(inp),zero)*vol(inp)/(te(inp)+small)
  enddo

  dif = vis
  call transport_scalar( te, dTEdxi, teo, teoo, 1, 'k_sgs', suv, spv, dif, wallDirichlet=.true. )

end subroutine


subroutine modify_mu_eff
!
! nut = Ck Delta sqrt(k_sgs); wall faces: molecular viscosity.
!
  implicit none
  integer :: inp, ib, i, iface, ijp, ijb, iWall
  real(dp) :: nxf, nyf, nzf, vtp, xtp, ytp, ztp, ut2, utau

  do inp=1,numCells
    vis(inp) = viscos + den(inp)*ck(inp)*delta(inp)*sqrt(max(te(inp),zero))
  enddo
  call updateBoundary( vis )

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

end subroutine


subroutine modify_viscosity_inlet_sgs_one_equation_dynamic
  implicit none
  integer :: i, ib, ijb

  do ib=1,numBoundaries
    if ( bctype(ib) == 'inlet' ) then
      do i=1,nfaces(ib)
        ijb = iBndValueStart(ib) + i
        vis(ijb) = viscos + den(ijb)*CK0*vol(owner(startFace(ib)+i))**(1.0_dp/3.0_dp)*sqrt(max(te(ijb),zero))
      end do
    endif
  enddo

end subroutine

end module
