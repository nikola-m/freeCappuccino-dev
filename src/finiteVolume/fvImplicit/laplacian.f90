subroutine laplacian(mu,phi)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of Laplacian operator: -div(mu*grad(phi)).
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix

  implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(numCells), intent(in) :: mu

  !
  ! Local variables
  !

  integer :: i, k, ijp, ijn
  integer :: ib,ijb,iface
  real(dp) :: cap, can
  real(dp) :: nxf,nyf,nzf
  real(dp) :: are,dfn,dcoef


  ! Initialize matrix array
  a = 0.0_dp

  ! > Assemble Laplacian system matrix

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxlaplacian(ijp, ijn, arx(i), ary(i), arz(i), facint(i), mu, cap, can)

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

  end do


  ! Modify matrix coefficients to reflect presence of Boundary Conditions in PDE problem.

  do ib=1,numBoundaries

    !if ( bctype(ib) == '...') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
        dfn = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        dcoef   = mu(ijp)*are/dfn

        a( diag(ijp) ) = a( diag(ijp) ) - dcoef  
        su(ijp) = su(ijp) - dcoef*phi(ijb)

      enddo

    !endif
    
  enddo


end subroutine





!***********************************************************************
!
subroutine faceFluxLaplacian(ijp, ijn, arx, ary, arz, lambda, mu, cap, can)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numCells,xc,yc,zc

  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numCells), intent(in) :: mu
  real(dp), intent(inout) :: cap, can

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: xpn,ypn,zpn
  real(dp) :: smdpn

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! _    _     _    _
  ! Sf . Sf / (Sf . dpn)
  smdpn = (arx*arx+ary*ary+arz*arz)/(arx*xpn+ary*ypn+arz*zpn)

  ! Coefficients of discretized Laplace equation
  cap = (fxp*mu(ijp)+fxn*mu(ijn))*smdpn
  can = cap

end subroutine
