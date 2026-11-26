module ebrsm
!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Elliptic-blending Reynolds stress model EB-RSM (Manceau & Hanjalic 2002, Manceau 2015), low-Re, y+ ~ 1 mesh.
!
!  Discussion:
!
!    Stress transport + epsilon equation of 'rsm_common' with the wall blending function alpha from
!
!      alpha - L^2 lapl(alpha) = 1,   alpha_w = 0,   L = CL max( k^1.5/eps, Ceta (nu^3/eps)^0.25 ),  CL = 0.161, Ceta = 80
!
!    n = grad alpha/|grad alpha| enters the wall pressure-strain Phi^w; the elliptic equation is solved by the scalar
!    engine without convection and time term (lconv=.false., ltrans=.false.), exactly as 'f' in k-epsilon-zeta-f.
!    T = max( k/eps, 6 sqrt(nu/eps) ); Ceps1 = 1.44 ( 1 + A1 (1 - alpha^3) P/eps ), A1 = 0.03.
!    Walls: R_ij = 0, alpha = 0 (Dirichlet), epsilon in the wall cell 2 nu k/y^2.
!    TurbModel%Scalar slots: 1..6 stresses, 7 epsilon, 8 alpha.
!
!  Reference:
!    [1] R. Manceau, K. Hanjalic, Phys. Fluids 14(2), 2002, pp. 744-754.
!    [2] R. Manceau, Int. J. Heat Fluid Flow 51, 2015, pp. 395-400.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use parameters
  use geometry
  use variables
  use TurbModelData, only: TurbModel
  use turb_transport
  use rsm_common

  implicit none

  real(dp), parameter :: CL = 0.161_dp, CETA = 80.0_dp
  integer, parameter :: iAlpha = 8
  real(dp), dimension(:), allocatable :: alphao, dif, suv, spv

  private

  public :: modify_viscosity_ebrsm
  public :: modify_viscosity_inlet_ebrsm

contains

subroutine modify_viscosity_ebrsm
  implicit none

  if ( .not.allocated(rs) ) then
    call rsm_allocate( .true. )
    allocate( alphao(numTotal), dif(numTotal), suv(numCells), spv(numCells) )
    alphao = alpha
  endif

  call rsm_shift_time
  call calcsc_alpha
  call rsm_sources( .true. )
  call rsm_solve_stresses( .true. )
  call rsm_update_k
  call rsm_solve_epsilon( .true. )
  call rsm_effective_viscosity( .true. )

end subroutine


subroutine calcsc_alpha
!
! Elliptic blending equation, alpha limited to [0,1].
!
  implicit none
  integer :: inp
  real(dp) :: nu, ll, lmax, szeta

  nu = viscos/densit
  lmax = maxval( wallDistance(1:numCells) )   ! length scale cannot exceed the size of the domain

  do inp=1,numCells
    szeta = sqrt(6.0_dp)*0.22_dp*twothirds*magStrain(inp) + small                      ! realizability bound as in zeta-f
    ll = CL*max( min( te(inp)**1.5_dp/(ed(inp)+small), sqrt(te(inp))/szeta ), CETA*sqrt(sqrt( nu**3/(ed(inp)+small) )) )
    ll = min( ll, lmax )
    dif(inp) = den(inp)*ll*ll
    suv(inp) = den(inp)*vol(inp)
    spv(inp) = den(inp)*vol(inp)
  enddo
  call copy_owner_to_boundary( dif )

  call transport_scalar( alpha, dAlphadxi, alphao, alphao, iAlpha, 'alpha', suv, spv, dif, &
                         lconv=.false., ltrans=.false., wallDirichlet=.true., positive=.false. )

  alpha(1:numCells) = max( 0.0_dp, min( 1.0_dp, alpha(1:numCells) ) )

end subroutine


subroutine modify_viscosity_inlet_ebrsm
  implicit none
  call rsm_inlet
end subroutine

end module
