module wall_distance

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing subroutines that compute cell center distance to the nearest wall.
!
!  Discussion:
!
!    Some turbulence models, such as k-omega SST, use wall distance field, where
!    by wall distance we consider distance to the nearest wall.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  ! use output, only: vtu_write_scalar_field

  implicit none


  private

  public :: wall_distance_poisson 

contains

subroutine wall_distance_poisson
!
!  Purpose: 
!
!    This subroutine calculates distance to the nearest wall using solution of Poisson equation.
!
!  Discussion:
!
!    The method where Poisson eqn is solved for wall distance is atributed to Brian Spalding.
!    If anyone knows the original reference, please let me know.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!   
!    xx/2015
!    12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    Referened on many places, but I'm not sure which one is original.
!
!  Parameters:
!
!    -
!
  use types
  use parameters
  use geometry
  use linear_solvers
  use gradients, only: grad_gauss

  implicit none

!
! Local variables
!
  integer :: i,ib,iface,ijp,ijb
  integer :: maxiter
  real(dp) :: resn,tolAbs,tolRel
  real(dp), dimension(:), allocatable :: q ! Source for the Poisson eqn.
  real(dp), dimension(:), allocatable :: mu
  real(dp), dimension(:), allocatable :: phi
  real(dp), dimension(:,:), allocatable :: dPdxi

  write(*,*) ' '
  write(*,*) ' Computing distance to the nearest wall:'
  write(*,*) ' '

  allocate( q(numCells), mu(numCells), phi(numTotal), dPdxi(3,numTotal) )

  ! Source term
  q = -Vol(1:numCells)

  ! Initialize solution
  phi = 0.0_dp

  !  Coefficient array for Laplacian
  mu = 1.0_dp       

  ! Laplacian operator and BCs         
  call laplacian(mu,phi) 

  ! We can hardcode these values related to linear solution of the discrete Poisson equation
  maxiter = 500
  tolAbs = 1e-12
  tolRel = 1e-10

  ! Solve linear system of equations
  call csrsolve('iccg', phi, q, resn, maxiter, tolAbs, tolRel, 'Wdis')

  ! Update values at constant gradient bc faces - we need these values for correct gradients

  do ib=1,numBoundaries

    if ( bctype(ib) /= 'wall' ) then
    ! All other boundary faces besides wall which has to be zero.

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        phi(ijb) = phi(ijp)

      enddo

    endif
  enddo


  ! Gradient of solution field stored in phi:
  call grad_gauss(phi,dPdxi)

  ! Wall distance computation from Poisson eq. solution stored in pp:
  wallDistance = -sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) ) + &
                  sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) + 2*phi(1:numCells) )


  deallocate(phi,dPdxi,mu,q)

  ! open(unit=8, file='wall_distance.vtu')
  ! rewind 8
  ! call vtu_write_scalar_field ( 8, 'wall_dist', wallDistance )
  ! close(8)

end subroutine

end module




