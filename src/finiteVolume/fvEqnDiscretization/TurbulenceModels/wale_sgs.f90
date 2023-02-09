module wale_sgs

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of the WALE sub-grid scale model
!
!  Discussion:
!
!    Reference:
!    F. Nicoud and F. Ducros, Subgrid-scale stress modelling based on the square of the velocity gradient tensor. 
!    Flow, Turbulence and Combustion, 62(3), 183-200, 1999.
!    
!    Also Fumiya Nozaki's blog was very useful for implementation details using tensor operators. Thanks Fumiya!
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use parameters
  use geometry
  use variables, only: u,v,w,den,vis,visw
  use TurbModelData, only: TurbModel

  implicit none

  real(dp), parameter :: Cw_ = 0.325

  private 

  public :: modify_viscosity_wale_sgs

contains


subroutine modify_viscosity_wale_sgs
!
!  Purpose: 
!
!    Main routine of the module - return mu_sgs for the WALE model.
!
!  Discussion:
!
!    We will do this with tensor field manipulation functionality.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    02 November 2020
!
!  Author:
!
!    Nikola Mirkov
!    Email: largeddysimulation@gmail.com / nikolamirkov@yahoo.com / nmirkov@vin.bg.ac.rs
!
!  Reference:
!
!    F. Nicoud and F. Ducros, Subgrid-scale stress modelling based on the square of the velocity gradient tensor. 
!    Flow, Turbulence and Combustion, 62(3), 183-200, 1999.
!
  use tensorFields  ! <--- Let us see these in action!
  use fvxGradient   ! <--

  implicit none

!
! Local variables
!
  integer :: i, ib, iface, ijp, ijn, ijb, ijbt, iWall, iPer

  real(dp) :: urf, rmin, rmax

  real(dp), parameter :: r13 = 1./3._dp

  type(volScalarField) :: magSqrSd

  type(volTensorField) :: D,Sd

  type(volVectorField) :: U_

  U_ = volVectorField("Velocity", u, v, w )

  ! Velocity gradient tensor.
  D = Grad( U_ )

  ! The traceless symmetric part of the square of the velocity gradient tensor.
  Sd = .dev.(.symm.(D*D))

  ! Its magnitude squared 
  magSqrSd = .magSq.Sd
   

  ! v---this should be 'musgs' but we reuse existing field
  magSqrSd = den * ( Cw_* Vol**r13 )**2 * &
                             power( magSqrSd, 3./2._dp ) &
           /  ( power( .magSq.(.symm.D), 5./2._dp ) + power( magSqrSd, 5./4._dp )  ) 


  urf = TurbModel%urfVis

  ! Update of the effective viscosity and underelaxation
  vis( 1:numCells ) = urf*( magSqrSd%mag(1:numCells) + viscos ) + (1.0-urf) * vis( 1:numCells )


  ! Boundary faces 

  iWall = 0
  iPer = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)

        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(viscos,zero)
        vis(ijb) = visw(iWall)

      enddo


    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1

      ! Faces trough periodic boundaries, Taiwo first
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        iface = startFaceTwin(iPer) + i
        ijn = owner(iface)

        vis(ijb) = half*( vis(ijp)+vis(ijn) )

        ! Now find where is the twin in field array
        ijbt = numCells + ( startFaceTwin(iPer) - numInnerFaces ) + i
        
        ! Twin takes the same values
        vis(ijbt) = vis(ijb)


      enddo


    else

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    endif 

  enddo

  write(*,*) " Updated effective viscosity - WALE model."

  rmin = minval(vis/viscos)
  rmax = maxval(vis/viscos)
 
  write(*,'(2x,es11.4,a,es11.4)') rmin,' <= Viscosity ratio <= ',rmax

end subroutine


end module wale_sgs