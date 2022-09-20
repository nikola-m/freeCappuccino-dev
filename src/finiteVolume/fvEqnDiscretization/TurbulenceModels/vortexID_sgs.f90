module vortexID_sgs

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of the Vortex ID sub-grid scale model
!
!  Discussion:
!
!    Reference:
!    XXXX
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

  real(dp), parameter :: CmQ_ = 3.4 
  real(dp), parameter :: Cml2_ = 2.0 
  ! real(dp), parameter :: CmQ_ = 0.325 
  ! real(dp), parameter :: CmQ_ = 0.325 

  private 

  public :: modify_viscosity_vortexID_sgs

contains


subroutine modify_viscosity_vortexID_sgs
!
!  Purpose: 
!
!    Main routine of the module - return mu_sgs for the vortex criteria SGS model.
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
!    20th September 2022
!
!  Author:
!
!    Nikola Mirkov
!    Email: largeddysimulation@gmail.com / nikolamirkov@yahoo.com / nmirkov@vin.bg.ac.rs
!
!  Reference:
!
!    
!
  use tensorFields  ! <--- Let us see these in action!
  use fvxGradient   ! <--

  implicit none

!
! Local variables
!
  integer :: i, ib, iface, ijp, ijn, ijb, ijbt, iWall, iPer

  real(dp) :: urf

  real(dp), parameter :: r23 = 2./3._dp

  type(volScalarField) :: Sij_, lambda2

  type(volTensorField) :: D

  type(volVectorField) :: U_

  U_ = volVectorField("Velocity", u, v, w )

  ! Velocity gradient tensor.
  D = Grad( U_ )

  ! Magnitude of the strain rate tensor (a symmetric part of velocity gradient tensor)
  Sij_ = .mag.(.symm.( D ) )  
   

  ! v---this should be 'musgs' but we reuse existing field
  Sij_ = den * Cm_* Vol**r23  * power( lambda2, 3./2._dp ) * Sij_ &
                           /  ( power( Sij_, 3 ) + power( lambda2, 3./2._dp )  ) 


  urf = TurbModel%urfVis

  ! Update of the effective viscosity and underelaxation
  vis( 1:numCells ) = urf*( Sij_%mag(1:numCells) + viscos ) + (1.0-urf) * vis( 1:numCells )



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

  write(*,*) " Updated effective viscosity - Vortex ID SGS model."

end subroutine


end module vortexID_sgs