module vortexID_sgs

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of the Vortex ID sub-grid scale model.
!
!  Discussion:
!
!    Reference:
!    X. Fang et al. Using vorted identifiers to build eddy-viscosity subgrid-scale models for large-eddy simulation, 
!    Phys Rev Fluids 4, 034606 (2019).
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
  use vortexIdentification, only: lambda2,setLambda2

  implicit none

  real(dp), parameter :: CmQsq = 3.4 
  real(dp), parameter :: Cml2sq = 2.0 

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
!    As in the module header.
!    
  use tensorFields 
  use fvxGradient   

  implicit none

!
! Local variables
!
  integer :: i, ib, iface, ijp, ijn, ijb, ijbt, iWall, iPer
  real(dp) :: urf, rmin, rmax
  real(dp), parameter :: r23 = 2./3._dp
  type(volScalarField) :: Sij_, l2p32
  type(volTensorField) :: D
  type(volVectorField) :: U_

  U_ = volVectorField("Velocity", u, v, w )

  ! Velocity gradient tensor.
  D = Grad( U_ )

  ! Magnitude of the strain rate tensor (a symmetric part of velocity gradient tensor)
  Sij_ = .mag.(.symm.( D ) )  
  
  ! Set value for vortex identification criteria
  call setLambda2 
  l2p32= new_volScalarField(numCells)
  l2p32%mag(1:numCells)  = max(lambda2**1.5, small)

  ! v---this should be 'musgs' but we reuse existing field
  Sij_ = den * Cml2sq * Vol**r23  * l2p32 * Sij_ &
                           /  ( power( Sij_, 3.0_dp ) + l2p32  ) 

  urf = TurbModel%urfVis

  ! Update of the effective viscosity
                          ! v---this should be 'musgs' but we reuse existing field
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

  write(*,*) " Updated effective viscosity - Vortex-ID SGS model."

  rmin = minval(vis/viscos)
  rmax = maxval(vis/viscos)
 
  write(*,'(2x,es11.4,a,es11.4)') rmin,' <= Viscosity ratio <= ',rmax

end subroutine


end module vortexID_sgs