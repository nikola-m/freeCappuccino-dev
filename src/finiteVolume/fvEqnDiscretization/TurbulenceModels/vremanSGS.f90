module vremanSGS

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of the Vreman (2004) sub-grid scale model for Large-Eddy Simulation.
!
!  Discussion:
!
!    Reference:
!    A.W. Vreman - An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications,
!    Physics of Fluids, Vol. 16, N0. 10, 2004.
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

  real(dp), parameter :: Cvsq = 0.0681_dp !(Cv = 2.5*Cs**2; Cs = 0.165 as Smagorinsky const.)


  private 

  public :: modify_viscosity_vreman_sgs

contains


subroutine modify_viscosity_vreman_sgs
!
!  Purpose: 
!
!    Main routine of the module - return mu_sgs for the Vreman SGS model.
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
!    4th February 2023
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
  real(dp), parameter :: r23 = 2./3.0_dp
  integer :: i, ib, iface, ijp, ijn, ijb, ijbt, iWall, iPer
  real(dp) :: urf, rmin, rmax
  type(volVectorField) :: U_
  type(volTensorField) :: D, G
  type(volScalarField) :: musgs


  U_ = volVectorField("Velocity", u, v, w )

  ! Velocity gradient tensor.
  D = Grad( U_ )

  ! Gij = Dik*Dkj
  G = .trans.D*D
  
  ! Second invariant of G in numerator, magnitude square of G in denominator, 
  ! than get square root of that.
  ! Here '/' is overloaded elementwise division of two scalar fields.
  ! and we have used overloaded sqrt operator, here applied to volScalarField object.
  musgs = sqrt(  half*( power(.tr.G, 2.0_dp) - .tr.(G*G) ) / .magSq.G )
  
  musgs%mag =  max(musgs%mag, small)

  musgs = den * Cvsq* Vol**r23 * musgs

  urf = TurbModel%urfVis


  ! Update of the effective viscosity
  vis( 1:numCells ) = urf*( musgs%mag(1:numCells) + viscos ) + (1.0-urf) * vis( 1:numCells )


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

  write(*,*) " Updated effective viscosity - Vreman SGS model."

  rmin = minval(vis/viscos)
  rmax = maxval(vis/viscos)
 
  write(*,'(2x,es11.4,a,es11.4)') rmin,' <= Viscosity ratio <= ',rmax

end subroutine


end module vremanSGS