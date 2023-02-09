module sigmaSGS

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing implementation of the Sigma sub-grid scale model for Large-Eddy Simulation.
!
!  Discussion:
!
!    Reference:
!    Nicoud F., Baya Toda H., Cabrit O., Bose S. and Lee J. Using singular values to build a subgrid-scale model 
!    for large eddy simulations. Physics of Fluids, 23 (8), pp. 085106, 2011.
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

  real(dp), parameter :: Cssq = 1.8225_dp 


  private 

  public :: modify_viscosity_sigma_sgs

contains


subroutine modify_viscosity_sigma_sgs
!
!  Purpose: 
!
!    Main routine of the module - return mu_sgs for the Sigma SGS model.
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
!    27th Janury 2023
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
  integer :: i, inp, ib, iface, ijp, ijn, ijb, ijbt, iWall, iPer
  real(dp) :: urf, rmin, rmax
  real(dp) :: a1,a2,a3,sigma1,sigma2,sigma3,musgs
  type(volVectorField) :: U_
  type(volTensorField) :: D, G
  type(volScalarField) :: I1,I2,I3


  U_ = volVectorField("Velocity", u, v, w )

  ! Velocity gradient tensor.
  D = Grad( U_ )

  ! Gij = Dik*Dkj
  G = .trans.D*D
  
  !  Invariants
  I1 = .tr.G
  I2 = half*( power(.tr.G, 2.0_dp) - .tr.(G*G) ) !..or: I2=half*( (.tr.G*.tr.G) - .tr.(G*G) )
  I3 = .det.G


  urf = TurbModel%urfVis

  do inp = 1,numCells

    a1 = I1%mag(inp)**2/9. - I2%mag(inp)/3.
    a2 = I1%mag(inp)**3/27. - I1%mag(inp)*I2%mag(inp)/6. + I3%mag(inp)/2.
    a3 = 1./3.*acos(a2/a1**1.5)

    sigma1 = sqrt( I1%mag(inp)/3. + 2*sqrt(a1)*cos(a3) )
    sigma2 = sqrt( I1%mag(inp)/3. - 2*sqrt(a1)*cos(Pi/3 + a3) )
    sigma3 = sqrt( I1%mag(inp)/3. - 2*sqrt(a1)*cos(Pi/3 - a3) )

    sigma1 = max(sigma1,small)
    sigma2 = max(sigma2,small)
    sigma3 = max(sigma3,small)

    musgs = den(inp) * Cssq * Vol(inp)**r23  * ( sigma3*(sigma1-sigma2)*(sigma2-sigma3) ) / ( sigma1*sigma1 )

    vis( inp ) = urf*( musgs + viscos ) + (1.0-urf)*vis( inp )

  enddo


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

  write(*,*) " Updated effective viscosity - Sigma SGS model."

  rmin = minval(vis/viscos)
  rmax = maxval(vis/viscos)
 
  write(*,'(2x,es11.4,a,es11.4)') rmin,' <= Viscosity ratio <= ',rmax

end subroutine


end module sigmaSGS