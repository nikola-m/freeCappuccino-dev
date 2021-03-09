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
  use variables!, only: den,vis,visw, dUdxi,dVdxi,dWdxi

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
  use tensor_fields   ! <--- Let us see these in action!
  ! use fvxGradient   ! <--

  implicit none

!
! Local variables
!

  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau

  integer :: i, ib, iface, ijp, ijb, iWall

  real(dp), parameter :: r13 = 1./3._dp

  type(volScalarField) :: magSqrSd

  type(volTensorField) :: D,Sd

  ! type(volVectorField) :: U_

  ! U_ = volVectorField("Velocity", u, v, w )

  ! ! Velocity gradient tensor.
  ! D = Grad( U_ )

  D = volTensorField("VelGrad", &
                     dUdxi(1,:), dUdxi(2,:), dUdxi(3,:), &
                     dVdxi(1,:), dVdxi(2,:), dVdxi(3,:), &
                     dWdxi(1,:), dWdxi(2,:), dWdxi(3,:) )

  ! The traceless symmetric part of the square of the velocity gradient tensor.
  Sd = .dev.(.symm.(D.o.D))

  ! Its magnitude squared 
  magSqrSd = .magSq.Sd
   

  ! v---this should be 'musgs' but we reuse existing field
  magSqrSd = den * ( Cw_* Vol**r13 )**2 * &
                             power( magSqrSd, 3./2._dp ) &
           /  ( power( .magSq.(.symm.D), 5./2._dp ) + power( magSqrSd, 5./4._dp )  ) 



  ! Update of the effective viscosity and underelaxation
  vis( 1:numCells ) = urf(ivis)*( magSqrSd%mag(1:numCells) + viscos ) &
               + (1.0-urf(ivis)) * vis( 1:numCells )


  ! Boundary faces 

  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' .or. bctype(ib) == 'symmetry' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1


        !
        ! Wall boundaries - update Visw and Ypl
        !

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Magnitude of a cell center velocity projected on boundary face normal
        Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

        ! Tangential velocity components 
        xtp = U(ijp)-Vnp*nxf
        ytp = V(ijp)-Vnp*nyf
        ztp = W(ijp)-Vnp*nzf

        ! Its magnitude
        Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)

        ! Tangent direction
        xtp = xtp/vtp
        ytp = ytp/vtp
        ztp = ztp/vtp

        ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        Tau(iWall) = viscos*Ut2/dnw(iWall)
        Utau = sqrt( Tau(iWall) / den(ijb) )
        ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos


        visw(iWall) = max(viscos,zero)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  ! MPI exchange
  call exchange(vis)

end subroutine


end module wale_sgs