module timescale

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module for pseudotransient simulation control.
!
!  Discussion:
!
!    Basic thing here is setting the timescale for pseudotransient timestepping.
!    The code written here is loosly modelled on CFX Theory guide.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!

  use types

  implicit none

  logical :: autotime = .false.
  character( len = 12 ) :: lenscale_option = 'conservative' ! 'conservative' or 'aggresive'

  private 

  public :: autotime, lenscale_option ! params
  public :: set_timescale

contains

subroutine set_timescale
!
!  Purpose: 
!
!    Set pseudotransient timestep - a fluid timescale estimate.
!
!  Discussion:
!
!    The routine produces a couple of different timescale estimates,
!    and chooses the minimal one.
!    Timescale is obtained as length scale over velocity scale.
!    Length scale is based on three different estimates, final one is chosen 
!    on two general criteria: 'conservative' or 'aggressive'.
! 
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    18th October 2022
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    ANSYS CFX Theory guide 12.1, 2009.
!
!  Parameters:
!
!    -
!
  use geometry, only: numCells,numBoundaries,iBndValueStart,bctype,nfaces,startFace,arx,ary,arz,x,y,z,vol
  use variables, only: u,v,w!,T,den
  use parameters, only: timestep,onethird!,lbuoy,gravx,gravy,gravz

  implicit none
!
! Local variables
!
integer :: i,ib,ijb,iface
real(dp) :: lenscale
real(dp) :: lenvol
real(dp) :: lenx,leny,lenz,lenext
real(dp) :: areabc,lenbc
real(dp) :: ubc,unew,uavg
real(dp) :: dtu
! real(dp) :: rhoavg,udp,presbcmax,presbcmin,csound,dtp,dtc
! real(dp) :: gacc,Tmax,Tmin,dtp,dtg,dtrot

!
! > Length scales
!

! > Lenth scale based on cube root of computational domain volume
lenvol = sum( Vol(1:numCells) )**onethird


! > Length scale based on spatial extent of the computation domain
lenx = maxval(x) - minval(x) 
leny = maxval(y) - minval(y) 
lenz = maxval(z) - minval(z) 
lenext = max(lenx, leny, lenz)

! > Length scale based on the spatial extent of the inlet area
areabc  = 0.0_dp
do ib=1,numBoundaries
  if ( bctype(ib) == 'inlet') then
    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      areabc = areabc + sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
    enddo
  endif
enddo
lenbc = sqrt(areabc)

! > Final lengthscale choice - a la CFX:
if(lenscale_option == 'conservative') then

  lenscale = min( lenvol,lenext )

elseif(lenscale_option == 'aggresive') then

  lenscale = max( lenvol,lenext )

else !(default = conservative)

 lenscale = min( lenvol,lenext )

endif

!
! > Velocity scales 
!

! > Max velocity at inlets
ubc = 0.0_dp
do ib=1,numBoundaries
  if ( bctype(ib) == 'inlet') then
    do i=1,nfaces(ib)
      ijb = iBndValueStart(ib) + i
      unew = sqrt(u(ijb)**2+v(ijb)**2+w(ijb)**2)
      ubc = max( ubc, unew )
    enddo
  endif
enddo

! Arithmetic average velocity in the domain interior 
uavg =  sum( sqrt(u(1:numCells)**2+v(1:numCells)**2+w(1:numCells)**2) ) / dble(numCells)

! ! > Velocity scale based on pressure difference over open boundaries like pressure inlets or outlets
! rhoavg = sum( den(1:numCells) ) / dble(numCells)
! presbcmax = maxval( p( numCells+1 : numTotal ) ) ! Boundary face values are stored in an array after numcells inner field values.
! presbcmin = minval( p( numCells+1 : numTotal ) )
! udp = sqrt( (presbcmax - presbcmin)  / rhoavg )

! ! > For compressible flow - hardcoded for ideal gasses - basically it is cell avg of (\partial rho / \partial p )^-1
! csound = sum( 1./( Rair*T(1:numCells) ) ) / dble(numCells)


!
! > Timescales
!
dtu = 0.3 * lenscale / max( ubc, uavg )

! ! > Timescale based on the pressure difference valocity scale.
! dtp = 0.3 * lenscale / udp

! ! > Buoyancy timescale.
! gacc = sqrt(gravx**2+gravy**2+gravz**2)
! Tmax = maxval( T(1:numCells) )
! Tmin = minval( T(1:numCells) )
! if ( lbuoy ) then
!   if(boussinesq) then
!     dtg = sqrt( lenscale / ( gacc*beta*(Tmax-Tmin) ) )
!   else
!     dtg = sqrt( lenscale / gacc )
!   endif
! endif

! ! > System rotation timescale.
! if (rotation) dtrot = 0.1_dp / OmegaRot

! ! > Compressible flow timescale.
! if ( compressible ) dtc = lenscale / max( uavg, ubc, udp, csound  )

! ! > Timescale estimate for pseudotransient simulation - overwrites global variable 'timestep'.
! timestep = min( dtu, dtp, dtg, dtrot, dtc )

! > Simplified:
timestep = dtu

write(*,'(a)') ' '
write(*,'(3a,e10.3)') '  Pseudotransient control (',trim(lenscale_option),'): dt = ', timestep

end subroutine

end module