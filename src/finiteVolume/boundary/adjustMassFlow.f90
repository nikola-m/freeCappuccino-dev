subroutine adjustMassFlow
!
! Purpose:
!   Correct mass flow at outflow boundaries to satisfy global mass conservation.
!
! Description:
!   Loop over all outlet boundaries and add up mass flows to get total mass flow - flowo.
!   Compare it to total mass flow that enters the computational doman - flomas. 
!   Find ratio - fac.
!   Adjust the massflow to outlet proportionally, also the velocity components at outlet boundary.
!
  use types
  use parameters, only: flomas, small
  use geometry, only: numBoundaries, owner, bctype, nfaces, startFace, iBndValueStart, arx, ary, arz
  use variables, only: u,v,w,flmass,den

  implicit none

  integer :: i,ib,ijb,ijp,iface
  real(dp) :: flowo,fac


  ! Mass flow trough all inlet boundaries (prescribed in routine 'bcin') adds up to give 'flmass'.

  ! Loop Outlet faces, first to get flowo, then again after calculating fac.

  ! Extrapolated velocity at outlet boundary, outlet mass fluxes
  flowo = 0.0_dp

  ! Loop over outlet boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        u(ijb) = u(ijp)
        v(ijb) = v(ijp)
        w(ijb) = w(ijp)

        flmass(iface) = den(ijp)*( u(ijb)*arx(iface)+v(ijb)*ary(iface)+w(ijb)*arz(iface) )        

        flowo = flowo + flmass(iface)

      end do

    endif 

  enddo


  ! Correct mass flux to satisfy global mass conservation & add to source
  fac = flomas/(flowo+small)


  ! Loop over outlet boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        flmass(iface) = flmass(iface)*fac

        u(ijb) = u(ijb)*fac
        v(ijb) = v(ijb)*fac
        w(ijb) = w(ijb)*fac
            
      end do

    endif 

  enddo


end subroutine