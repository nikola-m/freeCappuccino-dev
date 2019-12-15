subroutine adjustMassFlow
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables

  implicit none

  integer :: i,ib,ijb,ijp,iface
  real(dp) :: flowo,fac


  ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')

  ! Loop over inlet boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)

        ! Minus sign is there to make fmi(i) positive since it enters the cell.
        ! Check out comments in bcin.f90
        su(ijp) = su(ijp) - flmass(iface)

      end do

    endif 

  enddo  


  ! Loop Outlet faces, first to get flowo, then again after calculating fac.

  ! Extrapolated velocity at outlet boundary, outlet mass fluxes
  flowo=0.0_dp

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

  call global_sum( flowo )

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

        ! fmo is positive because of same direction of velocity and surface normal vectors
        ! but the mass flow is going out of the cell, therefore minus sign.
        su(ijp) = su(ijp) - flmass(iface)
            
      end do

    endif 

  enddo


end subroutine