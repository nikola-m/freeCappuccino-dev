subroutine updateVelocityAtBoundary
!  
!******************************************************************************
!
!     Updates values at symmetry boundaries 
! 
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

!
!     Local variables
!
  integer :: i,ijp,ijb,ib,iface
  real(dp) :: Unmag,flowo

  ! Update velocity components along outlet boundaries
  ! and correct mass flux to satisfy global mass conservation

  flowo=0.0_dp

  do ib=1,numBoundaries

    if ( bctype(ib) == 'outlet' ) then

      ! do i=1,nfaces(ib)

      !   iface = startFace(ib) + i
      !   ijp = owner(iface)
      !   ijb = iBndValueStart(ib) + i

      !   U(ijb) = U(ijp)
      !   V(ijb) = V(ijp)
      !   W(ijb) = W(ijp)

        ! flmass(iface) = den(ijp)*( u(ijb)*arx(iface)+v(ijb)*ary(iface)+w(ijb)*arz(iface) )
        
        ! flowo = flowo + flmass(iface)

      ! enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        ! Project velocity vector to face normal direction:
        Unmag = u(ijp)*arx(iface)+v(ijp)*ary(iface)+w(ijp)*arz(iface)

        U(ijb) = U(ijp)-Unmag*arx(iface)
        V(ijb) = V(ijp)-Unmag*ary(iface)
        W(ijb) = W(ijp)-Unmag*arz(iface)

      end do

    endif 

  enddo

  ! ! Ratio of inflow and outflow mass flux
  ! fac = flomas/(flowo+small)

  ! do ib=1,numBoundaries

  !   if ( bctype(ib) == 'outlet' ) then

  !     do i=1,nfaces(ib)

  !       iface = startFace(ib) + i
  !       ijb = iBndValueStart(ib) + i

  !       flmass(iface) = flmass(iface)*fac

  !       u(ijb) = u(ijb)*fac
  !       v(ijb) = v(ijb)*fac
  !       w(ijb) = w(ijb)*fac

  !     enddo

  !   endif 

  ! enddo


end subroutine