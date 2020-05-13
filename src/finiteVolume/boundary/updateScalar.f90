subroutine updateScalar(phi)
!  
!******************************************************************************
!
!     Updates values of a scalar field at boundaries 
! 
!******************************************************************************
!
  use types
  use parameters
  use geometry

  implicit none

  real(dp), dimension(numTotal), intent(inout) :: phi

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

      ! enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        ! iface = startFace(ib) + i
        ! ijp = owner(iface)
        ! ijb = iBndValueStart(ib) + i


      end do

    endif 

  enddo


end subroutine