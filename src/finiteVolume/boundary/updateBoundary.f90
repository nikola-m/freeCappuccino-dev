subroutine updateBoundary(phi)
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
  integer :: i,ijp,ijn,ijb,ib,iface,iPer,ijbt 


  !
  ! Update boundaries
  !

  iPer = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'outlet' .or. &
        bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' .or. &
        bctype(ib) == 'empty' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        phi(ijb) = phi(ijp)

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

        phi(ijb) = half*( phi(ijp)+phi(ijn) )

        ! Now find where is the twin in field array
        ijbt = numCells + ( startFaceTwin(iPer) - numInnerFaces ) + i
        
        ! Twin takes the same values
        phi(ijbt) = phi(ijb)

      enddo

    endif

  enddo


end subroutine