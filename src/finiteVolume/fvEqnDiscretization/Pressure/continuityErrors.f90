subroutine continuityErrors
!
! Purpose: 
!   Calculates and prints the continuity errors.
!
  use types
  use parameters, only: sumLocalContErr, globalContErr, cumulativeContErr, resor
  use geometry, only: numInnerFaces, owner, neighbour, numBoundaries, bctype, nfaces, startFace
  use sparse_matrix, only: res,apu
  use variables, only: flmass
  use fieldManipulation, only: volumeWeightedAverage

  implicit none

  integer :: i, ijp, ijn, ib, iface
  ! real(dp) :: sumLocalContErrPrev

  ! Initialize array with zero value.
  res = 0.0_dp

  ! Inner faces                                               
  do i=1,numInnerFaces                                                    
    ijp = owner(i)
    ijn = neighbour(i)
    res(ijp) = res(ijp)-flmass(ijp)
    res(ijn) = res(ijn)+flmass(ijp)
  enddo  


  ! Boundary faces
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)-flmass(iface)

      enddo

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)-flmass(iface)

      enddo

    endif

  enddo


  ! Used this so far..
  sumLocalContErr = volumeWeightedAverage( abs(res*apu) )

  ! Variant 2 - experimental, scaling as in Fluent theory guide
  ! sumLocalContErrPrev = sumLocalContErr
  ! sumLocalContErr = sum( abs( res ) )
  ! ! Scale continuity residual as in Fluent theory guide
  ! if ( (itime-itimes).le.5 .and. sumLocalContErrPrev < sumLocalContErr ) then
  !   res5Mass = sumLocalContErr
  ! endif    
  ! ! Scale residual:
  ! sumLocalContErr = sumLocalContErr / res5Mass

  ! Global mass conservation
  globalContErr = sum( res )

  cumulativeContErr = cumulativeContErr + globalContErr

  res = 0.0_dp

  ! For report of scaled residuals
  resor(4) = sumLocalContErr

  write(6,'(3(a,es10.3))') "  continuity errors : avg local = ", sumLocalContErr, &
 &                          ", global = ", globalContErr, &
 &                          ", cumulative = ", cumulativeContErr

end subroutine

