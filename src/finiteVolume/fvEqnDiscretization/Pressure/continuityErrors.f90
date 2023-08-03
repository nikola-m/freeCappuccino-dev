subroutine continuityErrors
!
! Purpose: 
!   Calculates and prints the continuity errors.
!
  use types
  use parameters, only: resor
  use geometry, only: numInnerFaces, owner, neighbour, numBoundaries, bctype, nfaces, startFace
  use sparse_matrix, only: su,apu
  use variables, only: flmass
  use fieldManipulation, only: volumeWeightedAverage

  implicit none

  integer :: i, ijp, ijn, ib, iface
  real(dp) :: sumLocalContErr, globalContErr
  real(dp), save :: cumulativeContErr = 0.0_dp

  ! Initialize array with zero value.
  su = 0.0_dp

  ! Inner faces                                               
  do i=1,numInnerFaces                                                    
    ijp = owner(i)
    ijn = neighbour(i)
    su(ijp) = su(ijp)-flmass(ijp)
    su(ijn) = su(ijn)+flmass(ijp)
  enddo  


  ! Boundary faces
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        su(ijp)=su(ijp)-flmass(iface)

      enddo

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        su(ijp)=su(ijp)-flmass(iface)

      enddo

    endif

  enddo


  sumLocalContErr = volumeWeightedAverage( abs(su*apu) )

  ! Global mass conservation
  globalContErr = sum( su )

  cumulativeContErr = cumulativeContErr + globalContErr

  su = 0.0_dp

  ! For report of scaled residuals
  resor(4) = sumLocalContErr

  write(6,'(3(a,es10.3))') "  continuity errors : avg local = ", sumLocalContErr, &
 &                          ", global = ", globalContErr, &
 &                          ", cumulative = ", cumulativeContErr

end subroutine

