subroutine continuityErrors
!
! Purpose: 
!   Calculates and prints the continuity errors.
!
  use types
  use parameters, only: sumLocalContErr, globalContErr, cumulativeContErr, resor, ip
  use geometry, only: numInnerFaces, owner, neighbour, numBoundaries, bctype, nfaces, startFace
  use sparse_matrix, only: res,apu
  use variables, only: flmass
  use fieldManipulation, only: volumeWeightedAverage

  implicit none

  integer :: i, ijp, ijn, ib, iface

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


  sumLocalContErr = volumeWeightedAverage( abs(res*apu) )

  ! globalContErr   = volumeWeightedAverage( res )

  ! sumLocalContErr = sum( abs( res ) ) 
  
  globalContErr = sum( res )

  cumulativeContErr = cumulativeContErr + globalContErr

  res = 0.0_dp

  ! For report of scaled residuals - scaled residual for continuity is this.
  resor(ip) = sumLocalContErr

  write(6,'(3(a,es10.3))') "  time step continuity errors : avg local = ", sumLocalContErr, &
 &                          ", global = ", globalContErr, &
 &                          ", cumulative = ", cumulativeContErr

end subroutine

