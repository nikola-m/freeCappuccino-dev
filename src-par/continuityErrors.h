!---------------------------------------------------------------------------
! Calculates and prints the continuity errors.
!---------------------------------------------------------------------------


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

    elseif ( bctype(ib) == 'process' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)-flmass(iface)

      enddo

    endif

  enddo


  sumLocalContErr = sum( abs( res ) ) 
  
  globalContErr = sum( res )

  call global_sum( sumLocalContErr )

  call global_sum( globalContErr )

  cumulativeContErr = cumulativeContErr + globalContErr

  write(*,*) 'After comm:',sumLocalContErr,globalContErr
  if ( myid .eq. 0 ) write(6,'(3(a,es10.3))') "  time step continuity errors : sum local = ", sumLocalContErr, &
 &                                            ", global = ", globalContErr, &
 &                                            ", cumulative = ", cumulativeContErr

