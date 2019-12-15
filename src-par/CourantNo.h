!
! Calculate and output the mean and maximum Courant Numbers.
!

 if (ltransient) then
 CoNum = 0.0_dp
 meanCoNum = 0.0_dp

 res = 0.0_dp

  !
  ! Suface sum of magnitude (i.e. absolute value) of mass flux phi, over inner faces only (which includes o- c- faces)
  !

  ! Inner faces                                               
  do i=1,numInnerFaces                                                                                                               
    ijp = owner(i)
    ijn = neighbour(i)
    res(ijp) = res(ijp)+abs(flmass(i))
    res(ijn) = res(ijn)+abs(flmass(i))                                                                                                               
  enddo     

  ! Faces on processor boundaries
  !do i=1,npro
  !  ijp = owner( iProcFacesStart + i )
  !  res(ijp) = res(ijp)+abs(fmpro(i))
  !end do

  ! Boundary faces
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then
      ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)+ abs(flmass(iface))
      enddo

    elseif ( bctype(ib) == 'outlet' ) then
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)+abs(flmass(iface))
      enddo

    elseif ( bctype(ib) == 'process' ) then
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)+abs(flmass(iface))
      enddo

    endif

  enddo

  ! Accumulate by looping trough cells
  suma = 0.0_dp

  do inp=1,numCells

    CoNum = max( CoNum , res(inp)/Vol(inp) )

    meanCoNum = meanCoNum + res(inp)
    
    suma = suma + Vol(inp)

  enddo

  res = 0.0_dp

  ! AllReduce using sum, values of suma and MeanCoNum.
  call global_sum(MeanCoNum)
  call global_sum(suma)

  CoNum = 0.5*CoNum*timestep

  ! Find global maximum Courant number in the whole field.
  call global_max(CoNum)


  ! Now use it to calculate mean Courant number
  meanCoNum = 0.5*meanCoNum/suma*timestep

  !// If we keep the value of Courant Number fixed
  if( CoNumFix .and. itime.ne.itimes ) then
      dt = timestep
      timestep = CoNumFixValue/CoNum * timestep

      CoNum = CoNum * timestep/dt
      meanCoNum  = meanCoNum * timestep/dt
  endif

  time = time + timestep

if(myid .eq. 0) then
 write(6,*)
 write(6,'(a,i0,a,es10.3,a,es10.3)') "  Time step no. : ",ITIME," dt : ",timestep," Time = ",time
 write(6,*)

 write(6,'(2(a,es10.3))') "  Courant Number mean: ", meanCoNum," max: ", CoNum
endif


endif
!// ************************************************************************* //
