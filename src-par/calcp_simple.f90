 !***********************************************************************
!
subroutine calcp_simple
!***********************************************************************
!
! Assemble and solve pressure correction equation in SIMPLE algorithm
! Correct mass fluxes, presure and velocity field
! Enables multiple pressure corrections for non-orthogonal meshes
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients
  use fieldManipulation
  use faceflux_mass
  use mpi

  implicit none

!
!***********************************************************************
!

  integer :: i, k, inp, iface, if, ijp, ijn, ib, ipro, istage
  real(dp) :: sum, suma, ppref, cap, can, fmcor
  ! real(dp) :: psum


  a = 0.0_dp
  apr = 0.0_dp
  su = 0.0_dp

  ! if( const_mflux ) call constant_mass_flow_forcing

  ! Tentative (!) velocity gradients used for velocity interpolation: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  ! > Assemble off diagonal entries of system matrix and find mass flux at faces using Rhie-Chow interpolation

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxmass(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i))

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources:

    su(ijp) = su(ijp) - flmass(i)
    su(ijn) = su(ijn) + flmass(i) 

  end do


    ! Contribution form boundaries

    iPro = 0

    do ib=1,numBoundaries

      if ( bctype(ib) == 'process' ) then

        ! Faces on process boundary

        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)
          ijn = iBndValueStart(ib) + i
          ipro = ipro+1

          call facefluxmass2(ijp, ijn, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), fpro(ipro), cap, can, flmass(if))

          ! > Off-diagonal elements:    
          apr(ipro) = can

          ! > Elements on main diagonal:

          ! (icell,icell) main diagonal element
          k = diag(ijp)
          a(k) = a(k) - can

          ! > Sources:

          su(ijp) = su(ijp) - flmass(if)   

        enddo

      endif 

    enddo ! Boundary loop


  if(.not.const_mflux) call adjustMassFlow

  ! Add contributions to source of the inlet and outlet boundaries.

  ! Loop over boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)

        ! Minus sign is there to make fmi(i) positive since it enters the cell.
        ! Check out comments in bcin.f90
        su(ijp) = su(ijp) - flmass(iface)

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)

        ! fmout is positive because of same direction of velocity and surface normal vectors
        ! but the mass flow is going out of the cell, therefore minus sign.
        su(ijp) = su(ijp) - flmass(iface)
            
      end do

    endif 

  enddo 

  ! Test continutity:
  if(ltest) then
    suma = sum(su)
    call global_sum(suma)
    if (myid .eq. 0) write(6,'(19x,a,1pe10.3)') ' Initial sum  =',suma
  endif

!=====Multiple pressure corrections=====================================
  do ipcorr=1,npcor

    ! Initialize pressure correction
    pp = 0.0_dp

    ! Solving pressure correction equation
    ! call bicgstab(pp,ip) 
    call iccg(pp,ip)
    ! call dpcg(pp,ip)
    ! call jacobi(pp,ip)
    ! call pmgmres_ilu ( numCells, nnz, ioffset, ja, a, diag, pp, ip, su, nsw(ip), 4, 1e-3, sor(ip) )
       
    ! SECOND STEP *** CORRECTOR STAGE
   
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad(pp,dPdxi)
  
    end do

    ! If simulation uses least-squares gradients call this to get conservative pressure correction gradients.
    if ( lstsq_qr .or. lstsq_dm .or. lstsq_qr ) call grad(pp,dPdxi,'gauss_corrected','no-limit')
    
    ! Reference pressure correction - p'
    if (myid .eq. iPrefProcess) then
      ppref = pp(pRefCell)
      call MPI_BCAST(ppref,1,MPI_DOUBLE_PRECISION,iPrefProcess,MPI_COMM_WORLD,IERR)
    endif

    ! psum = sum( pp(1:numCells) ) 
    ! call global_sum( psum )
    ! ppref = psum / dble(gloCells)

    !
    ! Correct mass fluxes at inner cv-faces only (only inner flux)
    !

    ! Inner faces:
    do iface=1,numInnerFaces

        ijp = owner(iface)
        ijn = neighbour(iface)

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_index(iface)

        flmass(iface) = flmass(iface) + a(k) * ( pp(ijn)-pp(ijp) )
  
    enddo


    ! Correct mass fluxes at processor boundaries
    iPro = 0

    do ib=1,numBoundaries

      if ( bctype(ib) == 'process' ) then

        ! Faces on process boundary

        do i=1,nfaces(ib)

          iface = startFace(ib) + i
          ijp = owner(iface)
          ijn = iBndValueStart(ib) + i
          ipro = ipro+1

          flmass(iface) = flmass(iface) + apr(ipro) * ( pp( ijn ) - pp( ijp ) )

        enddo

      endif 

    enddo ! Boundary loop


    !
    ! Correct velocities and pressure
    !      
    do inp=1,numCells

        u(inp) = u(inp) - dPdxi(1,inp) * vol(inp)*apu(inp)
        v(inp) = v(inp) - dPdxi(2,inp) * vol(inp)*apv(inp)
        w(inp) = w(inp) - dPdxi(3,inp) * vol(inp)*apw(inp)

        p(inp) = p(inp) + urf(ip)*(pp(inp)-ppref)

    enddo   

    call updateVelocityAtBoundary

    !.......................................................................................................!
    if(ipcorr.ne.npcor) then      
    !                                    
    ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    !

      ! Clean RHS vector
      su = 0.0d0

      do i=1,numInnerFaces                                                      
        ijp = owner(i)
        ijn = neighbour(i)

        call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)
        
        flmass(i) = flmass(i) + fmcor 

        su(ijp) = su(ijp) - fmcor
        su(ijn) = su(ijn) + fmcor   

      enddo                                                              


      ! Faces on processor boundary
      iPro = 0

      do ib=1,numBoundaries

        if ( bctype(ib) == 'process' ) then

          ! Faces on process boundary

          do i=1,nfaces(ib)

          iface = startFace(ib) + i
          ijp = owner(iface)
          ipro = ipro+1

          call fluxmc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fpro(ipro), fmcor)
          
          flmass(iface) = flmass(iface) + fmcor
          
          su(ijp) = su(ijp) - fmcor

          enddo

        endif 

      enddo ! Boundary loop                                                                                  

    endif                                                             
    !.......................................................................................................!


!=END: Multiple pressure corrections loop==============================
  enddo

  call exchange( u )
  call exchange( v )
  call exchange( w )
  call exchange( p )

  ! Write continuity error report:
  call continuityErrors
 
end subroutine
