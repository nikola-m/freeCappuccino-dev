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
  ! use LIS_linear_solver_library


  implicit none
!
!***********************************************************************
!

  ! character(5) :: maxno
  ! character(10) :: tol
  integer :: i, k, inp, ib, iface, istage
  integer :: ijp, ijn  
  real(dp) :: ppref, cap, can, fmcor,suma


  a = 0.0_dp
  su = 0.0_dp


  if( const_mflux ) call constant_mass_flow_forcing

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
    
  if(ltest) then  
    suma = sum(su)
    write(6,'(19x,a,1pe10.3)') ' Initial sum  =',suma
  endif  

!*Multiple pressure corrections loop *******************************************************************
  do ipcorr=1,npcor

    ! Initialize pressure correction
    pp=0.0d0

    ! Solving pressure correction equation
    ! call dpcg(pp,ip)
    ! call iccg(pp,ip)
    call bicgstab(pp,ip) 
    ! call pmgmres_ilu ( numCells, nnz, ioffset, ja, a, diag, pp(1:numCells), ip, su, 30, 4, 1e-7, sor(ip) )
    
    !  write(maxno,'(i5)') nsw(ip)
    !  write(tol,'(es9.2)') sor(ip)
    !  write(options,'(a)') "-i cg -p ilu -ilu_fill 1 -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    !  call solve_csr( numCells, nnz, ioffset, ja, a, su, pp )
       
    ! SECOND STEP *** CORRECTOR STAGE
   
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad(pp,dPdxi)

    end do

    ! If simulation uses least-squares gradients call this to get conservative pressure correction gradients.
    if ( lstsq_qr .or. lstsq_dm .or. lstsq_qr ) call grad(pp,dPdxi,'gauss_corrected','nolimit')


    ! Reference pressure correction - p'
    ppref = pp(pRefCell)

    !
    ! Correct mass fluxes at inner cv-faces only (only inner flux)
    !

    ! Inner faces:
    do iface=1,numInnerFaces

        ijp = owner(iface)
        ijn = neighbour(iface)

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_index(iface)

        flmass(iface) = flmass(iface) + a(k) * (pp(ijn)-pp(ijp))
  
    enddo

    !
    ! Correct velocities and pressure
    !      
    do inp=1,numCells
        u(inp) = u(inp) - dPdxi(1,inp) * vol(inp)*apu(inp)
        v(inp) = v(inp) - dPdxi(2,inp) * vol(inp)*apv(inp)
        w(inp) = w(inp) - dPdxi(3,inp) * vol(inp)*apw(inp)
        p(inp) = p(inp) + urf(ip)*(pp(inp)-ppref)
    enddo   

    ! Explicit correction of boundary conditions 
    call updateVelocityAtBoundary

    !.......................................................................................................!
    if(ipcorr.ne.npcor) then      
    !                                    
    ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    !

      ! Clean RHS vector
      su = 0.0_dp

      do i=1,numInnerFaces                                                      
        ijp = owner(i)
        ijn = neighbour(i)

        call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)

        flmass(i) = flmass(i)-fmcor

        su(ijp) = su(ijp)-fmcor
        su(ijn) = su(ijn)+fmcor 

      enddo                                                              

    endif                                                             
    !.......................................................................................................!


!*END: Multiple pressure corrections loop *******************************************************************
  enddo

  ! Write continuity error report:
  call continuityErrors

end subroutine
