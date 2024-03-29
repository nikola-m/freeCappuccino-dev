subroutine calcp_simple
!
!  Purpose: 
!
!    Contains discretization and solution of pressure correction equation 
!    as well as the correction step for velocities, pressure and mass fluxes 
!    as a part of the SIMPLE algorithm.
!
!  Discussion:
!
!    This rutine is called after momentum predictor step in SIMPLE iterative procedure 
!    to assemble and solve equation for pressure correction field
!    and to correct pressure, velocity and mass fluxes defined at faces.
!
!    The routine goes trough following steps:
!    Assemble and solve pressure correction equation for SIMPLE algorithm.
!    Correct mass fluxes, presure and velocity field.
!    Enables multiple pressure corrections for non-orthogonal meshes.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    Feb-2021: Capability for compressible flows using pressure-density-velocity coupling.
!
!  Author:
!
!    Nikola Mirkov
!    Email: largeddysimulation@gmail.com / nikolamirkov@yahoo.com
!
!  Reference:
!
!    Ferziger & Perić, 
!    Mirkov et al JCP 2015
!
!  Parameters:
!
!    None.
!

  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients
  use linear_solvers
  use faceflux_mass
  use nablap
  use velocity, only: updateVelocityAtBoundary
#ifdef LIS  
  use LIS_linear_solvers, only: lis_spsolve
#endif

  implicit none


  integer :: i, k, ib, iface, istage
  integer :: inp, ijp, ijn, ijb
  integer :: iPer, l, if, iftwin
  real(dp) :: ppref, cap, can, fmcor

  ! For compressible flows
  real(dp) :: Crho, suadd

  ! Clear matrix coefficient and default rhs vector arrays
  a = 0.0_dp
  su = 0.0_dp


  ! Tentative (!) velocity gradients used for velocity interpolation: 
  ! (comment it out if only CDS interpolation of Ui is used)
  ! call grad_gauss(U,dUdxi)
  ! call grad_gauss(V,dVdxi)
  ! call grad_gauss(W,dWdxi)

  ! > Assemble off diagonal entries of system matrix and find mass flux at faces using Rhie-Chow interpolation

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    ! The facefluxmass routine below may be used in cases of highly non-orthogonal grids.
    ! call facefluxmass( i, ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i) )
    ! One of these is the alternative - 2 is typical form from FVM papers, 3 is based on Fluent theory guide
    call facefluxmass2( i, ijp, ijn, arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i) )
    ! call facefluxmass3( i, ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i) )


    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = cap

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - cap

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources:
    su(ijp) = su(ijp) - flmass(i)
    su(ijn) = su(ijn) + flmass(i) 

  end do



  ! Add contributions to source of the boundaries.

  ! Scaled mass flow at outlet
  if(.not.const_mflux) call adjustMassFlow
 
  iPer = 0
  l = numInnerFaces

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


    elseif ( bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)
        ijb = iBndValueStart(ib) + i

        call facefluxmassPressBnd( ijp, ijb, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), cap, flmass(if) )

        ! > Elements on main diagonal:
        k = diag(ijp)
        a(k) = a(k) - cap

        ! > Sources:
        su(ijp) = su(ijp) - flmass(if)

        ! > Pressure Boundary value of pressure correction:
        pp(ijb) = 0.0_dp

      end do


    elseif (  bctype(ib) == 'periodic' ) then

      iPer = iPer + 1

      ! Faces trough periodic boundaries
      do i=1,nfaces(ib)

        if = startFace(ib) + i
        ijp = owner(if)

        iftwin = startFaceTwin(iPer) + i
        ijn = owner(iftwin)


        call facefluxmass2_periodic(i, ijp, ijn, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), half, cap, can, flmass(if))


        ! > Off-diagonal elements:

        ! l is in interval [numInnerFaces+1, numInnerFaces+numPeriodic]
        l = l + 1

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_index(l)
        a(k) = cap

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_index(l)
        a(k) = cap

        ! > Elements on main diagonal:

        ! (icell,icell) main diagonal element
        k = diag(ijp)
        a(k) = a(k) - cap

        ! (jcell,jcell) main diagonal element
        k = diag(ijn)
        a(k) = a(k) - cap

        ! > Sources:
        su(ijp) = su(ijp) - flmass(if)
        su(ijn) = su(ijn) + flmass(if) 


      end do 

    endif 

  enddo 


  !
  ! > Additional matrix and RHS terms for compressible flows.
  !

  if (compressible) then
 
    ! Volumetric sources time-stepping sources
    do inp = 1,numCells

      Crho = 1./(Rair*T(inp)+small)  ! rho = p/(RT) = Crho*p

      ! First order
      sp(inp) = sp(inp) + Crho*vol(inp)/timestep        ! <- goes to main diagonal
      su(inp) = su(inp) + deno(inp)*vol(inp)/timestep   ! <- rhs vector

      ! Second order
      ! sp(inp) = sp(inp) + Crho*1.5_dp*vol(inp)/timestep               
      ! su(inp) = ( 2*deno(inp) - 0.5_dp*denoo(inp) )*vol(inp)/timestep 

    end do

    ! Pressure correction gradient for high-resolution bounded interpolation in facefluxCompressible.
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad_gauss(pp,dPdxi)

    end do


    ! Internal faces:
    do i = 1,numInnerFaces

      ijp = owner(i)
      ijn = neighbour(i)

      call facefluxCompressible( ijp, ijn, xf(i), yf(i), zf(i), flmass(i), facint(i), pp, dpdxi, cap, can, suadd )

      ! > Off-diagonal elements:

      ! (icell,jcell) matrix element:
      k = icell_jcell_csr_index(i)
      ! a(k) already keeps coefficients of the pressure correction eqn. just add compressible part
      ! so that below when we perfom mass flux correction we have included also contribution due to
      ! compressibility.
      a(k) = a(k) + can  

      ! (jcell,icell) matrix element:
      k = jcell_icell_csr_index(i)

      a(k) = a(k) + cap

      ! > Elements on main diagonal:

      ! (icell,icell) main diagonal element
      k = diag(ijp)
      a(k) = a(k) - can

      ! (jcell,jcell) main diagonal element
      k = diag(ijn)
      a(k) = a(k) - cap

      ! > Sources:

      su(ijp) = su(ijp) + suadd
      su(ijn) = su(ijn) - suadd 

    end do

  endif

  ! > END: Additional matrix and RHS terms for compressible flows.

    
  if(ltest) write(6,'(19x,a,1pe10.3)') ' Initial sum  =',sum(su) 


  ! Multiple pressure corrections loop *******************************************************************
  do ipcorr=1,npcor

    ! Solving pressure correction equation
#ifdef LIS  
    call lis_spsolve( pp, su, 'p' )
#else
    call csrsolve( lSolverP, pp, su, resor(4), maxiterP, tolAbsP, tolRelP, 'p' ) 
#endif    
    !
    ! Correct mass fluxes
    !

    ! Inner faces:
    do iface=1,numInnerFaces

      ijp = owner(iface)
      ijn = neighbour(iface)

      ! (icell,jcell) matrix element:
      k = icell_jcell_csr_index(iface)

      flmass(iface) = flmass(iface) + a(k) * (pp(ijn)-pp(ijp)) 
  
    enddo


    ! Correct mass fluxes at periodic and pressure faces:
    iPer = 0
    l = numInnerFaces
            
    do ib=1,numBoundaries

      if (  bctype(ib) == 'periodic' ) then

        iPer = iPer + 1 ! found one periodic boundary pair

        ! Loop trough periodic boundary faces
        do i=1,nfaces(ib)

          iface = startFace(ib) + i
          ijp = owner(iface)

          iftwin = startFaceTwin(iPer) + i
          ijn = owner(iftwin)

          ! l is in interval [numInnerFaces+1, numInnerFaces+numPeriodic]
          l = l + 1

          ! (icell,jcell) matrix element - the same as (jcell,icel) 
          ! since pressure poisson equation is symmetric:
          k = icell_jcell_csr_index(l)

          ! Correct mass fluxes trough these faces mj = mj* + mj' (where j is the face index)
          flmass(iface)  = flmass(iface)  + a(k) * (pp(ijn)-pp(ijp)) 
          flmass(iftwin) = flmass(iface)

        end do 

      elseif ( bctype(ib) == 'pressure' ) then

        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)
          ijb = iBndValueStart(ib) + i

          ! Corrects mass fluxes at pressure boundary as well as velocity components
          call facefluxmassCorrPressBnd( ijp, ijb, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), flmass(if) )

        end do

      endif 

    enddo


    !
    ! Correct velocity components and pressure
    !   

    ! Reference pressure correction - p'
    ppref = pp(pRefCell)

    ! If there is a pressure boundary then set p'=0
    do ib=1,numBoundaries 
      if ( bctype(ib) == 'pressure' ) then 
        ppref = 0.0_dp
        exit
      endif
    enddo



    ! Use this to get: Source(u_i) = sum_f {-p'}_f * S_fi.
    call gradp_and_sources(pp)

    ! Correction is: u = u* + sum_f {-p'}_f * S_fi / ap.

    u(1:numCells) = u(1:numCells) + su*apu
    v(1:numCells) = v(1:numCells) + sv*apv
    w(1:numCells) = w(1:numCells) + sw*apw
    p(1:numCells) = p(1:numCells) + urfp*(pp(1:numCells)-ppref)


  
    ! > Recalculate density from equation of state:
    if ( compressible ) den = p / ( Rair*T + small ) 



    ! Explicit correction of velocity at symmetry boundaries
    call updateVelocityAtBoundary


    !........................................................................................................
    if(ipcorr.ne.npcor) then      
    !                                    
    ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    !

      ! Clean RHS vector
      su = 0.0_dp

      do i=1,numInnerFaces                                                      
        ijp = owner(i)
        ijn = neighbour(i)

        ! call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)
        call fluxmc2(ijp, ijn, arx(i), ary(i), arz(i), fmcor)

        flmass(i) = flmass(i)+fmcor

        su(ijp) = su(ijp)-fmcor
        su(ijn) = su(ijn)+fmcor 

      enddo                                                              

    endif                                                             
    !.........................................................................................................!


  ! END: Multiple pressure corrections loop *******************************************************************
  enddo


  ! Write continuity error report:
  call continuityErrors


  ! Correct driving force for a constant mass flow rate simulation:
  if(const_mflux) call constant_mass_flow_forcing

end subroutine

