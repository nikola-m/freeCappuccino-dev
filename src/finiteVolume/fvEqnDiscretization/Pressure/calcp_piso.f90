!***********************************************************************
!
subroutine calcp_piso
!
!***********************************************************************
!
! This implementation of PISO algorithm follows descripton given in
! Ferziger, Peric - Computational Methods for Fluid Dynamics, 2nd ed.
! It uses PRESSURE instead of PRESSURE CORRECTION as a variable.
! The same approach is also used in OpenFOAM library.
! The algorithm is summarised in following steps, with referenced eqns. 
! from the book given in braces: 
!
!  1) Ansemble and solve momentum eq. with pressure field from previous 
!     outer iteration or timestep (Eq. 7.30)
!     Obtained velocity doesn't satisfy continuity eqn.
!
!  2) Find "m* with tilde" velocities at each cell center (Eq. 7.31
!     but without the last term, i.e. the pressure gradient term)
!     ~ m*    ~ m*    ~ m*
!     u       v       w
!      P       P       P
!     "These are though as velocities from which the contribution of the 
!     pressure gradient has been removed"-Ferziger,Peric
!
!
!  3) Ansemble pressure equation (Eq. 7.35).
!     RHS is divergence of 
!         ~ m*        ~ m*         ~ m*
!     rho*u   ;   rho*v    ;   rho*w
!          P           P            P
!     Which means we interpolate these terms to cell face center, and 
!     dot them (perform scalar product) with cell face normal, i.e.
!     the face area vector.
!     Note, we don't need Rhie-Chow interpolation here. 
!     LHS matrix elements are divergence of unknown pressure gradient 
!     multiplied by rho/Ap, where Ap is diagonal term from momentum eq.
!     Note, divergence of gradient is Laplace operator, we can discretize
!     it in a routine that encapsulated implicit FVM calculation of  
!     Laplacian operator with rho/Ap as a coefficient.
!
!  4) Solve the system with tight tolerance (abs error below 1e-6 I guess)
!     to get new pressure field in cell centers.
!     Following that calculate new pressure gradients in cell centers.
!
!  5) Correct "m* with tilde" velocities to get velocities that satisfy
!     continuity equation (Eq. 7.34). 
!      m    ~ m*                             ~ m*      ~ m*
!     u   = u     - 1 / Apu * ( dp /dx)   ;  v   = ...  w = ...
!      P     P                                P          P 
!
!
!***********************************************************************
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
!
!***********************************************************************
!
  integer :: i, k, ib, if, istage
  integer :: ijp, ijn, ijb
  integer :: iPer, l, iftwin
  ! integer :: inp
  real(dp) :: cap, can
  real(dp) :: pavg, fmcor

  ! Before entering the corection loop backup a_nb coefficient arrays:
  h = a  

  ! if( const_mflux ) call constant_mass_flow_forcing

  
  do icorr=1,ncorr !== PISO Corrector loop ========================================================

    ! This is taken from cfd-online forum post:
    !// From the last solution of velocity, extract the diag. term from the matrix and store the reciprocal
    !// note that the matrix coefficients are functions of U due to the non-linearity of convection.
    !            volScalarField rUA = 1.0/UEqn.A();
    !// take a Jacobi pass and update U.  See Hrv Jasak's thesis eqn. 3.137 and Henrik Rusche's thesis, eqn. 2.43
    !// UEqn.H is the right-hand side of the UEqn minus the product of (the off-diagonal terms and U).
    !// Note that since the pressure gradient is not included in the UEqn. above, 
    !// this gives us U without the pressure gradient.  Also note that UEqn.H() is a function of U.
    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Posle ovoga imamo novo H(u)/ap, H(v)/ap ,i H(w)/ap A.K.A. "HbyA" smesteno u U,V, i W. To je polje brzina 
    ! bez uticaja gradijenta pritiska!

    ! First the right hand side of the momentum equations
    su = rU
    sv = rV
    sw = rW

    ! Assemble H(U) = - sum_j {a_j*U_pj}, j - runs trough neighbour indices
    ! Assemble H(V) = - sum_j {a_j*V_pj}, j - runs trough neighbour indices  
    ! Assemble H(W) = - sum_j {a_j*W_pj}, j - runs trough neighbour indices

    do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      k = icell_jcell_csr_index(i)
      su(ijp) = su(ijp) - h(k)*u(ijn)
      sv(ijp) = sv(ijp) - h(k)*v(ijn)
      sw(ijp) = sw(ijp) - h(k)*w(ijn)

      k = jcell_icell_csr_index(i)
      su(ijn) = su(ijn) - h(k)*u(ijp)
      sv(ijn) = sv(ijn) - h(k)*v(ijp)    
      sw(ijn) = sw(ijn) - h(k)*w(ijp)

    enddo
 
    ! HbyA
    u(1:numCells) = apu*su
    v(1:numCells) = apv*sv
    w(1:numCells) = apw*sw


    ! ! Tentative (!) velocity gradients used for velocity interpolation: 
    ! call updateVelocityAtBoundary
    ! call grad_gauss(U,dUdxi)
    ! call grad_gauss(V,dVdxi)
    ! call grad_gauss(W,dWdxi) 


    ! Initialize coefficient array and source:
    a = 0.0_dp
    su = 0.0_dp 

    ! > Assemble off diagonal entries of system matrix and find mass flux,
    !   accumulate diagonal entries of sysem matrix, and rhs vector stored in su array.

    ! Internal faces:
    do i = 1,numInnerFaces

      ijp = owner(i)
      ijn = neighbour(i)

      call facefluxmass_piso( i, ijp, ijn, &
                              arx(i), ary(i), arz(i), facint(i), &
                              cap, can, flmass(i) )
                              !, flmasso(i), flmassoo(i),flmassooo(i) )

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

      su(ijp) = su(ijp)-flmass(i)
      su(ijn) = su(ijn)+flmass(i) 

    end do


    !// "adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
    !// problem where a solution for pressure exists." - Comment in OF pisoFOAM code.
    if(.not.const_mflux) call adjustMassFlow

    ! Add contributions to source of the inlet and outlet boundaries.
 
    iPer = 0
    l = numInnerFaces

    ! Loop over boundaries
    do ib=1,numBoundaries
      
      if ( bctype(ib) == 'inlet' ) then

        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)

          ! Minus sign is there to make fmi(i) positive since it enters the cell.
          ! Check out comments in bcin.f90
          su(ijp) = su(ijp) - flmass(if)

        end do

      elseif ( bctype(ib) == 'pressure' ) then

        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)
          ijb = iBndValueStart(ib) + i

          call facefluxmassPressBnd( ijp, ijb, &
                                     xf(if), yf(if), zf(if), &
                                     arx(if), ary(if), arz(if), &
                                     cap, flmass(if) )

          ! > Elements on main diagonal:
          k = diag(ijp)
          a(k) = a(k) - cap

          ! > Sources:
          su(ijp) = su(ijp) - flmass(if)

          ! ! > Pressure Boundary value of pressure correction: - Do we need this? We work with pressure.
          ! pp(ijb) = 0.0_dp

        end do

      elseif ( bctype(ib) == 'outlet' ) then

        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)

          ! fmout is positive because of same direction of velocity and surface normal vectors
          ! but the mass flow is going out of the cell, therefore minus sign.
          su(ijp) = su(ijp) - flmass(if)
              
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


    !~~~~ Non orthogonal corrections loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! Multiple non-orthogonal passes amount to Discretizing Laplacian which includes adding 
    ! non-orthogonal contribution to RHS and forming RHS from mass fluxes formed using tentative velocity.
    !
    DO ipcorr=1,npcor

      ! !  "If you have a pressure equations with boundaries that do not fix pressure level, you have to fix a reference pressure." H.Jasak cfd-online forum
      ! !// In incompressible flow, only relative pressure matters.  Unless there is a pressure BC present,
      ! !// one cell's pressure has to be set to produce a unique pressure solution
      ! !    pEqn.setReference(pRefCell, pRefValue);
      ! !//
      ! ! So:
      ! a( ia(pRefCell):ia(pRefCell+1)-1 ) = 0.0_dp
      ! a( diag(pRefCell) ) = 1.0_dp
      ! ! Set rhs value for pRefCell to be reference pressure
      ! su(pRefCell) = p(pRefCell)
      
      !
      ! Solve pressure equation system
      !
#ifdef LIS  
      call lis_spsolve( pp, su, 'p' )
#else
      call csrsolve( lSolverP, pp, su, resor(4), maxiterP, tolAbsP, tolRelP, 'p' ) 
#endif

      ! We have pure Neumann problem - take out the average of the field as the additive constant
      pavg = sum(pp(1:numCells))/dble(numCells)

      ! Under-relaxation
      p(1:numCells) = (1.0_dp-urfP)*p(1:numCells) + urfP*(pp(1:numCells)-pavg)

      ! Pressure gradient
      do istage=1,nipgrad

        ! Pressure at boundaries (for correct calculation of p gradient)
        call bpres(p,istage)

        ! Calculate pressure gradient.
        call grad_gauss(p,dPdxi)

      end do

      !                                                                                  
      ! Laplacian source term modification due to non-orthogonality.
      !            
      if(ipcorr.ne.npcor) then                                            


        ! Add nonorthogonal terms from laplacian discretization to RHS of pressure eqn.
        do i=1,numInnerFaces

          ijp = owner(i)
          ijn = neighbour(i)

          call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)

          su(ijp) = su(ijp)-fmcor
          su(ijn) = su(ijn)+fmcor 

        enddo                                                                                                        

      !                                                                                        
      ! We have hit the last iteration of nonorthogonality correction: 
      !                     
      else ! that is if(ipcorr.eq.npcor)   

        !
        ! Correct mass fluxes at inner cv-faces only (only inner flux)
        !

        ! Inner faces:
        do if=1,numInnerFaces

            ijp = owner(if)
            ijn = neighbour(if)

            ! (icell,jcell) matrix element:
            k = icell_jcell_csr_index(if)

            flmass(if) = flmass(if) + a(k) * (p(ijn)-p(ijp))
      
        enddo
                                                                                                      

        ! Write continuity error report:
        call continuityErrors 

      endif 


    !~~~~ END: Non orthogonal corrections loop ~~~~~~~~~~~~~~~~~~~~~~~~~
    enddo


    !// Add pressure gradient to interior velocity and BC's.  Note that this pressure is not just a small
    !// correction to a previous pressure, but is the entire pressure field.  Contrast this to the use of p'
    !// in Ferziger & Peric, Eqn. 7.37.
    !// NOTE: This is whole pressure, opposite to what is done in SIMPLE: p(inp)+urf(ip)*(pp(inp)-ppref) !

    !
    ! Correct velocities
    !    

    !+++

    ! do inp=1,numCells
    !   u(inp) = u(inp) - apu(inp)*dPdxi(1,inp)*vol(inp)
    !   v(inp) = v(inp) - apv(inp)*dPdxi(2,inp)*vol(inp)
    !   w(inp) = w(inp) - apw(inp)*dPdxi(3,inp)*vol(inp)
    ! enddo 

    !+++

    ! ...or like this ...
    ! Do it in consistent way to how pressure was treated in calcuvw

    ! Use this to get: Source(u_i) = sum_f {-p}_f * S_fi.
    call gradp_and_sources(p)

    ! Correction is: u = u* + sum_f {-p}_f * S_fi / ap.

    u(1:numCells) = u(1:numCells) + su*apu
    v(1:numCells) = v(1:numCells) + sv*apv
    w(1:numCells) = w(1:numCells) + sw*apw

    !+++

    ! 
    ! Correct mass fluxes (and velocity) at pressure boundaries
    !

    ! Correct mass fluxes at periodic and pressure faces:
    iPer = 0
    l = numInnerFaces
            
    do ib=1,numBoundaries

      if (  bctype(ib) == 'periodic' ) then

        iPer = iPer + 1 ! found one periodic boundary pair

        ! Loop trough periodic boundary faces
        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)

          iftwin = startFaceTwin(iPer) + i
          ijn = owner(iftwin)

          ! l is in interval [numInnerFaces+1, numInnerFaces+numPeriodic]
          l = l + 1

          ! (icell,jcell) matrix element - the same as (jcell,icel) 
          ! since pressure poisson equation is symmetric:
          k = icell_jcell_csr_index(l)

          ! Correct mass fluxes trough these faces 
          flmass(if)  = flmass(if)  + a(k) * (p(ijn)-p(ijp)) 
          flmass(iftwin) = flmass(if)

        end do 

      
      elseif ( bctype(ib) == 'pressure' ) then

        do i=1,nfaces(ib)

          if = startFace(ib) + i
          ijp = owner(if)
          ijb = iBndValueStart(ib) + i

          call facefluxmassCorrPressBnd( ijp, ijb, xf(if), yf(if), zf(if), arx(if), ary(if), arz(if), flmass(if) )

        end do

      endif

    enddo


    ! Explicit correction of boundary conditions 
    call updateVelocityAtBoundary

  
  enddo !== END: PISO Corrector loop ==============================================================
               
  ! Correct driving force for a constant mass flow rate simulation:
  if(const_mflux) call constant_mass_flow_forcing                                                                                                        
      
end subroutine
