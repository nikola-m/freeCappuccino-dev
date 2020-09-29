module concentration
!
! Implementation of sclar transport equation for concentration.
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none

  ! Constants
  real(dp), parameter :: sigt = 0.9_dp  

  private 

  public :: calculate_concentration_field

contains


subroutine calculate_concentration_field()
!
! Main module routine to assemble and solve concentration field.
!
  use types
  use parameters
  use variables
  use gradients
  implicit none

  call calcsc(Con,dCondxi,icon) ! Assemble and solve concentration eq.

end subroutine



subroutine calcsc(Fi,dFidxi,ifi)
!
! Ansamble and solve transport eq. for a scalar field.
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use title_mod

  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi 

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iface, ipro
  real(dp) :: gam, prtr, apotime, urfrs, urfms
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: fimax,fimin


! Variable specific coefficients:
  gam = gds(ifi)      !< first to higher order convection scheme blending parameter gamma.
  prtr = 1.0d0/sigt   !< reciprocal value of Prandtl-Schmidt number

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
!=====================================
! VOLUME SOURCE TERMS 
!=====================================
  do inp=1,numCells
 
    !
    ! > UNSTEADY TERM
    !
    if( bdf .or. cn ) then
      apotime = den(inp)*vol(inp)/timestep
      su(inp) = su(inp) + apotime*cono(inp)
      sp(inp) = sp(inp) + apotime
    elseif( bdf2 ) then
      apotime=den(inp)*vol(inp)/timestep
      su(inp) = su(inp) + apotime*( 2*cono(inp) - 0.5_dp*conoo(inp) )
      sp(inp) = sp(inp) + 1.5_dp*apotime
    endif

  enddo



! Calculate terms integrated over surfaces

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, &
                     fi, dFidxi, prtr, cap, can, suadd )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! ! (icell,icell) main diagonal element
    ! k = diag(ijp)
    ! a(k) = a(k) - can
    sp(ijp) = sp(ijp) - can

    ! ! (jcell,jcell) main diagonal element
    ! k = diag(ijn)
    ! a(k) = a(k) - cap
    sp(ijn) = sp(ijn) - cap

    ! > Sources:

    su(ijp) = su(ijp) + suadd
    su(ijn) = su(ijn) - suadd 

  enddo



  !
  ! Boundary conditions
  !

  ipro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'process') then
    ! Faces on processor boundaries

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijn = iBndValueStart(ib) + i

        iPro = iPro + 1

        call facefluxsc( ijp, ijn, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), fpro(ipro), gam, &
                         fi, dfidxi, prtr, cap, can, suadd )

        ! > Off-diagonal elements:    
        apr(ipro) = can

        ! > Elements on main diagonal:
        sp(ijp) = sp(ijp) - can

        ! > Sources:
        su(ijp) = su(ijp) + suadd

      
      enddo

    elseif ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         Fi, dFidxi, prtr, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         FI, dFidxi, prtr, cap, can, suadd )

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

      end do

    elseif ( bctype(ib) == 'wall') then

      ! Isothermal wall boundaries

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i


      enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Adiabatic wall boundaries

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i


        fi(ijb)=fi(ijp)

      enddo


    endif
    
  enddo ! Boundary conditions loop  


  ! Modify coefficients for Crank-Nicolson
  if (cn) then

    a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

    do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      k = icell_jcell_csr_index(i)
      su(ijp) = su(ijp) - a(k)*cono(ijn)

      k = jcell_icell_csr_index(i)
      su(ijn) = su(ijn) - a(k)*cono(ijp)

    enddo

        ! Processor boundary
        
        iPro = 0

        do ib=1,numBoundaries

          if ( bctype(ib) == 'process') then
          ! Faces on processor boundaries

            do i=1,nfaces(ib)

              iface = startFace(ib) + i
              ijp = owner(iface)
              ijn = iBndValueStart(ib) + i

              iPro = iPro + 1

              ! This is relevant to previous loop over faces
              su(ijp) = su(ijp) - apr(ipro)*cono(ijn)

              ! This is relevant to next loop over cells
              su(ijp) = su(ijp) + apr(ipro)*cono(ijp)

            enddo

          endif
          
        enddo  

    do ijp=1,numCells
      apotime=den(ijp)*vol(ijp)/timestep
      off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) !- a(diag(ijp))
      su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*cono(ijp)
      sp(ijp) = sp(ijp)+apotime

    enddo

  endif

  ! Underrelaxation factors
  urfrs = urfr(ifi)
  urfms = urfm(ifi)

  ! Main diagonal term assembly:
  do inp = 1,numCells
    ! Main diagonal term assembly:
    ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
    ! we substract it from the sum, to eliminate it from the sum.
    off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) !- a(diag(ijp)) because = 0
    a(diag(inp)) = sp(inp) - off_diagonal_terms

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)

  enddo

  ! Solve linear system:
  call bicgstab(fi,ifi)

  !
  ! Update symmetry and outlet boundaries
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'outlet' .or. bctype(ib) == 'symmetry' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        fi(ijb)=fi(ijp)

      enddo

    endif

  enddo

! Report range of scalar values and clip if negative
  fimin = minval(fi(1:numCells))
  fimax = maxval(fi(1:numCells))

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

  call global_min(fimin)
  call global_max(fimax)

  if( myid .eq. 0 ) write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= ',chvar(ifi),' <= ',fimax

  ! MPI exchange
  call exchange(fi)

end subroutine calcsc

end module concentration