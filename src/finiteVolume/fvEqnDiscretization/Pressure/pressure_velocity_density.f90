module compressible
!
! Contains functions for compressible flow calculations using pressure based solver.
!  
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients
  use linear_solvers
  use fieldManipulation
  use faceflux_mass
  use nablap
  use velocity, only: updateVelocityAtBoundary
  implicit none

public :: pressure_velocity_density_coupling

contains


!***********************************************************************
!
subroutine pressure_velocity_density_coupling
!
!***********************************************************************
!
  use faceflux_mass

  implicit none

  integer :: i, k, inp, iface, ijp, ijn, istage
  real(dp) :: sum, ppref, capd, cand, capc, canc, fmcor
  real(prec) :: C_rho
  real(prec), parameter :: RVOZD = 287.058  ! J/(kg K)

!
!***********************************************************************
!  

  ! Reset soarse matrix coefficinet values and source
  a = 0.0_dp
  su = 0.0_dp
  
  ! Volume timestepping sources
  do inp = 1, numCells

    ! We consider only BDF2 timestepping here
    C_rho = 1./(RVOZD*T(inp))
    sp(inp) = C_rho*1.5_dp*vol(inp)/timestep
    su(inp) = ( 2*deno(inp) - 0.5_dp*denoo(inp) )*vol(inp)/timestep

  end do


 !**********************************************************************
 ! Contribution to fluxes and source as in incompressible case
 ! On the rhs we have -sum(m*_f). Mass flows at faces obtained
 ! using Rhie-Chow interpolation.
 !**********************************************************************

  ! Calculate gradients of tentative velocities, used for precise interpolation
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)


  ! > Assemble off diagonal entries of system matrix and find mass flux at faces using Rhie-Chow interpolation

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxmass(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), capd, cand, flmass(i))


    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, &
                     p, dpdxi, capc, canc, suadd )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_index(i)
    a(k) = cand + canc

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_index(i)
    a(k) = capd + capc

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - cand - canc

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - capd - capc

    ! > Sources:

    su(ijp) = su(ijp) - flmass(i) + suadd
    su(ijn) = su(ijn) + flmass(i) - suadd 

  end do


  if(.not.const_mflux) call adjustMassFlow


!*Multiple pressure corrections loop *******************************************************************
  do ipcorr=1,npcor

    ! Initialize pressure correction
    pp=0.0d0

    ! Solving pressure correction equation

    call bicgstab(pp,ip) 
 
    ! Using LIS solver
    ! write(maxno,'(i5)') nsw(ip)
    ! write(tol,'(es9.2)') sor(ip)
    ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! call solve_csr( numCells, nnz, ioffset, ja, a, su, pp )

   
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad(pp,dPdxi)

    end do

    ! If simulation uses least-squares gradinets call this to get conservative pressure correction gradients.
    if ( lstsq_qr .or. lstsq_dm .or. lstsq_qr ) call grad(pp,dPdxi,'gauss_corrected')

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

        flmass(iface) = flmass(iface) + a(k) * pp(ijn)

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_index(i)

        flmass(iface) = flmass(iface) + a(k) * ( -pp(ijp) )
  
    enddo

  
    ! Correct velocities and pressure
    do inp=1,numCells

        u(inp) = u(inp) - dPdxi(1,inp) * vol(inp)*apu(inp)
        v(inp) = v(inp) - dPdxi(2,inp) * vol(inp)*apv(inp)
        w(inp) = w(inp) - dPdxi(3,inp) * vol(inp)*apw(inp)

        p(inp) = p(inp) + urf(ip)*(pp(inp)-ppref)

    enddo   

    ! Explicit correction of boundary conditions 
    call correctBoundaryConditionsVelocity


    if(ipcorr.ne.npcor) then ! non-orthogonal corrections -> > >

    !                                    
    ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    !

      ! Clean RHS vector
      su = 0.0_dp

      do i=1,numInnerFaces                                                      
        ijp = owner(i)
        ijn = neighbour(i)

        call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)

        flmass(i) = flmass(i)+fmcor

        su(ijp) = su(ijp)-fmcor
        su(ijn) = su(ijn)+fmcor 

      enddo                                                              

    endif ! < < <- END: non-orthogonal corrections


  enddo ! END: Multiple pressure corrections loop 

  ! Write continuity error report:
  include 'continuityErrors.h'

return 
end



!***********************************************************************
!
subroutine facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                      fm, lambda, gam, FI, dFidxi, &
                      cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use variables
  use interpolation

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: fm
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: fii
  real(dp) :: fcfie,fcfii,ffic
  real(dp) :: fxp,fxn

!----------------------------------------------------------------------


  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda


  ! Coeff in convection-like term
  Crhof = 1./(RVOZD*T(inp)*den(inp))*fxp + 1./(RVOZD*T(ine)*den(ine))*fxe


  ! Convection fluxes - uds
 
  ! System matrix coefficients
  cap = min(fm*Crhof,zero)
  can = max(fm*Crhof,zero)


  ! Explicit higher order convection

  if( flmass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dFidxi)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dFidxi)
  endif

  fcfie = fm*fii

  ! Explicit first order convection
  fcfii = min(fm,zero)*fi(ijn)+max(fm,zero)*fi(ijp)
 

  ! Explicit part of fluxes - Deffered correction for convection = gama_blending*(high-low) 
  suadd = -gam*(fcfie-fcfii)

end subroutine




!***********************************************************************
!
subroutine facefluxsc_boundary(ijp, ijn, xf, yf, zf, arx, ary, arz, flmass, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use variables, only: vis

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numTotal), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: Cp,Ce
  real(dp) :: fm
  real(dp) :: fxp,fxn

!----------------------------------------------------------------------


  ! Face interpolation factor
  fxn=1.0_dp
  fxp=0.0_dp

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  cap = - max(fm,zero)
  can =   min(fm,zero)


end subroutine


end module