!***********************************************************************
!
subroutine asm_stress_terms
!
!***********************************************************************
!  Extra momentum source for turbulence models whose Reynolds stress is NOT the
!  Boussinesq one (explicit algebraic stress models, full Reynolds stress models).
!
!  The momentum equation of the code diffuses momentum with the effective viscosity
!
!        div( (mu + mut) (grad U + grad U^T) )
!
!  (implicit part + explicit transpose part, see facefluxuvw). The turbulence model gives a stress
!
!        T_ij = -rho <u_i' u_j'>  =  2 mut S_ij - 2/3 rho k delta_ij - rho k a_ij
!
!  where 'aij' (stored in variables, 6 components: xx,xy,xz,yy,yz,zz) is the anisotropy in EXCESS of the
!  Boussinesq one - 'a^(ex)' in the Turbulence Modeling Resource notation (https://turbmodels.larc.nasa.gov/easmko.html):
!  "Including nonlinear turbulence models like EASM is not simply a matter of computing mut alone. One must also
!  insure that the turbulent stress terms are computed appropriately to include the additional nonlinear components
!  in the Navier-Stokes equations."
!
!  The isotropic part -2/3 rho k delta_ij is absorbed in the pressure, so the difference we have to add as an
!  explicit source on the right hand side is
!
!        E_ij = T_ij - 2 mut S_ij + 2/3 rho k delta_ij = - rho k a_ij ,     su_i += sum_f ( E_ij )_f S_j .
!
!  For a full Reynolds stress model the 'mut' is only an artificial viscosity (stabilises the coupling of
!  stresses and velocities), and the same formula holds with  a_ij = ( <u_i u_j>/k - 2/3 delta_ij ) + 2 mut S_ij/(rho k).
!
!  Face values: linear interpolation of cell values with facint (weight of the neighbour cell). At inlet, outlet
!  and pressure boundaries the owner cell value is used; walls and symmetry planes are left to the wall treatment.
!
!***********************************************************************
  use types
  use geometry, only: numInnerFaces,numBoundaries,bctype,nfaces,startFace,facint,arx,ary,arz,owner,neighbour
  use sparse_matrix, only: su,sv,sw
  use variables, only: den,te,aij

  implicit none

  integer :: i, ib, iface, ijp, ijn
  real(dp) :: lam
  real(dp) :: exx, exy, exz, eyy, eyz, ezz
  real(dp) :: fx, fy, fz

  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    lam = facint(i)

    exx = -( (1.0_dp-lam)*den(ijp)*te(ijp)*aij(1,ijp) + lam*den(ijn)*te(ijn)*aij(1,ijn) )
    exy = -( (1.0_dp-lam)*den(ijp)*te(ijp)*aij(2,ijp) + lam*den(ijn)*te(ijn)*aij(2,ijn) )
    exz = -( (1.0_dp-lam)*den(ijp)*te(ijp)*aij(3,ijp) + lam*den(ijn)*te(ijn)*aij(3,ijn) )
    eyy = -( (1.0_dp-lam)*den(ijp)*te(ijp)*aij(4,ijp) + lam*den(ijn)*te(ijn)*aij(4,ijn) )
    eyz = -( (1.0_dp-lam)*den(ijp)*te(ijp)*aij(5,ijp) + lam*den(ijn)*te(ijn)*aij(5,ijn) )
    ezz = -( (1.0_dp-lam)*den(ijp)*te(ijp)*aij(6,ijp) + lam*den(ijn)*te(ijn)*aij(6,ijn) )

    fx = exx*arx(i) + exy*ary(i) + exz*arz(i)
    fy = exy*arx(i) + eyy*ary(i) + eyz*arz(i)
    fz = exz*arx(i) + eyz*ary(i) + ezz*arz(i)

    su(ijp) = su(ijp) + fx
    sv(ijp) = sv(ijp) + fy
    sw(ijp) = sw(ijp) + fz

    su(ijn) = su(ijn) - fx
    sv(ijn) = sv(ijn) - fy
    sw(ijn) = sw(ijn) - fz
  enddo

  do ib=1,numBoundaries

    if ( bctype(ib) /= 'inlet' .and. bctype(ib) /= 'outlet' .and. bctype(ib) /= 'pressure' ) cycle

    do i=1,nfaces(ib)
      iface = startFace(ib) + i
      ijp = owner(iface)

      exx = -den(ijp)*te(ijp)*aij(1,ijp)
      exy = -den(ijp)*te(ijp)*aij(2,ijp)
      exz = -den(ijp)*te(ijp)*aij(3,ijp)
      eyy = -den(ijp)*te(ijp)*aij(4,ijp)
      eyz = -den(ijp)*te(ijp)*aij(5,ijp)
      ezz = -den(ijp)*te(ijp)*aij(6,ijp)

      su(ijp) = su(ijp) + exx*arx(iface) + exy*ary(iface) + exz*arz(iface)
      sv(ijp) = sv(ijp) + exy*arx(iface) + eyy*ary(iface) + eyz*arz(iface)
      sw(ijp) = sw(ijp) + exz*arx(iface) + eyz*ary(iface) + ezz*arz(iface)
    enddo

  enddo

end subroutine
