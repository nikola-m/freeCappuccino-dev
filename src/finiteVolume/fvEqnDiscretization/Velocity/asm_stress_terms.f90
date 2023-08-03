!***********************************************************************
!
subroutine asm_stress_terms
!
!***********************************************************************
!  Calculates the additional terms for momentum equations
!  and ads them to SU, SV, SW rhs. vectors.
!
! "Including nonlinear turbulence models like EASM is not simply a matter
!  of computing \mu_t alone. One must also insure that the turbulent stress
!  terms \tau_{ij} are computed appropriately to include the additional 
!  nonlinear components in the Navier-Stokes equations."
!  Source: https://turbmodels.larc.nasa.gov/easmko.html
!
!  Ovo je dodatak koji uvodimo jer u fluksevima za momentum koristimo 
!  muEff*( gradU + (gradU)T ), a kada konstitutivna relacija
!  za Rejnoldsove stress tenzore nije EVM (Boussinesq-ova eddy viscosity hypotesis).
!  onda razliku moramo da dodamo kao source for the momentum eqn.
!
!
!***********************************************************************
  use types
  use parameters, only: viscos
  use geometry, only: numInnerFaces,facint,arx,ary,arz,owner,neighbour
  use sparse_matrix, only: su,sv,sw
  use variables, only: den,vis,uu,uv,uw,vv,vw,ww,dudxi,dvdxi,dwdxi

  implicit none
!
!***********************************************************************
!

!
! Local variables
!
  integer :: i, ijp, ijn
  real(dp) :: mutN, mutP
  real(dp) :: term1i, term2i,term3i
  real(dp) :: fxp, fxn

!-------------------------------------------------------
!     [U-VELOCITY COMPONENT: ]
!
!     TERM1 is \rho*tau_{1,1} + 2*dUdx*\mu_t
!     where: 
!       \mu_t = \mu_eff - \mu_molecular,
!       tau_{1,i} - first row of Reynolds stress tensor
!
!     TERM2 is \rho*tau_{1,2} + \mu_t*(dUdy+DVdx)
!
!     TERM3 is \rho*tau_{1,3} + \mu_t*(dUdz+DWdx)
!
!     Find divergence of (-rho*tau_ij-2*mut*S_ij)
!     that is the difference of the explicit Reynolds sress
!     and the one which we used implicitely in momentum diffussion.
!     Since we used Boussinesq assumption for momentum diffusion 
!     implicitely, we must add the differece to the RHS of momentum eqn.
!
!     Apply divergence usually by interpolating to face, dotting 
!     with face area normal vector and sum over faces.
!
!     Minus that belongs to Reynolds stresses term in momenutm
!     equation, - rho*tau_ij, goes in front of the bracket
!     and inside the '-' in front of eddy-viscosity terms turns into '+'.
!       _____________________________
!     -( \rho*tau_{1,j} + 2*mut*S_ij )_f*Sf_j; j=1,2,3.
!-------------------------------------------------------

do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-facint(ijp)

  mutN = vis(ijn)-viscos
  mutP = vis(ijp)-viscos

        ! Interpolate term1:
        term1i = ( den(ijn)*uu(ijn) + mutN*(dudxi(1,ijn)+dudxi(1,ijn)) )*fxn+ &
                 ( den(ijp)*uu(ijp) + mutP*(dudxi(1,ijp)+dudxi(1,ijp)) )*fxp

        ! Interpolate term2:
        term2i = ( den(ijn)*uv(ijn) + mutN*(dudxi(2,ijn)+dvdxi(1,ijn)) )*fxn+ &
                 ( den(ijp)*uv(ijp) + mutP*(dudxi(2,ijp)+dvdxi(1,ijp)) )*fxp

        ! Interpolate term3:
        term3i = ( den(ijn)*uw(ijn) + mutN*(dudxi(3,ijn)+dwdxi(1,ijn)) )*fxn+ &
                 ( den(ijp)*uw(ijp) + mutP*(dudxi(3,ijp)+dwdxi(1,ijp)) )*fxp

        su(ijp) = su(ijp) - ( term1i*arx(i) + term2i*ary(i) + term3i*arz(i) )
        su(ijn) = su(ijn) + ( term1i*arx(i) + term2i*ary(i) + term3i*arz(i) )
enddo


!...sada isto za V komponentu
do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-facint(ijp)

  mutN = vis(ijn)-viscos
  mutP = vis(ijp)-viscos

        ! Interpolate term1:
        term1i = ( den(ijn)*uv(ijn) + mutN*(dvdxi(1,ijn)+dudxi(2,ijn)) )*fxn+ &
                 ( den(ijp)*uv(ijp) + mutP*(dvdxi(1,ijp)+dudxi(2,ijp)) )*fxp

        ! Interpolate term2:
        term2i = ( den(ijn)*vv(ijn) + mutN*(dvdxi(2,ijn)+dvdxi(2,ijn)) )*fxn+ &
                 ( den(ijp)*vv(ijp) + mutP*(dvdxi(2,ijp)+dvdxi(2,ijp)) )*fxp

        ! Interpolate term3:
        term3i = ( den(ijn)*vw(ijn) + mutN*(dvdxi(3,ijn)+dwdxi(2,ijn)) )*fxn+ &
                 ( den(ijp)*vw(ijp) + mutP*(dudxi(3,ijp)+dwdxi(1,ijp)) )*fxp

        sv(ijp) = sv(ijp) - ( term1i*arx(i) + term2i*ary(i) + term3i*arz(i) )
        sv(ijn) = sv(ijn) + ( term1i*arx(i) + term2i*ary(i) + term3i*arz(i) )
enddo


!...sada isto za W komponentu
do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-facint(ijp)

  mutN = vis(ijn)-viscos
  mutP = vis(ijp)-viscos

        ! Interpolate term1:
        term1i = ( den(ijn)*uw(ijn) + mutN*(dwdxi(1,ijn)+dudxi(3,ijn)) )*fxn+ &
                 ( den(ijp)*uw(ijp) + mutP*(dwdxi(1,ijp)+dudxi(3,ijp)) )*fxp

        ! Interpolate term2:
        term2i = ( den(ijn)*vw(ijn) + mutN*(dwdxi(2,ijn)+dvdxi(3,ijn)) )*fxn+ &
                 ( den(ijp)*vw(ijp) + mutP*(dwdxi(2,ijp)+dvdxi(3,ijp)) )*fxp

        ! Interpolate term3:
        term3i = ( den(ijn)*ww(ijn) + mutN*(dwdxi(3,ijn)+dwdxi(3,ijn)) )*fxn+ &
                 ( den(ijp)*ww(ijp) + mutP*(dwdxi(3,ijp)+dwdxi(3,ijp)) )*fxp

        sw(ijp) = sw(ijp) - ( term1i*arx(i) + term2i*ary(i) + term3i*arz(i) )
        sw(ijn) = sw(ijn) + ( term1i*arx(i) + term2i*ary(i) + term3i*arz(i) )
enddo


end subroutine
