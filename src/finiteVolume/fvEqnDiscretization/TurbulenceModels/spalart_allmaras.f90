module Spalart_Allmaras
!
! Implementation of Spalart-Allmaras turbulence model.
! Reference:  
! * Spalart, P. R. and Allmaras, S. R., "A One-Equation Turbulence Model for Aerodynamic Flows",
!   Recherche Aerospatiale, No. 1, 1994, pp. 5-21.
! * https://turbmodels.larc.nasa.gov/spalart.html (for implementation)
! ***
! * NOTE: We use TKE array to store nu_tilda!
! ***
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none


  ! Turbulence model constants
  real(dp), parameter :: Cb1 = 0.1355_dp
  real(dp), parameter :: Cb2 = 0.622_dp
  real(dp), parameter :: sigma = 2.0_dp/3.0_dp
  real(dp), parameter :: Cw1 = 3.2390678167757287 ! = Cb1/cappa**2 + (1.0_dp+Cb2)/sigma
  real(dp), parameter :: Cw2 = 0.3_dp
  real(dp), parameter :: Cw3 = 2.0_dp
  real(dp), parameter :: Cv1 = 7.1_dp
  real(dp), parameter :: Ct3 = 1.2_dp
  real(dp), parameter :: Ct4 = 0.5_dp

  real(dp), parameter :: one_sixth = 1.0_dp/6.0_dp
  logical :: notf2 = .false. ! Set to .true. for SA model without f_t2 term. Also for DES and DDES.
  

  private 

  public :: correct_turbulence_spalart_allmaras
  public :: correct_turbulence_inlet_spalart_allmaras


contains


subroutine correct_turbulence_spalart_allmaras()
!
! Main module routine to solve turbulence model equations and update effective viscosity.
! ***
! * NOTE: We use TKE array to store nu_tilda here!
! ***
!
  use types
  use parameters
  use variables
  use gradients
  implicit none

  call calcsc(TE,dTEdxi,ite) ! Assemble and solve nu_tilda eq.
  call modify_mu_eff()

end subroutine



subroutine correct_turbulence_inlet_spalart_allmaras()
!
! Update effective viscosity at inlet
!
  implicit none

  call modify_mu_eff_inlet()

end subroutine



subroutine calcsc(Fi,dFidxi,ifi)
!
! Assemble and solve transport equation for scalar field.
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
  integer ::  i, k, inp, ijp, ijn, ijb, ib, iwall, iface
  real(dp) :: gam, prtr, apotime, urfrs, urfms, &
              genp, genn
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: fimax,fimin
  real(dp) :: nu_tilda,nu,xi,fv1,fv2,ft2,r,g,fw,strain_tilda
  real(dp) :: dnutdx,dnutdy,dnutdz


! Variable specific coefficients:
  gam=gds(ifi)

  prtr=1.0_dp

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
!

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================
  do inp=1,numCells

    nu_tilda = te(inp)   ! we store nu_tilda in tke field.
    nu = viscos/densit   ! kinematic viscosity
    xi = nu_tilda / nu   ! xi parameter
    fv1 = xi**3 / (xi**3 + Cv1**3) ! parameter
    fv2 = 1.0_dp - xi / (1.0_dp + xi*fv1) ! parameter
    strain_tilda = vorticity(inp) + (nu_tilda / (cappa*wallDistance(inp))**2) * fv2 ! Shat

    ! Limiting strain - this is very important, there are other ways besides clipping implemented here.
    ! Clipping implemented here originates from a note referencing private communication with P. Spalart on turbmodels.larc.nasa.gov website
    strain_tilda = max(strain_tilda, 0.3_dp*vorticity(inp)) 

    ft2 = Ct3 * exp(-Ct4 * xi**2)
    if (notf2) ft2 = 0.0_dp ! If we want to eliminate ft2 term, other method is setting ct3 = 0.

    gen(inp) = Cb1 * (1.0_dp - ft2) * strain_tilda * nu_tilda

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! > Add production term to the rhs:
    su(inp)=genp*vol(inp)  


    ! > Cross diffusion: 
    dnutdx=dTEdxi(1,inp)
    dnutdy=dTEdxi(2,inp)
    dnutdz=dTEdxi(3,inp)

    su(inp) = su(inp) + ( Cb2/sigma * (dnutdx**2+dnutdy**2+dnutdz**2) ) * vol(inp)


    ! > Add destruction term to the lhs:

    r = min(nu_tilda/(strain_tilda*(cappa*wallDistance(inp))**2 + small),10.0_dp) ! r parameter

    g = r + Cw2*(r**6-r) ! g parameter

    fw = g*((1.0_dp+Cw3**6)/(g**6+Cw3**6))**one_sixth ! fw parameter

    ! > estruction term
    sp(inp) = (Cw1*fw-(Cb1/cappa**2)*ft2)*(nu_tilda/wallDistance(inp))**2 * den(inp)*vol(inp)/(nu_tilda+small)

    ! > If production negative move it to the lhs:
    sp(inp) = sp(inp)-genn*vol(inp)/(nu_tilda+small)

    ! > UNSTEADY TERM
    if( bdf .or. cn ) then
      apotime = den(inp)*vol(inp)/timestep
      su(inp) = su(inp) + apotime*teo(inp)
      sp(inp) = sp(inp) + apotime
    elseif( bdf2 ) then
      apotime=den(inp)*vol(inp)/timestep
      su(inp) = su(inp) + apotime*( 2*teo(inp) - 0.5_dp*teoo(inp) )
      sp(inp) = sp(inp) + 1.5_dp*apotime
    endif

  ! End of volume source terms
  enddo



!
! CALCULATE TERMS INTEGRATED OVER FACES
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    nu_tilda = te(ijp)   ! we store nu_tilda in tke field.
    nu = viscos/densit   ! kinematic viscosity
    xi = nu_tilda / nu   ! xi parameter
    fv1 = xi**3 / (xi**3 + Cv1**3) ! parameter  
    ! expression which will give correct diff coef when inserted into routine below:  
    prtr = ( 3.0_dp + 3*xi -  2*den(ijp) )/( 2*xi*den(ijp)*fv1 + small) 

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
    k = diag(ijp)
    a(k) = a(k) - can

    ! ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources:

    su(ijp) = su(ijp) + suadd
    su(ijn) = su(ijn) - suadd 

  enddo

  !
  ! Boundary conditions
  !

  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        nu_tilda = te(ijp)   ! we store nu_tilda in tke field.
        nu = viscos/densit   ! kinematic viscosity
        xi = nu_tilda / nu   ! xi parameter
        fv1 = xi**3 / (xi**3 + Cv1**3) ! parameter  
        ! expression which will give correct diff coef when inserted into routine below:  
        prtr = ( 3.0_dp + 3*xi -  2*den(ijp) )/( 2*xi*den(ijp)*fv1 + small ) 

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

        nu_tilda = te(ijp)   ! we store nu_tilda in tke field.
        nu = viscos/densit   ! kinematic viscosity
        xi = nu_tilda / nu   ! xi parameter
        fv1 = xi**3 / (xi**3 + Cv1**3) ! parameter  
        ! expression which will give correct diff coef when inserted into routine below:  
        prtr = ( 3.0_dp + 3*xi -  2*den(ijp) )/( 2*xi*den(ijp)*fv1 + small ) 

        call facefluxsc( ijp, ijb, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         flmass(iface), &
                         FI, dFidxi, prtr, cap, can, suadd)

        Sp(ijp) = Sp(ijp)-can

        Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

      end do

    ! elseif ( bctype(ib) == 'wall') then

    !   do i=1,nfaces(ib)

    !     iface = startFace(ib) + i
    !     ijp = owner(iface)
    !     ijb = iBndValueStart(ib) + i
    !     iWall = iWall + 1

    !     nu_tilda = te(ijp)   ! we store nu_tilda in tke field.
    !     nu = viscos/densit   ! kinematic viscosity
    !     xi = nu_tilda / nu   ! xi parameter
    !     fv1 = xi**3 / (xi**3 + Cv1**3) ! parameter  
    !     ! expression which will give correct diff coef when inserted into routine below:  
    !     prtr = ( 3.0_dp + 3*xi -  2*den(ijp) )/( 2*xi*den(ijp)*fv1 + small ) 

    !     ! Condition for nu_tilda at wall
    !     te(ijb) = 0.0_dp

    !     call facefluxsc( ijp, ijb, &
    !                      xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
    !                      flmass(iface), &
    !                      FI, dFidxi, prtr, cap, can, suadd )

    !     Sp(ijp) = Sp(ijp)-can

    !     Su(ijp) = Su(ijp)-can*Fi(ijb) + suadd

    !   enddo

    endif
    
  enddo ! Boundary conditions loop  




  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_index(i)
            su(ijp) = su(ijp) - a(k)*teo(ijn)

            k = jcell_icell_csr_index(i)
            su(ijn) = su(ijn) - a(k)*teo(ijp)
        enddo

        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) - a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*teo(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo
        
  endif

  ! Underrelaxation factors
  urfrs=urfr(ifi)
  urfms=urfm(ifi)

  ! Main diagonal term assembly:
  do inp = 1,numCells

        ! Main diagonal term assembly:
        a(diag(inp)) = sp(inp) 
        do k = ioffset(inp),ioffset(inp+1)-1
          if (k.eq.diag(inp)) cycle
          a(diag(inp)) = a(diag(inp)) -  a(k)
        enddo

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
  fimin = minval(fi)
  fimax = maxval(fi)
  
  write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= ',chvar(ifi),' <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


subroutine modify_mu_eff()
!
! Update turbulent and effective viscosity.
!
  use types
  use parameters
  use geometry
  use variables
  implicit none

  integer :: i,ib,inp
  integer :: iface, iwall, ijp,ijb
  real(dp) :: visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau,viscw
  real(dp) :: nu_tilda,nu,xi,fv1

!
! Loop trough cells 
!

  do inp=1,numCells
        ! Store old value
        visold=vis(inp)

        nu_tilda = te(inp)
        nu = viscos/densit
        xi = nu_tilda / nu 
        fv1 = xi**3 / (xi**3 + Cv1**3)
        ! Update effective viscosity:
        vis(inp)=viscos + den(inp)*nu_tilda*fv1
        ! Underelaxation
        vis(inp)=urf(ivis)*vis(inp)+(1.0_dp-urf(ivis))*visold
  enddo

  !
  ! Boundary faces 
  !

  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'outlet' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      enddo

    elseif ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        !
        ! Wall boundaries - update Visw and Ypl
        !

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Magnitude of a cell center velocity projected on boundary face normal
        Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

        ! Tangential velocity components 
        xtp = U(ijp)-Vnp*nxf
        ytp = V(ijp)-Vnp*nyf
        ztp = W(ijp)-Vnp*nzf

        ! Its magnitude
        Vtp = xtp*xtp+ytp*ytp+ztp*ztp

        ! Tangent direction
        xtp = xtp/vtp
        ytp = ytp/vtp
        ztp = ztp/vtp

        ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        Tau(iWall) = viscos*Ut2/dnw(iWall)
        Utau = sqrt( Tau(iWall) / den(ijb) )
        ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

        ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
        ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
        ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

        viscw = zero

        if(ypl(iWall) > ctrans) then
          viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
        endif

        visw(iWall) = max(viscos,viscw)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo


end subroutine


subroutine modify_mu_eff_inlet()
!
! Update turbulent and effective viscosity at inlet.
!
  use types
  use parameters
  use geometry, only: numBoundaries,nfaces,iBndValueStart
  use variables

  implicit none

  integer :: i,ib,ijb
  real(dp) :: nu_tilda,nu,fv1,xi

  !
  ! Boundary faces - loop over inlet boundaries
  !
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then

      do i=1,nfaces(ib)

        ijb = iBndValueStart(ib) + i

        nu_tilda = te(ijb) ! nu_tilda stored in TKE array..
        nu = viscos/densit
        xi = nu_tilda / nu 
        fv1 = xi**3 / (xi**3 + Cv1**3)

        Vis(ijb) = viscos+den(ijb)*nu_tilda*fv1

      end do

    endif

  enddo    


end subroutine


end module Spalart_Allmaras