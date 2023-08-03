!***********************************************************************
!
subroutine calcstress
!***********************************************************************
!
! Calculates turbululent stresses,
!
! \tau_{ij} = 2*\mu_t*S_{ij} - 2/3*\rho*k*\delta_{ij}
!
! For Explicit Algebraic Reynolds Stress Models-EARSM we add anisotropy
!
! \tau_{ij} = 2 \mu_t S_{ij} - 2/3 \rho k \delta_{ij} + \rho k a_{ij}
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numCells
  use variables

  implicit none 
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: mut                    ! turbulent viscosity
  real(dp) :: dudx, dudy, dudz, &     ! dudxi - the velocity gradient
              dvdx, dvdy, dvdz, &     ! dvdxi - the velocity gradient
              dwdx, dwdy, dwdz        ! dwdxi - the velocity gradient

  do inp=1,numCells

    mut=(vis(inp)-viscos)

    dudx = dUdxi(1,inp)
    dudy = dUdxi(2,inp)
    dudz = dUdxi(3,inp)
    
    dvdx = dVdxi(1,inp)
    dvdy = dVdxi(2,inp)
    dvdz = dVdxi(3,inp)

    dwdx = dWdxi(1,inp)
    dwdy = dWdxi(2,inp)
    dwdz = dWdxi(3,inp)

  
    ! Reynolds stress tensor components
    uu(inp) = mut*(dudx+dudx) - twothirds*den(inp)*te(inp)
    vv(inp) = mut*(dvdy+dvdy) - twothirds*den(inp)*te(inp)
    ww(inp) = mut*(dwdz+dwdz) - twothirds*den(inp)*te(inp)

    uv(inp) = mut*(dudy+dvdx)
    uw(inp) = mut*(dudz+dwdx)
    vw(inp) = mut*(dvdz+dwdy)


    ! EARSM constitutive equation - add anisotropy tensor
    if (lasm) then

      uu(inp) = uu(inp) + den(inp)*te(inp)*aij(1,inp)
      vv(inp) = vv(inp) + den(inp)*te(inp)*aij(4,inp)
      ww(inp) = ww(inp) + den(inp)*te(inp)*aij(6,inp)

      uv(inp) = uv(inp) + den(inp)*te(inp)*aij(2,inp)
      uw(inp) = uw(inp) + den(inp)*te(inp)*aij(3,inp)
      vw(inp) = vw(inp) + den(inp)*te(inp)*aij(5,inp)

    endif

    ! Clip negative values
    uu(inp)=max(uu(inp),small)
    vv(inp)=max(vv(inp),small)
    ww(inp)=max(ww(inp),small)

  end do 


end subroutine
