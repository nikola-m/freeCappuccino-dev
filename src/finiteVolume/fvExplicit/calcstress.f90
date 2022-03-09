!***********************************************************************
!
subroutine calcstress
!***********************************************************************
!
! Calculates turbulent stresses
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

    mut=(vis(inp)-viscos) !/densit ! divide with density if tau is define with nu_t instead of mu_t

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


    ! ! Wallin-Johansson EARSM  - check this again!
    ! uu(inp) = mut*(dudx+dudx) - 2.*bij(1,inp)*te(inp) - twothirds*den(inp)*te(inp)
    ! vv(inp) = mut*(dvdy+dvdy) - 2.*bij(4,inp)*te(inp) - twothirds*den(inp)*te(inp)
    ! ! It seems that b(3,3)=-b(1,1)-b(2,2):
    ! ww(inp) = mut*(dwdz+dwdz) + 2.*(bij(1,inp)+bij(4,inp))*te(inp) - twothirds*den(inp)*te(inp)

    ! uv(inp) = mut*(dudy+dvdx) - 2.*bij(2,inp)*den(inp)*te(inp)
    ! uw(inp) = mut*(dudz+dwdx) - 2.*bij(3,inp)*den(inp)*te(inp)
    ! vw(inp) = mut*(dvdz+dwdy) - 2.*bij(5,inp)*den(inp)*te(inp)
 

    ! Clip negative values
    uu(inp)=max(uu(inp),small)
    vv(inp)=max(vv(inp),small)
    ww(inp)=max(ww(inp),small)

  end do 


end subroutine
