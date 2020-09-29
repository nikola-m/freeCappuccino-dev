!***********************************************************************
!
subroutine calc_statistics
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use variables
  use statistics

  implicit none
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: n_1n, nr

  integer :: iwall,ib,i,iface,ijp,ijb
  real(dp) :: are,nxf,nyf,nzf,Vnp,Vtp,xtp,ytp,ztp,Ut2

  n_sample = n_sample+1

  n_1n = dble(n_sample-1)/ dble(n_sample)

  nr = 1.0_dp / dble(n_sample) 

  do inp = 1,numCells

    ! Velocity field
    u_aver(inp) = u_aver(inp) * n_1n  + u(inp) * nr 
    v_aver(inp) = v_aver(inp) * n_1n  + v(inp) * nr 
    w_aver(inp) = w_aver(inp) * n_1n  + w(inp) * nr  

    if (lcal(ien)) t_aver(inp) = t_aver(inp) * n_1n + t(inp) * nr 
    if(lcal(icon)) con_aver(inp) = con_aver(inp) * n_1n + con(inp) * nr 


    if(lturb) then

      ! Reynolds stress components
      uu_aver(inp) = uu_aver(inp) * n_1n + (u(inp)-u_aver(inp))**2 * nr 
      vv_aver(inp) = vv_aver(inp) * n_1n + (v(inp)-v_aver(inp))**2 * nr 
      ww_aver(inp) = ww_aver(inp) * n_1n + (w(inp)-w_aver(inp))**2 * nr 

      te_aver(inp) = 0.5*(uu_aver(inp)+vv_aver(inp)+ww_aver(inp))

      uv_aver(inp) = uv_aver(inp) * n_1n + ((u(inp)-u_aver(inp))*(v(inp)-v_aver(inp))) * nr 
      uw_aver(inp) = uw_aver(inp) * n_1n + ((u(inp)-u_aver(inp))*(w(inp)-w_aver(inp))) * nr 
      vw_aver(inp) = vw_aver(inp) * n_1n + ((v(inp)-v_aver(inp))*(w(inp)-w_aver(inp))) * nr 
 
    endif

!    Concentration and turbulent fluxes
!    concon_aver(inp) = concon_aver(inp)+(con(inp)-con_nsample)**2 
!    ucon_aver(inp) = ucon_aver(inp)+((u(inp)-u_aver(inp))*(con(inp)-con_nsample))
!    vcon_aver(inp) = vcon_aver(inp)+((v(inp)-v_aver(inp))*(con(inp)-con_nsample))
!    wcon_aver(inp) = wcon_aver(inp)+((w(inp)-w_aver(inp))*(con(inp)-con_nsample))
  
  end do    

  ! Recompute wall shear stress in necessary
  iWall = 0

  do ib=1,numBoundaries
    if ( bctype(ib) == 'wall') then
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

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
        Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)

        ! Tangent direction
        xtp = xtp/vtp
        ytp = ytp/vtp
        ztp = ztp/vtp

        ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
        Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

        Tau(iWall) = viscos*Ut2/dnw(iWall)

      enddo
    endif 
  enddo

  ! Time averaged Wall Shear Stress at wall bundaries
  wss_aver(1:nwal) = wss_aver(1:nwal) * n_1n  + tau(1:nwal) * nr




end subroutine
