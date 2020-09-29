!***********************************************************************
!
subroutine calc_statistics
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numCells
  use variables
  use statistics

  implicit none
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: n_1n, nr


  n_sample = n_sample+1

  n_1n = dble(n_sample-1)/ dble(n_sample)

  nr = 1.0_dp / dble(n_sample) 

  do inp = 1,numCells

    ! Velocity field
    u_aver(inp) = u_aver(inp) * n_1n  + u(inp) * nr 
    v_aver(inp) = v_aver(inp) * n_1n  + v(inp) * nr 
    w_aver(inp) = w_aver(inp) * n_1n  + w(inp) * nr  
    ! con_aver(inp) = con_aver(inp) * n_1n + con(inp) * nr 


    if(lturb) then

      ! Reynolds stress components
      uu_aver(inp) = uu_aver(inp)+(u(inp)-u_aver(inp))**2
      vv_aver(inp) = vv_aver(inp)+(v(inp)-v_aver(inp))**2
      ww_aver(inp) = ww_aver(inp)+(w(inp)-w_aver(inp))**2

      te_aver(inp) = (uu_aver(inp)+vv_aver(inp)+ww_aver(inp)) / dble(2*n_sample)

      uv_aver(inp) = uv_aver(inp) + ((u(inp)-u_aver(inp))*(v(inp)-v_aver(inp)))
      uw_aver(inp) = uw_aver(inp) + ((u(inp)-u_aver(inp))*(w(inp)-w_aver(inp)))
      vw_aver(inp) = vw_aver(inp) + ((v(inp)-v_aver(inp))*(w(inp)-w_aver(inp)))

    endif

!  Other...
!    concon_aver(inp) = concon_aver(inp)+ &
!                 (con(inp)-con_nsample)**2 
!    ucon_aver(inp) = ucon_aver(inp)+ &
!                 ((u(inp)-u_aver(inp))*(con(inp)-con_nsample))
!    vcon_aver(inp) = vcon_aver(inp)+ &
!                 ((v(inp)-v_aver(inp))*(con(inp)-con_nsample))
!    wcon_aver(inp) = wcon_aver(inp)+ &
!                 ((w(inp)-w_aver(inp))*(con(inp)-con_nsample))
  
  end do    

end subroutine
