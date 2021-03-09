!***********************************************************************
!
subroutine write_restart_files
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use variables
  use statistics

  implicit none

  integer :: restart_unit, statistics_file

!
!***********************************************************************
!

!--------------------------------------------------------------
!    [ Writing restart file]
!--------------------------------------------------------------
  call get_unit ( restart_unit )

  open(unit=restart_unit,file=restart_file,form='unformatted')
  rewind restart_unit

    write(restart_unit) itime,time
      if(const_mflux) write(restart_unit) gradpcmf
    write(restart_unit) flmass
    write(restart_unit) u
    write(restart_unit) v
    write(restart_unit) w
    write(restart_unit) p
    write(restart_unit) te
    write(restart_unit) ed
    write(restart_unit) t
    write(restart_unit) vis
    write(restart_unit) visw    
    ! write(restart_unit) uu
    ! write(restart_unit) vv
    ! write(restart_unit) ww
    ! write(restart_unit) uv
    ! write(restart_unit) uw
    ! write(restart_unit) vw
    write(restart_unit) uo
    write(restart_unit) vo
    write(restart_unit) wo
    write(restart_unit) teo
    write(restart_unit) edo
    ! write(restart_unit) to
    ! write(restart_unit) varto
    ! write(restart_unit) con
    ! write(restart_unit) cono
    ! write(restart_unit) alph
    ! write(restart_unit) vart
    ! write(restart_unit) edd
    ! write(restart_unit) ret
    ! write(restart_unit) den
    ! write(restart_unit) utt
    ! write(restart_unit) vtt
    ! write(restart_unit) wtt

  rewind restart_unit
  close (restart_unit)

!--------------------------------------------------------------
!    [ Writing of the statistics restart file]
!--------------------------------------------------------------
  if (ltransient) then

    call get_unit ( statistics_file )

    open(unit=statistics_file,file='statistics')
    rewind statistics_file

    write(statistics_file,*) n_sample
    write(statistics_file,*) u_aver,v_aver,w_aver, &
                             uu_aver,vv_aver,ww_aver, &
                             uv_aver,uw_aver,vw_aver,te_aver, &
                             te_aver
    close (statistics_file)

  endif

  write(6,*)'=*=*= Simulation restart files have been written! =*=*='

end subroutine
