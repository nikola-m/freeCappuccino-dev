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
  use clsvof, only: twoPhase, clsvof_write_restart

  implicit none

  integer :: file_unit

!
!***********************************************************************
!

!--------------------------------------------------------------
!    [ Writing restart file]
!--------------------------------------------------------------
  call get_unit ( file_unit )

  open(unit=file_unit,file=restart_file,form='unformatted')
  rewind file_unit

    write(file_unit) itime,time
      if(const_mflux) write(file_unit) gradpcmf
    write(file_unit) flmass
    write(file_unit) u
    write(file_unit) v
    write(file_unit) w
    write(file_unit) p
    write(file_unit) te
    write(file_unit) ed
    write(file_unit) t
    ! write(file_unit) con
    write(file_unit) vis
    write(file_unit) visw    
    ! write(file_unit) uu
    ! write(file_unit) vv
    ! write(file_unit) ww
    ! write(file_unit) uv
    ! write(file_unit) uw
    ! write(file_unit) vw
    write(file_unit) uo
    write(file_unit) vo
    write(file_unit) wo
    write(file_unit) teo
    write(file_unit) edo
    write(file_unit) to
    ! write(file_unit) varto
    ! write(file_unit) cono
    ! write(file_unit) alph
    ! write(file_unit) vart
    ! write(file_unit) edd
    ! write(file_unit) ret
    ! write(file_unit) den
    ! write(file_unit) utt
    ! write(file_unit) vtt
    ! write(file_unit) wtt

  rewind file_unit
  close (file_unit)

!--------------------------------------------------------------
!    [ Writing of the statistics restart file]
!--------------------------------------------------------------
  if (ltransient) then

    call get_unit ( file_unit )

    open(unit=file_unit,file='statistics')
    rewind file_unit

    write(file_unit,*) n_sample
    write(file_unit,*) u_aver,v_aver,w_aver, &
                             uu_aver,vv_aver,ww_aver, &
                             uv_aver,uw_aver,vw_aver,te_aver, &
                             te_aver
    close (file_unit)

  endif

  if( twoPhase ) call clsvof_write_restart

  write(6,'(a)') ' '
  write(6,'(a)') ' =*=*= Simulation restart files have been written! =*=*='

end subroutine
