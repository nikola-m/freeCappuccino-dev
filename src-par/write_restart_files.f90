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
  use title_mod
  use statistics
  use utils

  implicit none

  integer :: restart_unit, statistics_file
  character( len = 5) :: nproc_char 

!
!***********************************************************************
!

  ! NOTE: nproc_char <- this (=myid + 1) written as left aligned string.
  call i4_to_s_left ( myid, nproc_char )

  call get_unit ( restart_unit )

  open ( unit = restart_unit,file=trim(restart_file)//'-'//trim(nproc_char),form='unformatted')
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
  write(restart_unit) to
  write(restart_unit) teo
  write(restart_unit) edo

  rewind restart_unit
  close (restart_unit)

!--------------------------------------------------------------
!    [ Writing of the statistics restart file]
!--------------------------------------------------------------
  if (ltransient) then

    call get_unit ( statistics_file )

    open( unit=statistics_file, file='statistics'//'-'//trim(nproc_char), form='unformatted' )  
    rewind statistics_file
    write(statistics_file) n_sample,u_aver,v_aver,w_aver,uu_aver,vv_aver,ww_aver,uv_aver,uw_aver,vw_aver,te_aver
    close (statistics_file)

  endif

  if( myid.eq.0 ) write(6,*)' **Created simulation restart files.'

end subroutine
