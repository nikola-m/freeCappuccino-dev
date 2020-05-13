!***********************************************************************
!
subroutine readfiles
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

  if( myid.eq.0 ) then
    write(6,'(a)') ' '    
    write(6,'(a)') ' **Reading simulation restart files.'
  endif

  ! NOTE: nproc_char <- this (=myid + 1) written as left aligned string.
  call i4_to_s_left ( myid, nproc_char )

  call get_unit ( restart_unit )

  open ( unit = restart_unit,file=trim(restart_file)//'-'//trim(nproc_char),form='unformatted') 

  rewind restart_unit

  read(restart_unit) itime,time
  if(const_mflux) read(restart_unit) gradpcmf
  read(restart_unit) flmass
  read(restart_unit) u
  read(restart_unit) v
  read(restart_unit) w
  read(restart_unit) p
  read(restart_unit) te
  read(restart_unit) ed
  read(restart_unit) t
  read(restart_unit) vis
  read(restart_unit) visw  
  ! read(restart_unit) uu
  ! read(restart_unit) vv
  ! read(restart_unit) ww
  ! read(restart_unit) uv
  ! read(restart_unit) uw
  ! read(restart_unit) vw
  read(restart_unit) uo
  read(restart_unit) vo
  read(restart_unit) wo
  read(restart_unit) to
  read(restart_unit) teo
  read(restart_unit) edo

  rewind restart_unit
  close (restart_unit)

  ! Initialize pressure correction with current pressure field.
  pp = p

  !------------------------------------------------
  ! [read statistics after first collection: ]
  !------------------------------------------------
  if (ltransient .and. lreadstat) then

    call get_unit ( statistics_file )

    open( unit=statistics_file, file='statistics'//'-'//trim(nproc_char), form='unformatted' )  
    rewind statistics_file
    read(statistics_file) n_sample,u_aver,v_aver,w_aver,uu_aver,vv_aver,ww_aver,uv_aver,uw_aver,vw_aver,te_aver
    close ( statistics_file )

  endif

end subroutine
