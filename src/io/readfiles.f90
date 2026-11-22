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
  use statistics 

  implicit none  

  integer :: file_unit

!
!***********************************************************************
!

  call get_unit( file_unit )
  open(unit=file_unit,file=restart_file,form='unformatted')
  rewind file_unit

    read(file_unit) itime,time
      if(const_mflux) read(file_unit) gradpcmf
    read(file_unit) flmass
    read(file_unit) u
    read(file_unit) v
    read(file_unit) w
    read(file_unit) p
    read(file_unit) te
    read(file_unit) ed
    read(file_unit) t
    ! read(file_unit) con
    read(file_unit) vis
    read(file_unit) visw    
    ! read(file_unit) uu
    ! read(file_unit) vv
    ! read(file_unit) ww
    ! read(file_unit) uv
    ! read(file_unit) uw
    ! read(file_unit) vw
    read(file_unit) uo
    read(file_unit) vo
    read(file_unit) wo
    read(file_unit) teo
    read(file_unit) edo
    read(file_unit) to
    ! read(file_unit) varto
    ! read(file_unit) cono
    ! read(file_unit) vart
    ! read(file_unit) edd
    ! read(file_unit) ret
    ! read(file_unit) den
    ! read(file_unit) utt
    ! read(file_unit) vtt
    ! read(file_unit) wtt

  rewind file_unit
  
  close (file_unit)

  !------------------------------------------------
  ! [read statistics after first collection: ]
  !------------------------------------------------
  ! if (ltransient) then

  !   call get_unit ( file_unit )

  !   open(unit=file_unit,file='statistics')

  !   rewind file_unit

  !   read(file_unit,*) n_sample
  !   read(file_unit,*) u_aver,v_aver,w_aver, &
  !                           uu_aver,vv_aver,ww_aver, &
  !                           uv_aver,uw_aver,vw_aver, &
  !                           te_aver
  !   close ( file_unit )

  ! endif

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Continuing simulaton from saved state - fields have been succesfully read! '
  write ( *, '(a)' ) ' '

end subroutine
