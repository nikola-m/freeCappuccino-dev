!***********************************************************************
!
subroutine openfiles
!
!***********************************************************************
!
  use title_mod
  use utils, only: show_logo
  
  implicit none
!
!***********************************************************************
!

  ! Open simulation log file
  open(unit=6,file=monitor_file)
  rewind 6

  ! Print Cappuccino logo in the header of monitor file
  call show_logo

  ! Open folder with data for postprocessing in Paraview
  call execute_command_line("mkdir VTK")


end subroutine
