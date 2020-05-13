!***********************************************************************
!
 subroutine writehistory
!
!***********************************************************************
!  Writes values at each monitoring point, at each time-step for 
!  transient simulation, to a specific monitor file.
!
!***********************************************************************
!
  use parameters, only: ltransient,time,mpoints,myid
  use variables, only: u,v,w,te,ed

  implicit none
!
!***********************************************************************
!
  integer :: i,inp,imon

  if(ltransient ) then

    read(89,*) mpoints

    do i=1,mpoints

      read(89,*) imon,inp

      write(91+imon,'(2x,1p7e14.5,2x)') time,u(inp),v(inp),w(inp),p(inp),te(inp),vis(inp)

    end do

    rewind 89

  end if

end subroutine
