module wallSurfaceFields
!
! Purpose:
!   Module contains functions for calculating flow quantities pertinent to wall boundaries.
! Description:
!
!  Author: Nikola Mirkov
!  Email: nikolamirkov@yahoo.com
!
!  Modified:
!    Apr 10, 2020.
!
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!
use types
use parameters
use geometry 
use variables
use gradients

implicit none


public 

contains


subroutine forcesAtWall(fshearx,fsheary,fshearz,fprx,fpry,fprz)
!
! Purpose:
!   Computes friction and pressure forces imposed on wall boundaries by fluid flow.
!

  implicit none

  ! Input
  real(dp), dimension(*), intent(inout) :: fshearx,fsheary,fshearz,fprx,fpry,fprz

  ! Local variables
  integer :: i,ib,iface,ijp,ijb,iWall
  real(dp) :: viss,are,nxf,nyf,nzf,vsol,upb,vpb,wpb


  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        viss=max(viscos,visw(iWall))

        ! Face area
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! Diffusion coef. 
        vsol = viss*srdw(iWall)

        ! Velocity difference vector components
        upb = u(ijp)-u(ijb)
        vpb = v(ijp)-v(ijb)
        wpb = w(ijp)-w(ijb)

        ! Shear forces at wall in x, y and z direction.
        fshearx(iWall) = vsol*( (u(ijb)-u(ijp))*(1.-nxf**2) + vpb*nyf*nxf                  + wpb*nzf*nxf )
        fsheary(iWall) = vsol*( upb*nxf*nyf                 + (v(ijb)-v(ijp))*(1.-nyf**2)  + wpb*nzf*nyf )
        fshearz(iWall) = vsol*( upb*nxf*nzf                 + vpb*nyf*nzf                  + (w(ijb)-w(ijp))*(1.-nzf**2) )

        ! Pressure forces ( NOTE: we assume that boundary face normals point outwards ,away from the fluid region)
        fprx(iWall) = p(ijb)*arx(iface)
        fpry(iWall) = p(ijb)*ary(iface)
        fprz(iWall) = p(ijb)*arz(iface)

      enddo

    endif 

  enddo

end subroutine