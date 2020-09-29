!***********************************************************************
! 
subroutine bpres(p,istage)
!
!***********************************************************************
! 
! Purpose:
!   Set pressures at boundaries.
!
! Discussion:
!   istage = 1 : At boundary faces, we set the values of owner cell,
!                we need those boundary values to be able to calculate
!                pressure or pressure correction gradients.
!   istage = 2 : We perform linear extrapolation from owner cell using 
!                previosuly calculated cell centered gradients.
!   istage = 3,etc. Same as istage=2.
!                Question is do we need more than 2 passes? 
!                Numerical experiments are needed.
!
!***********************************************************************
! 
  use types
  use parameters
  use variables, only: dPdxi
  use geometry
  implicit none
!
!***********************************************************************
! 
  real(dp), dimension(numTotal) :: p
  integer, intent(in) :: istage

  ! Locals:
  integer :: i, ijp, ijb, iface
  integer :: ib
  real(dp) :: xpb, ypb, zpb

  if ( istage.eq.1 ) then

    ! Loop Boundary faces:

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i

      ! Takes owner cell value
      p(ijb) = p(ijp)

    enddo

  else ! istage==2 and higher

    ! ! Loop Boundary faces:

    ! Wall faces
    do ib=1,numBoundaries
      
      if ( bctype(ib) == 'wall' ) then

        do i=1,nfaces(ib)

          iface = startFace(ib) + i
          ijp = owner(iface)
          ijb = iBndValueStart(ib) + i

          ! Distance vector
          xpb = xf(iface)-xc(ijp) 
          ypb = yf(iface)-yc(ijp)
          zpb = zf(iface)-zc(ijp)

          ! Linear extrapolation
          p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb

        end do

      endif 

    enddo


  endif ! istage

end subroutine
