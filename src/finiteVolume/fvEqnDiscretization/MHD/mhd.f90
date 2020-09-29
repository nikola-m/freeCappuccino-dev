module mhd
!
! Module for discretisation of equations for Magneto-Hydrodynamic flows.
!
! CURIX,CURIY,CURIZ : Induced current components.
! CURCX,CURCY,CURCZ : Conductive current components.
! CURX,CURY,CURZ    : Total current components.
! FLORX,FLORY,FLORZ : Lorentz force components.
!
use types
use parameters
use geometry
use variables

implicit none


real(dp), dimension(:), allocatable :: BMAGX,BMAGY,BMAGZ ! Magnetic inductions vector field components.
real(dp), dimension(:), allocatable :: EPOT              ! Electric potential.
real(dp), dimension(:), allocatable :: CURIX,CURIY,CURIZ ! Induced current vector field components.
real(dp), dimension(:), allocatable :: FLORX,FLORY,FLORZ ! Lorentz force vector field components.

real(dp), dimension(:,:), allocatable :: dEpotdxi
 
! Parameters defined in input file.       
real(dp), parameter :: SIGMA = 1.0_dp
real(dp) :: BNULA
real(dp) :: AMHD,BMHD,DMHD
real(dp) :: LBX,LBY,LBZ

public

contains

!***********************************************************************
!
subroutine calculate_Lorentz_force
!
!***********************************************************************
!
! Calculate Lorentz force used as a body force for the momentum equation.
!
! CURIX,CURIY,CURIZ : Induced current components.
! CURCX,CURCY,CURCZ : Conductive current components.
! CURX,CURY,CURZ    : Total current components.
! FLORX,FLORY,FLORZ : Lorentz force components.
!
!***********************************************************************
  use types
  use geometry
  use variables

  implicit none
!
!***********************************************************************

!
! Local variables
!
  integer :: inp,i,ib,ijb
  real(dp) :: curx,cury,curz       
  real(dp) :: curcx,curcy,curcz 

  ! Cells loop
  do inp=1,numCells

    !                             _ _
    ! > Calculate induced current uxB
    curix(inp)=v(inp)*bmagz(inp)-w(inp)*bmagy(inp)
    curiy(inp)=w(inp)*bmagx(inp)-u(inp)*bmagz(inp)
    curiz(inp)=u(inp)*bmagy(inp)-v(inp)*bmagx(inp)

    ! > Calculate conductive current
    ! Set the conductive current equal to minus the gradient of the electric potential.
    curcx=-dEpotdxi(1,inp)
    curcy=-dEpotdxi(2,inp)
    curcz=-dEpotdxi(3,inp)

    ! > Calculate total current
    curx=curix(inp)+curcx
    cury=curiy(inp)+curcy
    curz=curiz(inp)+curcz

    !                           _ _
    ! > Calculate Lorentz force JxB
    florx(inp)=cury*bmagz(inp)-curz*bmagy(inp)
    flory(inp)=curz*bmagx(inp)-curx*bmagz(inp)
    florz(inp)=curx*bmagy(inp)-cury*bmagx(inp)

  enddo

  ! Update induced current at boundaries (useful in fvxDiv( u x B ) calculation.)
  do ib=1,numBoundaries
    do i=1,nfaces(ib)      
      ijb = iBndValueStart(ib) + i
      curix(ijb)=v(ijb)*bmagz(ijb)-w(ijb)*bmagy(ijb)
      curiy(ijb)=w(ijb)*bmagx(ijb)-u(ijb)*bmagz(ijb)
      curiz(ijb)=u(ijb)*bmagy(ijb)-v(ijb)*bmagx(ijb)
    end do
  enddo

  ! write(6,*) "  Induced current: ", sum(curix), sum(curiy), sum(curiz)
  ! write(6,*) "  Force Lorentz: ", sum(florx)/numcells

end subroutine


!***********************************************************************
!
subroutine calculate_electric_potential
!
!*********************************************************************** 
!
! Assemble and solve Poisson equation for Electric potential field.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use title_mod
  use fieldManipulation, only: explDiv

  implicit none

!
! Local variables
!
  ! integer :: istage
  real(dp) :: ppref
  ! real(dp) :: fimax,fimin
  ! real(dp) :: epavg

  ! Initialize matrix and source arrays
  a  = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

  ! Nadji eksplicitnu divergenciju (u x B)
  su( 1:numCells ) = explDiv( curix, curiy, curiz )

  !  Coefficient array for Laplacian
  sv = 1.0_dp     

  ! Negative Laplacian operator     
  call updateEpotAtBoundaries(1)  
  call laplacian(sv,Epot) 

  ! Solve system
  pp = 0.0_dp
  call iccg(pp,iep) 
 
  ! First way:
  ppref = pp(pRefCell)
  Epot(1:numCells) = (1.0_dp-urf(iep))*Epot(1:numCells) + urf(iep)*(pp(1:numCells)-ppref)
  ! Second way:
  ! epavg = sum(pp(1:numCells))/dble(numCells)
  ! pp(1:numCells) = pp(1:numCells) - epavg
  ! Epot(1:numCells) = (1.0_dp-urf(iep))*Epot(1:numCells) + urf(iep)*pp(1:numCells)

  ! Update Epot at boundaries and calculate gradient right away, we'll need it.
  call updateEpotAtBoundaries(1)
  call grad(Epot,dEpotdxi)

  ! Report range of scalar values and clip if negative
  ! fimin = minval(Epot(1:numCells))
  ! fimax = maxval(Epot(1:numCells))
  
  ! write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= ',chvar(iep),' <= ',fimax

end subroutine


!***********************************************************************
! 
subroutine updateEpotAtBoundaries( istage )
!
!***********************************************************************
! 
! Purpose:
!   Set Electric potential at boundaries.
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
  use geometry

  implicit none
!
!***********************************************************************
! 
  integer, intent(in) :: istage

  ! Locals:
  integer :: i, ijp, ijb, iface
  integer :: ib
  real(dp) :: xpb, ypb, zpb

  if ( istage.eq.1 ) then

    ! Loop Boundary faces:
    do ib=1,numBoundaries
      if ( bctype(ib) == 'wall' ) then
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          ijp = owner(iface)     
          ijb = iBndValueStart(ib) + i
          Epot(ijb) = Epot(ijp)
        end do
      endif
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
          Epot(ijb) = Epot(ijp) + dEpotdxi(1,ijp)*xpb+dEpotdxi(2,ijp)*ypb+dEpotdxi(3,ijp)*zpb

        end do

      endif 

    enddo

  endif ! istage

end subroutine

end module

