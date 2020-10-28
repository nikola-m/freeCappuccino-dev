  program immdecomp
  !
  ! Purpose:
  ! An application implementing domain decomposition using Inertial Multi-Section Method (IMM)
  !
  ! References:
  ! Details have been taken from PhD thesis of Bojan Ničeno at TU Delft (2001) 
  ! (search online repository TU Delft library)
  !
  !+-----------------------------------------------------------------------------+

  use types
  use utils
  use geometry
  use mymodule
  use matrix_module, only: eig

  implicit none

  integer :: i
  integer :: LWORK, INFO
  real(dp) :: A( 3,3 ), VL( 3,3 ), VR( 3,3 ), WR( 3 ), WI( 3 ), WORK( 1000 )

  real(dp):: xm,ym,zm
  real(dp):: Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz
  ! real(dp), dimension(3,3) :: A,lam,v

!+-----------------------------------------------------------------------------+


  call show_logo( )


  ! Read mesh files and calculate mesh geometrical quantities

  call read_mesh_openfoam

  write(*,'(a)') ' '
  write(*,'(a)') 'Performing domain decomposition using Inertial Multi-Section Method'
  write(*,'(a)') ' '

  ! Centre of mass of the domain

  xm = volumeWeightedAverage(xc)
  ym = volumeWeightedAverage(yc)
  zm = volumeWeightedAverage(zc)

  write(*,'(a)') ' Centre of the mass of the domain:'
  write(*,'(3(e13.6,1x))') xm,ym,zm
  write(*,'(a)') ' '

  ! Moment of inertia matrix
  Ixx = sum( (yc(1:numCells)-ym)**2 + (zc(1:numCells)-zm)**2 )
  Ixy = sum( (xc(1:numCells)-xm)*(yc(1:numCells)-ym)**2 )
  Ixz = sum( (xc(1:numCells)-xm)*(zc(1:numCells)-zm)**2 )
  Iyx = sum( (yc(1:numCells)-ym)*(xc(1:numCells)-xm)**2 )
  Iyy = sum( (xc(1:numCells)-xm)**2 + (zc(1:numCells)-zm)**2  )
  Iyz = sum( (yc(1:numCells)-ym)*(zc(1:numCells)-zm)**2 )
  Izx = sum( (zc(1:numCells)-zm)*(xc(1:numCells)-xm)**2 )
  Izy = sum( (zc(1:numCells)-zm)*(yc(1:numCells)-ym)**2 )
  Izz = sum( (yc(1:numCells)-ym)**2 + (xc(1:numCells)-xm)**2  )


  write(*,'(a)') ' Inertia matrix of the domain:'
  write(*,'(3(e13.6,1x))') Ixx,Ixy,Ixz
  write(*,'(3(e13.6,1x))') Iyx,Iyy,Iyz
  write(*,'(3(e13.6,1x))') Izx,Izy,Izz
  write(*,'(a)') ' '

  A = reshape([Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz],[3,3])

  ! call eig(A,3,lam,v)

  ! write(*,'(a)') ' Eigenvalues of the inertia matrix:'
  ! write(*,'(3(e13.6,1x))') lam(1,1),lam(2,2),lam(3,3)
  ! write(*,'(a)') ' '


  ! All this but with Lapack
  LWORK = -1
  CALL DGEEV( 'Vectors', 'Vectors', 3, A, 3, WR, WI, VL, 3, &
              VR, 3, WORK, LWORK, INFO )
  LWORK = MIN( 1000, INT( WORK( 1 ) ) )
  ! Solve eigenproblem:
  CALL DGEEV( 'Vectors', 'Vectors', 3, A, 3, WR, WI, VL, 3, &
              VR, 3, WORK, LWORK, INFO )

  write(*,'(a)') ' Eigenvalues of the inertia matrix from LAPACK:'
  write(*,'(3(e13.6,1x))') (WR(i), i=1,3)
  write(*,'(3(e13.6,1x))') (WI(i), i=1,3)
  write(*,'(a)') ' '

  write(*,'(3(e13.6,1x))') (VL(1,i), i=1,3)
  write(*,'(3(e13.6,1x))') (VL(2,i), i=1,3)
  write(*,'(3(e13.6,1x))') (VL(3,i), i=1,3)
  write(*,'(a)') ' '

  call say_goodbye( )

end program