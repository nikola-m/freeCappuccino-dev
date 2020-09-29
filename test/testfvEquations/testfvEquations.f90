  program testFieldOperation
  use utils
  use geometry
  use tensorFields
  use gradients
  use fvxGradient
  use fvxInterpolation
  use fvxDivergence

  implicit none

  integer :: i

! Locals
  !integer :: ierr

  type(volVectorField) :: vec1, vec2
  type(volScalarField) :: phi,psi
  type(volVectorField) :: resVec
  type(volTensorField) :: T
  type(surfaceScalarField) :: sf


  ! Intro
!+-----------------------------------------------------------------------------+

  call show_logo( )


  ! Read mesh files and calculate mesh geometrical quantities
!+-----------------------------------------------------------------------------+

  call read_mesh

  ! Creation and manipulation fo tensor fields
!+-----------------------------------------------------------------------------+


!
! > Initialize vector fields
!

!vec1 = volVectorField("Vector1", xc**2+yc+zc, xc+yc**2+zc, xc+yc+zc**2 )

vec1 = volVectorField(             &
                      "Vector1",   &
                      xc**2+yc+zc, &
                      xc+yc**2+zc, &
                      xc+yc+zc**2  &
                     )

vec2 = volVectorField("Vector2", &
                      xc+yc+zc,  &
                      xc+yc+zc,  &
                      xc+yc+zc)

!
! Test dot product of two vector field
!
  phi = vec1*vec2

  write(*,'(a)') ' Dot product of two vectors:'
  write(*,'(e13.6)') phi%mag(1:5)
  write(*,'(a)') ' '

!
! Test cross product of two vector fields
!
  resVec = vec1.x.vec2

  write(*,'(a)') ' Cross product of two vectrors '
  write(*,'(3(e13.6,1x))') resVec%x(1:3),resVec%y(1:3),resVec%z(1:3)
  write(*,'(a)') ' '

!
! Tensor product of two vector fields
!
  T = vec1.o.vec2

  write(*,'(a)') ' Tensor product of two vectors: '
  write(*,'(3(e13.6,1x))') T%xx(5),T%xy(5),T%xz(5)
  write(*,'(3(e13.6,1x))') T%yx(5),T%yy(5),T%yz(5)
  write(*,'(3(e13.6,1x))') T%zx(5),T%zy(5),T%zz(5)
  write(*,'(a)') ' '

!
! Deviatoric part of a tensor
!
  T = .dev.T

  write(*,'(a)') ' Deviatoric part of that tensor (at one point): '
  write(*,'(3(e13.6,1x))') T%xx(5),T%xy(5),T%xz(5)
  write(*,'(3(e13.6,1x))') T%yx(5),T%yy(5),T%yz(5)
  write(*,'(3(e13.6,1x))') T%zx(5),T%zy(5),T%zz(5)
  write(*,'(a)') ' '

  !...curl of a derived tensor field

  resVec = .curl.(3.14_dp*T+T)

  ! ...set field name for vtk output file name:
  resVec % field_name = 'CurlField'

  write(*,'(a)') ' ...curl of a derived tensor field '
  write(*,'(1x,a)') resVec%field_name
  write(*,'(3(e15.8,1x))') resVec%x(1:3),resVec%y(1:3),resVec%z(1:3)
  write(*,'(a)') ' '
!
! Magnitude squared of a tensor field
!
  phi = .magSq.T

  write(*,'(a)') ' Magnitude squared of the same tensor field'
  write(*,'(e15.8)') phi%mag(1:5)
  write(*,'(a)') ' '


!
! Test gradient of a scalar field
!  

  ! First set up the linear field
  ! psi = volScalarField("linear_field", xc + yc + zc )
  psi = new_volScalarField( numTotal )
  ! psi%mag = 1.0
  psi%mag(1:numCells)  = xc(1:numCells) + yc(1:numCells) + zc(1:numCells)
  psi%mag(numCells+1:) = xf(numInnerFaces+1:)+yf(numInnerFaces+1:)+zf(numInnerFaces+1:)


  ! Now call the explicit gradient operation
  resVec = fvxGrad( psi )

  ! This is useful for book-keeping, lets give it a name:
  ! resVec%field_name = 'gradient_field'

  ! call grad_gauss(psi%mag, resVec%x, resVec%y, resVec%z)

  write(*,'(a)') ' Magnitude of a linear field'
  write(*,'(e15.8)') psi%mag(1:440)
  write(*,'(a)') ' '
  write(*,'(a)') ' Gradient of a scalar field - fvx operation: '
  write(*,'(1x,a)') resVec%field_name
  do i=1,440
    write(*,'(i0,1x,3(e15.8,1x))') i,resVec%x(i),resVec%y(i),resVec%z(i)
  enddo
  write(*,'(a)') ' '

!
! Interpolation of a scalar field
!
  sf = fvxInterpolate( psi )

  write(*,'(a)') ' Interpolated scalar field to cell faces'
  write(*,'(e15.8)') sf%mag(1:1280)
  write(*,'(a)') ' '

!
! Now test vector identities: divergence of the curl is zero
!

  Vec2%x = 1.0; Vec2%y = 1.0; Vec2%z = 1.0
  phi = fvxDiv( Vec2 )

  write(*,'(a)') ' Divergence of a curl is zero'
   do i=1,1280
    write(*,'(i0,1x,e15.8)') i,phi%mag(i)
   enddo
  write(*,'(a)') ' '

  ! 6) Write output files in Paraview .vtu format
!+-----------------------------------------------------------------------------+

!  ierr = write_volScalarField_field( phi )

!  ierr = write_volVectorField_field( resVec )


!
!  > Terminate.
!
   call say_goodbye( )

  end program

