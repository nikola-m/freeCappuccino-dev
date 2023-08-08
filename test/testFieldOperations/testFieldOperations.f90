  program testFieldOperation
  use types
  use utils
  use geometry
  use tensorFields
  use fvxGradient
  use fvxInterpolation
  use fvxDivergence

  implicit none

  integer :: i
  real(dp) :: alpha  = 3.14_dp

! Locals
  !integer :: ierr

  type(volVectorField) :: vec1, vec2, vec3
  type(volScalarField) :: phi,psi
  type(volVectorField) :: resVec
  type(volTensorField) :: T,D
  type(surfaceScalarField) :: sf


  ! Intro
!+-----------------------------------------------------------------------------+

  call show_logo( )


  ! Read mesh files and calculate mesh geometrical quantities
!+-----------------------------------------------------------------------------+

  call read_mesh_openFoam

  ! Creation and manipulation of tensor fields
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

  write(*,'(a)') ' Dot product of two vector fields:'
  write(*,'(e13.6)') phi%mag(1:5)
  write(*,'(a)') ' ...'

!
! Test cross product of two vector fields
!
  resVec = vec1.x.vec2

  write(*,'(a)') ' Cross product of two vector fields:'
  write(*,'(3(e13.6,1x))') resVec%x(1:5),resVec%y(1:5),resVec%z(1:5)
  write(*,'(a)') ' ...'

!
! Tensor product of two vector fields
!
  T = vec1.o.vec2

  write(*,'(a)') ' Tensor product of two vector fields at a field point: '
  write(*,'(3(e13.6,1x))') T%xx(5),T%xy(5),T%xz(5)
  write(*,'(3(e13.6,1x))') T%yx(5),T%yy(5),T%yz(5)
  write(*,'(3(e13.6,1x))') T%zx(5),T%zy(5),T%zz(5)
  write(*,'(a)') ' '

!
! Deviatoric part of a tensor
!
  T = .dev.T

  write(*,'(a)') ' Deviatoric part of that tensor: '
  write(*,'(3(e13.6,1x))') T%xx(5),T%xy(5),T%xz(5)
  write(*,'(3(e13.6,1x))') T%yx(5),T%yy(5),T%yz(5)
  write(*,'(3(e13.6,1x))') T%zx(5),T%zy(5),T%zz(5)
  write(*,'(a)') ' '

  !...curl of a derived tensor field

  resVec = .curl.(alpha*T+T)

  ! ...set field name for vtk output file name:
  resVec % field_name = 'CurlField'

  write(*,'(a)') 'curl using a Hodge dual - uses derived tensor field as input '
  write(*,'(1x,a)') resVec%field_name
  write(*,'(3(e15.8,1x))') resVec%x(1:2),resVec%y(1:20),resVec%z(1:20)
  write(*,'(a)') ' ...'

!
! Now test vector identities: divergence of the curl is zero
!

  phi = fvxDiv( resVec )

  write(*,'(a)') ' Divergence of a curl is zero'
   do i=1,40
    write(*,'(i0,1x,e15.8)') i,phi%mag(i)
   enddo
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

  psi%mag(1:numCells)  = xc(1:numCells) + yc(1:numCells) + zc(1:numCells)
  psi%mag(numCells+1: ) = xf(numInnerFaces+1: )+yf(numInnerFaces+1: )+zf(numInnerFaces+1: )


  ! Now call the explicit gradient operation
  resVec = Grad( psi )

  ! This is useful for book-keeping, lets give it a name:
  ! resVec%field_name = 'gradient_field'

  write(*,'(a)') ' Magnitude of a linear field'
  write(*,'(e15.8)') psi%mag(1:20)
  write(*,'(a)') ' ...\n'
  write(*,'(a)') ' Gradient of a scalar field - fvx operation: '
  write(*,'(1x,a)') resVec%field_name
  do i=1,400
    write(*,'(i0,1x,3(e15.8,1x))') i,resVec%x(i),resVec%y(i),resVec%z(i)
  enddo
  write(*,'(a)') ' '

!
! Interpolation of a scalar field
!
  sf = fvxInterpolate( psi )

  write(*,'(a)') ' Interpolated scalar field to cell faces'
  do i=1,400
    write(*,'(e15.8)') sf%mag(i)
  enddo
  write(*,'(a)') ' '

!
! The traceless symmetric part of the square of the velocity gradient tensor
!
  vec3 = new_volVectorField( numTotal )
  vec3%x(:) = 1.0_dp
  vec3%y(:) = 2.0_dp
  vec3%z(:) = 3.0_dp  

  vec3%x( 1:numCells ) = xc(1:numCells)
  vec3%y( 1:numCells ) = yc(1:numCells)
  vec3%z( 1:numCells ) = zc(1:numCells)

  D = Grad( vec3 )
  T = .dev.(.symm.(D*D))

  write(*,'(a)') ' Traceless symmetric part of the square of the velocity gradient tensor D (at one point): '
  do i=1,20
    write(*,'(i0,1x,3(e15.8,1x))') i
    write(*,'(3(e13.6,1x))') T%xx(i),T%xy(i),T%xz(i)
    write(*,'(3(e13.6,1x))') T%yx(i),T%yy(i),T%yz(i)
    write(*,'(3(e13.6,1x))') T%zx(i),T%zy(i),T%zz(i)
    write(*,'(a)') ' '
  enddo

  psi = .magSq.T
  write(*,'(a)') ' Magnitude square of that tensor field (this is a scalar field)'
  do i=1,20
    write(*,'(e15.8)') psi%mag(i)
  enddo
  write(*,'(a)') ' '

!
! This is how we get Q-criteria field for vortex identification,
! if we suppose that velocity gradient tensor is stored in D.
!
  ! D = Grad( vec1 )
  ! psi = 0.5_dp *( .sq.(.tr.D)  - .tr.( D*D ) )
  ! write(*,'(a)') ' Q criteria'
  ! do i=1,20
  !   write(*,'(e15.8)') psi%mag(i)
  ! enddo
  ! write(*,'(a)') ' '



!
!  > Terminate.
!
   call say_goodbye( )

  end program

