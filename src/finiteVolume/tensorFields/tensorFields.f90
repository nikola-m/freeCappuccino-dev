module tensorFields
!
! Definition of volume and surface tensor fields and operations on them.
!
use types

implicit none

!
! > Volume fields
!

type volScalarField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: mag
end type

type volVectorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: x, y, z
end type

type volSymmetricTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable ::     yy, yz
  real(dp), dimension(:), allocatable ::         zz
end type

type volTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable :: yx, yy, yz
  real(dp), dimension(:), allocatable :: zx, zy, zz
end type


!
! > Surface fields
!

type surfaceScalarField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: mag
end type

type surfaceVectorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: x, y, z
end type

type surfaceSymmetricTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable ::     yy, yz
  real(dp), dimension(:), allocatable ::         zz
end type

type surfaceTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable :: yx, yy, yz
  real(dp), dimension(:), allocatable :: zx, zy, zz
end type


!
! > Operations on vector and tensor fields
!

interface operator(*)
   module procedure calc_inner_product
   module procedure calc_inner_product_surface_vectors
   module procedure calc_inner_product_rank2_and_rank1_tensors
   module procedure calc_inner_product_rank1_and_rank2_tensors
end interface

interface operator(.x.)
   module procedure calc_cross_product
end interface

interface operator(.o.)
   module procedure calc_tensor_product
   module procedure calc_tensor_product_rank2_tensors
end interface

interface operator(**)
   module procedure calc_inner_product_rank2_tensor
   module procedure calc_inner_product_rank2_symmetric_tensor
end interface

interface operator(.trans.)
   module procedure transpose_rank2_tensor
end interface

interface operator(.tr.)
   module procedure trace_rank2_tensor
   module procedure trace_rank2_symmetric_tensor
end interface

interface operator(.Sq.)
   module procedure square_scalar_field
   module procedure square_rank2_tensor
   module procedure square_rank2_symmetric_tensor
end interface

interface operator(.det.)
   module procedure determinant_rank2_tensor
   module procedure determinant_rank2_symmetric_tensor
end interface

interface operator(.diagonal.)
   module procedure diagonal
end interface

interface operator(.hodge.)
   module procedure hodge_dual
end interface

interface operator(.curl.)
   module procedure curl
end interface

interface operator(.symm.)
   module procedure symm
end interface

interface operator(.skew.)
   module procedure skew
end interface

interface operator(.magSq.)
   module procedure magSqrTensorField
   module procedure magSqrSymmetricTensorField
end interface

interface operator(.mag.)
   module procedure magTensorField
   module procedure magSymmetricTensorField
end interface

interface operator(.dev.)
   module procedure deviatoric_part_rank2_tensor
end interface

interface operator(.devTwo.)
   module procedure deviatoric_part_rank2_tensor_23
end interface

interface operator(.hyd.)
   module procedure hydrostatic_part_rank2_tensor
end interface

interface power
   module procedure power_volScalarField
   module procedure power_volVectorField
   module procedure power_volTensorField
   module procedure power_volSymmetricTensorField
end interface


! > Operator overloading for vector and tensor fields

interface operator(+)
   module procedure add_tensors
   module procedure add_volScalarFields
end interface

interface operator(-)
   module procedure subtract_tensors
   module procedure subtract_volScalarFields
end interface


! Various ways of multiplying scalar and a vector/tensor
interface operator(*)
   module procedure scalar_volScalarField_multiply
   module procedure scalarField_volScalarField_multiply
   module procedure scalar_volVectorField_multiply
   module procedure scalarField_volVectorField_multiply
   module procedure volScalarField_volVectorField_multiply
   module procedure scalar_volTensorField_multiply
   module procedure scalarField_volTensorField_multiply
   module procedure volScalarField_volTensorField_multiply
   module procedure scalar_surfaceVectorField_multiply
   module procedure scalarField_surfaceVectorField_multiply
   module procedure surfaceScalarField_surfaceTensorField_multiply
end interface

interface operator(/)
   module procedure volScalarField_volScalarField_divide
end interface

public

contains

!
! > Create new fields
!

function new_volScalarField( num )
    implicit none
    integer, intent(in) :: num
    type(volScalarField) :: new_volScalarField

    new_volScalarField % field_name = 'unknownVolScalarField'

    allocate(new_volScalarField%mag ( num ))
end function new_volScalarField


function new_volVectorField( num )
    implicit none
    integer, intent(in) :: num
    type(volVectorField) :: new_volVectorField

    new_volVectorField % field_name = 'unknownVolVectorField'

    allocate(new_volVectorField%x ( num ))
    allocate(new_volVectorField%y ( num ))
    allocate(new_volVectorField%z ( num ))
end function new_volVectorField


function new_volSymmetricTensorField( num )
    implicit none
    integer, intent(in) :: num
    type(volSymmetricTensorField) :: new_volSymmetricTensorField

    new_volSymmetricTensorField % field_name = 'unknownVolSymmTensorField'

    allocate(new_volSymmetricTensorField%xx ( num ))
    allocate(new_volSymmetricTensorField%xy ( num ))
    allocate(new_volSymmetricTensorField%xz ( num ))

    allocate(new_volSymmetricTensorField%yy ( num ))
    allocate(new_volSymmetricTensorField%yz ( num ))

    allocate(new_volSymmetricTensorField%zz ( num ))
end function new_volSymmetricTensorField

function new_volTensorField( num )
    implicit none
    integer, intent(in) :: num
    type(volTensorField) :: new_volTensorField

    new_volTensorField % field_name = 'unknownVolTensorField'

    allocate(new_volTensorField%xx ( num ))
    allocate(new_volTensorField%xy ( num ))
    allocate(new_volTensorField%xz ( num ))

    allocate(new_volTensorField%yx ( num ))
    allocate(new_volTensorField%yy ( num ))
    allocate(new_volTensorField%yz ( num ))

    allocate(new_volTensorField%zx ( num ))
    allocate(new_volTensorField%zy ( num ))
    allocate(new_volTensorField%zz ( num ))
end function new_volTensorField

!

function new_surfaceScalarField( num )
    implicit none
    integer, intent(in) :: num
    type(surfaceScalarField) :: new_surfaceScalarField

    new_surfaceScalarField % field_name = 'unknownSurfaceScalarField'

    allocate(new_surfaceScalarField%mag ( num ))
end function new_surfaceScalarField

function new_surfaceVectorField(num)
    implicit none
    integer, intent(in) :: num
    type(surfaceVectorField) :: new_surfaceVectorField

    new_surfaceVectorField % field_name = 'unknownSurfaceVectorField'

    allocate(new_surfaceVectorField%x ( num ))
    allocate(new_surfaceVectorField%y ( num ))
    allocate(new_surfaceVectorField%z ( num ))
end function new_surfaceVectorField

function new_surfaceSymmetricTensorField( num )
    implicit none
    integer, intent(in) :: num
    type(surfaceSymmetricTensorField) :: new_surfaceSymmetricTensorField

    new_surfaceSymmetricTensorField % field_name = 'unknownSurfSymTensorField'

    allocate(new_surfaceSymmetricTensorField%xx ( num ))
    allocate(new_surfaceSymmetricTensorField%xy ( num ))
    allocate(new_surfaceSymmetricTensorField%xz ( num ))

    allocate(new_surfaceSymmetricTensorField%yy ( num ))
    allocate(new_surfaceSymmetricTensorField%yz ( num ))

    allocate(new_surfaceSymmetricTensorField%zz ( num ))
end function new_surfaceSymmetricTensorField

function new_surfaceTensorField(num)
    implicit none
    integer, intent(in) :: num
    type(surfaceTensorField) :: new_surfaceTensorField

    new_surfaceTensorField % field_name = 'unknownSurfaceTensorField'

    allocate(new_surfaceTensorField%xx ( num ))
    allocate(new_surfaceTensorField%xy ( num ))
    allocate(new_surfaceTensorField%xz ( num ))

    allocate(new_surfaceTensorField%yx ( num ))
    allocate(new_surfaceTensorField%yy ( num ))
    allocate(new_surfaceTensorField%yz ( num ))

    allocate(new_surfaceTensorField%xz ( num ))
    allocate(new_surfaceTensorField%yz ( num ))
    allocate(new_surfaceTensorField%zz ( num ))
end function new_surfaceTensorField

!
! > Operations over fields
!

! The .dot. operator defining scalar product between two vector fields and ...

function calc_inner_product(v1, v2)  result(inner_product)
    implicit none
    type(volVectorField), intent(in) :: v1, v2
    type(volScalarField)             :: inner_product
    integer :: i,num

    num = size( v1%x )

    ! Simple check
    if (num /= size(v2%x)) then
      write(*,*) "Vectors in dot product are not the same size!"
      stop
    endif

    inner_product = new_volScalarField(num)

    do i = 1,num
        inner_product%mag(i) = v1%x(i) * v2%x(i) + v1%y(i) * v2%y(i) + v1%z(i) * v2%z(i)
    enddo

end function


function calc_inner_product_surface_vectors(v1, v2)  result(inner_product)
    implicit none
    type(surfaceVectorField), intent(in) :: v1, v2
    type(surfaceScalarField)             :: inner_product
    integer :: i,num

    num = size( v1%x )

    ! Simple check
    if (num /= size(v2%x)) then
      write(*,*) "Vectors in dot product are not the same size!"
      stop
    endif

    inner_product = new_surfaceScalarField(num)

    do i = 1,num
        inner_product%mag(i) = v1%x(i) * v2%x(i) + v1%y(i) * v2%y(i) + v1%z(i) * v2%z(i)
    enddo
end function


! ... inner product between tensor and vector vi = Tij*vj, and...

function calc_inner_product_rank2_and_rank1_tensors(T, v1)  result(v2)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i,num

    num = size( v1%x )

    v2 = new_volVectorField(num)

    do i = 1,num
        v2%x (i) = T%xx(i) * v1%x(i) + T%xy(i) * v1%y(i) + T%xz(i) * v1%z(i)  
        v2%y (i) = T%yx(i) * v1%x(i) + T%yy(i) * v1%y(i) + T%yz(i) * v1%z(i)  
        v2%z (i) = T%zx(i) * v1%x(i) + T%zy(i) * v1%y(i) + T%zz(i) * v1%z(i)
    enddo
end function

! ... inner product between vector and tensor vi = vj*Tij

function calc_inner_product_rank1_and_rank2_tensors(v1,T)  result(v2)
    implicit none
    type(volVectorField), intent(in) :: v1
    type(volTensorField), intent(in) :: T
    type(volVectorField)             :: v2
    integer :: num

    num = size( v1%x )

    v2 = new_volVectorField(num)
!
!   bi = Tji*aj using derived operators:
!
!        transpose tensor Tij
!        |               inner product bi = Tij*aj
!        |               |
    v2 = .trans.T * v1

end function


! The cross product of two vector fields

function calc_cross_product(v1, v2)  result(v3)
    implicit none
    type(volVectorField), intent(in) :: v1, v2
    type(volVectorField)             :: v3
    integer :: i,num

    num = size(v1%x)

    v3 = new_volVectorField( num )

    do i = 1,num
        v3%x (i) = v1%y(i) * v2%z(i) - v1%z(i) * v2%y(i)  
        v3%y (i) = v1%z(i) * v2%x(i) - v1%x(i) * v2%z(i) 
        v3%z (i) = v1%x(i) * v2%y(i) - v1%y(i) * v2%x(i)
    enddo
end function


! ! The .o. operator defining tensor, or outer product between two column vectors, or ...

function calc_tensor_product(v1, v2)  result(T)
    implicit none
    type(volVectorField), intent(in) :: v1, v2
    type(volTensorField)             :: T
    integer :: i,num

    num = size(v1%x)

    T = new_volTensorField(num)

    do i = 1,num
        T%xx(i) = v1%x(i) * v2%x(i) 
        T%xy(i) = v1%x(i) * v2%y(i) 
        T%xz(i) = v1%x(i) * v2%z(i) 

        T%yx(i) = v1%y(i) * v2%x(i) 
        T%yy(i) = v1%y(i) * v2%y(i) 
        T%yz(i) = v1%y(i) * v2%z(i)

        T%zx(i) = v1%z(i) * v2%x(i) 
        T%zy(i) = v1%z(i) * v2%y(i) 
        T%zz(i) = v1%z(i) * v2%z(i)
    enddo
end function

! ... between two tensors

function calc_tensor_product_rank2_tensors(T1, T2)  result(T3)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volTensorField)             :: T3
    integer :: i,num

    num = size(T1%xx)

    T3 = new_volTensorField(num)

    do i = 1,num
        T3%xx(i) = T1%xx(i) * T2%xx(i) + T1%xy(i) * T2%yx(i) + T1%xz(i) * T2%zx(i)
        T3%xy(i) = T1%xx(i) * T2%xy(i) + T1%xy(i) * T2%yy(i) + T1%xz(i) * T2%zy(i) 
        T3%xz(i) = T1%xx(i) * T2%xz(i) + T1%xy(i) * T2%yz(i) + T1%xz(i) * T2%zz(i) 

        T3%yx(i) = T1%yx(i) * T2%xx(i) + T1%yy(i) * T2%yx(i) + T1%yz(i) * T2%zx(i) 
        T3%yy(i) = T1%yx(i) * T2%xy(i) + T1%yy(i) * T2%yy(i) + T1%yz(i) * T2%zy(i) 
        T3%yz(i) = T1%yx(i) * T2%xz(i) + T1%yy(i) * T2%yz(i) + T1%yz(i) * T2%zz(i)

        T3%zx(i) = T1%zx(i) * T2%xx(i) + T1%zx(i) * T2%yx(i) + T1%zz(i) * T2%zx(i)
        T3%zy(i) = T1%zx(i) * T2%xy(i) + T1%zy(i) * T2%yy(i) + T1%zz(i) * T2%zy(i) 
        T3%zz(i) = T1%zx(i) * T2%xz(i) + T1%zy(i) * T2%yz(i) + T1%zz(i) * T2%zz(i)
    enddo
end function


function calc_inner_product_rank2_tensor(T1, T2) result(inner_product)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volScalarField) :: inner_product
    integer :: i,num

    num = size(T1%xx)

    inner_product = new_volScalarField(num)

    do i = 1,num
        inner_product%mag(i) = T1%xx(i) * T2%xx(i) + T1%xy(i) * T2%xy(i) + T1%xz(i) * T2%xz(i)  &
                             + T1%yx(i) * T2%yx(i) + T1%yy(i) * T2%yy(i) + T1%yz(i) * T2%yz(i)  &
                             + T1%zx(i) * T2%zx(i) + T1%zy(i) * T2%zy(i) + T1%zz(i) * T2%zz(i)
    enddo
end function

function calc_inner_product_rank2_symmetric_tensor(T1, T2) result(inner_product)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T1, T2
    type(volScalarField) :: inner_product
    integer :: i,num

    num = size(T1%xx)

    inner_product = new_volScalarField(num)

    do i = 1,num
        inner_product%mag(i) = T1%xx(i) * T2%xx(i) + T1%yy(i) * T2%yy(i) + T1%zz(i) * T2%zz(i) &
                       + 2 * ( T1%xy(i) * T2%xy(i) + T1%xz(i) * T2%xz(i) + T1%yz(i) * T2%yz(i) )
                                                                         
    enddo
end function


function transpose_rank2_tensor(T1)  result(T2)
    implicit none
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i,num

    num = size(T1%xx)

    T2 = new_volTensorField(num)

    do i = 1,num
        T2%xx(i) = T1%xx(i)
        T2%xy(i) = T1%yx(i) 
        T2%xz(i) = T1%zx(i)

        T2%yx(i) = T1%xy(i)
        T2%yy(i) = T1%yy(i)
        T2%yz(i) = T1%zy(i)

        T2%zx(i) = T1%xz(i)
        T2%zy(i) = T1%yz(i)  
        T2%zz(i) = T1%zz(i)
    enddo
end function

function add_tensors(T1,T2)  result(T3)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volTensorField)             :: T3
    integer :: i,num

    num = size(T1%xx)

    T3 = new_volTensorField(num)

    do i = 1,num
        T3%xx(i) = T1%xx(i) + T2%xx(i)
        T3%xy(i) = T1%xy(i) + T2%xy(i)
        T3%xz(i) = T1%xz(i) + T2%xz(i)

        T3%yx(i) = T1%yx(i) + T2%yx(i)
        T3%yy(i) = T1%yy(i) + T2%yy(i)
        T3%yz(i) = T1%yz(i) + T2%yz(i)

        T3%zx(i) = T1%zx(i) + T2%zx(i)
        T3%zy(i) = T1%zy(i) + T2%zy(i)
        T3%zz(i) = T1%zz(i) + T2%zz(i)
    enddo
end function


function add_volScalarFields(s1,s2)  result(s3)
    implicit none
    type(volScalarField), intent(in) :: s1, s2
    type(volScalarField)             :: s3
    integer :: i,num

    num = size(s1%mag)

    s3 = new_volScalarField(num)

    do i = 1,num
        s3%mag(i) = s1%mag(i) + s2%mag(i)
    enddo
end function

function subtract_volScalarFields(s1,s2)  result(s3)
    implicit none
    type(volScalarField), intent(in) :: s1, s2
    type(volScalarField)             :: s3
    integer :: i,num

    num = size(s1%mag)

    s3 = new_volScalarField(num)

    do i = 1,num
        s3%mag(i) = s1%mag(i) - s2%mag(i)
    enddo
end function

function subtract_tensors(T1,T2)  result(T3)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volTensorField)             :: T3
    integer :: i,num

    num = size(T1%xx)

    T3 = new_volTensorField(num)

    do i = 1,num
        T3%xx(i) = T1%xx(i) - T2%xx(i)
        T3%xy(i) = T1%xy(i) - T2%xy(i)
        T3%xz(i) = T1%xz(i) - T2%xz(i)

        T3%yx(i) = T1%yx(i) - T2%yx(i)
        T3%yy(i) = T1%yy(i) - T2%yy(i)
        T3%yz(i) = T1%yz(i) - T2%yz(i)

        T3%zx(i) = T1%zx(i) - T2%zx(i)
        T3%zy(i) = T1%zy(i) - T2%zy(i)
        T3%zz(i) = T1%zz(i) - T2%zz(i)
    enddo
end function

function scalar_volScalarField_multiply(alpha,s1)  result(s2)
    implicit none
    real(dp), intent(in) :: alpha
    type(volScalarField), intent(in) :: s1
    type(volScalarField)             :: s2
    integer :: i,num

    num = size(s1%mag)

    s2 = new_volScalarField(num)

    do i = 1,num
        s2%mag(i) = alpha * s1%mag(i)  
    enddo
end function

function scalarField_volScalarField_multiply(alpha,s1)  result(s2)
    implicit none
    real(dp), dimension(*), intent(in) :: alpha
    type(volScalarField), intent(in) :: s1
    type(volScalarField)             :: s2
    integer :: i,num

    num = size(s1%mag)

    s2 = new_volScalarField(num)

    do i = 1,num
        s2%mag(i) = alpha(i) * s1%mag(i)  
    enddo
end function

function scalar_volVectorField_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), intent(in) :: alpha
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i,num

    num = size(v1%x)

    v2 = new_volVectorField(num)

    do i = 1,num
        v2%x (i) = alpha * v1%x(i)  
        v2%y (i) = alpha * v1%y(i) 
        v2%z (i) = alpha * v1%z(i)
    enddo
end function


function scalarField_volVectorField_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), dimension(*), intent(in) :: alpha
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i,num

    num = size(v1%x)

    v2 = new_volVectorField(num)

    do i = 1,num
        v2%x (i) = alpha(i) * v1%x(i)  
        v2%y (i) = alpha(i) * v1%y(i) 
        v2%z (i) = alpha(i) * v1%z(i)
    enddo
end function scalarField_volVectorField_multiply


function volScalarField_volVectorField_multiply(alpha,v1)  result(v2)
    implicit none
    type(volScalarField), intent(in) :: alpha
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i,num

    num = size(v1%x)

    v2 = new_volVectorField(num)

    do i = 1,num
        v2%x (i) = alpha%mag(i) * v1%x(i)  
        v2%y (i) = alpha%mag(i) * v1%y(i) 
        v2%z (i) = alpha%mag(i) * v1%z(i)
    enddo
end function


function scalar_volTensorField_multiply(alpha,T1)  result(T2)
    implicit none
    real(dp), intent(in) :: alpha
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i,num

    num = size(T1%xx)

    T2 = new_volTensorField(num)

    do i = 1,num
        T2%xx(i) = alpha * T1%xx(i)
        T2%xy(i) = alpha * T1%xy(i)
        T2%xz(i) = alpha * T1%xz(i)

        T2%yx(i) = alpha * T1%yx(i)
        T2%yy(i) = alpha * T1%yy(i)
        T2%yz(i) = alpha * T1%yz(i)

        T2%zx(i) = alpha * T1%zx(i)
        T2%zy(i) = alpha * T1%zy(i)
        T2%zz(i) = alpha * T1%zz(i)
    enddo
end function scalar_volTensorField_multiply


function scalarField_volTensorField_multiply(alpha,T1)  result(T2)
    implicit none
    real(dp), dimension(*), intent(in) :: alpha
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i,num

    num = size(T1%xx)

    T2 = new_volTensorField(num)

    do i = 1,num
        T2%xx(i) = alpha(i) * T1%xx(i)
        T2%xy(i) = alpha(i) * T1%xy(i)
        T2%xz(i) = alpha(i) * T1%xz(i)

        T2%yx(i) = alpha(i) * T1%yx(i)
        T2%yy(i) = alpha(i) * T1%yy(i)
        T2%yz(i) = alpha(i) * T1%yz(i)

        T2%zx(i) = alpha(i) * T1%zx(i)
        T2%zy(i) = alpha(i) * T1%zy(i)
        T2%zz(i) = alpha(i) * T1%zz(i)
    enddo
end function scalarField_volTensorField_multiply


function volScalarField_volTensorField_multiply(alpha,T1)  result(T2)
    implicit none
    type(volScalarField), intent(in) :: alpha
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i,num

    num = size(alpha%mag)

    T2 = new_volTensorField(num)

    do i = 1,num
        T2%xx(i) = alpha%mag(i) * T1%xx(i)
        T2%xy(i) = alpha%mag(i) * T1%xy(i)
        T2%xz(i) = alpha%mag(i) * T1%xz(i)

        T2%yx(i) = alpha%mag(i) * T1%yx(i)
        T2%yy(i) = alpha%mag(i) * T1%yy(i)
        T2%yz(i) = alpha%mag(i) * T1%yz(i)

        T2%zx(i) = alpha%mag(i) * T1%zx(i)
        T2%zy(i) = alpha%mag(i) * T1%zy(i)
        T2%zz(i) = alpha%mag(i) * T1%zz(i)
    enddo
end function volScalarField_volTensorField_multiply

function surfaceScalarField_surfaceTensorField_multiply(alpha,T1)  result(T2)
    implicit none
    type(surfaceScalarField), intent(in) :: alpha
    type(surfaceTensorField), intent(in) :: T1
    type(surfaceTensorField)             :: T2
    integer :: i,num

    num = size(alpha%mag)

    T2 = new_surfaceTensorField(num)

    do i = 1,num
        T2%xx(i) = alpha%mag(i) * T1%xx(i)
        T2%xy(i) = alpha%mag(i) * T1%xy(i)
        T2%xz(i) = alpha%mag(i) * T1%xz(i)

        T2%yx(i) = alpha%mag(i) * T1%yx(i)
        T2%yy(i) = alpha%mag(i) * T1%yy(i)
        T2%yz(i) = alpha%mag(i) * T1%yz(i)

        T2%zx(i) = alpha%mag(i) * T1%zx(i)
        T2%zy(i) = alpha%mag(i) * T1%zy(i)
        T2%zz(i) = alpha%mag(i) * T1%zz(i)
    enddo
end function surfaceScalarField_surfaceTensorField_multiply

function scalar_surfaceVectorField_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), intent(in) :: alpha
    type(surfaceVectorField), intent(in) :: v1
    type(surfaceVectorField)             :: v2
    integer :: i,num

    num = size(v1%x)

    v2 = new_surfaceVectorField(num)

    do i = 1,num
        v2%x (i) = alpha * v1%x(i)  
        v2%y (i) = alpha * v1%y(i) 
        v2%z (i) = alpha * v1%z(i)
    enddo
end function


function scalarField_surfaceVectorField_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), dimension(*), intent(in) :: alpha
    type(surfaceVectorField), intent(in) :: v1
    type(surfaceVectorField)             :: v2
    integer :: i,num

    num = size(v1%x)

    v2 = new_surfaceVectorField(num)

    do i = 1,num
        v2%x (i) = alpha(i) * v1%x(i)  
        v2%y (i) = alpha(i) * v1%y(i) 
        v2%z (i) = alpha(i) * v1%z(i)
    enddo
end function

function volScalarField_volScalarField_divide(s1,s2)  result(s3)
    implicit none
    type(volScalarField), intent(in) :: s1,s2
    type(volScalarField)             :: s3
    integer :: i,num

    num = size(s1%mag)

    s3 = new_volScalarField(num)

    do i = 1,num
        s3%mag (i) = s1%mag(i) / ( s2%mag(i) + 1e-30 )  
    enddo
end function

function trace_rank2_symmetric_tensor(T) result(trace)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T
    type(volScalarField)                      :: trace
    integer :: i,num

    num = size(T%xx)

    trace = new_volScalarField(num)

    do i = 1,num
        trace%mag(i) = T%xx(i) + T%yy(i) + T%zz(i)
    enddo
end function trace_rank2_symmetric_tensor


function trace_rank2_tensor(T) result(trace)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volScalarField)             :: trace
    integer :: i,num
   
    num = size(T%xx)

    trace = new_volScalarField(num)

    do i = 1,num
        trace%mag(i) = T%xx(i) + T%yy(i) + T%zz(i)
    enddo
end function trace_rank2_tensor

function square_scalar_field(s1)  result(s2)
    implicit none
    type(volScalarField), intent(in) :: s1
    type(volScalarField)             :: s2
    integer :: i,num

    num = size(s1%mag)

    s2 = new_volScalarField(num)

    do i = 1,num
        s2%mag(i) = s1%mag(i)**2
    enddo
end function

function square_rank2_tensor(T1)  result(T2)
    implicit none
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i,num

    num = size(T1%xx)

    T2 = new_volTensorField(num)

    do i = 1,num
        T2%xx(i) = T1%xx(i)**2
        T2%xy(i) = T1%xy(i)**2
        T2%xz(i) = T1%xz(i)**2

        T2%yx(i) = T1%yx(i)**2
        T2%yy(i) = T1%yy(i)**2
        T2%yz(i) = T1%yz(i)**2

        T2%zx(i) = T1%zx(i)**2
        T2%zy(i) = T1%zy(i)**2
        T2%zz(i) = T1%zz(i)**2
    enddo
end function square_rank2_tensor


function square_rank2_symmetric_tensor(T1)  result(T2)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T1
    type(volSymmetricTensorField)             :: T2
    integer :: i,num

    num = size(T1%xx)

    T2 = new_volSymmetricTensorField(num)

    do i = 1,num
        T2%xx(i) = T1%xx(i)**2
        T2%xy(i) = T1%xy(i)**2
        T2%xz(i) = T1%xz(i)**2

        T2%yy(i) = T1%yy(i)**2
        T2%yz(i) = T1%yz(i)**2

        T2%zz(i) = T1%zz(i)**2
    enddo
end function square_rank2_symmetric_tensor

function power_volScalarField(s,p) result(ss)
    implicit none
    type(volScalarField), intent(in) :: s
    real(dp), intent(in) :: p
    type(volScalarField) :: ss
    integer :: i,num

    num = size(s%mag)
 
    ss = new_volScalarField(num)

    do i = 1,num
        ss%mag(i) = s%mag(i)**p 
    enddo
end function

function power_volVectorField(v,p) result(vv)
    implicit none
    type(volVectorField), intent(in) :: v
    real(dp), intent(in) :: p
    type(volVectorField) :: vv
    integer :: i,num

    num = size(v%x)

    vv = new_volVectorField(num)

    do i = 1,num
        vv%x (i) = v%x(i)**p 
        vv%y (i) = v%y(i)**p
        vv%z (i) = v%z(i)**p
    enddo
end function

function power_volTensorField(T,p) result(T2)
    implicit none
    type(volTensorField), intent(in) :: T
    real(dp), intent(in) :: p
    type(volTensorField) :: T2
    integer :: i,num

    num = size(T%xx)

    T2 = new_volTensorField(num)

    do i = 1,num
        T2%xx(i) = T%xx(i)**p
        T2%xy(i) = T%xy(i)**p
        T2%xz(i) = T%xz(i)**p

        T2%yx(i) = T%yx(i)**p
        T2%yy(i) = T%yy(i)**p
        T2%yz(i) = T%yz(i)**p

        T2%zx(i) = T%zx(i)**p
        T2%zy(i) = T%zy(i)**p
        T2%zz(i) = T%zz(i)**p
    enddo
end function

function power_volSymmetricTensorField(T,p) result(T2)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T
    real(dp), intent(in) :: p
    type(volSymmetricTensorField) :: T2

    integer :: i,num

    num = size(T%xx)

    T2 = new_volSymmetricTensorField(num)

    do i = 1,num
        T2%xx(i) = T%xx(i)**p
        T2%xy(i) = T%xy(i)**p
        T2%xz(i) = T%xz(i)**p

        T2%yy(i) = T%yy(i)**p
        T2%yz(i) = T%yz(i)**p

        T2%zz(i) = T%zz(i)**p
    enddo
end function

function determinant_rank2_symmetric_tensor(T) result(determinant)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T
    type(volScalarField)                      :: determinant
    integer :: i,num

    num = size(T%xx)

    determinant = new_volScalarField(num)

    do i = 1,num
        determinant%mag(i) = ( T%xx(i) * ( T%yy(i) * T%zz(i) - T%yz(i)*T%yz(i) ) - &
                               T%xy(i) * ( T%xy(i) * T%zz(i) - T%yz(i)*T%xz(i) ) + &
                               T%xz(i) * ( T%xy(i) * T%yz(i) - T%yy(i)*T%xz(i) )   )
    enddo
end function


function determinant_rank2_tensor(T) result(determinant)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volScalarField)                      :: determinant
    integer :: i,num

    num = size(T%xx)

    determinant = new_volScalarField(num)

    do i = 1,num
        determinant%mag(i) = ( T%xx(i) * ( T%yy(i) * T%zz(i) - T%yz(i)*T%zy(i) ) - &
                               T%xy(i) * ( T%yx(i) * T%zz(i) - T%yz(i)*T%zx(i) ) + &
                               T%xz(i) * ( T%yx(i) * T%zy(i) - T%yy(i)*T%zx(i) )   )
    enddo
end function


function diagonal(T) result(v)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volVectorField)             :: v
    integer :: i,num

    num = size(T%xx)

    v = new_volVectorField(num)

    do i = 1,num
        v%x (i) = T%xx(i)  
        v%y (i) = T%yy(i) 
        v%z (i) = T%zz(i)
    enddo
end function


function hodge_dual(T) result(v)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volVectorField)             :: v
    integer :: i,num

    num = size(T%xx)

    v = new_volVectorField(num)

    do i = 1,num
        v%x (i) = T%yz(i)  
        v%y (i) =-T%xz(i) 
        v%z (i) = T%xy(i)
    enddo
end function


function eye(num) result(T)
    implicit none
    integer, intent(in) :: num
    type(volTensorField) :: T
    integer :: i

    T = new_volTensorField(num)

    do i = 1,num
        T%xx(i) = 1.0_dp  
        T%xy(i) = 0.0_dp  
        T%xz(i) = 0.0_dp  

        T%yx(i) = 0.0_dp 
        T%yy(i) = 1.0_dp 
        T%yz(i) = 0.0_dp 
 
        T%zx(i) = 0.0_dp 
        T%zy(i) = 0.0_dp 
        T%zz(i) = 1.0_dp
    enddo
end function

function symm(T)  result(D)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: D

!              overloaded operator - here '*' multiplies tensor fields by a constant scalar        
!              |     overloaded operator - here '+' adds two tensor fields
!              |     |
    D = 0.5_dp * ( T + .trans.T )

end function

function skew(T)  result(S)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: S

!              overloaded operator - here '*' multiplies tensor field by a constant scalar        
!              |     overloaded operator - here '-' subtracts two tensor fields
!              |     |
    S = 0.5_dp * ( T - .trans.T )

end function

function curl(D) result(v)
!
! Curl of a vector field is twice Hodge dual of skew-symmetric part of gradient tensor of that vector field.
! Make sure that D in this function call is gradient tensor, e.g. velocity gradient tensor, so it makes sense.
! To obtain it you will have to use gradient function from the 'fvx' module, and apply it to a vector field
! in question: [volTensorField] D = fvxGrad([volVectorField] u)
!
    implicit none
    type(volTensorField), intent(in) :: D
    type(volVectorField)             :: v

!              overloaded operator - here '*' multiplies vector field by a constant scalar              
!              |  derived operators - Hodge dual of skew-symmetric part of tensor T
!              |  |            
    v = 2.0_dp * (.hodge.(.skew.D))
    
end function


function deviatoric_part_rank2_tensor(T)  result(devT)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: devT
    type(volTensorField)             :: I
    integer :: num

    num = size(T%xx)

    I  = eye(num)

!            overloaded operator - here '-' subtracts two tensor fields
!            |             overloaded operator - here '*' multiplies tensor field by a constant scalar
!            |             |         overloaded operator - here '*' multiplies tensor fields by a real scalar array of size[1:numCells]
!            |             |         |
    devT = T - ( 1./3.0_dp * ( .tr.T * I) )

end function


function deviatoric_part_rank2_tensor_23(T)  result(devT)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: devT
    type(volTensorField)             :: I
    integer :: num

    num = size(T%xx)

    I  = eye(num)

!            overloaded operator - here '-' subtracts two tensor fields
!            |             overloaded operator - here '*' multiplies tensor field by a constant scalar
!            |             |         overloaded operator - here '*' multiplies tensor fields by a real scalar array of size[1:numCells]
!            |             |         |
    devT = T - ( 2./3.0_dp * ( .tr.T * I) )
                !^
                !!----2/3 here !

end function


function hydrostatic_part_rank2_tensor(T)  result(hydT)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: hydT
    type(volTensorField)             :: I
    integer :: num

    num = size(T%xx)

    I  = eye(num)
!                     overloaded operator - here '*' multiplies tensor field by a constant scalar
!                     |         overloaded operator - here '*' multiplies tensor fields by a real scalar array of size[1:numCells]
!                     |         |
    hydT =  1./3.0_dp * ( .tr.T * I ) 

end function

! .magSq.

function magSqrTensorField(T) result(scalar)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volScalarField)             :: scalar

        scalar = T**T

end function

function magSqrSymmetricTensorField(S) result(scalar)
    implicit none
    type(volSymmetricTensorField), intent(in) :: S
    type(volScalarField)             :: scalar

        scalar = S**S

end function

! .mag.

function magTensorField(T) result(scalar)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volScalarField)             :: scalar

        scalar = T**T
        scalar%mag = sqrt(scalar%mag)

end function

function magSymmetricTensorField(S) result(scalar)
    implicit none
    type(volSymmetricTensorField), intent(in) :: S
    type(volScalarField)             :: scalar

        scalar = S**S
        scalar%mag = sqrt(scalar%mag)

end function

end module tensorFields