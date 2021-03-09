module fvxGradient
!
! Purpose:
!  Module for explicit operations on discrete tensor fields - Gradient operator.
!
! Description:
!  Module contains procedures for explicit manipulation of discrete tensor fields based on 
!  finite volume computations and integral theorems (e.g. Gauss) of vector calculus.
!  Discrete tensor fields are defined on a given finite volume mesh.
!  That means we are expecting discrete volume/surface scalar/vector/tensor fields.
!  Included operations are:
!  fvxGrad         (vol<Scalar/Vector>Field ->vol<Vector/Tensor>Field)
!
!  Author: Nikola Mirkov
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!
use types
use geometry
use sparse_matrix, only: ioffset,ja,diag
use tensorFields

implicit none

! Gradient discretization approach 
logical :: lstsq    = .False.                         
logical :: lstsq_qr = .False.                      
logical :: lstsq_dm = .False.                   
logical :: gauss    = .True.   

logical, save :: createdLSQMat = .False. 

character(len=20) :: limiter = 'no-limit'          ! Gradient limiter. Options: none, Barth-Jespersen, Venkatakrishnan, MDL

real(dp),dimension(:,:), allocatable ::  Dmat      !  d(6,nxyz) - when using bn, or dm version of the subroutine
real(dp),dimension(:,:,:), allocatable ::  D       !  when using qr version of the subroutine size(3,6,nxyz)!

real(dp), parameter :: small = 1e-20
real(dp), parameter :: zero = 0.0_dp

interface grad
  module procedure grad_volScalarField
  module procedure grad_volVectorField    
end interface

interface sngrad
  module procedure sngrad_volScalarField
  module procedure sngrad_VolVectorField     
end interface


private

public :: lstsq, lstsq_qr, lstsq_dm, gauss, limiter, &
          grad, sngrad, grad_gauss, grad_gauss_corrected, grad_scalar_field_w_option, &
          create_lsq_grad_matrix


contains



function grad_volScalarField(phi) result(dPhi)
!
! Description:
!  Computes cell centered gradients of a scalar field.
! Usage:
!   [type(volVectorField)] dPhi = Grad( [type(volScalarField)] phi)
!

  implicit none

  type(volScalarField), intent(in) :: phi
!
! > Result
!
  type(volVectorField) :: dPhi

!+-----------------------------------------------------------------------------+

  dPhi = new_volVectorField( numCells )

  dPhi%field_name = 'gradient_field'

  call grad_scalar_field( phi%mag, dPhi%x, dPhi%y, dPhi%z )


end function grad_volScalarField


function grad_volVectorField(U) result(G)
!
! Description:
!  Computes cell centered gradients of a scalar field.
! Usage:
!   [type(volTensorField)] G = Grad( [type(volVectorField)] U)
!

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(volTensorField) :: G

!+-----------------------------------------------------------------------------+

  G = new_volTensorField( numCells )

  G%field_name = 'gradient_field'

  ! Go trough every vector component

  call grad_scalar_field( U%x, G%xx, G%xy, G%xz )
  call grad_scalar_field( U%y, G%yx, G%yy, G%yz )
  call grad_scalar_field( U%z, G%zx, G%zy, G%zz ) 

end function grad_volVectorField


function sngrad_volScalarField(phi) result(sndPhi)
!
! Description:
!  Computes cell centered gradients of a scalar field.
! Usage:
!   [type(volVectorField)] sndPhi = snGrad( [type(volScalarField)] phi)
!

  implicit none

  type(volScalarField), intent(in) :: phi
!
! > Result
!
  type(surfaceVectorField) :: sndPhi

!
! > Locals
!
  integer :: i,ijp,ijn,nrelax
  character(len=12) :: approach
  real(dp) :: dfixii, dfiyii, dfizii
  real(dp) :: dfixi, dfiyi, dfizi
  type(volVectorField) :: dPhi

!+-----------------------------------------------------------------------------+

  sndPhi = new_surfaceVectorField( numCells )

  sndPhi%field_name = 'gradient_field'

  dPhi = grad_volScalarField( Phi )

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    nrelax = 0
    approach  = 'skewness'

    call sngrad_scalar_field(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                             phi%mag, dPhi%x, dPhi%y, dPhi%z, nrelax, approach, &
                             dfixi, dfiyi, dfizi, &
                             dfixii, dfiyii, dfizii)

    sndPhi%x(i) = dfixii
    sndPhi%y(i) = dfiyii
    sndPhi%z(i) = dfizii

  enddo

  ! Update boundaries...


end function sngrad_volScalarField


function sngrad_volVectorField(U) result(G)
!
! Description:
!  Computes cell centered gradients of a scalar field.
! Usage:
!   [type(surfaceTensorField)] snG = snGrad( [type(volVectorField)] U)
!

  implicit none

  type(volVectorField), intent(in) :: U
!
! > Result
!
  type(surfaceTensorField) :: snG
!
! > Local
!
  integer :: i,ijp,ijn,nrelax
  character(len=12) :: approach
  real(dp) :: dfixii, dfiyii, dfizii
  real(dp) :: dfixi, dfiyi, dfizi  
  type(volTensorField) :: G
!+-----------------------------------------------------------------------------+

  snG = new_surfaceTensorField( numCells )

  snG%field_name = 'surface_grad_field'

  G = grad_volVectorField( U )

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    nrelax = 0
    approach  = 'skewness'

    call sngrad_scalar_field( ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                              U%x, G%xx, G%xy, G%xz, nrelax, approach, &
                              dfixi, dfiyi, dfizi, &
                              dfixii, dfiyii, dfizii )

    snG%xx(i) = dfixii
    snG%xy(i) = dfiyii
    snG%xz(i) = dfizii

    call sngrad_scalar_field( ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                              U%y, G%yx, G%yy, G%yz, nrelax, approach, &
                              dfixi, dfiyi, dfizi, &
                              dfixii, dfiyii, dfizii )

    snG%yx(i) = dfixii
    snG%yy(i) = dfiyii
    snG%yz(i) = dfizii

    call sngrad_scalar_field( ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                              U%z, G%zx, G%zy, G%zz, nrelax, approach, &
                              dfixi, dfiyi, dfizi, &
                              dfixii, dfiyii, dfizii )

    snG%zx(i) = dfixii
    snG%zy(i) = dfiyii
    snG%zz(i) = dfizii

  enddo

  ! Update boundaries...

end function sngrad_volVectorField


!***********************************************************************
!
subroutine create_lsq_grad_matrix
!
!***********************************************************************
!
!  Discussion:
!    Prepare System Matrix For Least-Squares Gradient Calculation.
!    It is done by setting this --v value to one.
!           call grad_lsq(U,dUdxi,1)
!
!***********************************************************************
!

implicit none


  call allocate_lsq_grad_matrix

  if (lstsq) then

    call grad_lsq_matrix

  elseif (lstsq_qr) then 

    call grad_lsq_qr_matrix

  elseif (lstsq_dm) then 

    call grad_lsq_dm_matrix

  endif 

  createdLSQMat = .True.

end subroutine


!***********************************************************************
!
subroutine allocate_lsq_grad_matrix
!
!***********************************************************************
!
implicit none
  
  integer :: ierr

  if( lstsq .or. lstsq_dm ) then

    allocate( Dmat(9,numCells), stat=ierr ) 
      if(ierr /= 0)write(*,*)"allocation error: Dmat"

  elseif( lstsq_qr ) then

    allocate( D(3,6,numCells), stat=ierr ) 
      if(ierr /= 0)write(*,*)"allocation error: D"
      
  endif

end subroutine


!***********************************************************************
!
subroutine grad_scalar_field(phi,dPhidx,dPhidy,dPhidz)
!
!***********************************************************************
!

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(numCells), intent(inout) :: dPhidx,dPhidy,dPhidz

  dPhidx = 0.0_dp
  dPhidy = 0.0_dp
  dPhidz = 0.0_dp    

  ! Before calling check this
  if ( .not.gauss .and. .not.createdLSQMat ) call create_lsq_grad_matrix
  
  if (lstsq) then 

    call grad_lsq(phi,dPhidx,dPhidy,dPhidz)

  elseif (lstsq_qr) then 

    call grad_lsq_qr(phi,dPhidx,dPhidy,dPhidz)

  elseif (lstsq_dm) then 

    call grad_lsq_dm(phi,dPhidx,dPhidy,dPhidz)

  elseif ( lstsq_dm .and. gauss  ) then

    call grad_lsq_dm(phi,dPhidx,dPhidy,dPhidz)

    call grad_gauss_corrected(phi,dPhidx,dPhidy,dPhidz) 

  else 

    call grad_gauss(phi,dPhidx,dPhidy,dPhidz)

  endif 

  !
  ! Gradient limiter:
  !
  if( limiter == 'Barth-Jespersen') then

    call slope_limiter_Barth_Jespersen(phi, dPhidx,dPhidy,dPhidz )

  elseif( limiter == 'Venkatakrishnan') then

    call slope_limiter_Venkatakrishnan(phi, dPhidx,dPhidy,dPhidz )

  elseif( limiter == 'MDL') then

    call slope_limiter_multidimensional(phi, dPhidx,dPhidy,dPhidz )

  else
    ! no-limit
  endif

end subroutine


!***********************************************************************
!
subroutine grad_scalar_field_w_option(phi,dPhidx,dPhidy,dPhidz,option,option_limiter)
!
!***********************************************************************
!
! The main reason why we write this subroutine is to correct velocities
! in SIMPLE algorithm with conservative gradients, which is possible
! with Gauss rule. 
! We noticed it is better for calculation precision.
! But calling gradients with option may be nice anyhow.
!
!***********************************************************************
!

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(numCells), intent(inout) :: dPhidx,dPhidy,dPhidz
  character( len=* ), intent(in) :: option
  character( len=* ), intent(in) :: option_limiter

  dPhidx = 0.0_dp
  dPhidy = 0.0_dp
  dPhidz = 0.0_dp  
 

  ! Before calling check this
  if ( .not.gauss .and. .not.createdLSQMat ) call create_lsq_grad_matrix 

  if ( option == 'lsq' ) then 

    call grad_lsq(phi, dPhidx,dPhidy,dPhidz )

  elseif ( option == 'lsq_qr' ) then 

    call grad_lsq_qr(phi, dPhidx,dPhidy,dPhidz )

  elseif ( option == 'wlsq' ) then 

    call grad_lsq_dm(phi, dPhidx,dPhidy,dPhidz )

  elseif ( option == 'gauss_corrected' ) then

    call grad_gauss_corrected(phi, dPhidx,dPhidy,dPhidz ) 

  elseif ( option == 'gauss' ) then

    call grad_gauss(phi, dPhidx,dPhidy,dPhidz )

  endif 

  !
  ! Gradient limiter:
  !
  if( option_limiter == 'Barth-Jespersen') then

    call slope_limiter_Barth_Jespersen(phi, dPhidx,dPhidy,dPhidz )

  elseif( option_limiter == 'Venkatakrishnan') then

    call slope_limiter_Venkatakrishnan(phi, dPhidx,dPhidy,dPhidz )

  elseif( option_limiter == 'MDL') then

    call slope_limiter_multidimensional(phi, dPhidx,dPhidy,dPhidz )

  else
    ! no-limit
  endif

end subroutine



!***********************************************************************
!
subroutine slope_limiter_Barth_Jespersen(phi, dPhidx,dPhidy,dPhidz )
!
!***********************************************************************
!
!     Calculates slope limiter and appiies to scalar gradient:
!     Barth and Jespersen slope limiter:
!
!     AIAA-89-0366, The design and application of upwind schemes
!     on unstructured meshes, T.J.Barth, D.C.Jespersen, 1989.
!
!***********************************************************************
!

  implicit none

  ! Input
  real(dp),dimension(numTotal), intent(in) :: phi
  real(dp),dimension(numCells), intent(inout) :: dPhidx,dPhidy,dPhidz


  ! Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  real(dp) :: fimax,fimin,deltamax,deltamin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ioffset(inp) ))
    phi_min = phi(ja( ioffset(inp) ))

    do k=ioffset(inp)+1, ioffset(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_max, phi(ja(k)) )      
    enddo


    deltamax = fimax - phi(inp)
    deltamin = fimin - phi(inp)

    slopelimit = 1.0_dp

    do k=ioffset(inp), ioffset(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      delta_face=dPhidx(ijp)*(xc(ijn)-xc(ijp))+dPhidy(ijp)*(yc(ijn)-yc(ijp))+dPhidz(ijp)*(zc(ijn)-zc(ijp)) 


     if( abs(delta_face) < 1.e-6 )then
       r = 1.0_dp
     else if( delta_face > 0.0 )then
       r = deltamax/delta_face
     else
       r = deltamin/delta_face
     endif

     slopelimit = min( slopelimit , r )

    enddo

    dPhidx(inp) = slopelimit*dPhidx(inp)
    dPhidy(inp) = slopelimit*dPhidy(inp)
    dPhidz(inp) = slopelimit*dPhidz(inp)

  enddo

end subroutine




!***********************************************************************
!
subroutine slope_limiter_Venkatakrishnan(phi, dPhidx,dPhidy,dPhidz )
!
!***********************************************************************
!
!     Calculates slope limiter and appiies to scalar gradient:
!     Venkatakrishnan slope limiter:
!
!    AIAA-93-0880, On the accuracy of limiters and convergence
!    to steady state solutions, V.Venkatakrishnan, 1993
!
!***********************************************************************
!

  implicit none

  ! Input
  real(dp),dimension(numTotal), intent(in) :: phi
  real(dp),dimension(numCells), intent(inout) :: dPhidx,dPhidy,dPhidz


  ! Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  real(dp) :: fimax,fimin,deltamax,deltamin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ioffset(inp) ))
    phi_min = phi(ja( ioffset(inp) ))

    do k=ioffset(inp)+1, ioffset(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_max, phi(ja(k)) )      
    enddo


    deltamax = fimax - phi(inp)
    deltamin = fimin - phi(inp)

    slopelimit = 1.0_dp

    do k=ioffset(inp), ioffset(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      delta_face=dPhidx(ijp)*(xc(ijn)-xc(ijp))+dPhidy(ijp)*(yc(ijn)-yc(ijp))+dPhidz(ijp)*(zc(ijn)-zc(ijp)) 


     if( abs(delta_face) < 1.e-6 )then
       r = 1.0_dp
     else if( delta_face > 0.0 )then
       r = deltamax/delta_face
     else
       r = deltamin/delta_face
     endif

     slopelimit = min( slopelimit , (r**2+2.0*r)/(r**2+r+2.0) )

    enddo

    dPhidx(inp) = slopelimit*dPhidx(inp)
    dPhidy(inp) = slopelimit*dPhidy(inp)
    dPhidz(inp) = slopelimit*dPhidz(inp)

  enddo

end subroutine




!***********************************************************************
!
subroutine slope_limiter_multidimensional(phi, dPhidx,dPhidy,dPhidz )
!
!***********************************************************************
!
!  Calculates slope limiter and applies to scalar gradient:
!  Multidimensional slope limiter
!  Ref.: SE Kim, B Makarov, D Caraeni - A Multi-Dimensional Linear 
!  Reconstruction Scheme for Arbitrary Unstructured Meshes, AIAA 2003-3990.
!  The same slope limiter is used in Fluent.
!
!***********************************************************************
!

  implicit none

  ! Input
  real(dp),dimension(numTotal), intent(in) :: phi
  real(dp),dimension(numCells), intent(inout) :: dPhidx,dPhidy,dPhidz 


  ! Locals
  integer :: inp,ijp,k,iface
  real(dp) :: dPhi,dPhimax,dPhimin
  real(dp) :: gtx,gty,gtz
  real(dp) :: gx,gy,gz,gn
  real(dp) :: xpn,ypn,zpn,dpn
  real(dp) :: nx,ny,nz
  real(dp), dimension(:), allocatable :: phimax,phimin

  allocate(phimax(numCells),phimin(numCells))

  ! Find phi_max and phi_min in  neighbours of Co, including itself.
  do inp = 1,numCells

    ! max and min values over current cell and neighbors
    phimax(inp) = phi(ja( ioffset(inp) ))
    phimin(inp) = phi(ja( ioffset(inp) ))

    do k=ioffset(inp)+1, ioffset(inp+1)-1
      phimax(inp) = max( phimax(inp), phi(ja(k)) )
      phimin(inp) = min( phimin(inp), phi(ja(k)) )      
    enddo

  enddo


  ! Loop over inner faces:
  do iface=1,numInnerFaces

    do k=1,2 ! Do the same for both owner and neighbour cell

      if (k==1) then
        ijp = owner(iface)
      else
        ijp = owner(iface)
      endif

      ! Initialize gradient vector with current unlimited value
      gx = dPhidx(ijp)
      gy = dPhidy(ijp)
      gz = dPhidz(ijp)

      ! Distance vector between cell center and face center
      xpn=xf(iface)-xc(ijp)
      ypn=yf(iface)-yc(ijp)
      zpn=zf(iface)-zc(ijp)

      dpn=sqrt(xpn**2+ypn**2+zpn**2)

      nx = xpn/dpn
      ny = ypn/dpn
      nz = zpn/dpn

      gn = gx*nx+gy*ny+gz*nz

      gtx = gx - gn*nx
      gty = gy - gn*ny
      gtz = gz - gn*nz

      ! Increment from P cell to face centroid f
      dPhi=gx*(xf(iface)-xc(ijp))+gy*(yf(iface)-yc(ijp))+gz*(zf(iface)-zc(ijp)) 

      dPhimax = phimax(ijp)-phi(ijp)
      dPhimin = phimin(ijp)-phi(ijp)

      ! Check for overshoots and undershoots and correct accordingly.

      if ( phimax(ijp) > phi(ijp) .and. dPhi > dPhimax ) then

        gx = gtx + nx*dPhimax
        gy = gty + ny*dPhimax
        gz = gtz + nz*dPhimax

      endif

      if ( phimin(ijp) < phi(ijp) .and. dPhi < dPhimin ) then

        gx = gtx + nx*dPhimin
        gy = gty + ny*dPhimin
        gz = gtz + nz*dPhimin

      endif

      dPhidx(ijp) = gx
      dPhidy(ijp) = gy
      dPhidz(ijp) = gz

    enddo

  enddo

  deallocate(phimax,phimin)

end subroutine



! Least square gradients
subroutine grad_lsq_matrix
!
!***********************************************************************
!
!  Purpose:
!  Calculates COEFFICIENT MATRIX for cell-centered gradients using UNWEIGHTED 
!  Least-Squares approach.
!  Coefficient matrix should be calculated only once and stored in memory.
!
!  Description:
!  Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.,
!  also in Muzaferija and Gossman JCP paper from 1995.
!
!
!***********************************************************************
!

  implicit none

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface
  real(dp) :: Dx,Dy,Dz
  real(dp) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,tmp
!
!***********************************************************************
!

  ! Initialize Dmat matrix:
  Dmat = 0.0_dp

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    Dx = xc(ijn)-xc(ijp)
    Dy = yc(ijn)-yc(ijp)
    Dz = zc(ijn)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + Dx*Dx
    Dmat(1,ijn) = Dmat(1,ijn) + Dx*Dx 

    Dmat(4,ijp) = Dmat(4,ijp) + Dy*Dy
    Dmat(4,ijn) = Dmat(4,ijn) + Dy*Dy

    Dmat(6,ijp) = Dmat(6,ijp) + Dz*Dz
    Dmat(6,ijn) = Dmat(6,ijn) + Dz*Dz  

    Dmat(2,ijp) = Dmat(2,ijp) + Dx*Dy
    Dmat(2,ijn) = Dmat(2,ijn) + Dx*Dy

    Dmat(3,ijp) = Dmat(3,ijp) + Dx*Dz
    Dmat(3,ijn) = Dmat(3,ijn) + Dx*Dz

    Dmat(5,ijp) = Dmat(5,ijp) + Dy*Dz
    Dmat(5,ijn) = Dmat(5,ijn) + Dy*Dz                                                                 
  enddo     
 


  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)

    Dx = xf(iface)-xc(ijp)
    Dy = yf(iface)-yc(ijp)
    Dz = zf(iface)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + Dx*Dx
    Dmat(4,ijp) = Dmat(4,ijp) + Dy*Dy
    Dmat(6,ijp) = Dmat(6,ijp) + Dz*Dz
    Dmat(2,ijp) = Dmat(2,ijp) + Dx*Dy 
    Dmat(3,ijp) = Dmat(3,ijp) + Dx*Dz 
    Dmat(5,ijp) = Dmat(5,ijp) + Dy*Dz 
  end do


  ! Prepare for storage:
  do inp=1,numCells 

    ! Copy from Coefficient matrix 
    D11 = Dmat(1,inp)
    D12 = Dmat(2,inp)
    D13 = Dmat(3,inp)

    D22 = Dmat(4,inp)
    D23 = Dmat(5,inp)
    D33 = Dmat(6,inp)

    ! Symmetric part
    D21 = D12
    D31 = D13
    D32 = D23

    ! Denominator used troughout
    tmp = 1./(d11*d22*d33 - d11*d23*d32 - d12*d21*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31 + small)

    Dmat(1,inp) = (d22*d33 - d23*d32) * tmp
    Dmat(2,inp) = (d21*d33 - d23*d31) * tmp
    Dmat(3,inp) = (d21*d32 - d22*d31) * tmp

    Dmat(4,inp) = (d11*d33 - d13*d31) * tmp
    Dmat(5,inp) = (d12*d33 - d13*d32) * tmp
    Dmat(6,inp) = (d11*d32 - d12*d31) * tmp

    Dmat(7,inp) = (d12*d23 - d13*d22) * tmp
    Dmat(8,inp) = (d11*d23 - d13*d21) * tmp
    Dmat(9,inp) = (d11*d22 - d12*d21) * tmp

  enddo  
  
end subroutine


subroutine grad_lsq( Phi, dPhidx,dPhidy,dPhidz )
!
!***********************************************************************
!
!  Purpose:
!  Calculates cell-centered gradients using UNWEIGHTED Least-Squares approach.
!
!  Description:
!  Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.,
!  also in Muzaferija and Gossman JCP paper from 1995.
!
!  Arguments:
!
!  Phi - field variable which gradient we look for.
!  dPhidx,dPhidy,dPhidz  - cell centered gradient - a three component gradient vector.
!
!***********************************************************************
!

  implicit none

  real(dp),dimension(numTotal), intent(in)   :: Phi
  real(dp),dimension(numCells), intent(inout) ::  dPhidx,dPhidy,dPhidz 

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface
  real(dp) :: b1,b2,b3 
  real(dp) :: Dx,Dy,Dz

!
!***********************************************************************
!

  !
  ! *** COMMENT: *** 
  ! We want to save space, so we will reuse dPhidx,dPhidy,dPhidz 
  ! to store rhs vectors. These were, in earlier version, denoted
  ! with 'b1','b2', 'b3' respectively.
  ! This way it is maybe hardrer to read code but introduces great savings,
  ! because each of 'b1','b2', 'b3' should be numCells long!
  !

  ! Initialize rhs vector
  dPhidx = 0.0_dp
  dPhidy = 0.0_dp
  dPhidz = 0.0_dp

  ! Inner faces:

  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    Dx = ( xc(ijn)-xc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
    Dy = ( yc(ijn)-yc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
    Dz = ( zc(ijn)-zc(ijp) ) * ( Phi(ijn)-Phi(ijp) )

    dPhidx(ijp) = dPhidx(ijp) + Dx 
    dPhidx(ijn) = dPhidx(ijn) + Dx

    dPhidy(ijp) = dPhidy(ijp) + Dy
    dPhidy(ijn) = dPhidy(ijn) + Dy 

    dPhidz(ijp) = dPhidz(ijp) + Dz 
    dPhidz(ijn) = dPhidz(ijn) + Dz    

  enddo     

  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells + i

    Dx = (Phi(ijn)-Phi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = (Phi(ijn)-Phi(ijp))*(yf(iface)-yc(ijp))
    Dz = (Phi(ijn)-Phi(ijp))*(zf(iface)-zc(ijp)) 

    dPhidx(ijp) = dPhidx(ijp) + Dx
    dPhidy(ijp) = dPhidy(ijp) + Dy
    dPhidz(ijp) = dPhidz(ijp) + Dz

  enddo


  !
  ! Solve the system A*X = B.
  ! 

  do inp=1,numCells 

    b1 = dPhidx(inp)
    b2 = dPhidy(inp)
    b3 = dPhidz(inp)
            
    dPhidx(inp) = b1*Dmat(1,inp) - b2*Dmat(2,inp) +  b3*Dmat(3,inp) 
    dPhidy(inp) = b1*Dmat(4,inp) - b2*Dmat(5,inp) -  b3*Dmat(6,inp) 
    dPhidz(inp) = b1*Dmat(7,inp) - b2*Dmat(8,inp) +  b3*Dmat(9,inp) 

  enddo 

end subroutine


! Least square gradients via QR decomposition
!***********************************************************************
!
subroutine grad_lsq_qr_matrix
!
!***********************************************************************
!
!      Purpose:
!      Calculates COEFFICIENT MATRIX for cell-centered gradients using 
!      Least-Squares approach.
!
!      Description:
!      Uses QR decomposition of system matrix via Householder or via
!      Gramm-Schmidt.
!      QR decomposition is precomputed and R^(-1)*Q^T is stored in 
!      D array for every cell.
!
!      Coefficient matrix - should be calculated only once!
!
!***********************************************************************
!
  use matrix_module

  implicit none

  integer, parameter :: n=3, m=6  ! m is the number of neighbours, e.g. for structured 3D mesh it's 6

  ! Locals
  integer ::  i,l,k,ijp,ijn,inp,iface

  integer, dimension(numCells) :: neighbour_index  

  real(dp), dimension(m,n) :: Dtmp
  real(dp), dimension(n,m) :: Dtmpt


  !REAL(dp), DIMENSION(m,n) :: R
  !REAL(dp), DIMENSION(m,m) :: Q
  !REAL(dp), DIMENSION(n,n) :: R1
  !REAL(dp), DIMENSION(n,m) :: Q1t

  INTEGER :: INFO
  REAL(dp), DIMENSION(n) :: TAU
  INTEGER, DIMENSION(n) :: WORK
  REAL(dp), DIMENSION(m) :: v1,v2,v3
  REAL(dp), DIMENSION(m,m) :: H1,H2,H3,Ieye
  REAL(dp), DIMENSION(n,n) :: R
  REAL(dp), DIMENSION(m,m) :: Q

 
!**************************************************************************************************
! Coefficient matrix - should be calculated only once 
!**************************************************************************************************
  Dtmp = 0.0d0
  neighbour_index = 0

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

      neighbour_index(ijp) = neighbour_index(ijp) + 1 
      l = neighbour_index(ijp)
      D(1,l,ijp) = xc(ijn)-xc(ijp)
      D(2,l,ijp) = yc(ijn)-yc(ijp)
      D(3,l,ijp) = zc(ijn)-zc(ijp)
  
      neighbour_index(ijn) = neighbour_index(ijn) + 1 
      l = neighbour_index(ijn)   
      D(1,l,ijn) = xc(ijp)-xc(ijn)
      D(2,l,ijn) = yc(ijp)-yc(ijn)
      D(3,l,ijn) = zc(ijp)-zc(ijn)
                                                            
  enddo     

  ! Boundary faces:

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ! ijb = numCells + i
        neighbour_index(ijp) = neighbour_index(ijp) + 1
        l = neighbour_index(ijp)
        D(1,l,ijp) = xf(iface)-xc(ijp)
        D(2,l,ijp) = yf(iface)-yc(ijp)
        D(3,l,ijp) = zf(iface)-zc(ijp)
  end do

  ! Form system matrix using QR decomposition:

  ! Cell loop

  do inp=1,numCells

  l = neighbour_index(inp)

  Dtmpt = D(:,:,inp)
  Dtmp = transpose(Dtmpt)

  !1 ...Decompose A=QR using Householder
  !      call householder_qr(Dtmp, m, n, Q, R)
  !2 ...Decompose A=QR using Gram-Schmidt
  !      call mgs_qr(Dtmp, m, n, Q, R)

  !      Q = transpose(Q)
  !      Q1t = Q(1:n,1:m)     ! NOTE: A=Q1R1 is so-called 'thin QR factorization' - see Golub & Van Loan
                              ! Here Q1 is actually Q1^T a transpose of Q1(thin Q - Q with m-n column stripped off)
  !      R1 = R(1:n,1:n)      ! our Q1 is thin transpose of original Q.
  !      R1 = inv(R1)         ! inv is a function in matrix_module, now works only for 3x3 matrices.
  !      Q1t  = matmul(R1,Q1t) ! this is actually R^(-1)*Q^T - a matrix of size n x m.
  !      D(:,:,INP) = Q1t     ! Store it for later.

  !3....LAPACK routine DGEQRF
  CALL DGEQRF( l, N, Dtmp, M, TAU, WORK, N, INFO )

  ! Upper triangular matrix R
  R(1:n,1:n)=Dtmp(1:n,1:n)

  ! Create reflectors
  !H(i) = I - TAU * v * v'
  Ieye=eye(l)
  !v(1:i-1) = 0. and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i)
  v1(1) = 1.; v1(2:l)=Dtmp(2:l,1)
  H1 = rank_one_update(Ieye,l,l,v1,v1,-TAU(1))
  v2(1) = 0.; v2(2) = 1.; v2(3:l)=Dtmp(3:l,2)
  H2 = rank_one_update(Ieye,l,l,v2,v2,-TAU(2))
  v3(1:2) = 0.; v3(3) = 1.; v3(4:l)=Dtmp(4:l,3)
  H3 = rank_one_update(Ieye,l,l,v3,v3,-TAU(3))
  ! The matrix Q is represented as a product of elementary reflectors H1, H2, ..., Hn
  Q=matmul(H1,H2)
  Q=matmul(Q,H3)

  ! Form R_1^(-1)*Q_1^T explicitely:
  do k=1,neighbour_index(inp)
    D(1,k,inp) = q(k,1)/r(1,1) - (r(1,2)*q(k,2))/(r(1,1)*r(2,2)) + (q(k,3)*(r(1,2)*r(2,3) - r(1,3)*r(2,2)))/(r(1,1)*r(2,2)*r(3,3))
    D(2,k,inp) = q(k,2)/r(2,2) - (r(2,3)*q(k,3))/(r(2,2)*r(3,3))
    D(3,k,inp) = q(k,3)/r(3,3)
  enddo

  enddo

end subroutine

!***********************************************************************
!
subroutine grad_lsq_qr( Phi, dPhidx,dPhidy,dPhidz )
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using Least-Squares approach.
!
!      Description:
!      Uses QR decomposition of system matrix via Householder or via
!      Gramm-Schmidt.
!      QR decomposition is precomputed and R^(-1)*Q^T is stored in 
!      D array for every cell.
!
!      Arguments:
!
!      Phi - field variable
!      dPhidx,dPhidy,dPhidz - cell centered gradient - a three component gradient vector.
!
!***********************************************************************
!

  implicit none

  real(dp), dimension(numTotal), intent(in)   :: phi
  real(dp), dimension(numCells), intent(inout) :: dPhidx,dPhidy,dPhidz

  ! Locals
  integer, parameter :: n=3, m=6  ! m is the number of neighbours, e.g. for structured 3D mesh it's 6
  integer ::  i,l,ijp,ijn,inp,iface

  integer, dimension(:), allocatable :: neighbour_index  
  real(dp), dimension(:,:), allocatable   :: b
 
  allocate( neighbour_index(numCells) )
  allocate( b(m,numCells) )

  ! RHS vector
  b = 0.0d0
  neighbour_index = 0

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    neighbour_index(ijp) = neighbour_index(ijp) + 1
    l = neighbour_index(ijp)
    b(l,ijp) = Phi(ijn)-Phi(ijp)

    neighbour_index(ijn) = neighbour_index(ijn) + 1        
    l = neighbour_index(ijn)   
    b(l,ijn) = Phi(ijp)-Phi(ijn)
                                                                                                               
  enddo     

  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells + i

    neighbour_index(ijp) = neighbour_index(ijp) + 1
    l = neighbour_index(ijp)
    b(l,ijp) = Phi(ijn)-Phi(ijp)

  enddo

! Solve overdetermined system in least-sqare sense

  ! Cell loop
  do inp=1,numCells

    l = neighbour_index(inp)     

    !  ...using precomputed QR factorization and storing R^(-1)*Q^T in D
    dPhidx(INP) = sum(D(1,1:l,inp)*b(1:l,inp))
    dPhidy(INP) = sum(D(2,1:l,inp)*b(1:l,inp))
    dPhidz(INP) = sum(D(3,1:l,inp)*b(1:l,inp))

  enddo

  deallocate( neighbour_index )
  deallocate( b )

end subroutine


! Weighted least square gradients
subroutine grad_lsq_dm_matrix
!
!***********************************************************************
!
!  Purpose:
!  Calculates COEFFICIENT MATRIX for cell-centered gradients using UNWEIGHTED 
!  Least-Squares approach.
!  Coefficient matrix should be calculated only once and stored in memory.
!
!  Description:
!  Approach taken from a paper:
!  Dimitry Mavriplis, "Revisiting the Least Squares procedure for Gradient Reconstruction on Unstructured Meshes." NASA/CR-2003-212683.
!  Weights based on inverse cell-centers distance is added to improve conditioning of the system matrix for skewed meshes.
!  Reduces storage requirements compared to QR subroutine.
!  System matrix is symmetric and can be solved efficiently using Cholesky decomposition or by matrix inversion. 
!
!
!***********************************************************************
!

  implicit none

  ! Locals
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: w
  real(dp) :: Dx,Dy,Dz
  real(dp) :: d11,d12,d13,d21,d22,d23,d31,d32,d33
  real(dp) :: tmp

!
!***********************************************************************
!

  ! Coefficient matrix - should be calculated only once 

  ! Initialize dmat matrix:
  Dmat = 0.0d0

  ! Inner faces:     

  do i=1,numInnerFaces                                                       
    ijp = owner(i)

    !***------------------------------------
    ijn = neighbour(i)

    w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

    Dx = xc(ijn)-xc(ijp)
    Dy = yc(ijn)-yc(ijp)
    Dz = zc(ijn)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
    Dmat(1,ijn) = Dmat(1,ijn) + w*Dx*Dx 

    Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
    Dmat(4,ijn) = Dmat(4,ijn) + w*Dy*Dy

    Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
    Dmat(6,ijn) = Dmat(6,ijn) + w*Dz*Dz  

    Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy
    Dmat(2,ijn) = Dmat(2,ijn) + w*Dx*Dy

    Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz
    Dmat(3,ijn) = Dmat(3,ijn) + w*Dx*Dz

    Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz
    Dmat(5,ijn) = Dmat(5,ijn) + w*Dy*Dz     
    !***------------------------------------

  ! !
  ! ! **Extended interpolation molecule: neighbours of neighbours**
  ! !
  !   nb_loop: do k = ioffset( neighbour(i) ), ioffset( neighbour(i)+1 )-1

  !     ijn = ja(k)

  !     if (ijn == ijp) cycle nb_loop ! nemoj ownera dirati

  !     ! Paznja kada je ja(k) = diag( neighbour(i) ) ijn ce uzeti samog cell_neighboura, 
  !     ! tako da ovaj loop moze uraditi isto sto i gornji, pa je u toj verziji
  !     ! racunanja gradijenata ovaj loop gore nepotreban. Eto simplifikacije koda...

  !     w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

  !     Dx = xc(ijn)-xc(ijp)
  !     Dy = yc(ijn)-yc(ijp)
  !     Dz = zc(ijn)-zc(ijp)

  !     Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
  !     Dmat(1,ijn) = Dmat(1,ijn) + w*Dx*Dx 

  !     Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
  !     Dmat(4,ijn) = Dmat(4,ijn) + w*Dy*Dy

  !     Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
  !     Dmat(6,ijn) = Dmat(6,ijn) + w*Dz*Dz  

  !     Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy
  !     Dmat(2,ijn) = Dmat(2,ijn) + w*Dx*Dy

  !     Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz
  !     Dmat(3,ijn) = Dmat(3,ijn) + w*Dx*Dz
 
  !     Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz
  !     Dmat(5,ijn) = Dmat(5,ijn) + w*Dy*Dz 

  !   enddo nb_loop
                                                           
  enddo     

 
  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)

    w = 1.0_dp/((xf(iface)-xc(ijp))**2+(yf(iface)-yc(ijp))**2+(zf(iface)-zc(ijp))**2)

    Dx = xf(iface)-xc(ijp)
    Dy = yf(iface)-yc(ijp)
    Dz = zf(iface)-zc(ijp)

    Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
    Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
    Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
    Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy 
    Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz 
    Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz 
  end do


  ! Prepare for storage:
  do inp=1,numCells 

    ! Copy from Coefficient matrix 
    D11 = Dmat(1,inp)
    D12 = Dmat(2,inp)
    D13 = Dmat(3,inp)

    D22 = Dmat(4,inp)
    D23 = Dmat(5,inp)
    D33 = Dmat(6,inp)

    ! Symmetric part
    D21 = D12
    D31 = D13
    D32 = D23

    ! Denominator used troughout
    tmp = 1./(d11*d22*d33 - d11*d23*d32 - d12*d21*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31 + small)

    Dmat(1,inp) = (d22*d33 - d23*d32) * tmp
    Dmat(2,inp) = (d21*d33 - d23*d31) * tmp
    Dmat(3,inp) = (d21*d32 - d22*d31) * tmp

    Dmat(4,inp) = (d11*d33 - d13*d31) * tmp
    Dmat(5,inp) = (d12*d33 - d13*d32) * tmp
    Dmat(6,inp) = (d11*d32 - d12*d31) * tmp

    Dmat(7,inp) = (d12*d23 - d13*d22) * tmp
    Dmat(8,inp) = (d11*d23 - d13*d21) * tmp
    Dmat(9,inp) = (d11*d22 - d12*d21) * tmp

  enddo  

end subroutine


subroutine grad_lsq_dm( Phi, dPhidx,dPhidy,dPhidz )
!
!***********************************************************************
!
!  Purpose:
!    Calculates cell-centered gradients using WEIGHTED Least-Squares approach.
!
!  Description:
!    Approach taken from a paper:
!    Dimitry Mavriplis, "Revisiting the Least Squares procedure for Gradient Reconstruction on Unstructured Meshes." NASA/CR-2003-212683.
!    Weights based on inverse cell-centers distance is added to improve conditioning of the system matrix for skewed meshes.
!    Reduces storage requirements compared to QR subroutine.
!    System matrix is symmetric and can be solved efficiently using Cholesky decomposition or by matrix inversion. 
!
!  Arguments:
!    Phi - field variable which gradient we look for.
!    dPhidx, dPhidy, dPhidz - cell centered gradient - a three component gradient vector.
!
!***********************************************************************
!


  implicit none

  real(dp),dimension(numTotal), intent(in)   :: Phi
  real(dp),dimension(numCells), intent(inout) ::  dPhidx,dPhidy,dPhidz 

  ! Locals
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: w
  real(dp) :: Dx,Dy,Dz
  real(dp) :: b1,b2,b3 
!
!***********************************************************************
!
  !
  ! *** COMMENT: *** 
  ! We want to save space, so we will reuse dPhidx,dPhidy,dPhidz 
  ! to store rhs vectors. These were, in earlier version, denoted
  ! with 'b1','b2', 'b3' respectively.
  ! This way it is maybe hardrer to read code but introduces great savings,
  ! because each of 'b1','b2', 'b3' should be numCells long!
  !

  ! Initialize rhs vector
  dPhidx = 0.0_dp
  dPhidy = 0.0_dp
  dPhidz = 0.0_dp

  ! Inner faces:

  do i=1,numInnerFaces                                                       
    ijp = owner(i)

    !***------------------------------------
    ijn = neighbour(i)

    w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

    Dx = w * ( xc(ijn)-xc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
    Dy = w * ( yc(ijn)-yc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
    Dz = w * ( zc(ijn)-zc(ijp) ) * ( Phi(ijn)-Phi(ijp) )

    dPhidx(ijp) = dPhidx(ijp) + Dx 
    dPhidx(ijn) = dPhidx(ijn) + Dx

    dPhidy(ijp) = dPhidy(ijp) + Dy
    dPhidy(ijn) = dPhidy(ijn) + Dy 

    dPhidz(ijp) = dPhidz(ijp) + Dz 
    dPhidz(ijn) = dPhidz(ijn) + Dz 

    !***------------------------------------
    
  ! !
  ! ! **Extended interpolation molecule: neighbours of neighbours**
  ! !
  !   nb_loop2: do k = ioffset( neighbour(i) ), ioffset( neighbour(i)+1 )-1

  !     ijn = ja(k)

  !     if (ijn == ijp) cycle nb_loop2 ! nemoj ownera dirati

  !     ! Paznja kada je ja(k) = diag( neighbour(i) ) ijn ce uzeti samo cell_neighboura, 
  !     ! tako da ovaj loop moze uraditi isto sto i gornji, pa je u toj verziji
  !     ! racunanja gradijenata ovaj loop gore nepotreban. Eto simplifikacije koda...

  !     w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

  !     Dx = w * ( xc(ijn)-xc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
  !     Dy = w * ( yc(ijn)-yc(ijp) ) * ( Phi(ijn)-Phi(ijp) )
  !     Dz = w * ( zc(ijn)-zc(ijp) ) * ( Phi(ijn)-Phi(ijp) )

  !     dPhidx(ijp) = dPhidx(ijp) + Dx 
  !     dPhidx(ijn) = dPhidx(ijn) + Dx

  !     b2(ijp) = b2(ijp) + Dy
  !     b2(ijn) = b2(ijn) + Dy 

  !     b3(ijp) = b3(ijp) + Dz 
  !     b3(ijn) = b3(ijn) + Dz  

  !   enddo nb_loop2
                                                            
  enddo     


  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells+i

    w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

    Dx = w*(Phi(ijn)-Phi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = w*(Phi(ijn)-Phi(ijp))*(yf(iface)-yc(ijp))
    Dz = w*(Phi(ijn)-Phi(ijp))*(zf(iface)-zc(ijp)) 

    dPhidx(ijp) = dPhidx(ijp) + Dx
    dPhidy(ijp) = dPhidy(ijp) + Dy
    dPhidz(ijp) = dPhidz(ijp) + Dz


  end do

  ! Calculate gradient

  do inp=1,numCells 

    b1 = dPhidx(inp)
    b2 = dPhidy(inp)
    b3 = dPhidz(inp)

    dPhidx(inp) = b1*Dmat(1,inp) - b2*Dmat(2,inp) +  b3*Dmat(3,inp) 
    dPhidy(inp) = b1*Dmat(4,inp) - b2*Dmat(5,inp) -  b3*Dmat(6,inp) 
    dPhidz(inp) = b1*Dmat(7,inp) - b2*Dmat(8,inp) +  b3*Dmat(9,inp) 

  enddo 

  
end subroutine


! Gauss gradients
subroutine grad_gauss(u,dUdx,dUdy,dUdz)
!
!***********************************************************************
!
!   Purpose:
!     Calculates cell centered gradient using Gauss theorem.
!   Parameters:
!     u - field, the gradient of which we are looking for
!     dUdx,dUdy,dUdz - arrays where the gradient components are stored
!
!     gauss gradient rule:
!     ------->                                 ->
!     grad(u) = 1/vol * sum_{i=1}^{i=nf} (u)_f*sf
!     where:
!     grad(u) - cell centered gradient vector
!     (u)_f   - face interpolated value of scalar u
!     vol     - cell volume
!     sf      - cell face area vector
!     nf      - number of faces in a cell
!
!***********************************************************************
!

  implicit none

  ! Arguments
  real(dp), dimension(numTotal), intent(in) :: u
  real(dp), dimension(numCells), intent(inout) :: dUdx,dUdy,dUdz

  ! Local
  integer :: i,ijp,ijn,ijb,lc,iface
  real(dp) :: volr,dfxe, dfye, dfze
  real(dp), dimension(:), allocatable :: dfxo,dfyo,dfzo


  allocate( dfxo(numCells) )
  allocate( dfyo(numCells) )
  allocate( dfzo(numCells) )

  ! Initialize
  dfxo = 0.0_dp
  dfyo = 0.0_dp
  dfzo = 0.0_dp


  ! Start iterative calculation of gradients
  do lc = 1,2
    
    ! Initialize new gradient
    dUdx = 0.0_dp
    dUdy = 0.0_dp
    dUdz = 0.0_dp

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      call gradco(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                  u, dfxo, dfyo, dfzo, dfxe, dfye, dfze)

      ! Accumulate contribution at cell center and neighbour
      dudx(ijp) = dudx(ijp) + dfxe
      dudy(ijp) = dudy(ijp) + dfye
      dudz(ijp) = dudz(ijp) + dfze
       
      dudx(ijn) = dudx(ijn) - dfxe
      dudy(ijn) = dudy(ijn) - dfye
      dudz(ijn) = dudz(ijn) - dfze

    enddo

    ! Contribution from boundaries

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i

      call gradbc( arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp) )

    enddo


    ! Calculate gradient components at cv-centers
    do ijp=1,numCells

      volr = 1.0_dp/vol(ijp)

      dudx(ijp) = dudx(ijp)*volr
      dudy(ijp) = dudy(ijp)*volr
      dudz(ijp) = dudz(ijp)*volr

    enddo

    ! Set old gradient = new gradient for the next iteration
    if(lc.lt.2) then
      dfxo = dudx(1:numCells)
      dfyo = dudy(1:numCells)
      dfzo = dudz(1:numCells)
    endif

  enddo ! lc-loop

  deallocate( dfxo )
  deallocate( dfyo )
  deallocate( dfzo )

end subroutine



! Corrected Gauss gradients
subroutine grad_gauss_corrected(u,dudx,dudy,dudz)
!
!***********************************************************************
!
!   Purpose:
!     Calculates cell centered gradient using Gauss theorem.
!   Parameters:
!     u - field, the gradient of which we are looking for
!     dudx - arrays where the gradient components are stored
!
!     Gauss gradient rule:
!     ------->                                 ->
!     grad(u) = 1/vol * sum_{i=1}^{i=nf} (u)_f*sf
!     where:
!     grad(u) - cell centered gradient vector
!     (u)_f   - face interpolated value of scalar u
!     vol     - cell volume
!     sf      - cell face area vector
!     nf      - number of faces in a cell
!
!***********************************************************************
!
  use types
  use geometry

  implicit none

  ! Arguments
  real(dp), dimension(numTotal), intent(in) :: u
  real(dp), dimension(numCells), intent(inout) :: dudx,dudy,dudz

  ! Local
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: volr,dfxe, dfye, dfze
  real(dp), dimension(:), allocatable :: dfxo,dfyo,dfzo


  allocate( dfxo(numCells) )
  allocate( dfyo(numCells) )
  allocate( dfzo(numCells) )

  ! Initialize gradient with lsq gradient
  dfxo = dudx(1:numCells)
  dfyo = dudy(1:numCells)
  dfzo = dudz(1:numCells)
  
  ! Initialize new gradient
  dUdx = 0.0_dp
  dUdy = 0.0_dp
  dUdz = 0.0_dp

  ! Calculate terms integrated over surfaces

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)

    call gradco(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                u, dfxo, dfyo, dfzo, dfxe, dfye, dfze)

    ! Accumulate contribution at cell center and neighbour
    dudx(ijp) = dudx(ijp) + dfxe
    dudy(ijp) = dudy(ijp) + dfye
    dudz(ijp) = dudz(ijp) + dfze
     
    dudx(ijn) = dudx(ijn) - dfxe
    dudy(ijn) = dudy(ijn) - dfye
    dudz(ijn) = dudz(ijn) - dfze

  enddo

  ! Contribution from boundaries
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijb = numCells + i

    call gradbc( arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp) )

  enddo


  ! Calculate gradient components at cv-centers
  do ijp=1,numCells

      volr = 1.0_dp/vol(ijp)

      dudx(ijp) = dudx(ijp)*volr
      dudy(ijp) = dudy(ijp)*volr
      dudz(ijp) = dudz(ijp)*volr

  enddo

end subroutine


subroutine gradco(ijp,ijn, &
                  xfc,yfc,zfc,sx,sy,sz,fif, &
                  fi,dfxo,dfyo,dfzo,dfxe,dfye,dfze)
!
!***********************************************************************
!
!     This routine calculates contribution to the gradient
!     vector of a scalar FI at the CV center, arising from
!     an inner cell face (cell-face value of FI times the 
!     corresponding component of the surface vector).
!
!***********************************************************************
!
  use types
  use parameters
  use geometry

  implicit none

  integer,    intent(in) :: ijp,ijn
  real(dp), intent(in) :: xfc,yfc,zfc
  real(dp), intent(in) :: sx,sy,sz
  real(dp), intent(in) :: fif
  real(dp), dimension(numTotal), intent(in) :: fi
  real(dp), dimension(numCells), intent(in) :: dfxo,dfyo,dfzo
  real(dp), intent(out)  :: dfxe,dfye,dfze


  real(dp) :: xi,yi,zi
  real(dp) :: dfxi,dfyi,dfzi
  real(dp) :: fie
  real(dp) :: fxn,fxp

  !
  ! Coordinates of point on the line connecting center and neighbor,
  ! old gradient vector components interpolated for this location.

  fxn = fif 
  fxp = 1.0d0-fxn

  xi = xc(ijp)*fxp+xc(ijn)*fxn
  yi = yc(ijp)*fxp+yc(ijn)*fxn
  zi = zc(ijp)*fxp+zc(ijn)*fxn

  dfxi = dfxo(ijp)*fxp+dfxo(ijn)*fxn
  dfyi = dfyo(ijp)*fxp+dfyo(ijn)*fxn
  dfzi = dfzo(ijp)*fxp+dfzo(ijn)*fxn

  ! Value of the variable at cell-face center
  fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

  ! (interpolated mid-face value)x(area)
  dfxe = fie*sx
  dfye = fie*sy
  dfze = fie*sz

end subroutine

subroutine gradbc(sx,sy,sz,fi,dfx,dfy,dfz)
!***********************************************************************
!
! This routine calculates the contribution of a 
! boundary cell face to the gradient at CV-center.
!
!***********************************************************************
  use types

  implicit none

  real(dp), intent(in) :: sx,sy,sz
  real(dp), intent(in) :: fi
  real(dp), intent(inout)  :: dfx,dfy,dfz

  dfx = dfx + fi*sx
  dfy = dfy + fi*sy
  dfz = dfz + fi*sz

end subroutine


!******************************************************************************
!
subroutine sngrad_scalar_field(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
                               Phi, dPhidx, dPhidy, dPhidz, nrelax, approach, &
                               dfixi, dfiyi, dfizi, &
                               dfixii, dfiyii, dfizii)
!
!******************************************************************************
!
!  Surface normal gradient with non-orthogonal correction done in two
!  possible ways - either by skewness correction of intersection point
!  offset.
!
!  Check out reference paper: 
!    Mirkov, Rasuo, Kenjeres, JCP, Vol. 287, 2015.
!
!******************************************************************************
!
  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numTotal), intent(in) :: Phi
  real(dp), dimension(numTotal), intent(in) :: dPhidx,dPhidy,dPhidz
  integer, intent(in) :: nrelax
  character(len=12) :: approach
  real(dp), intent(out) :: dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii
!
! Locals
!
  real(dp) :: are,vole
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ixi1,ixi2,ixi3
  real(dp) :: dpn,costheta,costn

  real(dp) :: d1x,d1y,d1z
  real(dp) :: d2x,d2y,d2z
  real(dp) :: fxp,fxn

  real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  real(dp) :: nablaFIxdnnp,nablaFIxdppp


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! Unit vectors of the face normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  costn = 1.0_dp

  if(nrelax == 1) then
    ! Minimal correction: nrelax = +1 :
    costn = costheta
  elseif(nrelax == 0) then
    ! Orthogonal correction: nrelax =  0 : 
    costn = 1.0_dp
  elseif(nrelax == -1) then
    ! Over-relaxed approach: nrelax = -1 :
    costn = 1.0_dp/costheta  
  endif

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Interpolate gradients defined at CV centers to faces
  dfixi = dPhidx(ijp)*fxp+dPhidx(ijn)*fxn
  dfiyi = dPhidy(ijp)*fxp+dPhidy(ijn)*fxn
  dfizi = dPhidz(ijp)*fxp+dPhidz(ijn)*fxn


  !-- Skewness correction -->
  if (adjustl(approach) == 'skewness') then

    ! Overrelaxed correction vector d2, where s=dpn+d2
    d1x = costn
    d1y = costn
    d1z = costn

    d2x = xpn*costn
    d2y = ypn*costn
    d2z = zpn*costn

    !.....du/dx_i interpolated at cell face:
    dfixii = dfixi*d1x + arx/vole*( Phi(ijn)-Phi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
    dfiyii = dfiyi*d1y + ary/vole*( Phi(ijn)-Phi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
    dfizii = dfizi*d1z + arz/vole*( Phi(ijn)-Phi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z )

  ! |-- Intersection point offset and skewness correction -->


  elseif (adjustl(approach) == 'offset') then

    ! Find points P' and Pj'
    xpp=xf-(xf-xc(ijp))*nxx
    ypp=yf-(yf-yc(ijp))*nyy 
    zpp=zf-(zf-zc(ijp))*nzz

    xep=xf-(xf-xc(ijn))*nxx 
    yep=yf-(yf-yc(ijn))*nyy 
    zep=zf-(zf-zc(ijn))*nzz     

    xpnp = xep-xpp 
    ypnp = yep-ypp 
    zpnp = zep-zpp

    volep = arx*xpnp+ary*ypnp+arz*zpnp

   ! Overrelaxed correction vector d2, where S=dpn+d2
    d1x = costn
    d1y = costn
    d1z = costn
    
    xpnp = xpnp*costn
    ypnp = ypnp*costn
    zpnp = zpnp*costn

    ! The cell face interpolated gradient (d phi / dx_i)_j:
    ! Nonorthogonal corrections:          ___
    ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
    ! And:                                ___
    ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
    nablaFIxdnnp = dPhidx(ijn)*(xep-xc(ijn))+dPhidy(ijn)*(yep-yc(ijn))+dPhidz(ijn)*(zep-zc(ijn))
    nablaFIxdppp = dPhidx(ijp)*(xpp-xc(ijp))+dPhidy(ijp)*(ypp-yc(ijp))+dPhidz(ijp)*(zpp-zc(ijp))

    dfixii = dfixi*d1x + arx/volep*( Phi(ijn)+nablaFIxdnnp-Phi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
    dfiyii = dfiyi*d1y + ary/volep*( Phi(ijn)+nablaFIxdnnp-Phi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
    dfizii = dfizi*d1z + arz/volep*( Phi(ijn)+nablaFIxdnnp-Phi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 
  !-- Uncorrected -->
  elseif (adjustl(approach) == 'uncorrected') then
    
    dfixii = dfixi
    dfiyii = dfiyi
    dfizii = dfizi

  endif


end subroutine


end module