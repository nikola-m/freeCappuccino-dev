module node_interpolation

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Functions for field interpolation from cell centers to nodes.
!
!  Discussion:
!
!    There are various ways to interpolate to a mesh node(vertex) from surrounding cells. 
!    One of the ways is inverse distace weighted interpolation, commonly used in CFD.
!    This is iportant for e.g. postprocessing, when values are interpolated to nodes
!    for visualisation in Paraview.
!
!    At the moment we have two options:
!      1) inverse distance weights, c2n = 'ids'
!      2) pseudo laplacian weights, c2n = 'psl' 
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use geometry, only: numTotal,numCells,numNodes,numFaces,numInnerFaces,owner,neighbour,xc,yc,zc,xf,yf,zf,x,y,z,node,nnodes,tiny

  implicit none



  integer :: nCellToNode, nFaceToNode

  real(dp), dimension(:), allocatable, save :: wCellToNode,wFaceToNode

  character( len=3 ) :: c2n = 'psl'  ! Cell-to-node method: 'ids', 'psl'

  private 

  public :: c2n ! param
  public :: set_interpolation_weights, interpolate_to_nodes 


contains


subroutine set_interpolation_weights
!
!  Purpose: 
!  
!   Driver subroutine for cell-to-node interpolation weights calculation.
!
!  Discussion:
!
!   At the moment we have two options:
!     1) inverse distance weights
!     2) pseudo laplacian weights
!
  ! integer :: i

  if (c2n == 'ids') then
    call weights_inverse_distance

  elseif (c2n == 'psl') then
    call weights_pseudolaplacian

  else

    write(*,*) " "
    write(*,*) "  Unknown method for cell to node interpolation!"
    write(*,*) "  => Setting default: inverse distance interpolation."
    call weights_inverse_distance

  endif 

  write(*,*) " "
  write(*,'(2a)') "  Cell-to-node interpolation weights are set. Method: ",c2n

  ! do i=1,nFaceToNode
  !   write(*,*) wFaceToNode(i)
  ! enddo
  ! stop


end subroutine



subroutine weights_inverse_distance
!
!  Purpose: 
!
!    Computes weights for interpolation from surrounding cells and boundary faces.
!    It allocates and sets values for module private arrays wCellToNode and wFaceToNode.
!
!  Discussion:
!
!    The implementation uses inverse distance weights for interpolation.
!    Weights are calculated as: wi = li/lsum 
!    (inverse distance from cell to node )/ (sum of inverse distances from all surrounding cells to this node)
!
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    Aug 10th, 2022.
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    -
!
!  Parameters:
!
!    -
!

  implicit none

!
! Local variables
!
integer :: i,k,ijp,ijn,ino,iface
real(dp) :: dcp
real(dp), dimension(:), allocatable :: dsum


! dsum vazi za svaki node pa ga ima numNodes;
allocate( dsum(numNodes) )

dsum = 0.0_dp ! Initialize sum of cell-to-node distances

nCellToNode = 0 ! integer, module private variable, initialized here

do iface=1,numInnerFaces

  ijp = owner(iface)
  ijn = neighbour(iface)

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    ! Owner cell center to node
    dcp = sqrt( ( xc(ijp)-x(ino) )**2 + ( yc(ijp)-y(ino) )**2 + ( zc(ijp)-z(ino) )**2 )

    dsum( ino ) =  dsum( ino ) + 1./dcp ! \\izvrti dok ne dobije sve

    ! usput mogu da vidim koliko ima ovih kuraca dcp pa da znam taj broj
    nCellToNode = nCellToNode + 1

    !---

    ! Neighbor cell center to node
    dcp = sqrt( ( xc(ijn)-x(ino) )**2 + ( yc(ijn)-y(ino) )**2 + ( zc(ijn)-z(ino) )**2 )

    dsum( ino ) =  dsum( ino ) + 1./dcp 

    nCellToNode = nCellToNode + 1

  enddo
  
enddo 

!
! Calculate weights now.
!

!\\ alociraj sad ovo sranje kada znas kolko je dugacko
allocate( wCellToNode( nCellToNode ) )

wCellToNode = 0.

i = 0

do iface=1,numInnerFaces

  ijp = owner(iface)
  ijn = neighbour(iface)

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    ! Owner cell center to node
    dcp = sqrt( ( xc(ijp)-x(ino) )**2 + ( yc(ijp)-y(ino) )**2 + ( zc(ijp)-z(ino) )**2 )

    i = i + 1

    wCellToNode(i) = 1. / ( dcp*dsum(ino) ) 

    ! Neighbor cell center to node
    dcp = sqrt( ( xc(ijn)-x(ino) )**2 + ( yc(ijn)-y(ino) )**2 + ( zc(ijn)-z(ino) )**2 )

    i = i + 1

    ! wi = li/lsum
    wCellToNode(i) = 1. / ( dcp*dsum(ino) ) 

  enddo
  
enddo 

!
! Boundary faces loop
!
nFaceToNode = 0

dsum = 0.0_dp ! clean arrays becasue we have contaminated some nodes that are on the boundary - from inner faces that touch boundary

do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    ! Face center to node
    dcp = sqrt( ( xf(iface)-x(ino) )**2 + ( yf(iface)-y(ino) )**2 + ( zf(iface)-z(ino) )**2 )

    dsum( ino ) =  dsum( ino ) + 1./dcp 

    ! usput mogu da vidim koliko ima ovih kuraca dcp pa da znam taj broj
    nFaceToNode = nFaceToNode + 1

  enddo

enddo

!\\ alociraj sad ovo sranje kada znas kolko je dugacko
allocate( wFaceToNode( nFaceToNode ) )

wFaceToNode = 0.

i = 0

do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    ! Face center to node
    dcp = sqrt( ( xf(iface)-x(ino) )**2 + ( yf(iface)-y(ino) )**2 + ( zf(iface)-z(ino) )**2 )

    i = i + 1

    wFaceToNode(i) = 1. / ( dcp*dsum(ino) ) 

  end do

enddo

deallocate(dsum)

end subroutine




subroutine weights_pseudolaplacian
!
! Purpose:
! Calculating weights for cell center to cell node weighted interpolation
! used for e.g. Gauss Node based gradients or for post processing.
!
! Description:
! Implementation based on pseudo_laplacian weighted interpolation method,
! proposed by Holmes and Connell, which is 2nd order accurate on unstructured meshes.
!
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    Aug 10th, 2022.
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    Holmes and Conell 
!
!  Parameters:
!
!    -
!

  implicit none

!
! Local variables
!
integer :: i,k,ijp,ijn,ino,iface
real(dp),dimension(:), allocatable :: rx,ry,rz,ixx,iyy,izz,ixy,ixz,iyz,wsum
real(dp) , parameter :: one = 1.0_dp
real(dp) :: D,lamx,lamy,lamz



allocate( rx(numNodes),ry(numNodes),rz(numNodes) )
allocate( ixx(numNodes),iyy(numNodes),izz(numNodes) )
allocate( ixy(numNodes),ixz(numNodes),iyz(numNodes) )

! Lambdas (lamx,lamy,lamz) are valid for one mesh node 
! and surrounding set of cell centers.
! We reinitialize these below because they participate in
! calculation of lambdas.
rx = 0.
ry = 0.
rz = 0.

ixx = 0.
iyy = 0.
izz = 0.

ixy = 0.
ixz = 0. 
iyz = 0.


nCellToNode = 0

do iface=1,numInnerFaces

  ijp = owner(iface)
  ijn = neighbour(iface)

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    rx( ino ) = rx( ino ) + ( xc( ijp ) - x( ino ) )
    ry( ino ) = ry( ino ) + ( yc( ijp ) - y( ino ) )
    rz( ino ) = rz( ino ) + ( zc( ijp ) - z( ino ) )

    ixx( ino ) = ixx( ino ) + ( xc( ijp ) - x( ino ) )**2
    iyy( ino ) = iyy( ino ) + ( yc( ijp ) - y( ino ) )**2
    izz( ino ) = izz( ino ) + ( zc( ijp ) - z( ino ) )**2

    ixy( ino ) = ixy( ino ) + ( xc( ijp ) - x( ino ) )*( yc( ijp ) - y( ino ) )
    ixz( ino ) = ixz( ino ) + ( xc( ijp ) - x( ino ) )*( zc( ijp ) - z( ino ) )
    iyz( ino ) = iyz( ino ) + ( yc( ijp ) - y( ino ) )*( zc( ijp ) - z( ino ) )


    nCellToNode = nCellToNode + 1

    !---

    rx( ino ) = rx( ino ) + ( xc( ijn ) - x( ino ) )
    ry( ino ) = ry( ino ) + ( yc( ijn ) - y( ino ) )
    rz( ino ) = rz( ino ) + ( zc( ijn ) - z( ino ) )

    ixx( ino ) = ixx( ino ) + ( xc( ijn ) - x( ino ) )**2
    iyy( ino ) = iyy( ino ) + ( yc( ijn ) - y( ino ) )**2
    izz( ino ) = izz( ino ) + ( zc( ijn ) - z( ino ) )**2

    ixy( ino ) = ixy( ino ) + ( xc( ijn ) - x( ino ) )*( yc( ijn ) - y( ino ) )
    ixz( ino ) = ixz( ino ) + ( xc( ijn ) - x( ino ) )*( zc( ijn ) - z( ino ) )
    iyz( ino ) = iyz( ino ) + ( yc( ijn ) - y( ino ) )*( zc( ijn ) - z( ino ) )

    nCellToNode = nCellToNode + 1

  enddo
  
enddo 

!
! Calculate weights now.
!

allocate( wCellToNode( nCellToNode ) )

i = 0

do iface=1,numInnerFaces

  ijp = owner(iface)
  ijn = neighbour(iface)

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node


    D = ixx(ino)*(iyy(ino)*izz(ino)-iyz(ino)**2) &
       -ixy(ino)*(ixy(ino)*izz(ino)-ixz(ino)*iyz(ino)) &
       +ixz(ino)*(ixy(ino)*iyz(ino)-iyy(ino)*ixz(ino)) + tiny
    
    lamx = (-rx(ino)*(ixx(ino)*izz(ino)-iyz(ino)**2) & 
           + ry(ino)*(ixy(ino)*izz(ino)-ixz(ino)*iyz(ino)) &
           - rz(ino)*(ixy(ino)*iyz(ino)-iyy(ino)*ixz(ino)) ) / D

    lamy = ( rx(ino)*(ixy(ino)*izz(ino)-ixz(ino)*iyz(ino)) &
           - ry(ino)*(ixx(ino)*izz(ino)-ixz(ino)**2)  &
           + rz(ino)*(ixx(ino)*iyz(ino)-ixy(ino)*ixz(ino)) ) / D
    
    lamz = (-rx(ino)*(ixy(ino)*iyz(ino)-iyy(ino)*ixz(ino)) &
           + ry(ino)*(ixx(ino)*iyz(ino)-ixy(ino)*ixz(ino)) &
           - rz(ino)*(ixx(ino)*iyy(ino)-ixy(ino)**2)  ) / D


    ! **Pseudo-Laplacian weight for each cell-to-node interpolation**
    ! NOTE: If you want to change mesh during simulation (moving mesh),
    ! (xn-xp), (yn-yp) and (zn-zp) changes over time, then
    ! save lamx,lamy,lamz, and update weight (wc2n) after each mesh update.

    i = i + 1

    wCellToNode(i) = 1.0_dp + lamx*( xc( ijp ) - x( ino ) ) &
                            + lamy*( yc( ijp ) - y( ino ) ) &
                            + lamz*( zc( ijp ) - z( ino ) )

    wCellToNode(i) = min( abs( wCellToNode(i) ), one) ! /// Diskutabilno, proveriti testovima


    i = i + 1

    wCellToNode(i) = 1.0_dp + lamx*( xc( ijn ) - x( ino ) ) &
                            + lamy*( yc( ijn ) - y( ino ) ) &
                            + lamz*( zc( ijn ) - z( ino ) )

    wCellToNode(i) = min( abs( wCellToNode(i) ), one) ! /// Diskutabilno, proveriti testovima

  enddo
  
enddo 

!
! Accumulate sum of weights for a node
! 
  allocate( wSum( numNodes ) ) 
  wSum = 0.

  i = 0

  do iface=1,numInnerFaces

    ijp = owner(iface)
    ijn = neighbour(iface)

    do k=1,nnodes(iface)

      ino = node(k,iface) !\\ true global index of the node

      i = i + 1

      wSum(ino) = wSum(ino) + wCellToNode(i) 

      i = i + 1

      wSum(ino) = wSum(ino) + wCellToNode(i) 

    enddo
    
  enddo 


!
! wi/sum_i weights as in Rausch Batina Yang, NASA TM 104039, Fig.1
!
  i = 0

  do iface=1,numInnerFaces

    ijp = owner(iface)
    ijn = neighbour(iface)

    do k=1,nnodes(iface)

      ino = node(k,iface) !\\ true global index of the node

      i = i + 1

      wCellToNode(i) =  wCellToNode(i) / wSum(ino)

      i = i + 1

      wCellToNode(i) = wCellToNode(i) / wSum(ino)

    enddo
    
  enddo 



!
! Boundary faces loop
!

rx = 0.
ry = 0.
rz = 0.

ixx = 0.
iyy = 0.
izz = 0.

ixy = 0.
ixz = 0. 
iyz = 0.

nFaceToNode = 0


do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ global index of the node

    ! Face center to node
    rx( ino ) = rx( ino ) + ( xf( iface ) - x( ino ) )
    ry( ino ) = ry( ino ) + ( yf( iface ) - y( ino ) )
    rz( ino ) = rz( ino ) + ( zf( iface ) - z( ino ) )

    ixx( ino ) = ixx( ino ) + ( xf( iface ) - x( ino ) )**2
    iyy( ino ) = iyy( ino ) + ( yf( iface ) - y( ino ) )**2
    izz( ino ) = izz( ino ) + ( zf( iface ) - z( ino ) )**2

    ixy( ino ) = ixy( ino ) + ( xf( iface ) - x( ino ) )*( yf( iface ) - y( ino ) )
    ixz( ino ) = ixz( ino ) + ( xf( iface ) - x( ino ) )*( zf( iface ) - z( ino ) )
    iyz( ino ) = iyz( ino ) + ( yf( iface ) - y( ino ) )*( zf( iface ) - z( ino ) )

    nFaceToNode = nFaceToNode + 1

  enddo

enddo


allocate( wFaceToNode( nFaceToNode ) )

wFaceToNode = 0.

i = 0

do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ global index of the node

    D = ixx(ino)*(iyy(ino)*izz(ino)-iyz(ino)**2) &
       -ixy(ino)*(ixy(ino)*izz(ino)-ixz(ino)*iyz(ino)) &
       +ixz(ino)*(ixy(ino)*iyz(ino)-iyy(ino)*ixz(ino)) + tiny
    
    lamx = (-rx(ino)*(ixx(ino)*izz(ino)-iyz(ino)**2) & 
           + ry(ino)*(ixy(ino)*izz(ino)-ixz(ino)*iyz(ino)) &
           - rz(ino)*(ixy(ino)*iyz(ino)-iyy(ino)*ixz(ino)) ) / D 

    lamy = ( rx(ino)*(ixy(ino)*izz(ino)-ixz(ino)*iyz(ino)) &
           - ry(ino)*(ixx(ino)*izz(ino)-ixz(ino)**2)  &
           + rz(ino)*(ixx(ino)*iyz(ino)-ixy(ino)*ixz(ino)) ) / D
    
    lamz = (-rx(ino)*(ixy(ino)*iyz(ino)-iyy(ino)*ixz(ino)) &
           + ry(ino)*(ixx(ino)*iyz(ino)-ixy(ino)*ixz(ino)) &
           - rz(ino)*(ixx(ino)*iyy(ino)-ixy(ino)**2)  ) / D
    
    ! Face center to node

    i = i + 1

    wFaceToNode(i) = 1.0_dp + lamx*( xf( iface ) - x( ino ) ) &
                            + lamy*( yf( iface ) - y( ino ) ) &
                            + lamz*( zf( iface ) - z( ino ) )

    wFaceToNode(i) = min( abs(wFaceToNode(i)), one) ! /// Diskutabilno, proveriti testovima

  end do

enddo


!
! Accumulate sum weights at boundary nodes
!

wSum = 0. ! clean

i = 0

do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface) ! Global index of the node

    i = i + 1

    wSum(ino) = wSum(ino) + wFaceToNode(i) 

  end do

enddo

!
! wi/sum_i weights as in Rausch Batina Yang, NASA TM 104039, Fig.1
!

i = 0

do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface)  ! Global index of the node

    i = i + 1

    wFaceToNode(i) = wFaceToNode(i) / wSum(ino)  

  end do

enddo


deallocate( rx,ry,rz )
deallocate( ixx,iyy,izz )
deallocate( ixy,ixz,iyz )
deallocate( wsum )

end subroutine



subroutine interpolate_to_nodes(phiC, phiN)
!
!  Purpose: 
!
!    Performs interpolation from cells and boundary faces to node.
!
!  Discussion:
!
!    This is general procedure in the sense that it doesn't matter how you set weights,
!    inverse distance weighted or pseudo laplacian.
!
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    Aug 10th, 2022.
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Reference:
!
!    -
!
!  Parameters:
!
!    Input, real(dp), phiC( numTotal ), Field variable which we interpolate to nodes.
!    Output, real(dp), phiN( numNodes ), nodal, interpolated, values of the field.
!

  implicit none

!
! Parameters
!
  real(dp), dimension(numTotal), intent(in) :: phiC
  real(dp), dimension(numNodes), intent(inout) :: phiN
!
! Local variables
!
  integer :: i,k,ib,ijp,ijn,ino,iface


phiN = 0.

i = 0

do iface=1,numInnerFaces

  ijp = owner(iface)
  ijn = neighbour(iface)

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    i = i + 1

    phiN(ino) = phiN(ino) + wCellToNode(i)*phiC(ijp)

    i = i + 1

    phiN(ino) = phiN(ino) + wCellToNode(i)*phiC(ijn)


  enddo
  
enddo 


!\\ We already interpolated to nodes at boundary from cells - what to do with that??

!
! Maybe first clean what we touched?
!
do iface=numInnerFaces+1,numFaces

  do k=1,nnodes(iface)

    ino = node(k,iface) ! global index of the node

    phiN(ino) = 0.

  end do

enddo

!
! Boundary faces loop and interpolate from faces
!
i = 0

do iface=numInnerFaces+1,numFaces

  ib = numCells+(iface-numInnerFaces) ! \\ boundary values are from [numCells+1 : numTotal]; (iface-numInnerFaces) is a boundary face.

  do k=1,nnodes(iface)

    ino = node(k,iface) !\\ true global index of the node

    i = i + 1

    phiN(ino) = phiN(ino) + wFaceToNode(i)*phiC(ib)

  end do

enddo

end subroutine



end module

