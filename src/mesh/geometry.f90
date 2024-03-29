
module geometry
!
! Purpose: Module for manipulation of unstructured meshes.
!
use types
use parameters, only: mesh_format
use utils, only: get_unit, file_row_count, r8vec_print_some, i4vec_print2

implicit none

! NOTE:
!***
! In variable arrays, field variables defined at cell centres are written in positions from 1 to numCells, after that we write
! variable values for boundary faces from numCells+1 to numTotal.  
!***

! General mesh data
integer :: numNodes         ! no. of nodes in mesh
integer :: numCells         ! no. of cells in the mesh
integer :: numFaces         ! no. of INNER+BOUNDARY faces in the mesh
integer :: numInnerFaces    ! no. of INNER cells faces in the mesh
integer :: numBoundaryFaces ! self explanatory
integer :: numPeriodic      ! number of periodic boundary faces (not counting them twice)
integer :: numTotal         ! number of volume field values + number of boundary field values numCells+numBoundaryFaces

! To define boundary regions
integer :: numBoundaries
integer, dimension(:), allocatable :: nfaces,startFace,iBndValueStart
character(len=30), dimension(:), allocatable :: bcname
character(len=15), dimension(:), allocatable :: bctype
real(dp), dimension(:,:), allocatable :: bcDFraction

integer :: nwal ! Total no. of boundary faces of type 'wall'
integer :: nsym ! Total no. of boundary faces of type 'symmetry'
integer :: ninl ! No. of inlet boundary faces
integer :: nout ! No. of outlet boundary faces
integer :: nper ! Stores the number of periodic boundary pairs, which will be treated as a single inner boundary               

integer, parameter :: nomax = 6 !!! Max number of nodes in face - determines size of some arrays, just change this if necessary. 
real(dp), parameter :: tiny = 1e-30

! Mesh file units
integer :: points_file, cells_file, faces_file, owner_file, neighbour_file, boundary_file, size_file

integer, parameter :: interpolation_coeff_variant = 2 ! (1,2,3) look at the code below.


! Mesh geometry

! Geometry parameters defined for mesh nodes [1:numNodes]
real(dp), dimension(:), allocatable :: x,y,z  ! Coordinates of mesh nodes 

! Geometry parameters defined cellwise [1:numCells]:
real(dp), dimension(:), allocatable :: xc,yc,zc     ! Coordinates of cell centers
real(dp), dimension(:), allocatable :: vol          ! Cell volume
real(dp), dimension(:), allocatable :: wallDistance ! Distance to the nearest wall - needed in some turb. models

! Geometry parameters defined for all (inner+boundary) cell faces [1:numFaces]
real(dp), dimension(:), allocatable :: arx, ary, arz   ! Cell face area x-, y- and z-component 
real(dp), dimension(:), allocatable :: xf, yf, zf      ! Coordinates of cell face center

! Geometry parameters defined for all inner cell faces [1:numInnerFaces]
real(dp), dimension(:), allocatable :: xpp, ypp, zpp   ! Coordinates of auxilliary points - owner cell 
real(dp), dimension(:), allocatable :: xnp, ynp, znp   ! Coordinates of auxilliary points - neighbour cell 
real(dp), dimension(:), allocatable :: facint          ! Interpolation factor 
real(dp), dimension(:), allocatable :: Df              ! Diffusion geometric coefficient

! Geometry parameters defined for boundary faces

real(dp), dimension(:), allocatable :: srds,dns    ! srds = |are|/|dns|, dns = normal distance to cell center from face |dpn*face_normal|
real(dp), dimension(:), allocatable :: srdw,dnw    ! srdw = |are|/|dnw|, dnw = normal distance to cell center from face |dpn*face_normal|

! Mesh topology information - connectivity of cells trough faces
integer, dimension(:), allocatable :: owner     ! Index of the face owner cell
integer, dimension(:), allocatable :: neighbour ! Index of the neighbour cell  - it shares the face with owner

integer, dimension(:,:), allocatable :: node  ! It will store global node numbers of face vertices
integer, dimension(:), allocatable :: nnodes  ! no. of nodes in specific face face size[numFaces].

! For periodic boundaries
integer :: startFaceTwin(3) ! we expect max triple periodic domains..


public 

contains


subroutine read_mesh
!
! Simple driver subroutine that calls one of the mesh reading and processing subrouines below, based on the specified format. 
!
  if (mesh_format == 'nativeMesh') then

    call read_mesh_native

  elseif(mesh_format == 'foamMesh') then

    call read_mesh_openfoam

  else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Fatal error!'
      write ( *, '(a)' ) '  Please specify mesh format type in the input file!'
      write ( *, '(a)' ) "  Supported options are 'nativeMesh' and 'foamMesh'"
      stop
  end if

  if (nPer > 0) then ! we have periodic boundaries: map faces at opposite boundaries.
    call face_mapping
  endif


end subroutine


subroutine read_mesh_native
!
!  Purpose: 
!    Reads mesh in native freeCappuccino format and provides basic geometrical mesh description.
!
!  Description:
!    Mesh is described in files of polyMesh format almost identical to the one used by OpenFoam.
!    The difference is that present one is adapted for Fortran numbering (starting from 1), etc.
!    Mesh files are 'points', 'faces', 'owner, 'neighbour', 'boundary'.
!    There is also a 'cells' file which is used only for Paraview postprocessing, i.e. when we write 
!    Praview unstructured .vtu files, when we need to give cell connectivity. Also, such information
!    is usually given by mesh generators, so we don't want to waste this information.
!    For the purposes of Finite Volume Computation, as it is implemented in freeCappuccino, the cell
!    connectivity is not necessary.
!
!  Date:
!    26/11/2015; Jan-2019; Dec-2019.
!
!  Author:
!    Nikola Mirkov (E-mail: largeddysimulation@gmail.com)
!
  implicit none

  ! Locals
  integer :: i,k
  integer :: ib,iwall,isym
  integer :: iface
  integer :: inp,inn,ijp

  character(len=80) :: line_string

  ! Parameters
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: third = 1./3._dp

  real(dp) :: px,py,pz, qx,qy,qz, nx,ny,nz, cx,cy,cz
  real(dp) :: riSi
  real(dp) :: are,areasum
  real(dp) :: xpn,ypn,zpn
  real(dp) :: xjp,yjp,zjp
  real(dp) :: dpn,djn,djp
  real(dp) :: nxf,nyf,nzf

  ! Array for temporary storage of doubles
  real(dp), dimension(:), allocatable :: r8tmp
 

!******************************************************************************
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'.
!..............................................................................

  call get_unit( points_file )
  open( unit = points_file,file = 'polyMesh/points' )
  rewind points_file

  ! call get_unit( cells_file )
  ! open( unit = cells_file,file = 'polyMesh/cells' )
  ! rewind cells_file

  call get_unit( faces_file )
  open( unit = faces_file, file='polyMesh/faces' )
  rewind faces_file

  call get_unit( owner_file )
  open( unit = owner_file, file='polyMesh/owner' )
  rewind owner_file

  call get_unit( neighbour_file )
  open( unit = neighbour_file, file='polyMesh/neighbour' )
  rewind neighbour_file

  call get_unit( boundary_file )
  open( unit = boundary_file, file='polyMesh/boundary' )
  rewind boundary_file

  call get_unit( size_file )
  open( unit = size_file, file='polyMesh/size' )
  rewind size_file
!
! > Read boundary conditions file. 
!

!   Discussion: 
!   Boundary conditions file consists of header and numBoundaries number of subsequent lines each containing:  
!   boundary condition given name, bc type, number of faces belonging to that bc and starting face for that bc.
!   Possible boundary contiions types are: inlet, outlet, symmtery, wall, pressure, etc. 
!   Please check are all of these implemented, because this is still in the development. Contact me if you have any questions.
!

  ! Number of rows in the file excluding #comment in header to find the number of prescribed boundaries.
  call file_row_count ( boundary_file, numBoundaries )

  read(boundary_file,'(a)') line_string ! First line is header.
  ! read(boundary_file,*) numBoundaries ! it doesn't need to read the number of boundaries because of the above

  ! Allocate 
  allocate ( bcname(numBoundaries) )
  allocate ( bctype(numBoundaries) )
  allocate ( nFaces(numBoundaries) )
  allocate ( startFace(numBoundaries) )
  allocate ( iBndValueStart(numBoundaries))

  ! Control how we update fields on boundaries using Dirichlet Fraction; equals one for Dirichlet boundary
  allocate ( bcDFraction(12,numBoundaries)) ! lets say for 12 fields and for every boundary

  ! Initialize all boundaries to Neumann, therefore Dirichlet fraction is zero.
  ! This is changed to one if we stumble upon "Dirichlet" keyword when reading initial fields in 0/ folder.
  ! checkout field_initialization module.
  bcDFraction = 0.0_dp 

  nwal = 0
  nsym = 0
  ninl = 0
  nper = 0        ! this one stores the number of periodic boundary pairs, which will be treated as a single inner boundary
  numPeriodic = 0 ! this one stores the number of cell faces at periodic boundaries

  do i=1,numBoundaries
    read(boundary_file,*) bcname(i), bctype(i), nfaces(i), startFace(i)

    ! We need total number of some bctype faces like wall and symmetry to make things easier for bc implementation
    ! so lets count them. 
    if( bctype(i) == 'wall') then
      nwal = nwal + nfaces(i)

    elseif( bctype(i) == 'symmetry') then
      nsym = nsym + nfaces(i)

    elseif( bctype(i) == 'inlet') then
      ninl = ninl + nfaces(i)

    elseif( bctype(i) == 'outlet') then
      nout = nout + nfaces(i)

    elseif( bctype(i) == 'periodic') then
      nper = nper + 1 ! number of periodic boundary pairs increases by one
      numPeriodic = numPeriodic + nfaces(i) ! number of periodic faces increases by nfaces(i)
      backspace( boundary_file ) ! Move one line back.
      ! Repeat reading the same line, but get additional information that should be available for perodic boundary,
      ! namely we have starting face for 'twin' boundary of this periodic boundary.  
      read(boundary_file,*) bcname(i), bctype(i), nfaces(i), startFace(i), startFaceTwin(nper)

    endif

  enddo  


!
! > Find out numCells, numNodes, numFaces, numInnerFaces, etc.
!

  read(size_file, *) numNodes
  read(size_file, *) numCells
  read(size_file, *) numInnerFaces
  read(size_file, *) numBoundaryFaces
  read(size_file, *) numFaces


  ! Number of boundary faces
  if (numBoundaryFaces /= numFaces - numInnerFaces) then
    write(*,*) "Something wrong with the mesh: numBoundaryFaces /= numFaces - numInnerFaces!"
    stop
  endif

  ! Size of arrays storing variables numCells+numBoundaryFaces
  numTotal = numCells + numBoundaryFaces

  ! NOTE:
  ! The variable values defined at boundary faces are stored in variable arrays after numCells. The length of variable arrays therefore becomes numTotal = numCells + numBoundaryFaces
  ! It is therefore important to do the following also:
  ! Define where are the boundary field values located in the variable array, for each boundary region.
  do i=1,numBoundaries
      iBndValueStart(i) = numCells + (startFace(i) - numInnerFaces)
  enddo 

!
! > Write report on mesh size into log file
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Mesh data: '

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of nodes, numNodes = ', numNodes

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of cells, numCells = ', numCells

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of cell-faces, numFaces = ', numFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of inner cell-faces, numInnerFaces = ', numInnerFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of cell-faces on boundary, numBoundaryFaces = ', numBoundaryFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary information (bcname, bctype, nFaces, startFace):'
  write ( *, '(a)' ) ' '
  
  do i=1,numBoundaries
    write(*,'(2x,a,1x,a,1x,2(i0,1x))') bcname(i), bctype(i), nfaces(i) ,startFace(i)
  enddo  
  write ( *, '(a)' ) ' '



!******************************************************************************
! > Allocate arrays for Mesh description
!..............................................................................

  ! Nodal coordinates
  allocate ( x(numNodes) )
  allocate ( y(numNodes) )
  allocate ( z(numNodes) )

  ! Coordinates of cell centers
  allocate ( xc(numCells) )
  allocate ( yc(numCells) )
  allocate ( zc(numCells) )

  ! Cell volumes
  allocate ( vol(numCells) )

  ! Face area vector components
  allocate ( arx(numFaces) )
  allocate ( ary(numFaces) )
  allocate ( arz(numFaces) )

  ! Coordinates of cell face centers
  allocate ( xf(numFaces) )
  allocate ( yf(numFaces) ) 
  allocate ( zf(numFaces) )

  ! Interpolation factor for inner faces
  allocate ( facint(numInnerFaces) ) 

  ! Indices of owner cell (inner+boundary faces) and indices of neighbours for every inner cell-face.
  allocate ( owner(numFaces) )
  allocate ( neighbour(numInnerFaces) )

  ! Array stroring the walue of distance to the nearest wall, needed only in some turb. models (maybe allocate only when needed)
  allocate ( wallDistance(numCells) )

  ! We need this only for wall and symmetry BCs, so it is not bad to have nsym and nwal parameters, which count home many of them 
  ! are there.
  allocate ( dns(nsym) )      
  allocate ( srds(nsym) )  

  allocate ( dnw(nwal) )      
  allocate ( srdw(nwal) )   

  ! allocate ( dnw(numBoundaryFaces) )      
  ! allocate ( srdw(numBoundaryFaces) )     

  allocate( Df(numInnerFaces) )
  
!******************************************************************************
! > Read and process Mesh files 
!..............................................................................

  ! The 'points' file
  do i=1,numNodes
    read(points_file,*) x(i),y(i),z(i) 
  end do


  ! The 'owner' file
  do i=1,numFaces
    read(owner_file,*) owner(i)
  end do


  ! The 'neighbour' file
  do i=1,numInnerFaces
    read(neighbour_file,*) neighbour(i)
  end do

  ! The 'faces' file

  ! read( faces_file, * )  line_string ! we don't need this info on first line

  ! Allocate tmp array of number of nodes for each face - nnodes, and node array
  allocate( nnodes(numFaces) )
  allocate( node(nomax,numFaces) ) ! We should put nomax instead of 4, because like this there can be max 4 nodes on a face.

  do iface=1,numFaces
    read( faces_file, * ) nnodes(iface),(node(k,iface), k=1,nnodes(iface) )
  enddo  

  ! Allocate tmp array of doubles 
  allocate( r8tmp(numCells))
  r8tmp = 0.0_dp


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Cell volumes, cell centers and cell face centers
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do iface=1,numFaces

    inp = owner(iface)

    ! Initialize total area of polygon.
    areasum = 0.0_dp

    ! We decompose a polygon to a series of triangles, all having first node in common.
    do i=1, nnodes(iface)-2 

      ! Vectors to vertices
      ! 2-1
      px = x( node(i+1,iface) )-x( node(1,iface) )
      py = y( node(i+1,iface) )-y( node(1,iface) )
      pz = z( node(i+1,iface) )-z( node(1,iface) )
      ! 3-1
      qx = x( node(i+2,iface) )-x( node(1,iface) )
      qy = y( node(i+2,iface) )-y( node(1,iface) )
      qz = z( node(i+2,iface) )-z( node(1,iface) )


      call triangle_area_vector( px,py,pz, qx,qy,qz, nx,ny,nz )

      !
      ! > Cell-face area vector components (Area vector lies in direction of face normal)
      ! 

      arx(iface) = arx(iface) + nx
      ary(iface) = ary(iface) + ny
      arz(iface) = arz(iface) + nz

      ! Face center for a triangle
      cx = third*( x( node(i+2,iface) ) + x( node(i+1,iface) ) + x( node(1,iface) ) )
      cy = third*( y( node(i+2,iface) ) + y( node(i+1,iface) ) + y( node(1,iface) ) )
      cz = third*( z( node(i+2,iface) ) + z( node(i+1,iface) ) + z( node(1,iface) ) ) 

      !
      ! > Cell-face centroid components - accumulation stage
      !

      are = sqrt(nx**2 + ny**2 + nz**2)

      xf(iface) = xf(iface) + (are*cx)
      yf(iface) = yf(iface) + (are*cy)
      zf(iface) = zf(iface) + (are*cz)

      ! Accumulate triangle areas to get total area of the polygon
      areasum = areasum + are


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! > Compute cell volumes and cell centers
      !
      ! We compute cell volumes by aplying divergence theorem to the position vector,
      ! see eq. (5) in [1].
      ! Cell center coordinates of an arbitrary polyhedron are computed using eq.(15) of ref. [1].
      !
      ! [1] Z.J. Wang - Improved Formulation for Geometric Properties of Arbitrary Polyhedra
      !     AIAA Journal, Vol. 37, No. 10, October 1999
      !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      riSi = ( cx*nx + cy*ny + cz*nz )

      vol(inp) = vol(inp) + third * riSi

      xc(inp) = xc(inp) + 0.75_dp * riSi * cx
      yc(inp) = yc(inp) + 0.75_dp * riSi * cy
      zc(inp) = zc(inp) + 0.75_dp * riSi * cz

      ! We use r8tmp array to store accumulated denominator
      r8tmp(inp) = r8tmp(inp) + riSi


      if ( iface <= numInnerFaces ) then 

        inn = neighbour(iface)

        riSi = -( cx*nx + cy*ny + cz*nz )

        vol(inn) = vol(inn) + third * riSi

        xc(inn) = xc(inn) + 0.75_dp * riSi * cx
        yc(inn) = yc(inn) + 0.75_dp * riSi * cy
        zc(inn) = zc(inn) + 0.75_dp * riSi * cz

        ! We use r8tmp array to store accumulated denominator
        r8tmp(inn) = r8tmp(inn) + riSi

      endif

    enddo

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! > Cell-face centroid components - final
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xf(iface) = xf(iface) / areasum
    yf(iface) = yf(iface) / areasum
    zf(iface) = zf(iface) / areasum
  
  enddo

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Cell centroid components - final
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Do one more loop over cell volumes to divide accumulated cell center values by
  ! denominator accumulated in wallDistance array for convenience.
  do inp=1,numCells
    xc(inp) = xc(inp) / r8tmp(inp)
    yc(inp) = yc(inp) / r8tmp(inp)
    zc(inp) = zc(inp) / r8tmp(inp)
  enddo

  ! Thank you.
  deallocate(r8tmp)


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Interpolation factor
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ( interpolation_coeff_variant == 1 ) then

    !
    ! > Interpolation factor > inner faces - Variant 1.
    !
    do iface=1,numInnerFaces

      inp = owner(iface)
      inn = neighbour(iface)

      xpn = xc(inn)-xc(inp)
      ypn = yc(inn)-yc(inp)
      zpn = zc(inn)-zc(inp)

      dpn = sqrt( xpn**2 + ypn**2 + zpn**2 ) 

      !
      ! > > Intersection point j' of line connecting centers with cell face, we are taking only three points assuming that other are co-planar
      !
      call find_intersection_point(&
                                   ! plane defined by three face vertices:
                                   x( node(1,iface) ), y( node(1,iface) ), z( node(1,iface) ), &
                                   x( node(2,iface) ), y( node(2,iface) ), z( node(2,iface) ), &
                                   x( node(3,iface) ), y( node(3,iface) ), z( node(3,iface) ), &
                                   ! line defined by cell center and neighbour center:
                                   xc(inp),yc(inp),zc(inp), &
                                   xc(inn),yc(inn),zc(inn), &
                                   ! intersection point (output):
                                   xjp,yjp,zjp &
                                  )
      xpn = xjp - xc(inp)
      ypn = yjp - yc(inp)
      zpn = zjp - zc(inp)

      djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

      ! Interpolation factor |P j'|/|P Pj| where P is cell center, Pj neighbour cell center and j' intersection point.
      facint(iface) = djn / dpn 
      
    enddo

  elseif ( interpolation_coeff_variant == 2 ) then

    !
    ! > Interpolation factor > inner faces - Variant 2.
    !
    do iface=1,numInnerFaces

      inp = owner(iface)
      inn = neighbour(iface)

      xpn = xf(iface) - xc(inp)
      ypn = yf(iface) - yc(inp)
      zpn = zf(iface) - zc(inp)

      djp = sqrt( xpn**2 + ypn**2 + zpn**2 )

      xpn = xf(iface) - xc(inn)
      ypn = yf(iface) - yc(inn)
      zpn = zf(iface) - zc(inn)

      djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

      ! Interpolation factor |jP|/(|jP|+|jN|)  where P is cell center, N neighbour cell center and j is face centroid.
      facint(iface) =  djp / ( djp + djn ) 
      
    enddo

  else ! (interpolation_coef_variant == 3)

    !
    ! > Interpolation factor > inner faces - Variant 3.
    !
    do iface=1,numInnerFaces

      inp = owner(iface)
      inn = neighbour(iface)

      xpn = xf(iface) - xc(inp)
      ypn = yf(iface) - yc(inp)
      zpn = zf(iface) - zc(inp)

      djp = abs(arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn)

      xpn = xc(inn)-xc(inp)
      ypn = yc(inn)-yc(inp)
      zpn = zc(inn)-zc(inp)

      dpn = abs(arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn)

      ! Interpolation factor = |Sf*djP|/|Sf*dPN|, where Sf face area (normal to face) vector,
      ! dPN is distance vector between cell centers, djP distance vector from face centroid j to 
      ! cell center P.
      facint(iface) =  djp / dpn 
      
    enddo  

  endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Calculate and store few recurrent geometric parameters:
  !   > Df - diffusion coefficient:
  !   _    _     _    _
  !   Sf . Sf / (Sf . dpn)
  !   It is recalculated frequently so it is best to keep it in memory
  !   to speed up the simulation.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  do iface=1,numInnerFaces

    inp = owner(iface)
    inn = neighbour(iface)

    ! Distance vector between cell centers
    xpn = xc(inn)-xc(inp)
    ypn = yc(inn)-yc(inp)
    zpn = zc(inn)-zc(inp)

    ! Instead to enable reuse and efficiency, first store dot product of vectors Sf*Sf:
    are = arx(iface)*arx(iface) + ary(iface)*ary(iface) + arz(iface)*arz(iface)

    ! Create Df coefficient used for implicit Lalplacian discretizaton, where 
    !      _    _     _    _
    ! Df = Sf . Sf / (Sf . dpn)
    Df(iface) = are /(arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn)

  !   ! Sqrt to get |Sf|
  !   are = sqrt( are  )

  !   ! Unit vectors of the normal
  !   nxx=arx(iface)/are
  !   nyy=ary(iface)/are
  !   nzz=arz(iface)/are

  !   ! Values at points p' and e' due to non-orthogonality. 
  !   xpp(inp)=xf-(xf-xc(inp))*nxx
  !   ypp(inp)=yf-(yf-yc(inp))*nyy
  !   zpp(inp)=zf-(zf-zc(inp))*nzz

  !   xnp(inn)=xf-(xf-xc(inn))*nxx
  !   ynp(inn)=yf-(yf-yc(inn))*nyy
  !   znp(inn)=zf-(zf-zc(inn))*nzz

  enddo



  ! Note: don't deallocate we use it for cell to node interpolation weights calculation (node_interpolation module) - only then deallocate.
  ! deallocate(nnodes)
  ! deallocate(node)  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Loop over boundary faces to calculate 
  ! normal distance from face center to cell center of the first cell layer - dnw,
  ! and cell face area divided by norma distance to the cell center - srdw.
  ! These are used in implementation of various boundary conditions.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iWall = 0
  iSym = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        iSym = iSym + 1

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
        dns(iSym) = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        ! Cell face area divided by distance to the cell center
        srds(iSym) = are/dns(iSym)

      end do

    elseif ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        iWall = iWall + 1

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
        dnw(iWall) = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        ! Cell face area divided by distance to the cell center
        srdw(iWall) = are/dnw(iWall)

      enddo

    endif 

  enddo

  !
  ! Same coefficients for all boundary faces - previous will be deprecated
  !

  ! do i=1,numBoundaryFaces
   
  !   iface=numInnerFaces+i
  !   inp = owner(iface)

  !   ! Face area 
  !   are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

  !   ! Face normals
  !   nxf = arx(iface)/are
  !   nyf = ary(iface)/are
  !   nzf = arz(iface)/are

  !   ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
  !   dnw(i) = (xf(iface)-xc(inp))*nxf + (yf(iface)-yc(inp))*nyf + (zf(iface)-zc(inp))*nzf

  !   ! Cell face area divided by distance to the cell center
  !   srdw(i) = are/dnw(i)


  ! enddo


!******************************************************************************
! > Report on geometrical quantities > I will leave this for debug purposes.
!..............................................................................

  ! write ( *, '(a)' ) ' '
  ! write ( *, '(a)' ) '  Cell data: '

  ! call r8vec_print_some ( numCells, vol, 1, numCells, &
  !     '  First 10 elements of cell volumes array:' )

  ! call r8vec_print_some ( numCells, xc, 1, 10, &
  !     '  First 10 elements of cell x-centers array:' )

  ! call r8vec_print_some ( numCells, yc, 1, 10, &
  !     '  First 10 elements of cell y-centers array:' )

  ! call r8vec_print_some ( numCells, zc, 1, 10, &
  !     '  First 10 elements of cell z-centers array:' )

  ! write ( *, '(a)' ) ' '
  ! write ( *, '(a)' ) '  Face data: '

  ! call i4vec_print2 ( 10, owner, neighbour, '  First 10 lines of owner and neighbour arrays:' )

  ! call r8vec_print_some ( numFaces, arx, 1, 10, &
  !     '  First 10 elements of Arx array:' )

  ! call r8vec_print_some ( numFaces, ary, 1, 10, &
  !     '  First 10 elements of Ary array:' )

  ! call r8vec_print_some ( numFaces, arz, 1, 10, &
  !     '  First 10 elements of Arz array:' )

  !   call r8vec_print_some ( numFaces, xf, 1, 10, &
  !     '  First 10 elements of xf array:' )

  ! call r8vec_print_some ( numFaces, yf, 1, 10, &
  !     '  First 10 elements of yf array:' )

  ! call r8vec_print_some ( numFaces, zf, 1, 10, &
  !     '  First 10 elements of zf array:' )

  ! call r8vec_print_some ( numInnerFaces, facint, 1, 10, &
  !     '  First 10 elements of interpolation factor (facint) array:' )

  write ( *, '(a)' ) ' '
  
!
!  > CLOSE polyMesh format file: 'points', 'faces', 'owner', 'neighbour', 'boundary'.
!
  close ( points_file )
  close ( faces_file )
  close ( owner_file )
  close ( neighbour_file)
  close ( boundary_file)
  close ( size_file )
!+-----------------------------------------------------------------------------+


end subroutine



subroutine read_mesh_openfoam
!
!  Description:
!    Calculates basic geometrical quantities of numerical mesh
!    defined in this module and needed for FVM computation.
!    Mesh is described in files of polyMesh format.
!    Mesh files are 'points', 'faces', 'owner, 'neighbour', 'boundary'
!
!  Date:
!    26/11/2015; Jan-2019
!
!  Author:
!    Nikola Mirkov nmirkov@vin.bg.ac.rs
!
  implicit none

  ! Locals
  integer :: i,j,k,l,itmp
  integer :: ib,iwall,isym
  integer :: iface
  integer :: inp,inn,ijp

  character(len=1) :: ch
  character(len=20) :: char_string,char_string2
  character(len=80) :: line_string

  integer, dimension(nomax) :: nodeOF  ! It will store global node numbers of cell vertexes
  integer :: nnodesOF              ! no. of nodes in face

  ! Parameters
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: third = 1./3._dp

  real(dp) :: px,py,pz, qx,qy,qz, nx,ny,nz, cx,cy,cz
  real(dp) :: riSi
  real(dp) :: are,areasum
  real(dp) :: xpn,ypn,zpn
  real(dp) :: xjp,yjp,zjp
  real(dp) :: dpn,djn,djp
  real(dp) :: nxf,nyf,nzf

  ! Array for temporary storage of doubles
  real(dp), dimension(:), allocatable :: r8tmp
 

!******************************************************************************
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'.
!..............................................................................

  call get_unit( points_file )
  open( unit = points_file,file = 'polyMesh/points' )
  rewind points_file

  call get_unit( faces_file )
  open( unit = faces_file, file='polyMesh/faces' )
  rewind faces_file

  call get_unit( owner_file )
  open( unit = owner_file, file='polyMesh/owner' )
  rewind owner_file

  call get_unit( neighbour_file )
  open( unit = neighbour_file, file='polyMesh/neighbour' )
  rewind neighbour_file

  call get_unit( boundary_file )
  open( unit = boundary_file, file='polyMesh/boundary' )
  rewind boundary_file



!
! > Read boundary conditions file. 
!

!
!   Boundary conditions file consists of header and numBoundaries number of subsequent lines each containing:  
!   boundary condition given name, bc type, number of faces belonging to that bc and starting face for that bc.
!   Possible boundary contiions types are: inlet, outlet, symmtery, wall, pressure, etc. 
!   Please check are all of these implemented, because this is still in the development. Contact me if you have any questions.
!
 

  ! Number of rows in the file excluding #comment in header to find the number of prescribed boundaries
  call file_row_count ( boundary_file, numBoundaries )

  read(boundary_file,'(a)') line_string ! Firts line is header.
  ! read(boundary_file,*) numBoundaries ! it doesn't need to read the number of boundaries because of the above

  ! Allocate 
  allocate ( bcname(numBoundaries) )
  allocate ( bctype(numBoundaries) )
  allocate ( nFaces(numBoundaries) )
  allocate ( startFace(numBoundaries) )
  allocate ( iBndValueStart(numBoundaries))

  ! Control how we update fields on boundaries using Dirichlet Fraction; equals one for Dirichlet boundary
  allocate ( bcDFraction(12,numBoundaries)) ! lets say for 12 fields and for every boundary

  ! Initialize all boundaries to Neumann, therefore Dirichlet fraction is zero.
  ! This is changed to one if we stumble upon "Dirichlet" keyword when reading initial fields in 0/ folder.
  ! checkout field_initialization module.
  bcDFraction = 0.0_dp 

  nwal = 0
  nsym = 0
  ninl = 0
  nper = 0        ! this one strores the number of periodic boundary pairs, which will be treated as a single inner boundary
  numPeriodic = 0 ! this one stores the number of cell faces at periodic boundaries
  
  do i=1,numBoundaries
    read(boundary_file,*) bcname(i), bctype(i), nfaces(i) ,startFace(i)

    ! We need total number of some bctype faces like wall and symmetry to make things easier for bc implementation
    ! so lets count them. 
    if( bctype(i) == 'wall') then
      nwal = nwal + nfaces(i)

    elseif( bctype(i) == 'symmetry') then
      nsym = nsym + nfaces(i)

    elseif( bctype(i) == 'inlet') then
      ninl = ninl + nfaces(i)

    elseif( bctype(i) == 'outlet') then
      nout = nout + nfaces(i)

    elseif( bctype(i) == 'periodic') then
      nper = nper + 1 ! number of periodic boundary pairs increases by one
      numPeriodic = numPeriodic + nfaces(i) ! number of periodic faces increases by nfaces(i)
      backspace( boundary_file ) ! Move one line back.
      ! Repeat reading but find some additional information is available for perodic boundary
      ! namely we have starting face for 'twin' boundary of this periodic boundary.  
      read(boundary_file,*) bcname(i), bctype(i), nfaces(i), startFace(i), startFaceTwin(nper)

    endif
  enddo  


!******************************************************************************
! > Find out numNodes, numFaces, numInnerFaces, etc.
!..............................................................................

  !
  ! Code here is tested for OpenFOAM version 4.0, if they don't change polyMesh format this should be OK.
  !       

  !
  ! The 'owner' file. After reading this we will have numNodes, numCells, numInnerFaces, numFaces
  !

  k=0
  l=0
  char_string = ' '

  owner_header_loop: do

    ! Trying to find the line with mesh size info, and read numCells
    read(owner_file,*) char_string,line_string
    ! write(*,*) "owner file ",char_string

    if (char_string == 'note') then

      ! Do probing for nPoints:
      do j=1,len_trim(line_string)-5
        if (line_string(j:j+6)=='nPoints') then
          k=j+8
        endif
        if (line_string(j:j+5)=='nCells') then
          l=j-2
          exit
        endif
      end do
      read(line_string(k:l),*) numNodes
      ! write(*,*) numNodes

      ! Do probing for nCells:
      do j=1,len_trim(line_string)-5
        if (line_string(j:j+5)=='nCells') then
          k=j+7
        endif
        if (line_string(j:j+5)=='nFaces') then
          l=j-2
          exit
        endif
      end do
      read(line_string(k:l),*) numCells
      ! write(*,*) numCells

      ! Do probing for nFaces:
      do j=1,len_trim(line_string)-5
        if (line_string(j:j+5)=='nFaces') then
          k=j+7
        endif
        if (line_string(j:j+13)=='nInternalFaces') then
          l=j-2
          exit
        endif
      end do
      read(line_string(k:l),*) numFaces
      ! write(*,*) numFaces

      ! Do probing for nInternalFaces:
      do j=1,len_trim(line_string)-5
        ! write(*,*)line_string(j:j+15)

        if (line_string(j:j+14)=='nInternalFaces:') then
          read(line_string(j+15:),*) numInnerFaces
          ! write(*,*) numInnerFaces
          exit
          
        endif
      end do

      exit owner_header_loop
    endif 

  end do owner_header_loop

  rewind owner_file

  ! NOTE: Trying to acces number of faces data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numFaces.

  ch = ' '
  owner_loop: do 
    read(owner_file,*) ch
    if (ch == '(') then
      ! Return two lines
      backspace(owner_file)
      backspace(owner_file)
      exit owner_loop
    endif
  end do owner_loop

  read(owner_file,*) itmp
  if (itmp /= numFaces ) then
    write(*,*) "Error reading polyMesh format. numFaces value is not confirmed in body of the 'owner' file."
    stop
  endif
  read(owner_file,*) char_string ! reads "(" again


  !
  ! The 'points' file
  !

  ! NOTE: Trying to acces number of points data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numNodes.

  ch = ' '
  point_header_loop: do
    read(points_file,*) ch
    if (ch == '(') then
      ! Return two lines
      backspace(points_file)
      backspace(points_file)
      exit point_header_loop
    endif
  end do point_header_loop

  read(points_file,*) itmp
  if (itmp /= numNodes ) then
    write(*,*) "Error reading polyMesh format. numNodes value is not confirmed in 'points' file."
    stop
  endif
  read(points_file,*) char_string ! reads "("

  !
  ! The 'neighbour' file
  !

  ! NOTE: Trying to acces number of inner faces data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numInnerFaces.

  ch = ' '
  neighbour_header_loop: do
    read(neighbour_file,*) ch
    if (ch == '(') then
      ! Return two lines
      backspace(neighbour_file)
      backspace(neighbour_file)
      exit neighbour_header_loop
    endif
  end do neighbour_header_loop

  read(neighbour_file,*) itmp
  if (itmp /= numInnerFaces ) then
    write(*,*) "Error reading polyMesh format. numInnerFaces value is not confirmed in 'neighbour' file."
    stop
  endif
  read(neighbour_file,*) char_string ! reads "("

  !
  ! The 'faces' file
  !

  ! NOTE: Trying to acces number of faces data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numFaces.

  ch = ' '
  faces_header_loop: do
    read(faces_file,*) ch
    if (ch == '(') then
      backspace(faces_file)
      backspace(faces_file)
      exit faces_header_loop
    endif
  end do faces_header_loop

  read(faces_file, *) itmp
  if (itmp /= numFaces ) then
    write(*,*) "Error reading polyMesh format. numFaces value is not confirmed in 'faces' file."
    stop
  endif
  read(faces_file,*) char_string

  ! Number of boundary faces
  numBoundaryFaces = numFaces - numInnerFaces

  ! Size of arrays storing variables numCells+numBoundaryFaces
  numTotal = numCells + numBoundaryFaces

  ! NOTE:
  ! The variable values defined at boundary faces are stored in variable arrays after numCells. The length of variable arrays therefore becomes numTotal = numCells + numBoundaryFaces
  ! It is therefore important to do the following also.
  ! Define where are the boundary field values located in the variable array, for each boundary region
  do i=1,numBoundaries
      iBndValueStart(i) = numCells + (startFace(i) - numInnerFaces)
  enddo 

!
! > Write report on mesh size into log file
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Mesh data: '

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of nodes, numNodes = ', numNodes

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of cells, numCells = ', numCells

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of cell-faces, numFaces = ', numFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of inner cell-faces, numInnerFaces = ', numInnerFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a,i0)' ) '  Number of cell-faces on boundary, numBoundaryFaces = ', numBoundaryFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary information (bcname, bctype, nFaces, startFace):'
  write ( *, '(a)' ) ' '
  
  do i=1,numBoundaries
    write(*,'(2x,a,1x,a,1x,2(i0,1x))') bcname(i), bctype(i), nfaces(i) ,startFace(i)
  enddo  
  write ( *, '(a)' ) ' '



!******************************************************************************
! > Allocate arrays for Mesh description
!..............................................................................

  ! Nodal coordinates
  allocate ( x(numNodes) )
  allocate ( y(numNodes) )
  allocate ( z(numNodes) )

  ! Coordinates of cell centers
  allocate ( xc(numCells) )
  allocate ( yc(numCells) )
  allocate ( zc(numCells) )

  ! Cell volumes
  allocate ( vol(numCells) )

  ! Face area vector components
  allocate ( arx(numFaces) )
  allocate ( ary(numFaces) )
  allocate ( arz(numFaces) )

  ! Coordinates of cell face centers
  allocate ( xf(numFaces) )
  allocate ( yf(numFaces) ) 
  allocate ( zf(numFaces) )

  ! Interpolation factor for inner faces
  allocate ( facint(numInnerFaces) ) 

  ! Indices of owner cell (inner+boundary faces) and indices of neighbours for every inner cell-face.
  allocate ( owner(numFaces) )
  allocate ( neighbour(numInnerFaces) )

  ! Array stroring the walue of distance to the nearest wall, needed only in some turb. models (maybe allocate only when needed)
  allocate ( wallDistance(numCells) )

  ! We need this only for wall and symmetry BCs, so it is not bad to have nsym and nwal parameters, which count home many of them 
  ! are there.
  allocate ( dns(nsym) )      
  allocate ( srds(nsym) )  

  allocate ( dnw(nwal) )      
  allocate ( srdw(nwal) )   

  ! allocate ( dnw(numBoundaryFaces) )      
  ! allocate ( srdw(numBoundaryFaces) )     

  allocate( Df(numInnerFaces) )                           


!******************************************************************************
! > Read and process Mesh files 
!..............................................................................

  ! The 'points' file
  do i=1,numNodes
    ! char_string reads (number, char_string reads number), we have to strip off the brackets later.
    read(points_file,*) char_string,y(i),char_string2 

    ! Char to double conversion + stripping off the brackets:
    read(char_string(2:),*) x(i)
    read(char_string2(1:len_trim(char_string2)-1),*) z(i)  
  end do


  ! The 'owner' file
  do i=1,numFaces
    read(owner_file,*) owner(i)
    owner(i) = owner(i) + 1 ! fortran starts from 1
  end do



  ! The 'neighbour' file
  do i=1,numInnerFaces
    read(neighbour_file,*) neighbour(i)
    neighbour(i) = neighbour(i) + 1 ! fortran starts from 1
  end do


  ! The 'faces' file

  ! Allocate tmp array of number of nodes for each face - nnodes, and node array
  allocate( nnodes(numFaces) )
  allocate( node(nomax,numFaces) )

  do iface=1,numFaces

    call read_line_faces_file_polyMesh(faces_file,nnodesOF,nodeOF,nomax)

    nnodes(iface) = nnodesOF

    do k=1,nnodesOF
      node(k,iface) = nodeOF(k)
    enddo

  enddo 

   
  ! Allocate tmp array of doubles 
  allocate( r8tmp(numCells))
  r8tmp = 0.0_dp


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Cell volumes, cell centers and cell face centers
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  do iface=1,numFaces

    inp = owner(iface)

    ! Initialize total area of polygon.
    areasum = 0.0_dp

    ! We decompose a polygon to a series of triangles, all having first node in common.
    do i=1, nnodes(iface)-2 

      ! Vectors to vertices
      ! 2-1
      px = x( node(i+1,iface) )-x( node(1,iface) )
      py = y( node(i+1,iface) )-y( node(1,iface) )
      pz = z( node(i+1,iface) )-z( node(1,iface) )
      ! 3-1
      qx = x( node(i+2,iface) )-x( node(1,iface) )
      qy = y( node(i+2,iface) )-y( node(1,iface) )
      qz = z( node(i+2,iface) )-z( node(1,iface) )


      call triangle_area_vector( px,py,pz, qx,qy,qz, nx,ny,nz )

      !
      ! > Cell-face area vector components (Area vector lies in direction of face normal)
      ! 

      arx(iface) = arx(iface) + nx
      ary(iface) = ary(iface) + ny
      arz(iface) = arz(iface) + nz

      ! Face center for a triangle
      cx = third*( x( node(i+2,iface) ) + x( node(i+1,iface) ) + x( node(1,iface) ) )
      cy = third*( y( node(i+2,iface) ) + y( node(i+1,iface) ) + y( node(1,iface) ) )
      cz = third*( z( node(i+2,iface) ) + z( node(i+1,iface) ) + z( node(1,iface) ) ) 

      !
      ! > Cell-face centroid components - accumulation stage
      !

      are = sqrt(nx**2 + ny**2 + nz**2)

      xf(iface) = xf(iface) + (are*cx)
      yf(iface) = yf(iface) + (are*cy)
      zf(iface) = zf(iface) + (are*cz)

      ! Accumulate triangle areas to get total area of the polygon
      areasum = areasum + are


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! > Compute cell volumes and cell centers
      !
      ! We compute cell volumes by aaplying divergence theorem to the position vector,
      ! see eq. (5) in [1].
      ! Cell center coordinates of an arbitrary polyhedron are computed using eq.(15) of ref. [1].
      !
      ! [1] Z.J. Wang - Improved Formulation for Geometric Properties of Arbitrary Polyhedra
      !     AIAA Journal, Vol. 37, No. 10, October 1999
      !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      riSi = ( cx*nx + cy*ny + cz*nz )

      vol(inp) = vol(inp) + third * riSi 

      xc(inp) = xc(inp) + 0.75_dp * riSi * cx
      yc(inp) = yc(inp) + 0.75_dp * riSi * cy
      zc(inp) = zc(inp) + 0.75_dp * riSi * cz

      ! We use r8tmp array to store accumulated denominator
      r8tmp(inp) = r8tmp(inp) + riSi


      if ( iface <= numInnerFaces ) then 
        inn = neighbour(iface)

        riSi = -( cx*nx + cy*ny + cz*nz )

        vol(inn) = vol(inn) + third * riSi

        xc(inn) = xc(inn) + 0.75_dp * riSi * cx
        yc(inn) = yc(inn) + 0.75_dp * riSi * cy
        zc(inn) = zc(inn) + 0.75_dp * riSi * cz

        ! We use r8tmp array to store accumulated denominator
        r8tmp(inn) = r8tmp(inn) + riSi

      endif

    enddo

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! > Cell-face centroid components - final
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xf(iface) = xf(iface) / areasum
    yf(iface) = yf(iface) / areasum
    zf(iface) = zf(iface) / areasum
  
  enddo

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Cell centroid components - final
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Do one more loop over cell volumes to divide accumulated cell center values by
  ! denominator accumulated in wallDistance array for convenience.
  do inp=1,numCells
    xc(inp) = xc(inp) / r8tmp(inp)
    yc(inp) = yc(inp) / r8tmp(inp)
    zc(inp) = zc(inp) / r8tmp(inp)
  enddo

  ! Thank you.
  deallocate(r8tmp)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Interpolation factor
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !
  ! > Interpolation factor > inner faces
  !

  if ( interpolation_coeff_variant == 1 ) then

    do iface=1,numInnerFaces

      inp = owner(iface)
      inn = neighbour(iface)

      xpn = xc(inn)-xc(inp)
      ypn = yc(inn)-yc(inp)
      zpn = zc(inn)-zc(inp)

      dpn = sqrt( xpn**2 + ypn**2 + zpn**2 ) 

      !
      ! > > Intersection point j' of line connecting centers with cell face, we are taking only three points assuming that other are co-planar
      !
      call find_intersection_point(&
                                   ! plane defined by three face vertices:
                                   x( node(1,iface) ), y( node(1,iface) ), z( node(1,iface) ), &
                                   x( node(2,iface) ), y( node(2,iface) ), z( node(2,iface) ), &
                                   x( node(3,iface) ), y( node(3,iface) ), z( node(3,iface) ), &
                                   ! line defined by cell center and neighbour center:
                                   xc(inp),yc(inp),zc(inp), &
                                   xc(inn),yc(inn),zc(inn), &
                                   ! intersection point (output):
                                   xjp,yjp,zjp &
                                  )
      xpn = xjp - xc(inp)
      ypn = yjp - yc(inp)
      zpn = zjp - zc(inp)

      djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

      ! Interpolation factor |P Pj'|/|P Pj| where P is cell center, Pj neighbour cell center and j' intersection point.
      facint(iface) = djn / dpn 
      
    enddo

  elseif ( interpolation_coeff_variant == 2 ) then

    !
    ! > Interpolation factor > inner faces - Variant 2.
    !
    do iface=1,numInnerFaces

      inp = owner(iface)
      inn = neighbour(iface)

      xpn = xf(iface) - xc(inp)
      ypn = yf(iface) - yc(inp)
      zpn = zf(iface) - zc(inp)

      djp = sqrt( xpn**2 + ypn**2 + zpn**2 )

      xpn = xf(iface) - xc(inn)
      ypn = yf(iface) - yc(inn)
      zpn = zf(iface) - zc(inn)

      djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

      ! Interpolation factor |PjF|/|PF + PjF|  where P is cell center, Pj neighbour cell center and F is face centroid.
      facint(iface) =  djp / ( djp + djn ) 
      
    enddo

  else ! (interpolation_coef_variant == 3)

    !
    ! > Interpolation factor > inner faces - Variant 3.
    !
    do iface=1,numInnerFaces

      inp = owner(iface)
      inn = neighbour(iface)

      xpn = xf(iface) - xc(inp)
      ypn = yf(iface) - yc(inp)
      zpn = zf(iface) - zc(inp)

      djp = abs(arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn)

      xpn = xc(inn)-xc(inp)
      ypn = yc(inn)-yc(inp)
      zpn = zc(inn)-zc(inp)

      dpn = abs(arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn)

      ! Interpolation factor = |Sf*djP|/|Sf*dPN|, where Sf face area (normal to face) vector,
      ! dPN is distance vector between cell centers, djP distance vector from face centroid j to 
      ! cell center P.
      facint(iface) =  djp / dpn 
      
    enddo  


  endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Calculate and store few recurrent geometric parameters:
  !   > Df - diffusion coefficient:
  !   _    _     _    _
  !   Sf . Sf / (Sf . dpn)
  !   It is recalculated frequently so it is best to keep it in memory
  !   to speed up the simulation.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  do iface=1,numInnerFaces

    inp = owner(iface)
    inn = neighbour(iface)

    ! Distance vector between cell centers
    xpn = xc(inn)-xc(inp)
    ypn = yc(inn)-yc(inp)
    zpn = zc(inn)-zc(inp)

    ! Instead to enable reuse and efficiency, first store dot product of vectors Sf*Sf:
    are = arx(iface)*arx(iface) + ary(iface)*ary(iface) + arz(iface)*arz(iface)

    ! Create Df coefficient used for implicit Lalplacian discretizaton, where 
    !      _    _     _    _
    ! Df = Sf . Sf / (Sf . dpn)
    Df(iface) = are /(arx(iface)*xpn + ary(iface)*ypn + arz(iface)*zpn)

  !   ! Sqrt to get |Sf|
  !   are = sqrt( are  )

  !   ! Unit vectors of the normal
  !   nxx=arx(iface)/are
  !   nyy=ary(iface)/are
  !   nzz=arz(iface)/are

  !   ! Values at points p' and e' due to non-orthogonality. 
  !   xpp(inp)=xf-(xf-xc(inp))*nxx
  !   ypp(inp)=yf-(yf-yc(inp))*nyy
  !   zpp(inp)=zf-(zf-zc(inp))*nzz

  !   xnp(inn)=xf-(xf-xc(inn))*nxx
  !   ynp(inn)=yf-(yf-yc(inn))*nyy
  !   znp(inn)=zf-(zf-zc(inn))*nzz

  enddo



  ! Note: don't deallocate we use it for cell to node interpolation weights calculation (node_interpolation module) - only then deallocate.
  ! deallocate(nnodes)
  ! deallocate(node)  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Loop over boundary faces to calculate 
  ! normal distance from face center to cell center of the first cell layer - dnw,
  ! and cell face area divided by norma distance to the cell center - srdw.
  ! These are used in implementation of various boundary conditions.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iWall = 0
  iSym = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'symmetry') then

      ! Symmetry
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        iSym = iSym + 1

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
        dns(iSym) = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        ! Cell face area divided by distance to the cell center
        srds(iSym) = are/dns(iSym)

      end do

    elseif ( bctype(ib) == 'wall') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        iWall = iWall + 1

        ! Face area 
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are

        ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
        dnw(iWall) = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

        ! Cell face area divided by distance to the cell center
        srdw(iWall) = are/dnw(iWall)

      enddo

    endif 

  enddo


!******************************************************************************
! > Report on geometrical quantities > I will leave this for debug purposes.
!..............................................................................

  ! write ( *, '(a)' ) ' '
  ! write ( *, '(a)' ) '  Cell data: '

  ! call r8vec_print_some ( numCells, vol, 1, 10, &
  !     '  First 10 elements of cell volumes array:' )

  ! call r8vec_print_some ( numCells, xc, 1, 10, &
  !     '  First 10 elements of cell x-centers array:' )

  ! call r8vec_print_some ( numCells, yc, 1, 10, &
  !     '  First 10 elements of cell y-centers array:' )

  ! call r8vec_print_some ( numCells, zc, 1, 10, &
  !     '  First 10 elements of cell z-centers array:' )

  ! write ( *, '(a)' ) ' '
  ! write ( *, '(a)' ) '  Face data: '

  ! call i4vec_print2 ( 10, owner, neighbour, '  First 10 lines of owner and neighbour arrays:' )

  ! call r8vec_print_some ( numFaces, arx, 1, 10, &
  !     '  First 10 elements of Arx array:' )

  ! call r8vec_print_some ( numFaces, ary, 1, 10, &
  !     '  First 10 elements of Ary array:' )

  ! call r8vec_print_some ( numFaces, arz, 1, 10, &
  !     '  First 10 elements of Arz array:' )

  !   call r8vec_print_some ( numFaces, xf, 1, 10, &
  !     '  First 10 elements of xf array:' )

  ! call r8vec_print_some ( numFaces, yf, 1, 10, &
  !     '  First 10 elements of yf array:' )

  ! call r8vec_print_some ( numFaces, zf, 1, 10, &
  !     '  First 10 elements of zf array:' )

  ! call r8vec_print_some ( numInnerFaces, facint, 1, 10, &
  !     '  First 10 elements of interpolation factor (facint) array:' )

  write ( *, '(a)' ) ' '
  
!
!  > CLOSE polyMesh format file: 'points', 'faces', 'owner', 'neighbour', 'boundary'.
!
  close ( points_file )
  close ( faces_file )
  close ( owner_file )
  close ( neighbour_file)
  close ( boundary_file)
!+-----------------------------------------------------------------------------+

end subroutine



subroutine triangle_area_vector(px,py,pz,qx,qy,qz,nx,ny,nz)
!
! Ai, i=1,n are n triangular faces enclosing polyhedron P.
! Vertices of Ai are (ai,bi,ci),
!                                _    _   _    _     _  _    _  _
! Unit normal to P on each Ai is ni = ni/|ni|, ni = (bi-ai)x(ci-ai)
! and 
! _    _  _    _    _  _
! p = (bi-ai); q = (ci-ai)
!
! Finally:              _
! Area A of Ai, A = 1/2|ni|
!
! Sources: 
! [1] Dr Robert Nurnberg, Imperial College London, www.ma.ic.ac.uk/~rn/centroid.pdf
! [2] paulbourke.net/geometry/polygonmesh
!
  implicit none
  real(dp), intent(in) :: px,py,pz,qx,qy,qz
  real(dp), intent(inout) :: nx,ny,nz

  real(dp), parameter :: half = 0.5_dp

  ! Cross products for triangle surface vectors
  nx = half * (py*qz-pz*qy)
  ny = half * (pz*qx-px*qz)
  nz = half * (px*qy-py*qx)

end subroutine


subroutine find_intersection_point( &
!                     plane defined by three face corners:
                       x1,y1,z1,&
                       x2,y2,z2, &
                       x3,y3,z3, &
!                      line defined by cell center and neighbour center:
                       x4,y4,z4, &
                       x5,y5,z5, &
!                      intersection point (output):
                       xjp,yjp,zjp &
                       )
!
!***********************************************************************
! Find intersection point (pjp={xjp,yjp,zjp}) of 
! plane (defined by points p1={x1,y1,z1}, p2={x2,y2,z2} and p3={x3,y3,z3}),
! and line (defined by points p4={x4,y4,z4} and p5={x5,y5,z5}).
! The intersection point j' is not the face center j on non-orthogonal meshes.
! There is an "intersection point offset" |jj'| which determines the level
! of nonorthogonality.
!
!
!       |1  1  1  1 |     |1  1  1  0    |
! t = - |x1 x2 x3 x4|  /  |x1 x2 x3 x5-x4|  (mind the minus sign!)
!       |y1 y2 y3 y4| /   |y1 y2 y3 y5-y4|
!       |z1 z2 z3 z4|     |z1 z2 z3 z5-z4|
!
! And intersection point is given by:
! xj = x4 +(x5-x4)*t
! yj = y4 +(y5-y4)*t
! zj = z4 +(z5-z4)*t
!
!
! Nikola Mirkov 2016.
!
!***********************************************************************
  implicit none 
!
!***********************************************************************
!

  real(dp), intent(in) :: x1,y1,z1,&
                          x2,y2,z2, &
                          x3,y3,z3, &
                          x4,y4,z4, &
                          x5,y5,z5
  real(dp), intent(inout) :: xjp,yjp,zjp

  real(dp) :: t

 ! Produced by MATLAB symbolic tool.
 t =-(x2*(y3*z4-y4*z3)-x1*(y3*z4-y4*z3)-x3*(y2*z4-y4*z2)+x1*(y2*z4-y4*z2)+x3*(y1*z4-y4*z1)-x2* &
     (y1*z4-y4*z1)+x4*(y2*z3-y3*z2)-x1*(y2*z3-y3*z2)-x4*(y1*z3-y3*z1)+x2*(y1*z3-y3*z1)+x4* &
     (y1*z2-y2*z1)-x3*(y1*z2-y2*z1)) &
    /(x2*(y3*(z5-z4)-(y5-y4)*z3)-x1*(y3*(z5-z4)-(y5-y4)*z3)-x3*(y2*(z5-z4)-(y5-y4)*z2)+x1* &
     (y2*(z5-z4)-(y5-y4)*z2)+x3*(y1*(z5-z4)-(y5-y4)*z1)-x2*(y1*(z5-z4)-(y5-y4)*z1)+(x5-x4)* &
     (y2*z3-y3*z2)-(x5-x4)*(y1*z3-y3*z1)+(x5-x4)*(y1*z2-y2*z1) + tiny)

  xjp = x4 +(x5-x4)*t
  yjp = y4 +(y5-y4)*t
  zjp = z4 +(z5-z4)*t

end subroutine


subroutine read_line_faces_file_polyMesh(faces_file,nn,nod,nmax)
  implicit none
  integer, intent(in) :: faces_file
  integer, intent(in) :: nmax
  integer, intent(out) :: nn
  integer, dimension(nmax), intent(out) :: nod
  integer :: j,m,n
  character(len=15) :: char_string,char_string2

    nn = 0
    nod = 0

    ! Read how many nodes in face
    read(faces_file,'(a)') char_string ! e.g. 4(1 22 463 442)
    read(char_string(1:1),*)j          ! in this example j=4
    backspace(faces_file)              ! go back so you are able to read this line again (what can I do...)

      read(faces_file,*) char_string,nod(2:j-1),char_string2

      ! Char to double conversion:
      read(char_string(1:1),*)j                          ! number of vertrices
      read(char_string(3:),*) m                          ! first vertex
      read(char_string2(1:len_trim(char_string2)-1),*) n ! last vertex

      nn = j
      nod(1) = m + 1              ! <
      nod(2:j-1) = nod(2:j-1) + 1 ! < We are going from zero to one based numbering because of Fortran
      nod(nn) = n + 1             ! <

end subroutine


subroutine face_mapping
!
!  Purpose: 
!
!    Used to map faces at periodic boundaries, that is 
!    to find pairs of periodic faces at two opposite periodic boundaries.
!
!  Discussion:
!
!    The order of the faces at two twin periodic boundaries is not the same in a way
!    that the 17th face of one periodic boundary region is not the match for the 17th
!    face in its twin periodic boundary.
!
!    At the moment we find twin periodic faces based on criteria that the distance
!    vecttor between real twins is minimized.
!
!    # NOTE: The concept used below is not universal!
!    In a case when two mapped periodic boundaries are not made by translation,
!    then an alternative has to be programmed.
!    I will put this warning in some user's guide, when describing periodic BC.
!    
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    circa April 2020.
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

  integer :: i,ib,ijp,ijn,iface
  integer :: iPer, itemp, indx
  real(dp) :: dist,mindist
  real(dp), parameter :: tol = 1e-10
  integer, dimension(:), allocatable :: tmp,itmp
  real(dp), dimension(:), allocatable :: vtmp

  iPer = 0

  do ib=1,numBoundaries

    if (  bctype(ib) == 'periodic' ) then


      allocate( tmp( nfaces(ib) ), itmp( nfaces(ib) ), vtmp( nfaces(ib) ) )

      iPer = iPer + 1 ! found one pair of periodic boundaries, increase counter by one

      ! Loop trough faces
      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)

        mindist = 1e+30
        indx = 0

        do itemp=1,nfaces(ib)   

          ijn = owner( startFaceTwin(iPer) + itemp )

          ! Maybe in the future in the 'boundary' file, on 'periodic' bctype line,
          ! give vector, eg. '4.0d0, 0.0D0, 0.0D0' which translates that boundary
          ! to it's periodic twin, and that's how we can map faces...
          ! For now we use the concept of minimizing distance between periodic faces to find twins.
          dist = abs(xf(iface)-xc(ijn)) + abs(yf(iface)-yc(ijn)) + abs(zf(iface)-zc(ijn))

          if ( dist < mindist ) then
            mindist = dist
            tmp(i)  = ijn
            itmp(i) = itemp 
          endif

        enddo

      enddo

      do i=1,nfaces(ib)
        owner( startFaceTwin(iPer) + i ) = tmp(i)
      enddo

      ! Permute face center coordinates
      do i=1,nfaces(ib)
        vtmp(i) = xf( startFaceTwin(iPer) + itmp(i) )
      enddo
      do i=1,nfaces(ib)
        xf( startFaceTwin(iPer) + i ) = vtmp(i)
      enddo

      do i=1,nfaces(ib)
        vtmp(i) = yf( startFaceTwin(iPer) + itmp(i) )
      enddo
      do i=1,nfaces(ib)
        yf( startFaceTwin(iPer) + i ) = vtmp(i)
      enddo

      do i=1,nfaces(ib)
        vtmp(i) = zf( startFaceTwin(iPer) + itmp(i) )
      enddo 
      do i=1,nfaces(ib)
        zf( startFaceTwin(iPer) + i ) = vtmp(i)
      enddo

      ! Permute face normal vector components
      do i=1,nfaces(ib)
        vtmp(i) = arx( startFaceTwin(iPer) + itmp(i) )
      enddo
      do i=1,nfaces(ib)
        arx( startFaceTwin(iPer) + i ) = vtmp(i)
      enddo

      do i=1,nfaces(ib)
        vtmp(i) = ary( startFaceTwin(iPer) + itmp(i) )
      enddo   
      do i=1,nfaces(ib)
        ary( startFaceTwin(iPer) + i ) = vtmp(i)
      enddo

      do i=1,nfaces(ib)
        vtmp(i) = arz( startFaceTwin(iPer) + itmp(i) )
      enddo
      do i=1,nfaces(ib)
        arz( startFaceTwin(iPer) + i ) = vtmp(i)
      enddo

      deallocate(tmp,itmp,vtmp)

    endif

  enddo 


  write(*,'(a)') '  **Mapped opposite faces at periodic boundaries-minimum distance criteria.'


end subroutine



end module