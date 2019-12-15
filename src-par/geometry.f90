!
! Module for geometry definition on unstructured meshes.
!
module geometry

use types
use parameters, only: myid
use utils, only: get_unit, file_row_count, r8vec_print_some, i4vec_print2, i4_to_s_left

implicit none

! NOTE:
! In variable arrays, field variables defined at cell centres are written in positions from 1 to numCells, after that we write
! variable values for boundary faces from numCells+1 to numTotal.  

! General mesh data
integer :: numNodes                           ! no. of nodes in mesh
integer :: numCells                           ! no. of cells in the mesh
integer :: numPCells                          ! no. of cells in the mesh + buffer cells
integer :: numFaces                           ! no. of INNER+BOUNDARY faces in the mesh
integer :: numInnerFaces                      ! no. of INNER cells faces in the mesh
integer :: numBoundaryFaces                   ! self explanatory
integer :: numTotal                           ! number of volume field values + number of boundary field values numCells+numBoundaryFaces

! To define boundary regions
integer :: numBoundaries
integer, dimension(:), allocatable :: nfaces,startFace,iBndValueStart
character(len=15), dimension(:), allocatable :: bcname, bctype

integer :: npro                               ! No. of boundary faces on 'process' boundary - connected to another computational domain
integer :: ninl                               ! No. of inlet boundary faces
integer :: nout                               ! No. of outlet boundary faces
integer :: nsym                               ! No. of symmetry boundary faces
integer :: nwal                               ! No. of wall boundary faces
          

integer, parameter :: nomax = 30              ! Max no. of nodes in face - determines size of some arrays, just change this if necessary.
real(dp), parameter :: tiny = 1e-30


! Mesh file units
integer :: points_file, faces_file, owner_file, neighbour_file, boundary_file, process_file

! Global number of nodes, cells, faces, inner faces, and boundary faces when summed from all domains
integer :: gloNodes,gloCells, gloFaces, gloIFaces, gloBFaces


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
real(dp), dimension(:), allocatable :: fpro            ! Interpolation factor for faces at process boundaries


! Geometry parameters defined for boundary faces
real(dp), dimension(:), allocatable :: srds,dns    ! srds = |are|/|dns|, dns = normal distance to cell center from face |dpn*face_normal_unit_vec|
real(dp), dimension(:), allocatable :: srdw,dnw    ! srdw = |are|/|dnw|, dnw = normal distance to cell center from face |dpn*face_normal_unit_vec|

! Mesh topology information - connectivity of cells trough faces
integer, dimension(:), allocatable :: owner     ! Index of the face owner cell
integer, dimension(:), allocatable :: neighbour ! Index of the neighbour cell  - it shares the face with owner


public 

contains


subroutine read_mesh
!
!  Description:
!    Calculates basic geometrical quantities of numerical mesh
!    defined in this module and needed for FVM computation.
!    Mesh is described in files similar to polyMesh format.
!    Mesh files are 'points', 'faces', 'owner, 'neighbour', 'boundary'
!
!  Date:
!    26/11/2015, 12/06/2019
!
!  Author:
!    Nikola Mirkov nmirkov@vinca.rs
!
 
  use my_mpi_module

  implicit none

  ! Locals
  integer :: i,j,k,l,itmp
  integer :: ib,iwall,isym,ipro
  integer :: iface
  integer :: inp,inn,ijp

  character(len=1) :: ch
  character(len=20) :: char_string,char_string2
  character(len=80) :: line_string
  character( len = 5) :: nproc_char

  integer, dimension(nomax) :: node  ! It will store global node numbers of cell vertexes
  integer :: nnodes                  ! no. of nodes in face

  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: third = 1./3._dp

  real(dp) :: px,py,pz, qx,qy,qz, nx,ny,nz, cx,cy,cz
  real(dp) :: riSi
  real(dp) :: are,areasum
  real(dp) :: xpn,ypn,zpn
  real(dp) :: xjp,yjp,zjp
  real(dp) :: dpn,djn
  real(dp) :: nxf,nyf,nzf
 
  ! Array for temporary storage of doubles
  real(dp), dimension(:), allocatable :: r8tmp

!******************************************************************************
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'.
!..............................................................................

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )


  call get_unit( points_file )
  open( unit = points_file,file = 'processor'//trim(nproc_char)//'/constant/polyMesh/points' )
  rewind points_file

  call get_unit( faces_file )
  open( unit = faces_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/faces' )
  rewind faces_file

  call get_unit( owner_file )
  open( unit = owner_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/owner' )
  rewind owner_file

  call get_unit( neighbour_file )
  open( unit = neighbour_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/neighbour' )
  rewind neighbour_file

  call get_unit( boundary_file )
  open( unit = boundary_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/boundary' )
  rewind boundary_file

  call get_unit( process_file )
  open( unit = process_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/process' )
  rewind process_file


!
! > Read 'process' file with domain connectivity information for MPI
!

  ! Number of rows in the file excluding #comment in header
  call file_row_count ( process_file, numConnections)

  read(process_file,'(a)') line_string ! header line
  ! read(process_file,*) numConnections
  
  ! Allocate domain connectivity arrays
  allocate ( neighbProcNo( numConnections ) )
  allocate ( neighbProcOffset( numConnections+1 ) )

  do i=1,numConnections
    read(process_file,*) neighbProcNo(i)
  enddo


!
! > Read boundary conditions file. 
!

!
!   Boundary conditions file consists of header and numBoundaries number of subsequent lines each containing:  
!   boundary condition given name, bc type, number of faces belonging to that bc and starting face for that bc.
!   Possible boundary contiions types are: process, inlet, outlet, symmetry, wall, wallIsoth, wallAdiab, wallQFlux, prOutlet, etc. 
!   Please check are all of these implemented, because this is still in the development. Contact me if you have any questions.
!
 

  ! Number of rows in the file excluding #comment in header to conclude the number of prescribed boundaries
  call file_row_count ( boundary_file, numBoundaries )

  read(boundary_file,'(a)') line_string ! Firts line is header.
  ! read(boundary_file,*) numBoundaries ! it doesn't need to read the number of boundaries because of the above

  ! Allocate 
  allocate ( bcname(numBoundaries) )
  allocate ( bctype(numBoundaries) )
  allocate ( nFaces(numBoundaries) )
  allocate ( startFace(numBoundaries) )
  allocate ( iBndValueStart(numBoundaries))

  nwal = 0
  nsym = 0
  ninl = 0
  npro = 0

  neighbProcOffset(1) = 1


  do i=1,numBoundaries
    read(boundary_file,*) bcname(i), bctype(i), nfaces(i) ,startFace(i) ! Now it reads four things btw

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

    elseif( bctype(i) == 'process') then
      npro = npro + nfaces(i)
      neighbProcOffset(i+1) = neighbProcOffset(i) + nfaces(i)

    endif
  enddo  

  ! This denotes position of the last plus one. 
  neighbProcOffset( numConnections+1 ) = npro + 1


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

  ! Size of array for gradients etc. which are numCells plus No. of buffer cells npro
  numPCells = numCells + npro

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

  gloNodes = numNodes
  call global_isum( gloNodes )

  gloCells = numCells
  call global_isum( gloCells )

  gloFaces = numFaces
  call global_isum( gloFaces )

  gloIFaces = numInnerFaces
  call global_isum( gloIFaces )
 
  gloBFaces = numBoundaryFaces
  call global_isum( gloBFaces )  

  if ( myid .eq. 0) then

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Mesh data: '



  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes, numNodes = ', gloNodes


  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cells, numCells = ', gloCells


  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cell-faces, numFaces = ', gloFaces


  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of inner cell-faces, numInnerFaces = ', gloIFaces

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary information:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cell-faces on boundary, numBoundaryFaces = ', gloBFaces
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary information (bcname, bctype, nFaces, startFace):'
  write ( *, '(a)' ) ' '
  
  do i=1,numBoundaries
    write(*,'(2x,a,1x,a,1x,2(i0,1x))') bcname(i), bctype(i), nfaces(i) ,startFace(i) ! Now it reads four things btw
  enddo  

  if( npro.gt.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of processor boundary faces  = ', npro
  endif

  write ( *, '(a)' ) ' '

endif


!******************************************************************************
! > Allocate arrays for Mesh description
!..............................................................................

  ! Nodal coordinates
  allocate ( x(numNodes) )
  allocate ( y(numNodes) )
  allocate ( z(numNodes) )

  ! Coordinates of cell centers
   
  ! For now we will use non efficient way to have c. centers and vol of size numTotal,
  ! this makes usage of flux routines much easier for 'process' boundaries,
  ! in the near future maybe we should change to smaller size of arrays (numPcells)
  ! since it is wastefull to have Volume and C. center for every boundary face.
  ! Let us recall, cell center values are stored in positions 1,numCells, boundary face values
  ! are stored in numCells+1, numTotal. If I want to exchange something I have only in the interior
  ! I have these things stored in positions numCells+1, numPCells and exchange them with
  ! subroutine exchange__ with double underscore. This one is inteded for thiese smaller arrays.

  ! allocate ( xc(numPCells) )
  ! allocate ( yc(numCPells) )
  ! allocate ( zc(numCPells) )

  allocate ( xc(numTotal) )
  allocate ( yc(numTotal) )
  allocate ( zc(numTotal) )

  ! Cell volumes
  ! allocate ( vol(numCells+Npro) )

  allocate ( vol(numTotal) )


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

  ! We need this only for wall and symmetry, so it is not bad to have nsym and nwal wich count home many of them there is.
  allocate ( dns(nsym) )      
  allocate ( srds(nsym) )  
  allocate ( dnw(nwal) )      
  allocate ( srdw(nwal) ) 

  ! Interpolation factor for process faces
  allocate ( fpro(npro) )                               

!
! > Allocate and initialize parameters for MPI communication, declared in mpi_module
!
  lenbuf = npro
  allocate ( bufind(lenbuf) )
  allocate ( buffer(lenbuf) )




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


  !
  ! > Array of buffer cell indices for MPI exchange
  !

  iPro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'process') then

      do i=1,nfaces(ib)

        iPro = iPro + 1

        iface = startFace(ib) + i
        bufind( ipro ) = owner( iface )

      enddo
    endif
  enddo
   



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Cell volumes, cell centers and cell face centers
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Allocate tmp array of doubles 
  allocate( r8tmp(numCells))
  r8tmp = 0.0_dp

  do iface=1,numFaces

    inp = owner(iface)

    ! Initialize array of node indexes that construct the face
    node = 0

    ! Read line in 'faces' file
    call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)

    ! Initialize total area of polygon.
    areasum = 0.0_dp

    ! We decompose a polygon to a series of triangles, all having first node in common.
    do i=1,nnodes-2 
      ! Vectors to vertices
      ! 2-1
      px = x(node(i+1))-x(node(1))
      py = y(node(i+1))-y(node(1))
      pz = z(node(i+1))-z(node(1))
      ! 3-1
      qx = x(node(i+2))-x(node(1))
      qy = y(node(i+2))-y(node(1))
      qz = z(node(i+2))-z(node(1))


      call triangle_area_vector( px,py,pz, qx,qy,qz, nx,ny,nz )

      !
      ! > Cell-face area vector components (Area vector lies in direction of face normal)
      ! 

      arx(iface) = arx(iface) + nx
      ary(iface) = ary(iface) + ny
      arz(iface) = arz(iface) + nz

      ! Face center for a triangle
      cx = third*( x(node(i+2)) + x(node(i+1)) + x(node(1)) )
      cy = third*( y(node(i+2)) + y(node(i+1)) + y(node(1)) )
      cz = third*( z(node(i+2)) + z(node(i+1)) + z(node(1)) ) 

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

  ! We need some geometry in the buffer cells, we fill these by exchanging info with other processes
  ! call exchange__( xc )
  ! call exchange__( yc )
  ! call exchange__( zc )
  ! call exchange__( Vol )


  call exchange( xc ) 
  call exchange( yc )
  call exchange( zc )
  call exchange( Vol )


  ! Rewind 'faces' file for one more sweep
  rewind( faces_file )
  ch = ' '
  do
    read(faces_file,*) ch
    if (ch == "(") then
      exit
    endif
  end do 


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! > Interpolation factor
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  do iface=1,numFaces

    !
    ! > Interpolation factor > inner faces
    !

    if(iface <= numInnerFaces) then

      inp = owner(iface)
      inn = neighbour(iface)

      node = 0

      ! Read line in 'faces' file
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)

      xpn = xc(inn)-xc(inp)
      ypn = yc(inn)-yc(inp)
      zpn = zc(inn)-zc(inp)

      dpn = sqrt( xpn**2 + ypn**2 + zpn**2 ) 

      !
      ! > > Intersection point j' of line connecting centers with cell face, we are taking only three points assuming that other are co-planar
      !
      call find_intersection_point(&
                                   ! plane defined by three face vertices:
                                   x(node(1)),y(node(1)),z(node(1)),&
                                   x(node(2)),y(node(2)),z(node(2)), &
                                   x(node(3)),y(node(3)),z(node(3)), &
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

    endif
    
  enddo


  !
  ! > Interpolation factor > faces on process boundaries
  !

  iPro = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'process') then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        inp = owner(iface)
        inn = iBndValueStart(ib) + i

        iPro = iPro + 1


        node = 0

        ! Read line in 'faces' file
        call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)


        xpn = xc(inn)-xc(inp)
        ypn = yc(inn)-yc(inp)
        zpn = zc(inn)-zc(inp)

        dpn = sqrt( xpn**2 + ypn**2 + zpn**2 )

        ! > > Intersection point j' of line connecting centers with cell face, we are taking only three points assuming that other are co-planar
        call find_intersection_point( &
                                     ! plane defined by three face vertices:
                                     x(node(1)),y(node(1)),z(node(1)),&
                                     x(node(2)),y(node(2)),z(node(2)), &
                                     x(node(3)),y(node(3)),z(node(3)), &
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
        fpro(iPro) = djn/dpn

      enddo

    else
    ! Not a proces boundary

      do i=1,nfaces(ib)

        ! Read line in 'faces' file
        call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)

      enddo

    endif

  enddo ! interpolation factor loop

  ! Loop over wall boundaries to calculate normal distance from cell center of the first cell layer - dnw, and
  ! loop over symmetry boundaries to calculate normal distance from cell center of the first cell layer - dns.

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

  ! if (myid .eq. 0 ) then

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

  ! endif


!
!  > CLOSE polyMesh format file: 'points', 'faces', 'owner', 'neighbour', 'boundary', 'process'.
!
  close ( points_file )
  close ( faces_file )
  close ( owner_file )
  close ( neighbour_file)
  close ( boundary_file)
  close ( process_file)
!+-----------------------------------------------------------------------------+


end subroutine read_mesh


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


pure function cell_volume_part_polymesh(ax,ay,az,nx,ny,nz) result(volume_part)
!
! Ai, i=1,n are n triangular faces enclosing polyhedron P.
! Vertices of Ai are (ai,bi,ci),
!                                _    _   _    _
! Unit normal to P on each Ai is ni = ni/|ni|, ni = (bi-ai)x(ci-ai)
! Volume of P is given by:
!                         _  _                     _  _ 
!    V = 1/6 * sum_i=1^n (ai.ni) or sum_i=1^n 1/6*(ai.ni)
!
! Sources: 
! [1] Dr Robert Nurnberg, Imperial College London, www.ma.ic.ac.uk/~rn/centroid.pdf
! [2] paulbourke.net/geometry/polygonmesh

  implicit none
  real(dp), intent(in) :: ax,ay,az
  real(dp), intent(in) :: nx,ny,nz
  real(dp) :: volume_part

  volume_part = 1./6._dp*( ax*nx + ay*ny + az*nz )
end function


pure function centroid_component_part_polymesh(ax,bx,cx,nx,vol) result(cc_part)
!
! Ai, i=1,n are n triangular faces enclosing polyhedron P.
! Vertices of Ai are (ai,bi,ci),
!                                _    _   _    _
! Unit normal to P on each Ai is ni = ni/|ni|, ni = (bi-ai)x(ci-ai)
! Cell center (cc) component of P is given by:
!                                  
!    cc.ed = 1/2V * sum_i=1^n (x.ed)^2*(ni.ed), where ed is unit basis of R^3, d=1,2,3.
!
! Thios function calculates part under summation which is accumulated in outer loop. 
!
! Sources: 
! [1] Dr Robert Nurnberg, Imperial College London, www.ma.ic.ac.uk/~rn/centroid.pdf
! [2] paulbourke.net/geometry/polygonmesh

  implicit none
  real(dp), intent(in) :: ax,bx,cx
  real(dp), intent(in) :: vol,nx
  real(dp) :: cc_part

  cc_part = 1./(2*vol) * 1./24._dp*nx * ( (ax+bx)**2 + (bx+cx)**2 + (cx+ax)**2 )

end function  


  function determinant(a1,a2,a3,b1,b2,b3,q1,q2,q3)
!
! Calculates determinant of 3x3 matrix
!
  implicit none
  real(dp) :: a1,a2,a3,b1,b2,b3,q1,q2,q3
  real(dp) :: determinant

  determinant = (a2*b3-b2*a3)*q1 &
               +(b1*a3-a1*b3)*q2 &
               +(a1*b2-a2*b1)*q3
  end function



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
! Nikola Mirkov @2016
!
! example usage: 
! call find_intersection_point( &
!!                            plane defined by face corners:
!                             x(inp),y(inp),z(inp),&
!                             x(inp-idns),y(inp-idns),z(inp-idns), &
!                             x(inp-idtb),y(inp-idtb),z(inp-idtb), &
!!                            line defined by cell center and neighbour center:
!                             xc(inp),yc(inp),zc(inp), &
!                             xc(inp+idew),yc(inp+idew),zc(inp+idew), &
!!                            intersection point:
!                             xj,yj,zj &
!                             )
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



end module
