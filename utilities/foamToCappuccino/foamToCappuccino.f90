program foamToCappuccino
!
!    _____                     _________________                                          .__               
!  _/ ____\_________    _____  \_____  \_   ___ \_____  ______ ______  __ __   ____  ____ |__| ____   ____  
!  \   __\/  _ \__  \  /     \  /  ____/    \  \/\__  \ \____ \\____ \|  |  \_/ ___\/ ___\|  |/    \ /  _ \ 
!   |  | (  <_> ) __ \|  Y Y  \/       \     \____/ __ \|  |_> >  |_> >  |  /\  \__\  \___|  |   |  (  <_> )
!   |__|  \____(____  /__|_|  /\_______ \______  (____  /   __/|   __/|____/  \___  >___  >__|___|  /\____/ 
!                   \/      \/         \/      \/     \/|__|   |__|               \/    \/        \/  
!  Purpose:
!    Converts polyMesh file format to polyMesh format compatible with freeCappuccino code.
!  
!  Description:
!    OpenFOAM's native polyMesh format consists of several files.
!    These are: 'points', 'faces', 'owner, 'neighbour', 'boundary'.
!    freeCappuccino uses very similar format except a few changes.
!    Numbering starts from 1, headers are different, 
!    no braces '()', no semicolons ';'
!    Boudnary file is structured differently, one line for every boundary region
!    describing BC name, BC type, number of faces and start face. 
!    For postprocessing we need also 'cells' file with cell connectivity
!    forwarded to Paraview or TecPlot. For that we use 'cellConnectivity' utility,
!    see 'utilities' folder.
!
!  Date:
!    23/12/2019
!
!  Author:
!    Nikola Mirkov (Email: largeddysimulation@gmail.com)
!
  use utils, only: get_unit, timestamp!, file_row_count, r8vec_print_some, i4vec_print2

  implicit none

  ! Double precision type
  integer, parameter :: dp = kind(1.0d0)

  ! General mesh data
  integer :: numNodes         ! no. of nodes in mesh
  integer :: numCells         ! no. of cells in the mesh
  integer :: numFaces         ! no. of INNER+BOUNDARY faces in the mesh
  integer :: numInnerFaces    ! no. of INNER cells faces in the mesh
  integer :: numBoundaryFaces ! self explanatory
  ! integer :: numTotal         ! number of volume field values + number of boundary field values numCells+numBoundaryFaces

  ! ! To define boundary regions
  ! integer :: numBoundaries
  ! integer, dimension(:), allocatable :: nfaces,startFace,iBndValueStart
  ! character(len=15), dimension(:), allocatable :: bcname, bctype

  ! Mesh file units
  integer :: points_file, faces_file, owner_file, neighbour_file, size_file
  integer :: points_file_of, faces_file_of, owner_file_of, neighbour_file_of

  real(dp) :: x,y,z  ! Coordinates of mesh nodes 

  integer :: i,j,k,l,itmp
  integer :: owner, neighbour
  integer :: iface

  character(len=1) :: ch
  character(len=20) :: char_string,char_string2
  character(len=80) :: line_string

  integer :: node(4)  ! It will store global node numbers of cell vertexes
  integer :: nnodes                  ! no. of nodes in face

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    _____                     _________________                                          .__               '
  write ( *, '(a)' ) '  _/ ____\_________    _____  \_____  \_   ___ \_____  ______ ______  __ __   ____  ____ |__| ____   ____  '
  write ( *, '(a)' ) '  \   __\/  _ \__  \  /     \  /  ____/    \  \/\__  \ \____ \\____ \|  |  \_/ ___\/ ___\|  |/    \ /  _ \ '
  write ( *, '(a)' ) '   |  | (  <_> ) __ \|  Y Y  \/       \     \____/ __ \|  |_> >  |_> >  |  /\  \__\  \___|  |   |  (  <_> )'
  write ( *, '(a)' ) '   |__|  \____(____  /__|_|  /\_______ \______  (____  /   __/|   __/|____/  \___  >___  >__|___|  /\____/ '
  write ( *, '(a)' ) '                   \/      \/         \/      \/     \/|__|   |__|               \/    \/        \/        '
  write ( *, '(a)' ) ' '
 

!******************************************************************************
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'.
!..............................................................................

  call get_unit( points_file_of )
  open( unit = points_file_of,file = 'polyMesh-OF/points' )
  rewind points_file_of

  call get_unit( faces_file_of )
  open( unit = faces_file_of, file='polyMesh-OF/faces' )
  rewind faces_file_of

  call get_unit( owner_file_of )
  open( unit = owner_file_of, file='polyMesh-OF/owner' )
  rewind owner_file_of

  call get_unit( neighbour_file_of )
  open( unit = neighbour_file_of, file='polyMesh-OF/neighbour' )
  rewind neighbour_file_of

  ! call get_unit( boundary_file_of )
  ! open( unit = boundary_file_of, file='polyMesh-OF/boundary' )
  ! rewind boundary_file_of

!******************************************************************************
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'
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

  ! call get_unit( boundary_file )
  ! open( unit = boundary_file, file='polyMesh/boundary' )
  ! rewind boundary_file

  call get_unit( size_file )
  open( unit = size_file, file='polyMesh/size' )
  rewind size_file

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
    read(owner_file_of,*) char_string,line_string
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

  rewind owner_file_of

  ! NOTE: Trying to acces number of faces data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numFaces.

  ch = ' '
  owner_loop: do 
    read(owner_file_of,*) ch
    if (ch == '(') then
      ! Return two lines
      backspace(owner_file_of)
      backspace(owner_file_of)
      exit owner_loop
    endif
  end do owner_loop

  read(owner_file_of,*) itmp
  if (itmp /= numFaces ) then
    write(*,*) "Error reading polyMesh format. numFaces value is not confirmed in body of the 'owner' file."
    stop
  endif
  read(owner_file_of,*) char_string ! reads "(" again


  !
  ! The 'points' file
  !

  ! NOTE: Trying to acces number of points data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numNodes.

  ch = ' '
  point_header_loop: do
    read(points_file_of,*) ch
    if (ch == '(') then
      ! Return two lines
      backspace(points_file_of)
      backspace(points_file_of)
      exit point_header_loop
    endif
  end do point_header_loop

  read(points_file_of,*) itmp
  if (itmp /= numNodes ) then
    write(*,*) "Error reading polyMesh format. numNodes value is not confirmed in 'points' file."
    stop
  endif
  read(points_file_of,*) char_string ! reads "("

  !
  ! The 'neighbour' file
  !

  ! NOTE: Trying to acces number of inner faces data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numInnerFaces.

  ch = ' '
  neighbour_header_loop: do
    read(neighbour_file_of,*) ch
    if (ch == '(') then
      ! Return two lines
      backspace(neighbour_file_of)
      backspace(neighbour_file_of)
      exit neighbour_header_loop
    endif
  end do neighbour_header_loop

  read(neighbour_file_of,*) itmp
  if (itmp /= numInnerFaces ) then
    write(*,*) "Error reading polyMesh format. numInnerFaces value is not confirmed in 'neighbour' file."
    stop
  endif
  read(neighbour_file_of,*) char_string ! reads "("

  !
  ! The 'faces' file
  !

  ! NOTE: Trying to acces number of faces data. So we go check line by line, 
  ! when we get to "(" we go two lines back and read numFaces.

  ch = ' '
  faces_header_loop: do
    read(faces_file_of,*) ch
    if (ch == '(') then
      backspace(faces_file_of)
      backspace(faces_file_of)
      exit faces_header_loop
    endif
  end do faces_header_loop

  read(faces_file_of, *) itmp
  if (itmp /= numFaces ) then
    write(*,*) "Error reading polyMesh format. numFaces value is not confirmed in 'faces' file."
    stop
  endif
  read(faces_file_of,*) char_string

  ! Number of boundary faces
  numBoundaryFaces = numFaces - numInnerFaces

  ! Size of arrays storing variables numCells+numBoundaryFaces
  ! numTotal = numFaces + numBoundaryFaces


!******************************************************************************
! > Write headers 
!..............................................................................

  ! write( points_file, '(i0,1x,a)' ) numNodes, 'numNodes'
  ! write( owner_file, '(i0,1x,a)' ) numFaces, 'numFaces'
  ! write( neighbour_file, '(i0,1x,a)' ) numInnerFaces, 'numInnerFaces'
  ! write( faces_file, '(i0,1x,a)') numFaces, 'numFaces'

  write(size_file,'(I0,a)') numNodes,    " Vertices"
  write(size_file,'(I0,a)') numCells,    " Cells"
  write(size_file,'(I0,a)') numInnerFaces, " Inner faces"
  write(size_file,'(I0,a)') numBoundaryFaces, " Boundary faces"
  write(size_file,'(I0,a)') numFaces, " Faces Total"

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


!******************************************************************************
! > Read Mesh files 
!..............................................................................

  ! The 'points' file
  do i=1,numNodes
    ! char_string reads (number, char_string reads number), we have to strip off the brackets later.
    read(points_file_of,*) char_string,y,char_string2 

    ! Char to double conversion + stripping off the brackets:
    read(char_string(2:),*) x
    read(char_string2(1:len_trim(char_string2)-1),*) z 

    write(points_file,*) x,y,z

  end do


  ! The 'owner' file
  do i=1,numFaces
    read(owner_file_of,*) owner
    owner = owner + 1 ! fortran starts from 1
    write(owner_file,'(i0)') owner
  end do



  ! The 'neighbour' file
  do i=1,numInnerFaces
    read(neighbour_file_of,*) neighbour
    neighbour = neighbour + 1 ! fortran starts from 1
    write(neighbour_file, '(i0)' ) neighbour
  end do

  ! The 'faces' file
  do iface=1,numFaces

    ! Initialize array of node indexes that construct the face
    node = 0

    ! Read line in 'faces' file
    call read_line_faces_file_polyMesh(faces_file_of,nnodes,node,4)

    write(faces_file, '( 5(i0,1x) )') nnodes, (node(k), k=1,nnodes)

  end do


  write ( *, '(a)' ) ' Completed succesfully!'
  
!
!  > CLOSE polyMesh format file: 'points', 'faces', 'owner', 'neighbour', 'boundary'.
!
  close ( points_file_of )
  close ( faces_file_of )
  close ( owner_file_of )
  close ( neighbour_file_of )
  ! close ( boundary_file_of )

  close ( points_file )
  close ( faces_file )
  close ( owner_file )
  close ( neighbour_file )
  ! close ( boundary_file )
  close (size_file )

!+-----------------------------------------------------------------------------+

end program



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