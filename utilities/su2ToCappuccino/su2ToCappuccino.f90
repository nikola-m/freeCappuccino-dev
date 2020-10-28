program su2ToCappuccino
!
! 
!             _____ _____     _____                                  _             
!            / __  \_   _|   /  __ \                                (_)            
!  ___ _   _ `' / /' | | ___ | /  \/ __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___  
! / __| | | |  / /   | |/ _ \| |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \ 
! \__ \ |_| |./ /___ | | (_) | \__/\ (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |
! |___/\__,_|\_____/ \_/\___/ \____/\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/ 
!                                        | |   | |                                 
!                                        |_|   |_|                                 
!
! Description:
!
! Converts Stanford Unstructured su2 file to mesh format compatible with freeCappuccino solver.
!
! Author:
!   Nikola Mirkov (largeddysimulation@gmail.com)
!
! Date:
!   10th October 2020

use utils

implicit none

 integer :: iarg
 integer :: iargc
 integer :: ios
 integer :: num_arg
 character ( len = 255 ) prefix
 character ( len = 255 ) input_filename

 integer :: su2_file = 4

 integer, parameter :: dp = kind(1.0d0)
 integer, parameter :: nonoel = 8 ! no. of nodes in element-here Hex
 integer, parameter :: nonofa = 4 ! no, of nodes in element face-here Hex
 integer, parameter :: nofaelmax     = 6 ! no. of faces in element
 integer :: nel    ! no. of elements in mesh
 integer :: ndim   ! dimension of the problem. 3 for 3D.
 integer :: i,iel,k,nm,l,np
 integer :: iface,jface
 integer :: numhex ! No. of Hex elements 
 integer :: numpri ! No. of Prism elements
 integer :: numtet ! No. of Tet elements 
 integer :: numpyr ! No. of Pyr elements 
 integer :: numPoints ! no. of nodes in mesh
 integer :: numCells
 integer :: nfacesTotal
 integer :: nBoundFaces 
 integer :: nInnerFaces
 integer :: nInnerFacePairs
 integer :: listLength
 integer :: NE, NTYPE, NODE(8), fVert(4) 
 real(dp) :: XCOO(3)
 integer :: NMARK ! no. of boundary regions
 integer, dimension(:,:), allocatable :: fV
 integer, dimension(:), allocatable :: owner,neighbour
 integer, dimension(:), allocatable :: own,nbr,prmt
 integer, dimension(:,:), allocatable :: fV1
 integer, dimension(:), allocatable :: bcSize
 character( len = 30 ), dimension(:), allocatable :: bcName
 character( len = 10 ) :: key


!+-----------------------------------------------------------------------------+
  call timestamp ( )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) "              _____ _____     _____                                  _             "          
  write ( *, '(a)' ) "             / __  \_   _|   /  __ \                                (_)            "              
  write ( *, '(a)' ) "   ___ _   _ `' / /' | | ___ | /  \/ __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___  " 
  write ( *, '(a)' ) "  / __| | | |  / /   | |/ _ \| |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \ "
  write ( *, '(a)' ) "  \__ \ |_| |./ /___ | | (_) | \__/\ (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |"
  write ( *, '(a)' ) "  |___/\__,_|\_____/ \_/\___/ \____/\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/ "
  write ( *, '(a)' ) "                                         | |   | |                                 "
  write ( *, '(a)' ) "                                         |_|   |_|                                 "
  write ( *, '(a)' ) " "           
  write ( *, '(a)' ) 'su2ToCappuccino'
  write ( *, '(a)' ) '  A preprocessor program for freeCappuccino code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reads mesh in .su2 format and produces  '
  write ( *, '(a)' ) '  mesh files for freeCappuccino solver.          '
  write ( *, '(a)' ) ' '
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the filename prefix:'
    read ( *, '(a)', iostat = ios ) prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, prefix )

  end if
!
!  Create the filenames.
!
  input_filename = trim ( prefix ) // '.su2'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is "' // trim ( input_filename ) // '".'
  write ( *, '(a)' ) ' '
!+-----------------------------------------------------------------------------+


! OPEN fcp polyMesh format file: 'boundary','faces', 'owner', 'neighbour'.
  open(unit=8,file='polyMesh/boundary')
  open(unit=9,file='polyMesh/faces')
  open(unit=10,file='polyMesh/owner')
  open(unit=11,file='polyMesh/neighbour')
  open(unit=12,file='polyMesh/cells')
  open(unit=13,file='polyMesh/size')

  rewind 8  
  rewind 9
  rewind 10
  rewind 11
  rewind 12
  rewind 13

!
! > Read input from su2 mesh file (.su2)
!
  open(unit=su2_file,file=input_filename,status='old')
  rewind su2_file

  read(su2_file,*) key,ndim
  read(su2_file,*) key,nel  

  ! List length for the hash map
  listLength = nofaelmax*nel*2

  allocate( owner(listLength), neighbour(listLength),fV(nonofa,listLength) )

  ! Initialize
  fVert = 0
  jface = 0 ! counter

  ! Initialize total number of hex cells, prism cells, etc.
  numhex = 0
  numpri = 0
  numtet = 0
  numpyr = 0

  ! Initialize number of faces.
  nfacesTotal = 0
  nBoundFaces = 0
  nInnerFaces = 0

  owner = 0
  neighbour = 0
  fV = 0

  !
  ! > Read elements
  !
  element_loop: do iel=1,nel

  read(su2_file,*) NTYPE, (NODE(k), k=1,noel(NTYPE)), NE

  np = noel(NTYPE)

  nfacesTotal = nfacesTotal + nofael(NTYPE)

  ! > Write into 'cells' polyMesh file
  write(12,'(I0,1X,8(I0,1x))') NTYPE,(NODE(k), k=1,np)


! write(*,*)iel
!  Element Type  Identifier
!
!  Line           3
!  Triangle       5
!  Quadrilateral  9
!  Tetrahedral   10
!  Hexahedral    12
!  Prism         13
!  Pyramid       14

! From su2format specs:
! Note again that the global index values for the nodes and elements stored 
! within SU2 are zero-based...
! The ordering of the nodes given in the connectivity list for a specific element 
! is important, and the user is referred to the VTK format guide for the correct 
! ordering for each supported element type (page 9).


  if (NTYPE.eq.10) then
  ! TETRAHEDRON

    do iface = 1,nofael(NTYPE)
      fVert = face_vertices_NTYPE10(iface,NODE,np)
      call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
    enddo

    numtet = numtet + 1 


  elseif (NTYPE.eq.12) then
  ! HEXAHEDRON

    do iface = 1,nofael(NTYPE)
      fVert = face_vertices_NTYPE12(iface,NODE,np)
      call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
    enddo

    numhex = numhex + 1 

  elseif (NTYPE.eq.13) then
  ! PRISM

    do iface = 1,nofael(NTYPE)
      fVert = face_vertices_NTYPE13(iface,NODE,np)
      call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
    enddo

    numpri = numpri + 1 


  elseif (NTYPE.eq.14) then
  ! PYRAMID

    do iface = 1,nofael(NTYPE)
      fVert = face_vertices_NTYPE14(iface,NODE,np)
      call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
    enddo

    numpyr = numpyr + 1 


  endif 

  jface = jface + nofael(NTYPE)

  end do element_loop

!
! > Report after reading al the meesh elements
!

  numCells = numhex+numtet+numpyr+numpri

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) ' Gambit file:'
  write ( *, '(4x,a,I0)' ) 'Total no. of HEX cells: ', numhex
  write ( *, '(4x,a,I0)' ) 'Total no. of TET cells: ', numtet
  write ( *, '(4x,a,I0)' ) 'Total no. of PYR cells: ', numpyr
  write ( *, '(4x,a,I0)' ) 'Total no. of PRI cells: ', numpri
  write ( *, '(3x,a)' )    '+--------------------------------='
  write ( *, '(4x,a,I0)' ) 'Total no. of cells: ', numCells
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) 'Normal end of reading .su2 file.'
  write ( *, '(a)' ) ' '

  if(nel.ne.numCells) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  Unexpected error: Did not read all cells in the file!'
    stop  
  endif  


  write(*,*)" Found",nInnerFaces,"inner faces!"

  allocate( own(nInnerFaces), nbr(nInnerFaces), prmt(nInnerFaces), fV1(nonofa,nInnerFaces) )

  write(*,'(a)') " Sorting 'owner' and 'neighbour' arrays..."

  ! Fill arrays
  l = 0
  do i=1,listLength
    if ( owner(i).ne.0 .and. neighbour(i).ne.0 ) then
      l = l+1
      own(l) = owner(i)
      nbr(l) = neighbour(i)
      fV1(:,l) = fV(:,i)
      prmt(l) = l
    endif
  enddo

  ! Sort in ascending order
  call i4vec3_sort_a ( nInnerFaces, own, nbr, prmt )

  write(*,'(a)') " Writting polyMesh files..."

  ! Write to files
  do i=1,nInnerFaces

    np = 4
    if( fV1(4,prmt(i)) == 0 ) np=3

    ! Add one bcs zero numbering in su2, reverse order in faces (np:1:-1, instead 1:np), 
    ! and pick face based on permutation done in sorting, which is written in 'prmt' array.
    write(9, '(5(I0,1x))') np, fV1( np:1:-1, prmt(i) )+1 
    write(10,'(I0)') own(i)
    write(11,'(I0)') nbr(i)

    ! Note: I had to reverse the order of indices in faces because I have read them in CCW order when 
    ! looking from inside cell. Reversing order while writing to file is easier than changing everything 
    ! in 'get_face_vertices_XXX' functions in utils module. Nikola

  enddo

  deallocate(own,nbr,prmt,fV1)

  !
  ! > Continue reading mesh file, first read node cooradinates and write them to a file
  !

  ! OPEN text file: 'points'
  open(unit=7,file='polyMesh/points')
  rewind 7

  read(su2_file,*) key,numPoints ! NPOIN in su2 file


  ! Read nodal coordinates
  do i=1,numPoints

    read(su2_file,*) (XCOO(k),k=1,ndim)

    ! Write to file: 'points'
    write(7,'(3E20.11)') (XCOO(k), k=1,ndim)

  enddo

  ! CLOSE file: 'points'
  close (7)


! Read boudary conditions one by one

  nBoundFaces = 0

  nFacesTotal = nInnerFaces

  iel = 0  ! Need to reset it so the neighbour inside hashmap_insert remains 0 - a hack.  

  write(8,'(a)') '# bcName bcType nFaces startFace'

  read(su2_file,*) key, NMARK ! NMARK will hold how many boundary regions there is.

  allocate( bcName(NMARK), bcSize(NMARK) )

  do nm = 1,NMARK

    read(su2_file,*) key,bcName(nm)
    read(su2_file,*) key,bcSize(nm)

    ! Write down data in 'boundary' file
    write(8,'(a,1x,a,2(1x,I0))') trim(bcName(nm)), 'patch', bcSize(nm), nInnerFaces+nBoundFaces

    ! Now update total number of boundary faces
    nBoundFaces = nBoundFaces + bcSize(nm)


    boundary_loop: do i=1,bcSize(nm)

      read(su2_file,*) NTYPE, (fVert(k), k=1,noel(NTYPE))

      np = noel(NTYPE)

      ! > Write into 'cells' polyMesh file
      write(12,'(I0,1X,4(I0,1X))') NTYPE,(NODE(k), k=1,np)

      call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nFacesTotal,ne)

      write(10,'(I0)') ne ! owner cell on return form hashmap_insert
      write(9,'(5(I0,1x))') np,(fVert(k)+1, k=np,1,-1) ! face vertices of boundary face, reversed.

      ! Note: I had to reverse the order of indices in faces because I have read them in CCW order when 
      ! looking from inside cell. Reversing order while writing to file is easier than changing everything 
      ! in 'get_face_vertices_XXX' functions in utils module. Nikola

    enddo boundary_loop

  enddo

  ! End reading and processing data from .su2 mesh file.
  close(su2_file)

  deallocate( owner, neighbour,fV )

  ! Write 'size' file containing mesh size parameters.
  nInnerFacePairs = jface-nBoundFaces
  nInnerFaces = nInnerFacePairs/2
  nFacesTotal = nInnerFaces + nBoundFaces
 
  write(*,*) "Inner faces: ",nInnerFaces,"Boundary faces: ",nBoundFaces,"Total: ",nfacesTotal

  write(13,'(I0,a)') numPoints,   " Vertices"
  write(13,'(I0,a)') nel,         " Cells"
  write(13,'(I0,a)') nInnerFaces, " Inner faces"
  write(13,'(I0,a)') nBoundFaces, " Boundary faces"
  write(13,'(I0,a)') nfacesTotal, " Faces Total"

  ! Close polyMesh format files 
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)

  !
  ! Other releant task for definition of simulation case...
  !

  ! Create template files for vector and scalar fields in 0/ folder
  call create_field_init_files(bcName,NMARK)

  ! Create mesh file for Paraview in .vtm format, which will include all regions, inner and boundaries.
  write(*,'(a)') " Writing mesh in Paraview .vtu/.vtp format files..."

  call create_paraview_files( numPoints, nel, nInnerFaces, bcSize, NMARK)
 
  deallocate( bcName, bcSize )


end program