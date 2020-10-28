 program gmshToCappuccino
! Description:
!
! Converts Gmsh .msh files to one for use in Cappuccino solver.
!
! Author:
!   Nikola Mirkov
! 
! Date:
!   5. November 2015.
!

 use utils

 implicit none

 integer :: iarg
 integer :: iargc
 integer :: ios
 integer :: num_arg
 character ( len = 255 ) prefix
 character ( len = 255 ) input_filename

 integer :: gmsh_file = 4

 integer, parameter :: dp = kind(1.0d0)

 integer, parameter :: nonoel = 8 ! no. of nodes in element-here Hex
 integer, parameter :: nonofa = 4 ! no, of nodes in element face-here Hex
 integer, parameter :: nofaelmax     = 6 ! no. of faces in element
 integer, parameter :: nofael_NTYPE4 = 4 ! no. of faces in element-here Tet
 integer, parameter :: nofael_NTYPE5 = 6 ! no. of faces in element-here Hex
 integer, parameter :: nofael_NTYPE6 = 5 ! no. of faces in element-here Prism
 integer, parameter :: nofael_NTYPE7 = 5 ! no. of faces in element-here Pyramid
 integer, parameter :: ndim = 3          ! dimension of the problem. 3 for 3D.

 integer, dimension(:,:), allocatable :: fV
 integer, dimension(:), allocatable :: owner
 integer, dimension(:), allocatable :: neighbour

 integer :: nonome ! no. of nodes in mesh
 integer :: nel ! no. of elements in mesh
 integer :: numhex ! No. of Hex elements 
 integer :: numpri ! No. of Prism elements
 integer :: numtet ! No. of Tet elements 
 integer :: numpyr ! No. of Pyr elements 
 integer :: nInnerFaces
 integer :: nBoundFaces
 integer :: nFacesTotal
 integer :: nInnerFacePairs
 integer :: listLength

 ! For boundary conditions
 integer :: nbc
 integer, dimension(:), allocatable :: iBc, iZone
 character( len = 30 ), dimension(:), allocatable :: bcName

 integer :: NELEM
 integer :: NP
 integer :: NE, NTYPE, NODE(8), NTAGS, TAGS(4),fVert(4)
 real(dp) :: XCOO (3)

 integer :: i,k
 integer :: iel
 integer :: iface,jface
 integer :: njump
 character( len = 80 ) :: aLine

! 1) Intro
!+-----------------------------------------------------------------------------+
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write ( *, '(a)' ) 'gmshToCappucino'
  write ( *, '(a)' ) '  A preprocessor program for the Cappuccino code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reads mesh in .msh format and produces  '
  write ( *, '(a)' ) '  mesh files for freeCappuccino solver.          '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
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
  input_filename = trim ( prefix ) // '.msh'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Input file is "' // trim ( input_filename ) // '".'
  write ( *, '(a)' ) ' '
!+-----------------------------------------------------------------------------+


! OPEN fCp polyMesh format files:
  open(unit=7,file='points')
  open(unit=8,file='boundary')
  open(unit=9,file='faces')
  open(unit=10,file='owner')
  open(unit=11,file='neighbour')
  open(unit=12,file='cells')

  rewind 7
  rewind 8  
  rewind 9
  rewind 10
  rewind 11
  rewind 12


! > Read input from Gmsh mesh file
!+-----------------------------------------------------------------------------+
  open(unit=gmsh_file,file=input_filename,status='old')
  rewind gmsh_file


! Skip header, version info etc..
  do i=1,4
    read(gmsh_file,'(a)') aLine
  enddo

! You should have approached $PhysicalNames now
  read(gmsh_file,*) nbc
 
  allocate( iZone(nbc), iBc(nbc), bcName(nbc) )

  do i=1,nbc
    read(gmsh_file,*) iZone(i), iBc(i), bcName(i)
  enddo

  !$EndPhysicalNames
  !$Nodes
  do i=1,2
    read(gmsh_file,'(a)') aLine
  enddo

  ! Read number of nodes in the mesh
  read(gmsh_file,*) nonome

  write( *, '(6x,a)') 'NumPoints'
  write( *, '(1(1X,I9))') nonome
  write( *, '(a)' ) ' '

  ! Read nodal coordinates
  do i=1,nonome

    read(gmsh_file,*) NP,(XCOO(k),k=1,ndim)

    ! Write to file: 'points'
    write(7,'(3E20.11)') (xcoo(k), k=1,ndim)

  enddo

  ! CLOSE fCp format file: 'points'
  close (7)

  ! Skip rows '$EndNodes' and '$Elements'
  do i=1,2
    read(gmsh_file,*)
  enddo

  ! NELEM  - nuber of cells in the mesh
  read(gmsh_file,*) NELEM

  njump=0

  ! List length for the hash map
  listLength = nofaelmax*NELEM*2


  ! Read elements
  boundary_face_loop: do i=1,NELEM

!$Elements
! Format:
!number-of-elements
!elm-number elm-type number-of-tags < tag > … node-number-list
!…
!$EndElements
!
!Variable 	Description
!NE 	Global element number (not required to be sequential or continuous)
!NTYPE 	Element geometry type:
!	1 = 2-node line.
!	2 = 3-node triangle.
!	3 = 4-node quadrangle.
!	4 = 4-node tetrahedron.
!	5 = 8-node hexahedron.
!	6 = 6-node prism.
!	7 = 5-node pyramid.
!       ... (Whole list and descritpion is in README section: gmsh-mah-format.)
!NTAGS 	Number of tags for the element
!TAGS   List of tags for the element
!NODE 	List of nodes that define the element

    read(gmsh_file,*) NE, NTYPE, NTAGS, (TAGS(k), k=1,NTAGS), (NODE(k), k=1,noel(NTYPE))


    ! > Maybe we have face describing boundary conditions...
    if (NTYPE.eq.2 .or. NTYPE.eq.3) then ! Triangular or Quad boundary face:
      njump = njump+1
      write(8,*) TAGS(2),NTYPE,(NODE(k), k=1,noel(NTYPE))
    elseif(NTYPE.eq.1) then
      print*, 'We have hit an edge in the elements list, we need faces or cells!'
      stop
    else
      if (i.eq.1) njump = 0 
      exit boundary_face_loop
    endif

  end do boundary_face_loop

! Rewind file
  rewind gmsh_file
! Skip header
! Aproach the line where $Elements are:
  do
    read(gmsh_file,'(a)') aLine
    if (aLine(1:9).ne.'$Elements') then
      cycle
    else
      exit
    endif
  enddo

  ! > Go to the place in the file where cells start:
  do i=1,njump+1
    read(gmsh_file,*)
  enddo

  nel = NELEM-njump   ! The true number of cells
  nBoundFaces = njump ! Number of boundary faces just counted.

! > We can read cell now

  write( *, '(6x,a)') 'NumCells'
  write( *, '(1(1X,I9))') nel
  write( *, '(a)' ) ' '

  allocate( owner(listLength), neighbour(listLength),fV(nonofa,listLength) )

  ! Initialize arrays
  owner = 0
  neighbour = 0
  fV = 0

  fVert = 0

  numhex = 0
  numpri = 0
  numtet = 0
  numpyr = 0

  ! counter
  jface = 0 

  ! Initialize number of faces.
  nfacesTotal = 0
  nInnerFaces = 0

! Read elements
  element_loop: do iel=1,nel

!$Elements
! Format:
!$Elements
!number-of-elements
!elm-number elm-type number-of-tags < tag > … node-number-list
!…
!$EndElements
!
!Variable 	Description
!NE 	Global element number (not required to be sequential or continuous)
!NTYPE 	Element geometry type:
!	1 = 2-node line.
!	2 = 3-node triangle.
!	3 = 4-node quadrangle.
!	4 = 4-node tetrahedron.
!	5 = 8-node hexahedron.
!	6 = 6-node prism.
!	7 = 5-node pyramid.
!       ... (Whole list and descritpion is in README section: gmsh-mah-format.)
!NTAGS 	Number of tags for the element
!TAGS   List of tags for the element
!NODE 	List of nodes that define the element

  read(gmsh_file,*) NE, NTYPE, NTAGS, (TAGS(k), k=1,NTAGS), (NODE(k), k=1,noel(NTYPE))

  np = noel(NTYPE)

  node(1:np) = node(1:np) + 1 ! 1 based indexing of nodes, gmsh is zero based.

  nfacesTotal = nfacesTotal + nofael(NTYPE)

!
! > Write into 'cells' polyMesh file
!
  !if (NTYPE.eq.5) then
  !  write(12,'(I1,1X,8I8)') ntypeCappuccino(NTYPE),NODE(5),NODE(6),NODE(2),NODE(1),NODE(8),NODE(7),NODE(3),NODE(4) ! Order of nodes is important!
  !else
    write(12,'(I1,1X,4I8:(4I8:))') ntypeCappuccino(NTYPE),(NODE(k), k=1,noel(NTYPE))
  !endif


 if (NTYPE.eq.4) then
! ELEMENT is 4-node TETRAHEDRON
! Tetrahedron:                          Tetrahedron10:
! 
!                    v
!                  .
!                ,/
!               /
!            2                                     2                              
!          ,/|`\                                 ,/|`\                          
!        ,/  |  `\                             ,/  |  `\       
!      ,/    '.   `\                         ,6    '.   `5     
!    ,/       |     `\                     ,/       8     `\   
!  ,/         |       `\                 ,/         |       `\ 
! 0-----------'.--------1 --> u         0--------4--'.--------1
!  `\.         |      ,/                 `\.         |      ,/ 
!     `\.      |    ,/                      `\.      |    ,9   
!        `\.   '. ,/                           `7.   '. ,/     
!           `\. |/                                `\. |/       
!              `3                                    `3        
!                 `\.
!                    ` w

   do iface = 1,nofael_NTYPE4
     fVert = get_face_vertices_NTYPE4(iface,NODE,noel(NTYPE))
     call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
   enddo

   numtet = numtet + 1

 elseif (NTYPE.eq.5) then
! ELEMENT is 8-node HEX
! Hexahedron:             Hexahedron20:          Hexahedron27:
!
!        v
! 3----------2            3----13----2           3----13----2     
! |\     ^   |\           |\         |\          |\         |\    
! | \    |   | \          | 15       | 14        |15    24  | 14  
! |  \   |   |  \         9  \       11 \        9  \ 20    11 \  
! |   7------+---6        |   7----19+---6       |   7----19+---6 
! |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23| 
! 0---+---\--1   |        0---+-8----1   |       0---+-8----1   | 
!  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
!   \ |     \  \ |         10 |        12|        10 |  21    12| 
!    \|      w  \|           \|         \|          \|         \| 
!     4----------5            4----16----5           4----16----5 
!

   do iface = 1,nofael_NTYPE5
    fVert = get_face_vertices_NTYPE5(iface,NODE,noel(NTYPE))
    call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
   enddo

   numhex = numhex + 1

 elseif (NTYPE.eq.6) then
! ELEMENT is 6-node PRISM
!  Prism:                      Prism15:               Prism18:
!
!            w
!            ^
!            |
!            3                       3                      3        
!          ,/|`\                   ,/|`\                  ,/|`\      
!        ,/  |  `\               12  |  13              12  |  13    
!      ,/    |    `\           ,/    |    `\          ,/    |    `\  
!     4------+------5         4------14-----5        4------14-----5 
!     |      |      |         |      8      |        |      8      | 
!     |    ,/|`\    |         |      |      |        |    ,/|`\    | 
!     |  ,/  |  `\  |         |      |      |        |  15  |  16  | 
!     |,/    |    `\|         |      |      |        |,/    |    `\| 
!    ,|      |      `\        10     |      11       10-----17-----11
!  ,/ |      0      | `\      |      0      |        |      0      | 
! u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    | 
!     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  | 
!     |,/         `\|         |,/         `\|        |,/         `\| 
!     1-------------2         1------9------2        1------9------2 
!

   do iface = 1,nofael_NTYPE6
     fVert = get_face_vertices_NTYPE6(iface,NODE,noel(NTYPE))
     call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
   enddo

   numpri = numpri + 1


 elseif (NTYPE.eq.7) then
! ELEMENT is 5-node PYRAMID
!  Pyramid:                     Pyramid13:                   Pyramid14:
! 
!                4                            4                            4
!              ,/|\                         ,/|\                         ,/|\
!            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
!          ,/   | | \                   ,/   | | \                   ,/   | | \
!        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
!      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
!    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
!  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
! 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
!  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
!    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
!      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \ 
!        `\.'      `\     `\`          `\.'             `\`         `\.'             `\` 
!           1----------------2            1--------8-------2           1--------8-------2
!                     `\
!                        u
!

   do iface = 1,nofael_NTYPE7
      fVert = get_face_vertices_NTYPE7(iface,NODE,noel(NTYPE))
      call hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nInnerFaces,ne)
   enddo

   numpyr = numpyr + 1

 elseif (NTYPE.eq.2) then
    cycle ! We have hit the triangular boundary face.
 elseif (NTYPE.eq.3) then
     cycle ! We have hit the quad boundary face.
 else
      write(*,*) 'Fatal error: Non-existing cell type!'
      stop
 endif


 jface = jface + nofael(NTYPE)

end do element_loop

!
! Report
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Gmsh file:'
  write ( *, '(a,I7)' ) '  Total no. of HEX cells: ', numhex
  write ( *, '(a,I7)' ) '  Total no. of TET cells: ', numtet
  write ( *, '(a,I7)' ) '  Total no. of PYR cells: ', numpyr
  write ( *, '(a,I7)' ) '  Total no. of PRI cells: ', numpri
  write ( *, '(a)' )    ' +--------------------------------='
  write ( *, '(a,I7)' ) '  Total no. of cells: ', numhex+numtet+numpyr+numpri
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Normal end of reading .msh file.'
  write ( *, '(a)' ) ' '


!+-----------------------------------------------------------------------------+
! > END: Read input from Gmsh mesh file.


  ! End reading and processing data from .su2 mesh file.
  close(gmsh_file)

  nInnerFacePairs = jface-nBoundFaces
  nInnerFaces = nInnerFacePairs/2
  nFacesTotal = nInnerFaces + nBoundFaces
 
  write(*,*) "Inner faces: ",nInnerFaces,"Boundary faces: ",nBoundFaces,"Total: ",nfacesTotal

  write(12,'(I0,a)') nonome,      " Vertices"
  write(12,'(I0,a)') nel,         " Cells"
  write(12,'(I0,a)') nInnerFaces, " Inner faces"
  write(12,'(I0,a)') nBoundFaces, " Boundary faces"
  write(12,'(I0,a)') nfacesTotal, " Faces Total"

  ! Close polyMesh format files 
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)

!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Preprocessor:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

 end program