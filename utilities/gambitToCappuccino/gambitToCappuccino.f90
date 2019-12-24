 program gambitToCappuccino
!
!                       _     _ _   ___   _____                                 _              
!                      | |   (_) | |__ \ / ____|                               (_)             
!  __ _  __ _ _ __ ___ | |__  _| |_   ) | |     __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___   
! / _` |/ _` | '_ ` _ \| '_ \| | __| / /| |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \  
!| (_| | (_| | | | | | | |_) | | |_ / /_| |___| (_| | |_) | |_) | |_| | (_| (__| | | | | (_) | 
! \__, |\__,_|_| |_| |_|_.__/|_|\__|____|\_____\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/  
!  __/ |                                            | |   | |                                                               
! |___/                                             |_|   |_|                                  
!
! Description:
!
! Converts Gambit .neu file to those compatible with freeCappuccino solver.
!
! Author:
!   Nikola Mirkov (largeddysimulation@gmail.com)
! 
! Date:
!   5. November 2015.
! Modification:
!   23. November 2019. (November is coding time!)
!

 use utils
 use qsort_c_module

 implicit none

 integer, parameter :: dp = kind(1.0d0)
 integer :: nonome ! no. of nodes in mesh
 integer, parameter :: nonoel = 8 ! no. of nodes in element-here Hex
 integer, parameter :: nonofa = 4 ! no, of nodes in element face-here Hex
 integer, parameter :: nofaelmax     = 6 ! no. of faces in element
 integer, parameter :: nofael_NTYPE4 = 6 ! no. of faces in element-here Hex
 integer, parameter :: nofael_NTYPE5 = 5 ! no. of faces in element-here Hex
 integer, parameter :: nofael_NTYPE6 = 4 ! no. of faces in element-here Tetrahedron
 integer, parameter :: nofael_NTYPE7 = 5 ! no. of faces in element-here Pyramid
 integer :: nel    ! no. of elements in mesh
 integer :: ndim   ! dimension of the problem. 3 for 3D.
 integer, dimension(:), allocatable :: cface_indx_start ! Where faces list start for each element
 integer, dimension(:,:), allocatable :: fVert ! size[6,value_of(cface(6,nel))] or [6, 6*nel]
 integer, dimension(:,:), allocatable :: fVUnsrt ! faces list that will not we sorted, see code below.
 integer, dimension(:,:), allocatable :: fV         ! faces list that will be sorted, see code below.
 character(len=32), dimension(:), allocatable :: bcName
 integer, dimension(:), allocatable :: bcType  
 integer, dimension(:), allocatable :: bcSize
 integer, dimension(:), allocatable :: owner
 integer, dimension(:), allocatable :: neighbour

! Locals
 integer :: i,k
 integer :: iel, jel
 integer :: indx
 integer :: iface, jface, imatch, ifound
 integer :: numhex ! No. of Hex elements 
 integer :: numpri ! No. of Prism elements
 integer :: numtet ! No. of Tet elements 
 integer :: numpyr ! No. of Pyr elements 
 integer :: numCells

 integer :: iarg
 integer :: iargc
 integer :: ios
 integer :: num_arg
 character ( len = 255 ) prefix
 character ( len = 255 ) input_filename

 integer :: nfacesTotal
 integer :: nBoundFaces 
 integer :: nInnerFaces
 integer :: nInnerFacePairs
 integer :: nBoundary
 integer :: numVrtx 

! Gambit related
 integer :: NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL
 integer :: NP
 integer :: NE, NTYPE, NDP, NODE(8), elem(4)
 integer :: ITYPE, NENTRY, NVALUES, IBCODE1
 integer :: NGP, NELGP, MTYP, NFLAGS
 character( len = 82 ) :: inLine
 character (len = 32) :: string(4)
 real(dp), dimension(3) :: XCOO


! 1) Intro
!+-----------------------------------------------------------------------------+
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                       _     _ _   ___   _____                                 _              '
  write ( *, '(a)' ) '                      | |   (_) | |__ \ / ____|                               (_)             '
  write ( *, '(a)' ) '  __ _  __ _ _ __ ___ | |__  _| |_   ) | |     __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___   '
  write ( *, '(a)' ) " / _` |/ _` | '_ ` _ \| '_ \| | __| / /| |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \  "
  write ( *, '(a)' ) '| (_| | (_| | | | | | | |_) | | |_ / /_| |___| (_| | |_) | |_) | |_| | (_| (__| | | | | (_) | '
  write ( *, '(a)' ) ' \__, |\__,_|_| |_| |_|_.__/|_|\__|____|\_____\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/  '
  write ( *, '(a)' ) '  __/ |                                            | |   | |                                  '                             
  write ( *, '(a)' ) ' |___/                                             |_|   |_|                                  '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write ( *, '(a)' ) 'gambitToCappuccino'
  write ( *, '(a)' ) '  A preprocessor program for the Cappuccino code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reads mesh in Gambit .neu format and produces  '
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
  input_filename = trim ( prefix ) // '.neu'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is "' // trim ( input_filename ) // '".'
  write ( *, '(a)' ) ' '
!+-----------------------------------------------------------------------------+



! > Read input from Gambit mesh file (.neu)
!+-----------------------------------------------------------------------------+
  open(unit=4,file=input_filename,status='old')
  rewind 4

! Skip header
  do i=1,6
    read(4,*)
  enddo

! First thing we read is below line:
! NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
  read(4,'(6(1X,I9))') NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL

  write( *, '(a)') '      NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL'
  write( *, '(6(1X,I9))') NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL
  write( *, '(a)' ) ' '

  nonome = NUMNP
  nel = NELEM
  ndim = NDFCD
  if( ndim .ne. 3 ) then
      write(*,*) 'Fatal error: Only 3D meshes are accepted!'
      stop
  end if 

  ! fVert stores four (or three) vertices that make a face, 
  ! an element number the face belongs to, and 0/1 depending 
  ! if face turns out to be a boundary face.
  allocate ( fVert(nonofa+2,nofaelmax*nel) )
  allocate ( cface_indx_start(nel+1) )

! Initialize
  cface_indx_start(1) = 0

! Skip rows:
! ENDOFSECTION
!    NODAL COORDINATES
  do i=1,2
    read(4,*)
  enddo

! OPEN text file: 'points'
  open(unit=7,file='points')
  rewind 7

  ! Write size of cells arrays
  write(7, '(i8,a)') nonome, ' numPoints'

! Read nodal coordinates
  do i=1,nonome

    read(4,'(I10,3E20.11)') NP,(XCOO(k),k=1,ndim)

    ! Write to file: 'points'
    write(7,'(3E20.11)') (XCOO(k), k=1,ndim)

  enddo

! CLOSE file: 'points'
  close (7)

! Skip rows
! ENDOFSECTION
!       ELEMENTS/CELLS
  do i=1,2
    read(4,*)
  enddo

  ! Initialize total number of hex cells, prism cells, etc.
  numhex = 0
  numpri = 0
  numtet = 0
  numpyr = 0

  ! Initialize number of faces.
  nfacesTotal = 0
  nBoundFaces = 0
  nInnerFaces = 0

! OPEN file: 'cells',
  open(unit=12,file='cells')
  rewind 12

  ! Write size of cells arrays
  write(12, '(i0,a)') nel, ' numCells'

! Read elements
  element_loop: do i=1,nel

!ELEMENTS/CELLS
!Variable 	Description
!NE 	Global element number (not required to be sequential or continuous)
!NTYPE 	Element geometry type:
!	1 = Edge
!	2 = Quadrilateral
!	3 = Triangle
!	4 = Brick
!	5 = Wedge (Prism)
!	6 = Tetrahedron
!	7 = Pyramid
!NDP 	Number of nodes that define the element
!NODE 	List of nodes that define the element

  read(4,'(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))') NE, NTYPE, NDP, (NODE(k), k=1,NDP)

  cface_indx_start(i+1) = cface_indx_start(i) + nofael(NTYPE)

  nfacesTotal = nfacesTotal + nofael(NTYPE)

!
! > Write into 'cells' polyMesh file; 
!

!
!   NOTE: We write this for Paraview which starts counting from 0, therefore we substract 1 from a node number
!

  !write(12,'(7I8:(7I8:))') (NODE(k), k=1,NDP)

  if (NTYPE.eq.4) then
  
  write(12,'(I2,1X,8I8)') paraview_ntype(NTYPE), NODE(5)-1,NODE(6)-1,NODE(2)-1,NODE(1)-1,NODE(7)-1,NODE(8)-1,NODE(4)-1,NODE(3)-1 ! Order of nodes is important!
  
  elseif (NTYPE.eq.7) then
  
  write(12,'(I2,1X,5I8)') paraview_ntype(NTYPE),NODE(1)-1,NODE(2)-1,NODE(4)-1,NODE(3)-1,NODE(5)-1 ! Order of nodes is important!
  
  else
  
  write(12,'(I2,1X,4I8:(4I8:))') paraview_ntype(NTYPE),(NODE(k)-1, k=1,NDP)
  
  endif



 if (NTYPE.eq.4) then
! NTYPE4 is 8-node HEX

!Element Type and Node-Numbering Conventions
!	   2----------------3
!	  /|               /|
!	 / |              / |
!	6----------------7  |
!	|  |             |  |
!	|  |             |  |
!	|  |             |  |
!	|  0----------------1
!	| /              | /
!	|/               |/
!	4----------------5
!
!Brick, 8-Node
!Edge 	Nodes 	Face 	Nodes
!1 	0,4 	1 	0,1,5,4
!2 	0,1 	2 	1,3,7,5
!3 	1,5 	3 	3,2,6,7
!4 	4,5 	4 	2,0,4,6
!5 	1,3 	5 	1,0,2,3
!6 	3,7 	6 	4,5,7,6
!7 	5,7 		
!8 	2,3 		
!9 	2,6 		
!10 	6,7 		
!11 	0,2 		
!12 	4,6

! All this is packed in function:

   ! How many faces before
   jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   do iface = 1,nofael_NTYPE4
     fVert(:, jface+iface ) = get_face_vertices_NTYPE4(iface,NODE,NDP)
     fVert(5,jface+iface) = i ! = element number
   enddo

   numhex = numhex + 1 

elseif (NTYPE.eq.5) then
! NTYPE6 is 6-node PRISM

   ! How many faces before
   jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   do iface = 1,nofael_NTYPE5
     fVert(:, jface+iface ) = get_face_vertices_NTYPE5(iface,NODE,NDP)
     fVert(5,jface+iface) = i ! = element number
   enddo

   numpri = numpri + 1


 elseif (NTYPE.eq.6) then
! NTYPE6 is 4-node TETRAHEDRON
! Tetrahedron:                 
! 
!                    v
!                  .
!                ,/
!               /
!            1                                                                  
!          ,/|`\                                                          
!        ,/  |  `\                                  
!      ,/    '.   `\                            
!    ,/       |     `\                   
!  ,/         |       `\                 
! 0-----------'.--------3 --> u         
!  `\.         |      ,/               
!     `\.      |    ,/                      
!        `\.   '. ,/                              
!           `\. |/                            
!              `2                                       
!                 `\.
!                    ` w
! Ede Nodes Face Nodes
! 1 0,1    1   1,0,2
! 2 1,2    2   0,1,3
! 3 2,0    3   1,2,3
! 4 0,3    4   2,0,3
! 5 1,3 
! 6 2,3
 
   ! How many faces before
   jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   do iface = 1,nofael_NTYPE6
     fVert(:, jface+iface ) = get_face_vertices_NTYPE6(iface,NODE,NDP)
     fVert(5,jface+iface) = i ! = element number
   enddo

   numtet = numtet + 1

 elseif (NTYPE.eq.7) then
! NTYPE7 is 5-node PYRAMID

   ! How many faces before
   jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   do iface = 1,nofael_NTYPE7
     fVert(:, jface+iface ) = get_face_vertices_NTYPE7(iface,NODE,NDP)
     fVert(5,jface+iface) = i ! = element number
   enddo

   numpyr = numpyr + 1

 else
      write(*,*) 'Fatal error: Non-existing cell type!'
      stop
 endif



end do element_loop

!
! > CLOSE polyMesh format file: 'cells',
!
  close(12)

!
! Report
!

  numCells = numhex+numtet+numpyr+numpri

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) ' Gambit file:'
  write ( *, '(4x,a,I7)' ) 'Total no. of HEX cells: ', numhex
  write ( *, '(4x,a,I7)' ) 'Total no. of TET cells: ', numtet
  write ( *, '(4x,a,I7)' ) 'Total no. of PYR cells: ', numpyr
  write ( *, '(4x,a,I7)' ) 'Total no. of PRI cells: ', numpri
  write ( *, '(3x,a)' )    '+--------------------------------='
  write ( *, '(4x,a,I7)' ) 'Total no. of cells: ', numCells
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) 'Normal end of reading .neu file.'
  write ( *, '(a)' ) ' '

!
! > OPEN polyMesh format file: 'boundary','faces', 'owner', 'neighbour'.
!
  open(unit=8,file='boundary')
  open(unit=9,file='faces')
  open(unit=10,file='owner')
  open(unit=11,file='neighbour')

  rewind 8  
  rewind 9
  rewind 10
  rewind 11

!
! Element groupds and Boundary conditions
!

  ! Aproach the line where BOUNDARY CONDITIONS and ELEMENT GROUP sections are.
  bc_loop: do

    read(4,'(a)', iostat = ios ) inLine
    if(ios /= 0) then
        exit bc_loop  
    endif

    if ( adjustl(inLine(1:26)) == 'BOUNDARY CONDITIONS') then


      ! Now read the info to enter bc_loop 
      read(4,'(A32, 4I10)') inLine, ITYPE, NENTRY, NVALUES, IBCODE1
      write(8,'(A32,4I10)') inLine, ITYPE, NENTRY, NVALUES, IBCODE1

      ! bc_read_loop: do

        ! add to total number of boundary faces
        nBoundFaces = nBoundFaces + NENTRY    

        ! if(ibcode1.eq.6) then 
        ! ELEMENT_SIDE
        do i=1,NENTRY
          read(4,*) iel, NTYPE, iface
          indx = cface_indx_start(iel) + iface 
          
          ! Write to boundary file
          if(fVert(4,indx) == 0) then ! triangle face
            write(8,'(i0,1x,a,1x,3(i0,1x))') iel,' 3',fVert(1:3,indx)
          else                        ! quad face
            write(8,'(i0,1x,a,1x,4(i0,1x))') iel,' 4', fVert(1:4,indx)
          endif

          ! Mark face as boundary face: on position 6 put 1.
          fVert(6,indx) = 1

        enddo

        read(4,'(a)') inLine ! ENDOFSECTION string

      cycle bc_loop

    elseif ( adjustl(inLine(1:20)) == 'ELEMENT GROUP') then ! We cut last characters which hold the file version...

      ! Read element group data
      read(4, '(a)' ) inLine
      read( inLine(8:17), * ) NGP
      read( inLine(29:38), * ) NELGP
      read( inLine(50:50), * ) MTYP    
      read( inLine(60:69), * ) NFLAGS  

      write(string(1),*) NGP
      open(unit=13,file='cellgroup-'//adjustl(string(1)) )
      rewind 13

      write(13,'(2I8,a)') NELGP,MTYP,' NELGP,MTYP'

      ! For next two lines I cannot still find use.
      read(4,'(a)') inLine
      read(4,'(a)') inLine

      ! NELGP/10 lines with 10 numbers
      do i=1, NELGP/10
        read(4,'(10I8)') (elem(k), k=1,10)
        do k=1,10
          write(13,'(I8)') elem(k)
        enddo
      enddo

      ! One last line with less then 10 elements
      indx = NELGP - 10*(NELGP/10) 
      if (indx /= 0) then
        read(4,'(10I8)') (elem(k), k=1,indx)
        do k=1,indx
          write(13,'(I8)') elem(k)
        enddo
      endif

      read(4,'(a)') inLine ! ENDOFSECTION string

      close(13)

      cycle bc_loop

    else

      cycle bc_loop

    endif

  enddo bc_loop


!+-----------------------------------------------------------------------------+
! > END: Read input from Gambit mesh file.

  ! Number of inner face pairs.
  nInnerFacePairs = nFacesTotal-nBoundFaces 

  ! Backup face vertices in two arrays, one will be sorted - fV, other unsorted-fVUnsrt
  ! Note: 6 = nonofa + 2  
  allocate ( fV(6,nInnerFacePairs) )
  allocate ( fVUnsrt(5,nInnerFacePairs) )

  ! Populate these backup arrays with only inner faces.
  iface = 1
  do i = 1,nfacesTotal
    if(iface > nInnerFacePairs ) exit   
    if ( fVert(6,i) == 0 ) then ! it's inner face
      ! Backup face vertex indices and owner cell index.
      fV(1:5,iface) = fVert(1:5,i)    
      iface = iface + 1        
    endif
  enddo


  ! Backup unsorted
  fVUnsrt = fV(1:5,:)

  ! Free big array with a lot of zeros
  deallocate(fVert)

  write ( *, '(2x,a)' ) 'Sort fVert: Begin'
    do i = 1,nInnerFacePairs
      call sortIntArray(fV(1:4,i),4)
    enddo
  write ( *, '(2x,a)' ) 'Sort fVert: End'
  write ( *, '(a)' ) ' '


!
! Sortiram po drugom, trecem, i cetvrtom indeksu. 
! Tri su dovoljna da bi se utvrdilo poklapanje dve celije po jednom licu. 
!

! > QUICKSORT over vertices:

  write ( *, '(2x,a)' ) 'Quicksort fVert: Begin'

    ! Do quicksort on a first vertex  - should speed things up considerably
    call QsortC(fV(2,:), &
                fV(3,:),fV(4,:),fV(1,:),fV(5,:), &
                fVUnsrt(1,:),fVUnsrt(2,:),fVUnsrt(3,:),fVUnsrt(4,:))

    ! Sada sortiramo samo po drugom indeksu po grupama gde imaju isti prvi indeks
    ! kolika ce grupa biti zavisi, zato povecavamo k koliko mozemo

    i=1
    k=i+1
    do

      if( i > nInnerFacePairs ) exit

      ! List of faces that have the same index 2 grows by one 
      if(fV(2,i)==fV(2,k)) then    
        k=k+1
     
      ! Or if it doesn't grow anymore, we can now quick sort it.
      else
      
        ! Ready to Go:
        call QsortC(fV(3,i:k-1), &
                    fV(4,i:k-1),fV(1,i:k-1),fV(2,i:k-1),fV(5,i:k-1), &
                    fVUnsrt(1,i:k-1),fVUnsrt(2,i:k-1),fVUnsrt(3,i:k-1),fVUnsrt(4,i:k-1))
        i=k
        k=i+1
      
        cycle
     
      endif

    enddo


    ! Sada sortiramo samo po trecem po redu indeksu po grupama gde imaju isti prvi i drugi  indeks
    ! kolika ce grupa biti zavisi, zato povecavamo k koliko mozemo
    i=1
    k=i+1
    do

      if( i > nInnerFacePairs) exit

      ! List of faces that have the same index 3 grows by one 
      if(fV(3,i)==fV(3,k)) then
      k=k+1

      ! Or if it doesn't grow anymore, we can now quick sort it.
      else

        ! Ready to Go:
        call QsortC(fV(4,i:k-1), &
                  fV(1,i:k-1),fV(2,i:k-1),fV(3,i:k-1),fV(5,i:k-1), &
                  fVUnsrt(1,i:k-1),fVUnsrt(2,i:k-1),fVUnsrt(3,i:k-1),fVUnsrt(4,i:k-1))
        
        i=k
        k=i+1

        cycle 

      endif

    enddo

! > END: QUICKSORT over vertices.

  write ( *, '(2x,a)' ) 'Quicksort fVert: End'

!=============================================================================

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) 'Match faces, find owner neighbour: Start'


  ! Write number of elements to be written in every file.
  nFacesTotal = nInnerFacePairs/2 + nBoundFaces
  nInnerFaces = nInnerFacePairs/2

  allocate(owner(nInnerFaces))
  allocate(neighbour(nInnerFaces))

  ifound = 0

  iface_loop: do iface=1,nInnerFacePairs-1

      ! Lets check if this face has already been matched with its pair 
      if (fV(6,iface) == 1) cycle iface_loop

      imatch_loop: do imatch = iface+1,nInnerFacePairs

        if ( fV(4,iface) == fV(4,imatch) ) then
          if ( fV(3,iface) == fV(3,imatch) ) then
            if ( fV(2,iface) == fV(2,imatch) ) then
             
              ! We found match - so lets put one on place 6
              fV(6,iface) = 1
              fV(6,imatch) = 1

              iel = fV(5,iface)
              jel = fV(5,imatch)

              ifound = ifound + 1

              if (iel < jel) then

                owner(ifound) = iel
                neighbour(ifound) = jel
                fV(1:4,ifound) = fVUnsrt(1:4,iface)

              else

                owner(ifound) = jel
                neighbour(ifound) = iel 
                fV(1:4,ifound) = fVUnsrt(1:4,imatch)    

              endif
     
              exit imatch_loop

            endif
          endif
        endif

      enddo imatch_loop 

  enddo iface_loop
!=============================================================================!

  write ( *, '(2x,a)' ) 'Match faces, find owner neighbour: End'

  call QsortC2( owner(:), neighbour(:), &
                fV(1,1:nInnerFaces), fv(2,1:nInnerFaces), fv(3,1:nInnerFaces), fV(4,1:nInnerFaces))

  !
  ! Write to files
  !

  ! write owner size:
  write(10,'(i0,a)') nFacesTotal, ' numFaces'

  ! write neighbour szie
  write(11,'(i0,a)') nInnerFaces, ' numInnerFaces'

  ! write face face file size:
  write(9,'(i0,a)') nFacesTotal, ' numFaces'

  do iface=1,nInnerFaces

    numVrtx = 4
    if ( fV(4,iface) == 0 ) numVrtx = 3 ! tj. kad je face trougao.

    ! write owner:
    write(10,'(i0)') owner(iface)
    ! write neighbour
    write(11,'(i0)') neighbour(iface)
    ! write face into the face file:
    write(9,'(5(i0,1x))')  numVrtx, fV(1:numVrtx,iface)
  enddo


!
! This part is not very efficient because we have slow process of rewritting
! from one file to the other - but it does what we need. Maybe in the future
! if anyone wants, we can improve this part..
!

! Move data from boundary file to owner and faces files, just keep BC info
  rewind 8

  nBoundary = NBSETS

  allocate(bcName(nBoundary))
  allocate(bcType(nBoundary))  
  allocate(bcSize(nBoundary)) 

  do k=1,nBoundary

  ! Now read the infor to enter bc_loop 
  read(8,'(A32, 2I10)') bcName(k), bcType(k), bcSize(k)

    do i=1,bcSize(k)

      read(8,*) iel,numVrtx,fV(1:numVrtx,1)

        ! write owner:
        write(10,'(i0)') iel

        ! write face into the face file:
        write(9,'(5(i0,1x))')  numVrtx, fV(1:numVrtx,1)


    enddo

  enddo

! Rewrite boundary file

  rewind 8

  write(8,'(a)') '# bcName bcType nFaces startFace'

  ! NOTE: Below bcType is an integer - if correct BCtype is designated in
  ! mesh generator you can use table of integer codes of BC types which
  ! is documented in GAMBIT documentation. We may also just use a function
  ! which transforms integer value to proper string, eg. 'wall', 'inlet', etc.
  ! that can be a return value from a function like bcTypeToString( bcType(k) )
  do k=1,nBoundary

    write(8,'(a,4(i0,1x))') adjustl(bcName(k)), bcType(k), bcSize(k), nInnerFaces

    nInnerFaces = nInnerFaces+bcSize(k)

  enddo

!
!  > CLOSE polyMesh format file: 'boundary', 'faces', 'owner', 'neighbour'.
!
  close(8)
  close(9)
  close(10)
  close(11)
  close(4)


!
!  > Administrative tasks.
!
  deallocate ( fV )
  deallocate ( fVUnsrt )
  deallocate ( cface_indx_start )
  deallocate( owner )
  deallocate( neighbour )
  deallocate(bcName)
  deallocate(bcType)  
  deallocate(bcSize) 

!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'gambitToCappuccino:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

 end program gambitToCappuccino

