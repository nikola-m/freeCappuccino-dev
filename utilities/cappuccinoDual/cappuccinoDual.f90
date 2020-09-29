 program cappuccinoDual
!                            
!
! Description:
!
! Reads Gambit .neu file and prepares geometry data for Mean-Meadian Dual mesh
! discretization with freeCappuccino library.
!
! Author:
!   Nikola Mirkov (largeddysimulation@gmail.com)
! 
! Date:
!   July 2020
!

 use utils
 use qsort_c_module

 implicit none

 integer, parameter :: dp = kind(1.0d0)
 integer :: nonome ! no. of nodes in mesh
 integer, parameter :: nonoel = 8 ! no. of nodes in element-here Hex
 integer, parameter :: nonofa = 4 ! no, of nodes in element face-here Hex
 ! integer, parameter :: nofael_NTYPE4 = 6 ! no. of faces in element-here Hex
 ! integer, parameter :: nofael_NTYPE5 = 5 ! no. of faces in element-here Hex
 ! integer, parameter :: nofael_NTYPE6 = 4 ! no. of faces in element-here Tetrahedron
 ! integer, parameter :: nofael_NTYPE7 = 5 ! no. of faces in element-here Pyramid
 integer :: nel    ! no. of elements in mesh
 integer :: ndim   ! dimension of the problem. 3 for 3D.

 real(dp), parameter :: 1./3.
 real(dp), dimension(:,:), allocatable :: XCOO ! node coordinates
 real(dp), dimension(3) :: cc,e1,e2,e3,e4,e5,e6,f1,f2,f3,f4

 ! integer, dimension(:), allocatable :: cface_indx_start ! Where faces list start for each element
 ! integer, dimension(:,:), allocatable :: fVert ! size[6,value_of(cface(6,nel))] or [6, 6*nel]
 ! integer, dimension(:,:), allocatable :: fVUnsrt ! faces list that will not we sorted, see code below.
 ! integer, dimension(:,:), allocatable :: fV         ! faces list that will be sorted, see code below.
 ! character(len=32), dimension(:), allocatable :: bcName
 ! integer, dimension(:), allocatable :: bcType  
 ! integer, dimension(:), allocatable :: bcSize
 ! integer, dimension(:), allocatable :: owner
 ! integer, dimension(:), allocatable :: neighbour

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


! Gambit related
 integer :: NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL
 integer :: NP
 integer :: NE, NTYPE, NDP, NODE(8), elem(4)
 integer :: ITYPE, NENTRY, NVALUES, IBCODE1
 integer :: NGP, NELGP, MTYP, NFLAGS
 character( len = 82 ) :: inLine
 character (len = 32) :: string(4)



! 1) Intro
!+-----------------------------------------------------------------------------+
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                      _            ______             _  '
  write ( *, '(a)' ) '                                     (_)           |  _  \           | | '
  write ( *, '(a)' ) '  ___ __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___ | | | |_   _  __ _| | '
  write ( *, '(a)' ) " / __/ _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \| | | | | | |/ _` | | "
  write ( *, '(a)' ) '| (_| (_| | |_) | |_) | |_| | (_| (__| | | | | (_) | |/ /| |_| | (_| | | '
  write ( *, '(a)' ) ' \___\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/|___/  \__,_|\__,_|_| '
  write ( *, '(a)' ) '          | |   | |                                                      '
  write ( *, '(a)' ) '          |_|   |_|                                                      '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write ( *, '(a)' ) 'cappuccinoDual'
  write ( *, '(a)' ) '  A preprocessor program for the Cappuccino code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reads mesh in Gambit .neu format and produces  '
  write ( *, '(a)' ) '  mean median dual mesh files for freeCappuccino '
  write ( *, '(a)' ) '  solver.'
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

  ! Allocate arrays for nodal coordinates
  call allocate( xcoo(nonome, 3) )

! Read nodal coordinates
  do i=1,nonome

    read(4,'(I10,3E20.11)') NP,(XCOO(i,k),k=1,ndim)

    ! Write to file: 'points'
    write(7,'(3E20.11)') (XCOO(i,k), k=1,ndim)

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

  ! ! Initialize number of faces.
  ! nfacesTotal = 0
  ! nBoundFaces = 0
  ! nInnerFaces = 0

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

  ! nfacesTotal = nfacesTotal + nofael(NTYPE)

!
! > Write into 'cells' polyMesh file; 
!

! !
! !   NOTE: We write this for Paraview which starts counting from 0, therefore we substract 1 from a node number
! !

!   !write(12,'(7I8:(7I8:))') (NODE(k), k=1,NDP)

!   if (NTYPE.eq.4) then
  
!   write(12,'(I2,1X,8I8)') paraview_ntype(NTYPE), NODE(5)-1,NODE(6)-1,NODE(2)-1,NODE(1)-1,NODE(7)-1,NODE(8)-1,NODE(4)-1,NODE(3)-1 ! Order of nodes is important!
  
!   elseif (NTYPE.eq.7) then
  
!   write(12,'(I2,1X,5I8)') paraview_ntype(NTYPE),NODE(1)-1,NODE(2)-1,NODE(4)-1,NODE(3)-1,NODE(5)-1 ! Order of nodes is important!
  
!   else
  
!   write(12,'(I2,1X,4I8:(4I8:))') paraview_ntype(NTYPE),(NODE(k)-1, k=1,NDP)
  
!   endif



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

   ! ! How many faces before
   ! jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   ! do iface = 1,nofael_NTYPE4
   !   fVert(:, jface+iface ) = get_face_vertices_NTYPE4(iface,NODE,NDP)
   !   fVert(5,jface+iface) = i ! = element number
   ! enddo

   numhex = numhex + 1 

elseif (NTYPE.eq.5) then
! NTYPE6 is 6-node PRISM

   ! ! How many faces before
   ! jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   ! do iface = 1,nofael_NTYPE5
   !   fVert(:, jface+iface ) = get_face_vertices_NTYPE5(iface,NODE,NDP)
   !   fVert(5,jface+iface) = i ! = element number
   ! enddo

   numpri = numpri + 1


 elseif (NTYPE.eq.6) then
!
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
!      e1    'e2  `e5                            
!    ,/   f1  |     `\                   
!  ,/        f2    f3 `\                 
! 0-------e4--'.--------3 --> u         
!  `\.         |      ,/               
!     `\.   f4 |    ,e6                      
!        e3    '. ,/                              
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
 

  subcelvol = hexvolume( v1, e4, e1, f2, e3, f4, )

  ! How many faces before
  ! jface  = numhex*nofael_NTYPE4 + numpri*nofael_NTYPE5 + numtet*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

  ! do iface = 1,nofael_NTYPE6
  !  fVert(:, jface+iface ) = get_face_vertices_NTYPE6(iface,NODE,NDP)

  !  ! prikupi sve vertexe po njihovom redosledu, npr u vrtx nizu...
  !  ! xcoo(1:3,fVert(1:4)) = >

  ! enddo

  ! Primary cell centroid
  cc(1:3) = 0.25*sum( XCOO( NODE(1:4), 1:3), 1 )

  ! Edge centers
  e12(1:3) = 0.5*( XCOO( NODE(1), 1:3) + XCOO( NODE(2), 1:3) ) !(0,1)    
  e23(1:3) = 0.5*( XCOO( NODE(2), 1:3) + XCOO( NODE(3), 1:3) ) !(1,2)  
  e31(1:3) = 0.5*( XCOO( NODE(3), 1:3) + XCOO( NODE(1), 1:3) ) !(2,0)
  e14(1:3) = 0.5*( XCOO( NODE(1), 1:3) + XCOO( NODE(4), 1:3) ) !(0,3)   
  e24(1:3) = 0.5*( XCOO( NODE(2), 1:3) + XCOO( NODE(4), 1:3) ) !(1,3)
  e34(1:3) = 0.5*( XCOO( NODE(3), 1:3) + XCOO( NODE(4), 1:3) ) !(2,3)

  ! Face centers
  f213(1:3) = third*( XCOO( NODE(2), 1:3) + XCOO( NODE(1), 1:3) + XCOO( NODE(3), 1:3) ) !1,0,2
  f124(1:3) = third*( XCOO( NODE(1), 1:3) + XCOO( NODE(2), 1:3) + XCOO( NODE(4), 1:3) ) !0,1,3
  f234(1:3) = third*( XCOO( NODE(2), 1:3) + XCOO( NODE(3), 1:3) + XCOO( NODE(4), 1:3) ) !1,2,3
  f314(1:3) = third*( XCOO( NODE(3), 1:3) + XCOO( NODE(1), 1:3) + XCOO( NODE(4), 1:3) ) !2,0,3

  ! Zatim subface center (sfc) na osnovu cc, edge centers i adjecent face centers - subface je quad i moras imati
  !  rutinu za centroid quada sa datim vertexima.

  ! Sa istim vertexima pravis is sfn - sub-face normal vektor, Ima ih 6  - za svaki edge koji spaja
  ! da temena. Subface se sastoji id edge centra, cc, i da susjedna face centra koje si sve vec nasao

  ! Imas i 4 sektorska cc, za doprinose zapr integralu za svaki vertex
  ! scelcc(1:3) (sector cell centre x,y and z coodinates)
  scelcc(1:3) = 17/96*(45/17*ax+bx+cx+dx, ... )

  ! Volume of hexaheral subcel - subcvol(1:4)
  ! Njih radis preko formule centroida za hex.
  v0= XCOO( NODE(1), 1:3)

  subcelvol = hexvolume( v0, )


   ! rezime
   ! imas niz duzine broj elemenata
   ! tu ima sest subfejsova
   ! svaki subfejs ima svoj sub-face center sfc(1:3) koordinate i ima sub-face-normal sfn(1:3) koponente povrsine izrazene preko normale

   ! na kraju mi trebaju sledeci podaci
   ! sfn(1:6,1:3) - vezan za edge koji spaja dva noda
   ! sfc(1:6,1:3) - vezan za edge koji spaja dva noda
   ! scelcc(1:4, 1:3) - vezan za cell node iliti vertex
   ! scelvol(1:4) - vezan za svaki vertex

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
  ! open(unit=8,file='boundary')
  ! open(unit=9,file='faces')
  ! open(unit=10,file='owner')
  ! open(unit=11,file='neighbour')

  ! rewind 8  
  ! rewind 9
  ! rewind 10
  ! rewind 11

!
! Element groupds and Boundary conditions
!

  ! ! Aproach the line where BOUNDARY CONDITIONS and ELEMENT GROUP sections are.
  ! bc_loop: do

  !   read(4,'(a)', iostat = ios ) inLine
  !   if(ios /= 0) then
  !       exit bc_loop  
  !   endif

  !   if ( adjustl(inLine(1:26)) == 'BOUNDARY CONDITIONS') then


  !     ! Now read the info to enter bc_loop 
  !     read(4,'(A32, 4I10)') inLine, ITYPE, NENTRY, NVALUES, IBCODE1
  !     write(8,'(A32,4I10)') inLine, ITYPE, NENTRY, NVALUES, IBCODE1

  !     ! bc_read_loop: do

  !       ! add to total number of boundary faces
  !       nBoundFaces = nBoundFaces + NENTRY    

  !       ! if(ibcode1.eq.6) then 
  !       ! ELEMENT_SIDE
  !       do i=1,NENTRY
  !         read(4,*) iel, NTYPE, iface
  !         indx = cface_indx_start(iel) + iface 
          
  !         ! Write to boundary file
  !         if(fVert(4,indx) == 0) then ! triangle face
  !           write(8,'(i0,1x,a,1x,3(i0,1x))') iel,' 3',fVert(1:3,indx)
  !         else                        ! quad face
  !           write(8,'(i0,1x,a,1x,4(i0,1x))') iel,' 4', fVert(1:4,indx)
  !         endif

  !         ! Mark face as boundary face: on position 6 put 1.
  !         fVert(6,indx) = 1

  !       enddo

  !       read(4,'(a)') inLine ! ENDOFSECTION string

  !     cycle bc_loop

  !   elseif ( adjustl(inLine(1:20)) == 'ELEMENT GROUP') then ! We cut last characters which hold the file version...

  !     ! Read element group data
  !     read(4, '(a)' ) inLine
  !     read( inLine(8:17), * ) NGP
  !     read( inLine(29:38), * ) NELGP
  !     read( inLine(50:50), * ) MTYP    
  !     read( inLine(60:69), * ) NFLAGS  

  !     write(string(1),*) NGP
  !     open(unit=13,file='cellgroup-'//adjustl(string(1)) )
  !     rewind 13

  !     write(13,'(2I8,a)') NELGP,MTYP,' NELGP,MTYP'

  !     ! For next two lines I cannot still find use.
  !     read(4,'(a)') inLine
  !     read(4,'(a)') inLine

  !     ! NELGP/10 lines with 10 numbers
  !     do i=1, NELGP/10
  !       read(4,'(10I8)') (elem(k), k=1,10)
  !       do k=1,10
  !         write(13,'(I8)') elem(k)
  !       enddo
  !     enddo

  !     ! One last line with less then 10 elements
  !     indx = NELGP - 10*(NELGP/10) 
  !     if (indx /= 0) then
  !       read(4,'(10I8)') (elem(k), k=1,indx)
  !       do k=1,indx
  !         write(13,'(I8)') elem(k)
  !       enddo
  !     endif

  !     read(4,'(a)') inLine ! ENDOFSECTION string

  !     close(13)

  !     cycle bc_loop

  !   else

  !     cycle bc_loop

  !   endif

  ! enddo bc_loop


!+-----------------------------------------------------------------------------+
! > END: Read input from Gambit mesh file.


!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cappuccinoDual:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

 end program


contains

function volpart(a1,a2,a3,b1,b2,b3,q1,q2,q3) result(volpart)

  implicit none

  real(dp), intent(in) :: a1,a2,a3,b1,b2,b3,q1,q2,q3
  real(dp) :: volpart

    volpart = 1./6.*(a2*b3-b2*a3)*q1+(b1*a3-a1*b3)*q2+(a1*b2-a2*b1)*q3 

end function

function hexvolume() result (volum)
!
! Cell subvolumes are hex cell, therefoe we use a function which
! calculates the volume of a hex.
!
!   v2----------------v3
!   /|               /|
!  / |              / |
!v6----------------v7 |
! |  |             |  |
! |  |             |  |
! |  |             |  |
! | v0----------------v1
! | /              | /
! |/               |/
!v4----------------v5
!
!   imjk-------------ijk
!   /|               /|
!  / |              / |
!imjmk-------------ijmk 
! |  |             |  |
! |  |             |  |
! |  |             |  |
! | ijkm-nj----------ijkm
! | /              | /
! |/               |/
!ijka--------------ijkm-1
!
!
! IJK=LIK+J
! IMJK=IJK-NJ
! IMJMK=IMJK-1
! IJMK=IJK-1
! IJKM=IJK-NIJ
! IJKA=IJKM-NJ-1
!

  ! XA=X(IJKA) ! v4
  ! YA=Y(IJKA)
  ! ZA=Z(IJKA)
  ! DXAD=X(IJMK)-XA ! v7-v4
  ! DYAD=Y(IJMK)-YA
  ! DZAD=Z(IJMK)-ZA
  ! DXAB=X(IJKM-1)-XA ! v5-v4
  ! DYAB=Y(IJKM-1)-YA
  ! DZAB=Z(IJKM-1)-ZA
  ! DXAC=X(IMJMK)-XA ! v6-v4
  ! DYAC=Y(IMJMK)-YA
  ! DZAC=Z(IMJMK)-ZA
  ! DXAE=X(IJKM-NJ)-XA ! v0 - v4
  ! DYAE=Y(IJKM-NJ)-YA
  ! DZAE=Z(IJKM-NJ)-ZA
  ! DXAF=X(IJKM)-XA ! v1-v4
  ! DYAF=Y(IJKM)-YA
  ! DZAF=Z(IJKM)-ZA
  ! DXAG=X(IMJK)-XA !v2-v4
  ! DYAG=Y(IMJK)-YA
  ! DZAG=Z(IMJK)-ZA
  ! DXAH=X(IJK)-XA ! v3-v4
  ! DYAH=Y(IJK)-YA
  ! DZAH=Z(IJK)-ZA


  XA=X(IJKA) ! v4
  YA=Y(IJKA)
  ZA=Z(IJKA)
  DXAD=X(IJMK)-XA ! v7-v4
  DYAD=Y(IJMK)-YA
  DZAD=Z(IJMK)-ZA
  DXAB=X(IJKM-1)-XA ! v5-v4
  DYAB=Y(IJKM-1)-YA
  DZAB=Z(IJKM-1)-ZA
  DXAC=X(IMJMK)-XA ! v6-v4
  DYAC=Y(IMJMK)-YA
  DZAC=Z(IMJMK)-ZA
  DXAE=X(IJKM-NJ)-XA ! v0 - v4
  DYAE=Y(IJKM-NJ)-YA
  DZAE=Z(IJKM-NJ)-ZA
  DXAF=X(IJKM)-XA ! v1-v4
  DYAF=Y(IJKM)-YA
  DZAF=Z(IJKM)-ZA
  DXAG=X(IMJK)-XA !v2-v4
  DYAG=Y(IMJK)-YA
  DZAG=Z(IMJK)-ZA
  DXAH=X(IJK)-XA ! v3-v4
  DYAH=Y(IJK)-YA
  DZAH=Z(IJK)-ZA
  VOLUM=volpart(DXAD,DYAD,DZAD,DXAB,DYAB,DZAB,DXAH,DYAH,DZAH)+ &
        volpart(DXAC,DYAC,DZAC,DXAD,DYAD,DZAD,DXAH,DYAH,DZAH)+ &
        volpart(DXAG,DYAG,DZAG,DXAC,DYAC,DZAC,DXAH,DYAH,DZAH)+ &
        volpart(DXAB,DYAB,DZAB,DXAF,DYAF,DZAF,DXAH,DYAH,DZAH)+ &
        volpart(DXAF,DYAF,DZAF,DXAE,DYAE,DZAE,DXAH,DYAH,DZAH)+ &
        volpart(DXAE,DYAE,DZAE,DXAG,DYAG,DZAG,DXAH,DYAH,DZAH)

  end function

