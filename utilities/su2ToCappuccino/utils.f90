module utils
contains



 pure integer function noel(NTYPE)
!
! Number of nodes in element of NTYPE su2 type of element
!
 implicit none
 integer, intent(in) :: NTYPE
 

 if(NTYPE.eq.3) then      ! -> su2 Tri
   noel = 2
 elseif(NTYPE.eq.5) then  ! -> su2 Tri
   noel = 3 
 elseif(NTYPE.eq.9) then  ! -> su2 Quad
   noel = 4 
 elseif(NTYPE.eq.10) then ! -> su2 Tet
   noel = 4
 elseif(NTYPE.eq.12) then ! -> su2 Hex
   noel = 8
 elseif(NTYPE.eq.13) then ! -> su2 Prism
   noel = 6
 elseif(NTYPE.eq.14) then ! -> su2 Pyramid
   noel = 5
 else
   noel = 0
 endif

 end function


 pure integer function nofael(NTYPE)
!
! Number of faces in element of NTYPE Gambit type of element
!
 implicit none
 integer, intent(in) :: NTYPE

 if(NTYPE.eq.10) then ! -> su2 Tet
   nofael = 4
 elseif(NTYPE.eq.12) then ! -> su2 Hex
   nofael = 6
 elseif(NTYPE.eq.13) then ! -> su2 Prism
   nofael = 5
 elseif(NTYPE.eq.14) then ! -> su2 Pyramid
   nofael = 5
 else
   nofael = 0
 endif

 end function


 function face_vertices_NTYPE12(iface,NODE,NDP) result(vertices)
!
! SU2 follows rules from Paraveiw format file, page 9.
! NTYPE12 is 8-node HEX
!
! NOTE: SU2 starts counting from 0, we here from 1!
!
 implicit none
 integer, intent(in) :: iface
 integer, intent(in) :: NDP
 integer, dimension(NDP), intent(in) :: NODE
!
! Result
!
 integer, dimension(4) :: vertices

 vertices = (/ 0,0,0,0 /)

 select case (iface)
    case (1)
       vertices(:) = (/ NODE(1), NODE(2), NODE(3), NODE(4) /)
    case (2)
       vertices(:) = (/ NODE(6), NODE(7), NODE(3), NODE(2) /)
    case (3)
       vertices(:) = (/ NODE(6), NODE(5), NODE(8), NODE(7) /)
    case (4)
       vertices(:) = (/ NODE(8), NODE(5), NODE(1), NODE(4) /)
    case (5)
       vertices(:) = (/ NODE(1), NODE(5), NODE(6), NODE(2) /)
    case (6)
       vertices(:) = (/ NODE(3), NODE(7), NODE(8), NODE(4) /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function


 function face_vertices_NTYPE14(iface,NODE,NDP) result(vertices)
!
! NTYPE14 is 5-node PYRAMID
!
! So the SU2 convention is:
! 	Face 	Nodes
! 	1 	0,1,2,3
! 	2 	0,4,1
! 	3 	1,4,2
! 	4 	2,4,3
! 	5 	3,4,0
! NOTE: SuU2 starts counting from 0, we here from 1!
!
 implicit none
 integer, intent(in) :: iface
 integer, intent(in) :: NDP
 integer, dimension(NDP), intent(in) :: NODE
!
! Result
!
 integer, dimension(4) :: vertices

 vertices = (/ 0,0,0,0 /)

 select case (iface)
    case (1)
       vertices(:) = (/ NODE(1), NODE(2), NODE(3), NODE(4) /)
    case (2)
       vertices(:) = (/ NODE(1), NODE(5), NODE(2), 0/)
    case (3)
       vertices(:) = (/ NODE(2), NODE(5), NODE(3), 0 /)
    case (4)
       vertices(:) = (/ NODE(3), NODE(5), NODE(4), 0 /)
    case (5)
       vertices(:) = (/ NODE(4), NODE(5), NODE(1), 0 /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function


 function face_vertices_NTYPE10(iface,NODE,NDP) result(vertices)
!
! NTYPE10 is 4-node TETRAHEDRON
!
! So the SU2 convention is:
! 	Face 	Nodes
! 	1 	0,1,2
! 	2 	0,3,1
! 	3 	1,3,2
! 	4 	0,2,3
! NOTE: SU2 starts counting from 0, we here from 1!
!
 implicit none
 integer, intent(in) :: iface
 integer, intent(in) :: NDP
 integer, dimension(NDP), intent(in) :: NODE
!
! Result
!
 integer, dimension(4) :: vertices

 vertices = (/ 0,0,0,0 /)

 select case (iface)
    case (1)
       vertices(:) = (/ NODE(1), NODE(2), NODE(3), 0 /)
    case (2)
       vertices(:) = (/ NODE(1), NODE(4), NODE(2), 0/)
    case (3)
       vertices(:) = (/ NODE(2), NODE(4), NODE(3), 0 /)
    case (4)
       vertices(:) = (/ NODE(1), NODE(3), NODE(4), 0 /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function

 function face_vertices_NTYPE13(iface,NODE,NDP) result(vertices)
!
! NTYPE13 is 6-node PRISM
!
! So the SU2 convention is:
! 	Face 	Nodes
! 	1 	0,1,4,3
! 	2 	1,2,5,4
! 	3 	2,0,3,5
! 	4 	0,2,1
! 	5 	3,4,5
! NOTE: SU2 starts counting from 0, we here from 1!
!
 implicit none
 integer, intent(in) :: iface
 integer, intent(in) :: NDP
 integer, dimension(NDP), intent(in) :: NODE
!
! Result
!
 integer, dimension(4) :: vertices

 vertices = (/ 0,0,0,0 /)

 select case (iface)
    case (1)
       vertices(:) = (/ NODE(1), NODE(2), NODE(5), NODE(4) /)
    case (2)
       vertices(:) = (/ NODE(2), NODE(3), NODE(6), NODE(5) /)
    case (3)
       vertices(:) = (/ NODE(3), NODE(1), NODE(4), NODE(6) /)
    case (4)
       vertices(:) = (/ NODE(1), NODE(3), NODE(2), 0 /)
    case (5)
       vertices(:) = (/ NODE(4), NODE(5), NODE(6), 0 /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function


 subroutine sortIntArray(IntArray,length) !result(newIntArray)
! Ascending sort of an integer array of size[length]
 implicit none
 integer, intent(in) :: length
 integer , dimension(length), intent(inout) :: IntArray
!
! Result
!
 !integer , dimension(length) :: newIntArray
!
! Locals
!
 integer :: icount, indx

! newIntArray = 0.

 do icount = 1, length
     indx = minloc( IntArray(icount:length), dim = 1 )
     indx = indx + icount - 1 
     call swapInt( IntArray(icount), IntArray( indx ) )
 enddo

 end subroutine sortIntArray


 subroutine swapInt(a, b)
 integer, intent(inout) :: a,b
 integer :: c

   c = a
   a = b
   b = c
 end subroutine


subroutine hashmap_insert(iel,fVert,owner,neighbour,fV,listLength,nFaces,ne)

  implicit none

  integer, intent(in) :: iel
  integer , intent(in):: listLength
  integer, dimension(4), intent(in) :: fVert
  integer, dimension(listLength), intent(inout) :: owner, neighbour
  integer, dimension(4,listLength), intent(inout) :: fV
  integer , intent(inout) :: nFaces
  integer, intent(out) :: ne

! integer, parameter :: MyLongIntType = selected_int_kind(18) !9223372036854775807
! integer (MyLongIntType) :: LongInt
  integer :: LongInt
  integer :: lp
  integer :: hash
  integer, dimension(4) :: fV1, fV2


  ! ## Hashing ## 
  ! We get a large number by summing the face indices and we spread it a little bit
  ! more by adding the max and min vertex index, so we decrease multiple collisions,
  ! Many combinations of indices may otherwise give same sum.
  LongInt = sum( fVert(:) )+maxval(fVert(:))+minval(fVert(:))
  ! Printing hash function for debugging purposes only.
  ! print*,mod(listLength+LongInt-(minInt-1),listLength)

  ! Linear probing to solve multiple collisions
  ! Note: We can also use quadratic probing if hash below is: hash = mod(LongInt+lp*lp,listLength) 
  linear_probing_loop: do lp=0,listLength

    ! Before creating a hash we add listLength and subtract (minInt-1) to
    ! get well behaving number, ie. biger than listLength, bigger by exaclty 1 at least in once case.
    ! Hash function:
    hash = mod(LongInt+lp*lp,listLength) 
    ! print*, hash,fVert ! To check which combination of indices gave which hash.

    ! First check is owner at that location free, if it is free, write it there
    if (owner(hash) == 0) then
      
      owner(hash) = iel
      
      fV(:, hash ) = fVert
      
      exit linear_probing_loop

    ! Otherwise check if the neighbour(hash) is free - if yes perform checking if these two are the same face.
    elseif(neighbour(hash) == 0) then

      ! Perform check - have we found a matching pair.
      ! A check is performed first by sorting array of face idices for both faces,
      ! and if three vertices coincide than it is enough to claim the faces are matching.
      ! We have found two neighbouring cells
      fV1 = fVert
      call sortIntArray( fV1, 4 ) 

      fV2 = fV(:, hash)
      call sortIntArray( fV2, 4 ) 

      if ( fV2(2)==fV1(2) .and. fV2(3)==fV1(3) .and. fV2(4)==fV1(4) ) then !( three is enough for matching)
        
        !!* We have a matching pair *!!
        neighbour(hash) = iel

        ne = owner(hash) ! a hack to return owner when going trough boundary faces, does nothing in element loop.
        
        nFaces = nFaces+1
        
        exit linear_probing_loop

      else ! Not matched, lets look for the alternative

        cycle linear_probing_loop

      endif

    ! If cell is full and they don't match with cell vertices we have to do linear probing and
    ! find next free location for this face!
    else ! if( (neighbour(hash) /= 0).and.owner(hash)/=0 ) then

      cycle linear_probing_loop

    endif

  enddo linear_probing_loop

end subroutine


subroutine create_field_init_files(bcName,nbc)
!
! Create two template files in 0/ folder to help setting 
! initial conditions for scalar and vector fields.
!
 implicit none

 integer, intent(in) :: nbc ! no of boundaries
 character( len = 30 ), dimension(nbc) :: bcName

 integer :: i
 integer :: os

 ! Write vector field template
 call get_unit( os )

 open(unit=os, file='0/U')
 rewind os

 write(os,'(a)') 'internalField'
 write(os,'(2x,a)') 'uniform'
 write(os,'(4x,a)') '0.0 0.0 0.0'
 write(os,'(a)') 'boundaryField'
 do i=1,nbc
   write(os,'(a)') trim( bcName(i) )
   write(os,'(2x,a)') 'Dirichlet/Neumann'
   write(os,'(4x,a)') 'uniform/zeroGradient'
   write(os,'(6x,a)') '0.0 0.0 0.0'
 enddo

 close(os)

 ! Write scalar field template
 call get_unit( os )

 open(unit=os, file='0/T.template')
 rewind os

 write(os,'(a)') 'internalField'
 write(os,'(2x,a)') 'uniform'
 write(os,'(4x,a)') '0.0'
 write(os,'(a)') 'boundaryField'
 do i=1,nbc
   write(os,'(a)') trim( bcName(i) )
   write(os,'(2x,a)') 'Dirichlet/Neumann'
   write(os,'(4x,a)') 'uniform/zeroGradient'
   write(os,'(6x,a)') '0.0'
 enddo

 close(os)

end subroutine


subroutine create_paraview_files(numPoints,numCells,numInnerFaces,bcSize,nbc)
!
! Write Paraview files for all regions which will serve 
! in the freeCappuccino solver later as templates for creating
! output files. Also it may serve to view the mesh we have created.
!

 implicit none

 integer, intent(in) :: numPoints,numCells,numInnerFaces
 integer, intent(in) :: nbc
 integer, dimension(nbc) :: bcSize

 integer :: i,k
 integer :: os,is
 integer :: off, ntype
 integer :: nBoundFaces 
 character( len = 2 ) :: ch2
 ! character( len = 10 ) :: ch
 real(kind = 8), dimension(:,:), allocatable :: xcoo 
 integer, dimension(:), allocatable :: types
 integer, dimension(8) :: node

 !
 ! Create .vtu file for interior of the domain.
 !

 call get_unit( os )

 open(unit=os, file='vtk/mesh/mesh_0_0.vtu')
 rewind os

 write(os,'(a)') '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
 write(os,'(2x,a)') '<UnstructuredGrid>'
 write(os,'(4x,2(a,i0),a)') '<Piece NumberOfPoints="',numPoints,'" NumberOfCells="',numCells,'">'
 write(os,'(6x,a)') '<PointData>'
 write(os,'(6x,a)') '</PointData>'
 write(os,'(6x,a)') '<CellData>'
 write(os,'(6x,a)') '</CellData>'
 write(os,'(6x,a)') '<Points>'
 write(os,'(8x,a)') '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">'

 call get_unit( is )
 open(unit=is, file='polyMesh/points')
 rewind is

 allocate(xcoo(3,numPoints))

 do k=1,numPoints
   read( is,*) xcoo(1:3,k)
   write(os,'(10x,3E20.11)') xcoo(1:3,k)
 enddo

 close(is)

 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Points>'
 write(os,'(6x,a)') '<Cells>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="connectivity" format="ascii">'

 call get_unit( is )
 open(unit=is, file='polyMesh/cells')
 rewind is

 allocate( types(numCells) )

 do k=1,numCells
   read( is,*) ntype,node(1:noel(ntype))
   types(k) = ntype
   write(os,'(10x,8(i0,1x))') node(1:noel(ntype))
 enddo

 close(is)

 write(os,'(8x,a)') '</DataArray>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="offsets" format="ascii">'
 off = 0
 do k=1,numCells
   ntype = types(k)
   off = off + noel(ntype)
   write(os,'(10x,i0)') off
 enddo

 write(os,'(8x,a)') '</DataArray>'
 write(os,'(8x,a)') '<DataArray type="UInt8" Name="types" format="ascii">'
 do k=1,numCells
   write(os,'(10x,i0)') types(k)
 enddo
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Cells>'
 write(os,'(4x,a)') '</Piece>'
 write(os,'(2x,a)') '</UnstructuredGrid>'
 write(os,'(a)')    '</VTKFile>'

 close(os)

 deallocate(types)

 !
 ! Write .vtp files for each boundary region.
 !

 ! allocate array for face types
 nBoundFaces = maxval( bcSize(:) )
 allocate( types( nBoundFaces ) )

 ! Prepare faces file for reading of boundaries - fast-forward it to
 ! the place where boundarry faces start
 call get_unit( is )
 open(unit=is, file='polyMesh/faces')
 rewind is
 do k=1,numInnerFaces
  read(is,*)
 enddo

 do i=1,nbc

 nBoundFaces = bcSize(i)

 call get_unit( os )

 write(ch2,'(i0)') i
 open(unit=os, file='vtk/mesh/mesh_'//trim(ch2)//'_0.vtp')
 rewind os

 write(os,'(a)') '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
 write(os,'(2x,a)') '<PolyData>'
 write(os,'(4x,2(a,i0),a)') '<Piece NumberOfPoints="',numPoints,'" NumberOfVerts="0" NumberOfLines="0"&
 & NumberOfStrips="0" NumberOfPolys="',nBoundFaces,'">'
 write(os,'(6x,a)') '<PointData>'
 write(os,'(6x,a)') '</PointData>'
 write(os,'(6x,a)') '<CellData>'
 write(os,'(6x,a)') '</CellData>'
 write(os,'(6x,a)') '<Points>'
 write(os,'(8x,a)') '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">'
 do k=1,numPoints
   write(os,'(10x,3E20.11)') xcoo(1:3,k)
 enddo
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Points>'
 write(os,'(6x,a)') '<Verts>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Verts>'
 write(os,'(6x,a)') '<Lines>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Lines>'
 write(os,'(6x,a)') '<Strips>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Strips>'
 write(os,'(6x,a)') '<Polys>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="connectivity" format="ascii">'

 do k=1,nBoundFaces
   read(is,*) ntype, node(1:ntype)
   types(k) = ntype
   write(os,'(10x,4(i0,1x))') node(ntype:1:-1)-1
 enddo

 write(os,'(8x,a)') '</DataArray>'
 write(os,'(8x,a)') '<DataArray type="Int64" IdType="1" Name="offsets" format="ascii">'

 off = 0
 do k=1,nBoundFaces
   ntype = types(k)
   off = off + ntype
   write(os,'(10x,i0)') off
 enddo

 write(os,'(8x,a)') '</DataArray>'
 write(os,'(6x,a)') '</Polys>'
 write(os,'(4x,a)') '</Piece>'
 write(os,'(2x,a)') '</PolyData>'
 write(os,'(a)')    '</VTKFile>'

 close(os)
 
 enddo

 close(is)

 deallocate( xcoo,types )

end subroutine

!
! Useful routines by John Burkardt
!


subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) temp

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      temp  = a1(i)
      a1(i) = a1(j)
      a1(j) = temp

      temp  = a2(i)
      a2(i) = a2(j)
      a2(j) = temp
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end


subroutine i4vec3_sort_a ( n, a1, a2, a3 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers 
!! plus a permutation vector A3.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!    Nikola Mirkov (modified to i4vec3_sort_a fro i4vec2_sort_a)
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) a3(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) temp

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      temp  = a1(i)
      a1(i) = a1(j)
      a1(j) = temp

      temp  = a2(i)
      a2(i) = a2(j)
      a2(j) = temp

      temp  = a3(i)
      a3(i) = a3(j)
      a3(j) = temp      
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end


subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end

subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end

subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end

subroutine i4vec_print2 ( n, a, b, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), first vector to be printed.
!
!    Input, integer ( kind = 4 ) B(N), second vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n),b(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i4,2x,i4)' ) i, ':', a(i), b(i)
  end do

  return
end


subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end


end module utils
