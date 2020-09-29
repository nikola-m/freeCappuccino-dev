module output
use types
use parameters
use geometry
use utils, only: get_unit, i4_to_s_left

public 

contains

subroutine vtu_write_XML_header ( output_unit )
!
! Writes header data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit

  integer :: icell

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )


  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'

! > Open XML tag for recording data in cell centers
  write ( output_unit, '(6x,a)' ) '<CellData>'

  call i4_to_s_left ( numCells-1, cells_num_string )

  write ( output_unit, '(8x,a)' ) &
  & '<DataArray type="Int32" Name="cellID" format="ascii" RangeMin="0" RangeMax="'//trim( cells_num_string )//'">'

    do icell=1,numCells
      write( output_unit, '(10x,i8)') icell
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

end subroutine


subroutine vtu_write_XML_meshdata ( output_unit )
!
! Writes unstructured mesh data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit

  integer :: i,k
  integer :: icell
  integer :: ntype
  integer :: offset
  integer :: cells_file
  integer, dimension(8) :: node


  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file


! > End of recording data in cell centers
  write ( output_unit, '(6x,a)' ) '</CellData>'

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3e15.7)' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype,(node(k), k=1,noel(ntype))
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file

  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine



subroutine vtu_write_XML_scalar_field ( output_unit, scalar_name, scalar_field )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numCells), intent(in) :: scalar_field

  integer :: icell

! <Scalars in cell-centers>
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers and boundary faces > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,e15.7)') scalar_field(icell) 
    enddo


  write ( output_unit, '(8x,a)' ) '</DataArray>'
! </Scalars in cell-centers>

end subroutine vtu_write_XML_scalar_field


subroutine vtu_write_XML_scalar_field_boundary ( output_unit, scalar_name, scalar_field, istart, iend )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit, istart, iend
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numTotal), intent(in) :: scalar_field

  integer :: icell

! <Scalars in cell-centers>
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers and boundary faces > write scalar data
!
    do icell=istart,iend
      write( output_unit, '(10x,e15.7)') scalar_field(icell) 
    enddo


  write ( output_unit, '(8x,a)' ) '</DataArray>'
! </Scalars in cell-centers>

end subroutine vtu_write_XML_scalar_field_boundary


subroutine vtu_write_XML_vector_field ( output_unit, field_name, u, v, w )
!
! Writes vector field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: field_name
  real(dp), dimension(numCells), intent(in) :: u, v, w

  integer :: icell

! <Vectors in cell-centers>
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',trim( field_name ),'" NumberOfComponents="3" Format="ascii">'

!
! > Scalars in cell-centers and boundary faces > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,3(1x,e15.7))') u(icell), v(icell), w(icell)
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
! </Vectors in cell-centers>

end subroutine vtu_write_XML_vector_field


subroutine vtu_write_XML_vector_field_boundary ( output_unit, field_name, u, v, w, istart, iend )
!
! Writes vector field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit, istart, iend
  character ( len = * ), intent(in) :: field_name
  real(dp), dimension(numTotal), intent(in) :: u, v, w

  integer :: icell

! <Vectors in cell-centers>
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',trim( field_name ),'" NumberOfComponents="3" Format="ascii">'

!
! > Scalars in cell-centers and boundary faces > write scalar data
!
    do icell=istart,iend
      write( output_unit, '(10x,3(1x,e15.7))') u(icell), v(icell), w(icell)
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
! </Vectors in cell-centers>

end subroutine vtu_write_XML_vector_field_boundary

subroutine vtm_write_scalar_field ( scalar_name, scalar_field, timechar )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  character( len = * ), intent(in) :: scalar_name, timechar
  real(dp), dimension(numTotal), intent(in) :: scalar_field

  character( len = 20 ) :: node_num_string
  character( len = 20 ) :: cells_num_string
  character( len = 2 ) :: ch2
  character( len = 1 ) :: ch

  integer :: i,k,ib,ijb
  integer :: icell, iface
  integer, dimension(numCells) :: ntype(numCells)
  integer :: offset
  integer :: output_unit, cells_file, faces_file
  integer :: nnodes
  integer, dimension(8) :: node

  ! Open folder with data for postprocessing in Paraview
  ! call execute_command_line("mkdir VTK")
  call execute_command_line("mkdir VTK/"//trim( timechar ) )
  call execute_command_line("mkdir VTK/"//trim( timechar )//"/boundary")  

!
! > Open and write a .vtm file for multi-block datasets for interior + boundary regions data.
!
  call get_unit( output_unit )

  open(unit=output_unit,file='VTK/scalar_name_'//trim( timechar )//'.vtm')

  write ( output_unit, '(a)' )    '<?xml version="1.0"?>'
  write ( output_unit, '(2x,a)' ) '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
  write ( output_unit, '(4x,a)' ) '<vtkMultiBlockDataSet>'
  write ( output_unit, '(6x,a)' ) '<Block index="0" name="Cells">'
  write ( output_unit, '(8x,a)' ) '<DataSet index="0" name="interior" file="'//trim(timechar)//'/interior.vtu" format="appended">'
  write ( output_unit, '(8x,a)' ) '</DataSet>'
  write ( output_unit, '(6x,a)' ) '</Block>'
  write ( output_unit, '(6x,a)' ) '<Block index="1" name="Boundaries">'

  do i=1,numBoundaries

    call i4_to_s_left ( i-1, ch2 )

    write ( output_unit, '(8x,a)' ) '<DataSet index="'//trim(ch2)//'" name="'//trim( bcname(i) )//&
    &'" file="'//trim(timechar)//'/boundary/'//trim( bcname(i) )//'.vtu" format="appended">'
    
    write ( output_unit, '(8x,a)' ) '</DataSet>'

  enddo ! Boundary loop

  write ( output_unit, '(6x,a)' ) '</Block>'
  write ( output_unit, '(4x,a)' ) '</vtkMultiBlockDataSet>'
  write ( output_unit, '(2x,a)' ) '</VTKFile>'

  close(output_unit)

!
! > Open and write a unstructuctured .vtu file for the interior cells
!
  call get_unit( output_unit )

  open(unit=output_unit,file='VTK/'//trim(timechar)//'/interior.vtu')

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  read(cells_file, *) ch ! A line where numCells is.

  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'

!
! > Scalars in cell-centers > open xml tags
!
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,e15.7)') scalar_field(icell) 
    enddo
!
! > Scalars in cell-centers > close xml tags
!
  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</CellData>'

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3(e15.7,1x))' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype(icell),(node(k), k=1,noel( ntype(icell) ))
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel( ntype(icell) )) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells
      offset = offset+noel( ntype(icell) )
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells
      write( output_unit, '(10x,i2)') ntype(icell) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' )    '</VTKFile>'

  close( output_unit )


!
! > Open and write a .vtm file for multi-block datasets ONLY for boundary regions data.
!
  call get_unit( output_unit )

  open(unit=output_unit,file='VTK/boundary_'//trim(timechar)//'.vtm')

  write ( output_unit, '(a)' )    '<?xml version="1.0"?>'
  write ( output_unit, '(2x,a)' ) '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
  write ( output_unit, '(4x,a)' ) '<vtkMultiBlockDataSet>'
  write ( output_unit, '(6x,a)' ) '<Block index="0" name="Boundaries">'

  do i=1,numBoundaries

    call i4_to_s_left ( i-1, ch2 )

    write ( output_unit, '(8x,a)' )'<DataSet index="'//trim(ch2)//'" name="'//trim( bcname(i) )//&
    &'" file="'//trim(timechar)//'/boundary/'//trim( bcname(i) )//'.vtu" format="appended">'
    
    write ( output_unit, '(8x,a)' ) '</DataSet>'

  enddo ! Boundary loop

  write ( output_unit, '(6x,a)' ) '</Block>'
  write ( output_unit, '(4x,a)' ) '</vtkMultiBlockDataSet>'
  write ( output_unit, '(2x,a)' ) '</VTKFile>'

  close(output_unit)


! 
! > Open and rewind faces file to the point where THIS (ib) BOUNDARY faces start
!
  call get_unit( faces_file )
  open( unit = faces_file, file='polyMesh/faces' )
  rewind faces_file

  ! The thing below is ONLY FOR OpenFOAM polyMesh; change it later for native format!
  rewind( faces_file )
  ch = ' '
  do
    read(faces_file,*) ch ! Read anything... 
    if (ch == "(") then
      exit
    endif
  end do 

  do iface=1,startFace(1)  ! <--- go to place where boundary faces start
    read(faces_file,*) ch  ! Read anything... 
  enddo

!
! > Open and write a unstructuctured .vtu file for EACH boundary region.
!  
  do ib=1,numBoundaries

    call get_unit( output_unit )

    open(unit=output_unit,file='VTK/'//trim(timechar)//'/boundary/'//trim( bcname(ib) )//'.vtu')

    !
    ! > Header
    !
    call i4_to_s_left ( numNodes, node_num_string )
    call i4_to_s_left ( nFaces(ib), cells_num_string )


    write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
    write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
    write ( output_unit, '(4x,5a)' ) &
    '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'

    !
    ! > Scalars in face-centers > open xml tags
    !
    write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
    write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

    !
    ! > Scalars in face-centers > write scalar data
    !
      do i=1,nFaces(ib)

        ijb = iBndValueStart(ib) + i
        write( output_unit, '(10x,e15.7)') scalar_field(ijb) 

      enddo

    !
    ! > Scalars in face-centers > close xml tags
    !
    write ( output_unit, '(8x,a)' ) '</DataArray>'
    write ( output_unit, '(6x,a)' ) '</CellData>'

    !
    ! > Mesh data
    !
    write ( output_unit, '(6x,a)' ) '<Points>'
    write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

    !
    ! > Mesh data > Nodal coordinates
    !
    do i=1,numNodes
      write ( output_unit, '(10x,3(e15.7,1x))' ) x(i),y(i),z(i)
    enddo

    write ( output_unit, '(8x,a)' ) '</DataArray>'
    write ( output_unit, '(6x,a)' ) '</Points>'

    !
    ! > Mesh data > Connectivity
    !
    write ( output_unit, '(6x,a)' ) '<Cells>'

    write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'

      do i=1,nFaces(ib)

        call read_line_faces_file_polyMesh(faces_file,nnodes,node,4) 

        !read( faces_file, * ) nnodes,(node(k), k=1,nnodes )

        write( output_unit, '(10x,4i8)') (node(k)-1, k=1,nnodes )  
       
        ! Accumulate for offsets
        ntype( i ) = nnodes

      enddo

    write ( output_unit, '(8x,a)' ) '</DataArray>'

    !
    ! > Mesh data > Offsets
    !
    offset = 0

    write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

      do i=1,nFaces(ib)
        offset = offset + ntype(i) 
        write( output_unit, '(10x,i12)') offset 
      enddo

    write ( output_unit, '(8x,a)' ) '</DataArray>'

    !
    ! > Mesh data > Types
    !
    write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

      do i=1,nFaces(ib)

        ! A hack for how to write element type for each face
        if ( ntype( i ) == 3 ) then   ! i.e. it is a triangular face
          write(ch2,'(a)') '5 '       ! '5' - paraview identificator for triagnular element       
        else                          ! ntype(icell) == 4 i.e. it is a quad face
          write(ch2,'(a)') '9 '       ! '9' - paraview identificator for quad element
        endif

        write( output_unit, '(10x,a)') trim( ch2 )
      enddo

    write ( output_unit, '(8x,a)' ) '</DataArray>'

    write ( output_unit, '(6x,a)' ) '</Cells>'
    write ( output_unit, '(4x,a)' ) '</Piece>'
    write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
    write ( output_unit, '(a)' )    '</VTKFile>'

    close( output_unit )

  enddo ! Boundary loop

  close( faces_file )

end subroutine


subroutine vtu_write_scalar_field ( output_unit, scalar_name, scalar_field )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numCells), intent(in) :: scalar_field

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string

  integer :: i,k
  integer :: icell
  integer, dimension(numCells) :: ntype(numCells)
  integer :: offset
  integer :: cells_file
  integer, dimension(8) :: node

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file


  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'


! <Scalars in cell-centers>
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,e15.7)') scalar_field(icell) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</CellData>'
! </Scalars in cell-centers>

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3(1x,e15.7))' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype(icell),(node(k), k=1,noel( ntype(icell) ))
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel( ntype(icell) )) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  !rewind cells_file

  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells
      !read( cells_file, * ) ntype
      offset = offset+noel( ntype(icell) )
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  !rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells
      !read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype(icell) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_scalar_field




subroutine vtu_write_vector_field ( output_unit, field_name, u, v, w )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!

  implicit none
  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: field_name
  real(dp), dimension(numCells), intent(in) :: u, v, w

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string

  integer :: i,k
  integer :: icell
  integer :: ntype
  integer :: offset
  integer :: cells_file
  integer, dimension(8) :: node

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'


! <Scalars in cell-centers>
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',trim( field_name ),'" NumberOfComponents="3" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,3(1x,e15.7))') u(icell), v(icell), w(icell)
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</CellData>'
! </Scalars in cell-centers>

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3(1x,e15.7))' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype,(node(k), k=1,noel(ntype))
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_vector_field



subroutine vtu_write_mesh ( output_unit )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string

  integer :: i,k
  integer :: icell
  integer :: ntype
  integer :: offset
  integer :: cells_file

  integer, dimension(8) :: node

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3(1x,e15.7))' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1,1x,4i8:(4i8:))' ) ntype,(node(k), k=1,noel(ntype))
      write( output_unit, '(4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do icell=1,numCells
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'
    do icell=1,numCells
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_mesh


 pure integer function paraview_ntype(NTYPE)
!
! Element type in Paraview corresponding to Cappuccino type of element given by NTYPE.
!
 implicit none
 integer, intent(in) :: NTYPE

 if(NTYPE.eq.1) then ! -> Cappuccino line
   paraview_ntype = 3 
 elseif(NTYPE.eq.2) then ! -> Cappuccino Tri
   paraview_ntype = 5 
 elseif(NTYPE.eq.3) then ! -> Cappuccino Quad
   paraview_ntype = 9 
 elseif(NTYPE.eq.4) then ! -> Cappuccino Tet
   paraview_ntype = 10   
 elseif(NTYPE.eq.5) then ! -> Cappuccino Hex
   paraview_ntype = 12 
 elseif(NTYPE.eq.6) then ! -> Cappuccino Prism
   paraview_ntype = 13 
 elseif(NTYPE.eq.7) then ! -> Cappuccino Pyramid
   paraview_ntype = 14 
 else
   paraview_ntype = 0
 endif

 end function

 pure integer function noel(NTYPE)
!
! Number of nodes in element of NTYPE type of element
!
 implicit none
 integer, intent(in) :: NTYPE


 if(NTYPE.eq.12) then ! -> Paraview Hex
   noel = 8  
 elseif(NTYPE.eq.13) then ! -> Paraview Prism
   noel = 6 
 elseif(NTYPE.eq.14) then ! -> Paraview Pyramid
   noel = 5
 elseif(NTYPE.eq.10) then ! -> Paraview Tet
   noel = 4
 elseif(NTYPE.eq.9) then ! -> Paraview Quad
   noel = 4
 elseif(NTYPE.eq.5) then ! -> Paraview Tri
   noel = 3
 else
   noel = 0
 endif

 end function

end module