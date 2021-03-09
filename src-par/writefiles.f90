!***********************************************************************
!
subroutine writefiles
!
! Write output files in Paraview .vtm format
!
!***********************************************************************
!
  use types
  use parameters
  use title_mod
  use geometry
  use variables
  use statistics
  use sparse_matrix
  use output
  use utils, only: i4_to_s_left
  use vortexIdentification
  ! use mhd

  implicit none
!
!***********************************************************************
!
  character( len = 200 ) :: line
  character( len = 4 ) :: nproc_char
  character( len = 2 ) :: ch2

  integer :: i,ib
  integer :: istart, iend
  integer :: stat
  integer :: iWall
  integer :: output_unit, mesh_file

  ! Write in a char variable current timestep number and create a folder with this name
  call i4_to_s_left ( itime, timechar )  
  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

  !+-----------------------------------------------------------------------------+  

  ! Open folder with data for postprocessing in Paraview
  call execute_command_line('mkdir '//trim( timechar ) )
  call execute_command_line('mkdir '//trim( timechar )//'/processor'//trim(nproc_char) )
  call execute_command_line('mkdir '//trim( timechar )//'/processor'//trim(nproc_char)//'/boundary')  

  
!
! > Open and write a .vtm file for multi-block datasets for interior + boundary regions data.
!
  call get_unit( output_unit )

  open(unit=output_unit,file=trim( timechar )//'/solution_fields_proc'//trim( nproc_char )//'.vtm')

  write ( output_unit, '(a)' )    '<?xml version="1.0"?>'
  write ( output_unit, '(2x,a)' ) '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
  write ( output_unit, '(4x,a)' ) '<vtkMultiBlockDataSet>'
  write ( output_unit, '(6x,a)' ) '<Block index="0" name="Cells">'
  write ( output_unit, '(8x,a)' ) '<DataSet index="0" name="interior" file="processor'&
                                  &//trim(nproc_char)//'/interior.vtu" format="appended">'
  write ( output_unit, '(8x,a)' ) '</DataSet>'
  write ( output_unit, '(6x,a)' ) '</Block>'
  write ( output_unit, '(6x,a)' ) '<Block index="1" name="Boundaries">'

  do i=1,numBoundaries

    call i4_to_s_left ( i-1, ch2 )

    write ( output_unit, '(8x,a)' ) '<DataSet index="'//trim(ch2)//'" name="'//trim( bcname(i) )//&
    &'" file="processor'//trim(nproc_char)//'/boundary/'//trim( bcname(i) )//'.vtp" format="appended">'
    
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

  open(unit=output_unit,file=trim(timechar)//'/processor'//trim(nproc_char)//'/interior.vtu')

!
! > Header
!
  call get_unit( mesh_file )
  open( unit = mesh_file, file='processor'//trim(nproc_char)//'/vtk/mesh/mesh_0_0.vtu' )
  rewind mesh_file

  do i=1,6
    read( mesh_file, '(a)' ) line! A line where numCells is.
    write ( output_unit, '(a)' ) adjustl( line )
  enddo  

!
! > Vectors in cell-centers 
!

  call vtu_write_XML_vector_field( output_unit, 'U', u, v, w )

  ! if( lcal(iep) ) call vtu_write_XML_vector_field( output_unit, 'uxB', curix, curiy, curiz )


!
! > Scalars in cell-centers 
!

  call vtu_write_XML_scalar_field ( output_unit, 'p', p )

  if( lturb ) call vtu_write_XML_scalar_field ( output_unit, 'mueff', vis )

  if(solveTKE) call vtu_write_XML_scalar_field ( output_unit, 'k', te )

  if( solveEpsilon ) call vtu_write_XML_scalar_field ( output_unit, 'epsilon', ed )

  if( solveOmega ) call vtu_write_XML_scalar_field ( output_unit, 'omega', ed )

  if( lcal(ien) ) call vtu_write_XML_scalar_field ( output_unit, 'T', t )

  ! if( lcal(iep) ) call vtu_write_XML_scalar_field ( output_unit, 'Epot', Epot )

  !
  ! > Time averaged fields in the interior
  !
  if(ltransient) then

    call vtu_write_XML_vector_field( output_unit, 'Uavg', u_aver,v_aver,w_aver )

    if(lturb) then
      call vtu_write_XML_scalar_field ( output_unit, '05uiui', te_aver )
      call vtu_write_XML_scalar_field ( output_unit, 'uu', uu_aver )
      call vtu_write_XML_scalar_field ( output_unit, 'vv', vv_aver )
      call vtu_write_XML_scalar_field ( output_unit, 'ww', ww_aver )
      call vtu_write_XML_scalar_field ( output_unit, 'uv', uv_aver )
      call vtu_write_XML_scalar_field ( output_unit, 'uw', uw_aver )
      call vtu_write_XML_scalar_field ( output_unit, 'vw', vw_aver )                            
    endif

    if(lcal(ien) ) call vtu_write_XML_scalar_field ( output_unit, 'Tavg', t_aver )

  endif


  !
  ! > Identification of vortices
  !
  ! if( lturb ) then
  !   call setQvortex
  !   call vtu_write_XML_scalar_field ( output_unit, 'Qcrit', Qvortex )
  ! endif


!
! > Mesh data
!
  mesh_data_loop: do

    read(mesh_file, '(a)', IOSTAT=stat) line

    if ( stat < 0 ) then
      exit mesh_data_loop
    else
      write ( output_unit, '(a)' ) adjustl( line )
      cycle mesh_data_loop
    endif
 
  enddo mesh_data_loop

  close( mesh_file )
  close( output_unit )


!
! > Open and write a .vtm file for multi-block datasets ONLY for boundary regions data.
!
  call get_unit( output_unit )

  open(unit=output_unit,file=trim(timechar)//'/boundary_proc'//trim(nproc_char)//'.vtm')

  write ( output_unit, '(a)' )    '<?xml version="1.0"?>'
  write ( output_unit, '(2x,a)' ) '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
  write ( output_unit, '(4x,a)' ) '<vtkMultiBlockDataSet>'
  write ( output_unit, '(6x,a)' ) '<Block index="0" name="Boundaries">'

  do i=1,numBoundaries

    call i4_to_s_left ( i-1, ch2 )

    write ( output_unit, '(8x,a)' )'<DataSet index="'//trim(ch2)//'" name="'//trim( bcname(i) )//&
    &'" file="processor'//trim(nproc_char)//'/boundary/'//trim( bcname(i) )//'.vtp" format="appended">'
    
    write ( output_unit, '(8x,a)' ) '</DataSet>'

  enddo ! Boundary loop

  write ( output_unit, '(6x,a)' ) '</Block>'
  write ( output_unit, '(4x,a)' ) '</vtkMultiBlockDataSet>'
  write ( output_unit, '(2x,a)' ) '</VTKFile>'

  close(output_unit)

!
! > Open and write a unstructuctured .vtp file for EACH boundary region.
!  
  iWall = 1

  do ib=1,numBoundaries

    call get_unit( output_unit )

    open(unit=output_unit,file=trim(timechar)//'/processor'//trim(nproc_char)//'/boundary/'//trim( bcname(ib) )//'.vtp')

    !
    ! > Header
    !

    ! Serial number of the specific boundary region - with it's own file
    call i4_to_s_left ( ib, ch2 )

    call get_unit( mesh_file )
    open( unit = mesh_file, file='processor'//trim(nproc_char)//'/vtk/mesh/mesh_'//trim(ch2)//'_0.vtp' )
    rewind mesh_file

    do i=1,6
      read(mesh_file,  '(a)' ) line
      write ( output_unit, '(a)' ) adjustl( line )
    enddo  


    istart = iBndValueStart(ib) + 1
    iend = iBndValueStart(ib) + nFaces(ib)

    !
    ! > Vectors in face-centers 
    !

    call vtu_write_XML_vector_field_boundary( output_unit, 'U', u, v, w, istart, iend )

    ! if( lcal(iep) ) call vtu_write_XML_vector_field_boundary( output_unit, 'uxB', curix, curiy, curiz, istart, iend )

    !
    ! > Scalars in face-centers 
    !

    call vtu_write_XML_scalar_field_boundary ( output_unit, 'p', p, istart, iend )

    if( lturb ) call vtu_write_XML_scalar_field_boundary ( output_unit, 'mueff', vis, istart, iend )

    if(solveTKE) call vtu_write_XML_scalar_field_boundary ( output_unit, 'k', te, istart, iend )

    if( solveEpsilon ) call vtu_write_XML_scalar_field_boundary ( output_unit, 'epsilon', ed, istart, iend )

    if( solveOmega ) call vtu_write_XML_scalar_field_boundary ( output_unit, 'omega', ed, istart, iend )

    if( lcal(ien) ) call vtu_write_XML_scalar_field_boundary ( output_unit, 'T', t, istart, iend )

    ! if( lcal(iep) ) call vtu_write_XML_scalar_field_boundary ( output_unit, 'Epot', Epot, istart, iend )


    !
    ! > Time averaged fields in the boundary
    !
    if(ltransient) then

      call vtu_write_XML_vector_field_boundary ( output_unit, 'Uavg', u_aver,v_aver,w_aver, istart, iend )

      if(lturb) then
        call vtu_write_XML_scalar_field_boundary ( output_unit, '05uiui', te_aver, istart, iend )
        call vtu_write_XML_scalar_field_boundary  ( output_unit, 'uu', uu_aver, istart, iend )
        call vtu_write_XML_scalar_field_boundary  ( output_unit, 'vv', vv_aver, istart, iend )
        call vtu_write_XML_scalar_field_boundary  ( output_unit, 'ww', ww_aver, istart, iend )
        call vtu_write_XML_scalar_field_boundary  ( output_unit, 'uv', uv_aver, istart, iend )
        call vtu_write_XML_scalar_field_boundary  ( output_unit, 'uw', uw_aver, istart, iend )
        call vtu_write_XML_scalar_field_boundary  ( output_unit, 'vw', vw_aver, istart, iend )  
      endif      

      if(lcal(ien) ) call vtu_write_XML_scalar_field_boundary ( output_unit, 'Tavg', t_aver, istart, iend )

    endif


    ! Write y+ and shear force at wall boundary regions
    if ( bctype(ib) == 'wall' ) then

      istart = iWall
      iend = iWall + nFaces(ib)

      call vtu_write_XML_scalar_field_boundary ( output_unit, 'ypl', ypl, iWall, iend )

      call vtu_write_XML_scalar_field_boundary ( output_unit, 'tau', tau, iWall, iend )

      if(ltransient) call vtu_write_XML_scalar_field_boundary ( output_unit, 'tauAvg', wss_aver, iWall, iend )

      iWall = iWall + nFaces(ib)

    endif


    !
    ! > Mesh data
    !
    mesh_data_loop2: do

      read(mesh_file, '(a)', IOSTAT=stat) line

      if ( stat < 0 ) then
        exit mesh_data_loop2
      else
        write ( output_unit, '(a)' ) adjustl( line )
        cycle mesh_data_loop2
      endif
   
    enddo mesh_data_loop2

    close( mesh_file )
    close( output_unit )

  enddo ! Boundary loop




end subroutine


