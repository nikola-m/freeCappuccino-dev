module boundary_module

 implicit none

  character(len=20) :: key(10), value(10)


 public


 contains
 

 subroutine process_boundary_file(boundary_of,boundary_file)

  implicit none

  integer, intent(in) :: boundary_of
  integer, intent(in) :: boundary_file

  integer :: numBoundaries
  integer :: nFaces, startFace
  character(len=1) :: ch
  character(len=20) :: line_string
  character(len=30) :: char_string,bc_type,bc_name


    do
      read(boundary_of,*) ch
      if (ch == "(") then
        ! Return two lines
        backspace(boundary_of)
        backspace(boundary_of)
        read(boundary_of,*) numBoundaries
        exit
      endif
    end do

    outer_loop: do 

       !**************************************
       !
       read(boundary_of,*) line_string
       !
       !**************************************

        if (line_string(1:1) == "{") then
                backspace(boundary_of)
                backspace(boundary_of)
                read(boundary_of,*) bc_name              ! Read BC name
                read(boundary_of,*) ch                   ! This reads '{' again
                read(boundary_of,*) char_string, bc_type ! Read BC type
                call check_dict(bc_type)                 ! Check in Dictionary if there is a substitute word for this type          
                call check_dict(bc_name)                 ! Check in Dictionary if there is a substitute word for this name

        inner_loop: do
        read(boundary_of,*) char_string, line_string


            if (char_string == "nFaces") then
                read(line_string,*) nFaces

            elseif (char_string == "startFace" ) then 
                read(line_string,*) startFace
                if(nFaces > 0) then
                  write(boundary_file,'(a,1x,a,2(1x,i0))')  trim(bc_name),trim(bc_type),nFaces,startFace
                endif
                exit inner_loop

            endif

        enddo inner_loop

        elseif (line_string(1:1) == "}") then
            cycle outer_loop
        endif

        if (line_string(1:1) == ")") then
            exit outer_loop
        endif

    end do outer_loop


end subroutine



subroutine process_boundary_file_par(boundary_of,boundary_file, process_file)

  implicit none

  integer, intent(in) :: boundary_of
  integer, intent(in) :: boundary_file
  integer, intent(in) :: process_file

  integer :: numBoundaries
  integer :: nFaces, startFace
  character(len=1) :: ch
  character(len=20) :: line_string
  character(len=30) :: char_string,bc_type,bc_name
  logical :: processor_boundary


    do
      read(boundary_of,*) ch
      if (ch == "(") then
        ! Return two lines
        backspace(boundary_of)
        backspace(boundary_of)
        read(boundary_of,*) numBoundaries
        exit
      endif
    end do

    outer_loop: do 

       processor_boundary = .false.

       !**************************************
       !
       read(boundary_of,*) line_string
       !
       !**************************************

        if (line_string(1:1) == "{") then
                backspace(boundary_of)
                backspace(boundary_of)
                read(boundary_of,*) bc_name              ! Read BC name
                read(boundary_of,*) ch                   ! This reads '{' again
                read(boundary_of,*) char_string, bc_type ! Read BC type
                call check_dict(bc_type)                 ! Check in Dictionary if there is a substitute word for this type   
                call check_dict(bc_name)                 ! Check in Dictionary if there is a substitute word for this name

        inner_loop: do
        read(boundary_of,*) char_string, line_string

            if ( line_string(1:9) == 'processor' ) then
                processor_boundary = .true.
            else
                processor_boundary = .false.
            endif

            if (char_string == "nFaces") then
                read(line_string,*) nFaces             

            elseif (char_string == "startFace" .and. .not.processor_boundary ) then 
                read(line_string,*) startFace
                if(nFaces > 0) then
                  write(boundary_file,'(2a,2(1x,i0))')  bc_name,bc_type,nFaces,startFace
                endif
                exit inner_loop

            elseif (char_string == "startFace" .and. processor_boundary) then 
                read(line_string,*) startFace

            elseif (char_string == "neighbProcNo") then 
                ! read(line_string,*) neighbProcNo
                write(process_file,*) line_string !neighbProcNo
                write(boundary_file,'(2a,2(1x,i0))')  trim(bc_name),trim(bc_type),nFaces,startFace
                exit inner_loop

            endif

        enddo inner_loop

        elseif (line_string(1:1) == "}") then
            cycle outer_loop
        endif

        if (line_string(1:1) == ")") then
            exit outer_loop
        endif

    end do outer_loop

end subroutine


subroutine check_dict(bc_type)
!
! We give a possibility of program to 'translate' boundary types
! To do this specify key value pairs for dictionary at command line
! when we call this program, eg. 
! > readBoundaryFile -d patch zeroGrad
! The in our 'boudnary' file - instead of 'key' string, the 'value' string will be written.
!

 implicit none

  character(len=20), intent(inout) :: bc_type

  integer :: i,len

  len = size(key)

  do i=1,len

    if( trim(bc_type) == trim(key(i)) ) then
      bc_type = trim( value(i) )
    endif

  enddo

 return
end subroutine


subroutine create_field_init_files
!
! Create two template files in 0/ folder to help setting 
! initial conditions for scalar and vector fields.
!
 implicit none

 integer :: boundary_file
 integer :: nbc ! no of boundaries
 character( len = 30 ), dimension(:), allocatable :: bcName

 integer :: i
 integer :: os

 call get_unit( boundary_file )
 open(unit=boundary_file, file='polyMesh/boundary')
 rewind boundary_file

 call file_row_count ( boundary_file, nbc )

 allocate( bcName(nbc) )

 ! Rewind boundary file and skip header line
 read( boundary_file, * ) 
 do i=1,nbc
   read( boundary_file, * ) bcName(i)
 end do

 close( boundary_file )

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

subroutine file_row_count ( input_unit, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!    11 October 2017
!
!  Author:
!
!    John Burkardt
!
!  Modified:
!
!    Nikola Mirkov
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  ! character ( len = * ) input_file_name
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  ! call get_unit ( input_unit )

  ! open ( unit = input_unit, file = input_file_name, status = 'old', &
  !   iostat = input_status )

  ! if ( input_status /= 0 ) then
  !   row_num = -1;
  !   ierror = 1
  !   write ( *, '(a)' ) ' '
  !   write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
  !   write ( *, '(a,i8)' ) '  Could not open the input file "' // &
  !     trim ( input_file_name ) // '" on unit ', input_unit
  !   stop
  ! end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  ! close ( unit = input_unit )
  rewind input_unit

  return
end

end module boundary_module