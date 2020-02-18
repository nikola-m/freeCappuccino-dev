!***********************************************************************
!
subroutine writefiles
!
! Write output files in Paraview .vtu format
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

  implicit none
!
!***********************************************************************
!
  integer :: output_unit
  character( len = 5) :: nproc_char

  ! Write in a char variable current timestep number
  call i4_to_s_left ( itime, timechar )

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

  !+-----------------------------------------------------------------------------+

  call get_unit( output_unit )

  open(unit=output_unit,file='vtk/solution_fields_'//'proc_'//trim(nproc_char)//'_'//trim(timechar),form='unformatted')

  write(output_unit) numCells
  write(output_unit) numTotal
  write(output_unit) u
  write(output_unit) v
  write(output_unit) w
  write(output_unit) p
  write(output_unit) te
  write(output_unit) ed
  write(output_unit) vis

  ! Header
  ! call vtu_write_XML_header ( output_unit )

  ! !+-----------------------------------------------------------------------------+

  ! ! Write fields to VTU file

  ! call vtu_write_XML_vector_field( output_unit, 'U', u, v, w )


  ! call vtu_write_XML_scalar_field ( output_unit, 'p', p )


  ! if( lturb ) then

  !   call vtu_write_XML_scalar_field ( output_unit, 'mueff', vis )

  ! endif


  ! if(solveTKE) then

  !   call vtu_write_XML_scalar_field ( output_unit, 'k', te )

  ! endif


  ! if( solveEpsilon ) then

  !   call vtu_write_XML_scalar_field ( output_unit, 'epsilon', ed )

  ! endif


  ! if( solveOmega ) then

  !   call vtu_write_XML_scalar_field ( output_unit, 'omega', ed )

  ! endif      

  ! !+-----------------------------------------------------------------------------+

  ! Mesh data
  ! call vtu_write_XML_meshdata ( output_unit )

  ! close( output_unit )


end subroutine