module simulation_control

 implicit none


subroutine read_input_file
!
! Open & Read and Process Input File
!

  implicit none



  open( unit=inpt, file=input_file ) 
  rewind inpt

  do ! until EOF

  read (inpt, *) key
   if (key == '#' .or. key == ' ' ) cycle ! # for a comment in input file or empty string

   if (key == 'Temperature') then
     read (inpt,*) lcal(7)

   if (key == 'limiter') then
     read (inpt,*) limiter



  enddo

end subroutine

end module