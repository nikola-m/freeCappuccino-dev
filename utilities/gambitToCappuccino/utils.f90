module utils
contains

! Gambit element types NTYPE:
!	1 = Edge
!	2 = Quadrilateral
!	3 = Triangle
!	4 = Brick
!	5 = Wedge (Prism)
!	6 = Tetrahedron
!	7 = Pyramid 

 pure integer function noel(NTYPE)
!
! Number of nodes in element of NTYPE Gmsh type of element
!
 implicit none
 integer, intent(in) :: NTYPE
 integer, parameter :: noel_NTYPE1 = 2 ! no. of nodes in element-here Line
 integer, parameter :: noel_NTYPE2 = 4 ! no. of nodes in element-here Quad
 integer, parameter :: noel_NTYPE3 = 3 ! no. of nodes in element-here Tri
 integer, parameter :: noel_NTYPE4 = 8 ! no. of nodes in element-here Hex
 integer, parameter :: noel_NTYPE5 = 6 ! no. of nodes in element-here Prism
 integer, parameter :: noel_NTYPE6 = 4 ! no. of nodes in element-here Tet
 integer, parameter :: noel_NTYPE7 = 5 ! no. of nodes in element-here Pyramid
 

 if(NTYPE.eq.1) then
   noel = noel_NTYPE1 
 elseif(NTYPE.eq.2) then
   noel = noel_NTYPE2 
 elseif(NTYPE.eq.3) then
   noel = noel_NTYPE3 
 elseif(NTYPE.eq.4) then ! -> Gambit Hex
   noel = noel_NTYPE4 
 elseif(NTYPE.eq.5) then! -> Gambit Prism
   noel = noel_NTYPE5 
 elseif(NTYPE.eq.6) then! -> Gambit Tet
   noel = noel_NTYPE6 
 elseif(NTYPE.eq.7) then! -> Gambit Pyramid
   noel = noel_NTYPE7 
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
 integer, parameter :: nofael_NTYPE1 = 1
 integer, parameter :: nofael_NTYPE2 = 4 
 integer, parameter :: nofael_NTYPE3 = 3 
 integer, parameter :: nofael_NTYPE4 = 6 
 integer, parameter :: nofael_NTYPE5 = 5 
 integer, parameter :: nofael_NTYPE6 = 4 
 integer, parameter :: nofael_NTYPE7 = 5 
! Because...
! 1 = Edge
! 2 = Quadrilateral
! 3 = Triangle
! 4 = Brick
! 5 = Wedge (Prism)
! 6 = Tetrahedron
! 7 = Pyramid

 if(NTYPE.eq.1) then
   nofael = nofael_NTYPE1 
 elseif(NTYPE.eq.2) then
   nofael = nofael_NTYPE2 
 elseif(NTYPE.eq.3) then
   nofael = nofael_NTYPE3 
 elseif(NTYPE.eq.4) then
   nofael = nofael_NTYPE4 
 elseif(NTYPE.eq.5) then
   nofael = nofael_NTYPE5 
 elseif(NTYPE.eq.6) then
   nofael = nofael_NTYPE6 
 elseif(NTYPE.eq.7) then
   nofael = nofael_NTYPE7 
 else
   nofael = 0
 endif

 end function

 pure integer function paraview_ntype(NTYPE)
!
! Element type in Paraview corresponding to Gmbit type of element given by NTYPE.
!
 implicit none
 integer, intent(in) :: NTYPE

 if(NTYPE.eq.1) then ! -> Gambit line
   paraview_ntype = 3 
 elseif(NTYPE.eq.2) then ! -> Gambit Tri
   paraview_ntype = 5 
 elseif(NTYPE.eq.3) then ! -> Gambit Quad
   paraview_ntype = 9 
 elseif(NTYPE.eq.6) then ! -> Gambit Tet
   paraview_ntype = 10   
 elseif(NTYPE.eq.4) then ! -> Gambit Hex
   paraview_ntype = 12 
 elseif(NTYPE.eq.5) then ! -> Gambit Prism
   paraview_ntype = 13 
 elseif(NTYPE.eq.7) then ! -> Gambit Pyramid
   paraview_ntype = 14 
 else
   paraview_ntype = 0
 endif

 end function

 pure integer function ntypeCappuccino(NTYPE)
!
! Element type in Cappuccino corresponding to NTYPE Gmsh type of element
! NOTE: this is not neccessary because Cappuccino element types are the
! same as in GMSH, but we leave it here, if we change our mind in future.
!
 implicit none
 integer, intent(in) :: NTYPE
 integer, parameter :: cappuccinoNTYPE1 = 1 !  Line
 integer, parameter :: cappuccinoNTYPE2 = 2 !  Tri
 integer, parameter :: cappuccinoNTYPE3 = 3 !  Quad
 integer, parameter :: cappuccinoNTYPE4 = 4 !  Tet
 integer, parameter :: cappuccinoNTYPE5 = 5 !  Hex
 integer, parameter :: cappuccinoNTYPE6 = 6 !  Prism
 integer, parameter :: cappuccinoNTYPE7 = 7 !  Pyramid
 

 if(NTYPE.eq.1) then
   ntypeCappuccino = 1 
 elseif(NTYPE.eq.2) then
   ntypeCappuccino = 2 
 elseif(NTYPE.eq.3) then
   ntypeCappuccino = 3 
 elseif(NTYPE.eq.4) then ! -> Gambit Hex
   ntypeCappuccino = 5 
 elseif(NTYPE.eq.5) then ! -> Gambit Prism
   ntypeCappuccino = 6 
 elseif(NTYPE.eq.6) then ! -> Gambit Tetrahedron
   ntypeCappuccino = 4 
 elseif(NTYPE.eq.7) then ! -> Gambit Pyramid
   ntypeCappuccino = 7
 else
   ntypeCappuccino = 0
 endif

 end function

 function get_face_vertices_NTYPE4(iface,NODE,NDP) result(vertices)
! Gambit gives a convetion of how to generate faces from element(cell) vertex data
! So in Gambit faces are not given explicitely, but we have a rule how to find them. 
! Something similar happens in other mesh generation codes  but the rule is different
! so to make Cappuccino preprocessor for any other mesh file format you have to consult
! their documentation on how to generate faces.
!
! NTYPE4 is 8-node HEX
!
! So the Gambit convention is:
! 	Face 	Nodes
!       1 	0,1,5,4
! 	2 	1,3,7,5
! 	3 	3,2,6,7
!	4 	2,0,4,6
! 	5 	1,0,2,3
! 	6 	4,5,7,6
! NOTE: Gambit starts counting from 0, we here from 1!
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
       vertices(:) = (/ NODE(1), NODE(2), NODE(6), NODE(5) /)
    case (2)
       vertices(:) = (/ NODE(2), NODE(4), NODE(8), NODE(6) /)
    case (3)
       vertices(:) = (/ NODE(4), NODE(3), NODE(7), NODE(8) /)
    case (4)
       vertices(:) = (/ NODE(3), NODE(1), NODE(5), NODE(7) /)
    case (5)
       vertices(:) = (/ NODE(2), NODE(1), NODE(3), NODE(4) /)
    case (6)
       vertices(:) = (/ NODE(5), NODE(6), NODE(8), NODE(7) /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function


 function get_face_vertices_NTYPE7(iface,NODE,NDP) result(vertices)
!
! NTYPE7 is 5-node PYRAMID
!
! So the Gambit convention is:
! 	Face 	Nodes
! 	1 	0,2,3,1
! 	2 	0,1,4
! 	3 	1,3,4
! 	4 	3,2,4
! 	5 	2,0,4
! NOTE: Gambit starts counting from 0, we here from 1!
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
       vertices(:) = (/ NODE(1), NODE(3), NODE(4), NODE(2) /)
    case (2)
       vertices(:) = (/ NODE(1), NODE(2), NODE(5), 0/)
    case (3)
       vertices(:) = (/ NODE(2), NODE(4), NODE(5), 0 /)
    case (4)
       vertices(:) = (/ NODE(4), NODE(3), NODE(5), 0 /)
    case (5)
       vertices(:) = (/ NODE(3), NODE(1), NODE(5), 0 /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function


 function get_face_vertices_NTYPE6(iface,NODE,NDP) result(vertices)
!
! NTYPE6 is 4-node TETRAHEDRON
!
! So the Gambit convention is:
! 	Face 	Nodes
! 	1 	1,0,2
! 	2 	0,1,3
! 	3 	1,2,3
! 	4 	2,0,3
! NOTE: Gambit starts counting from 0, we here from 1!
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
       vertices(:) = (/ NODE(2), NODE(1), NODE(3), 0 /)
    case (2)
       vertices(:) = (/ NODE(1), NODE(2), NODE(4), 0/)
    case (3)
       vertices(:) = (/ NODE(2), NODE(3), NODE(4), 0 /)
    case (4)
       vertices(:) = (/ NODE(3), NODE(1), NODE(4), 0 /)
    case default
       vertices = (/ 0,0,0,0 /)
 end select

 end function

 function get_face_vertices_NTYPE5(iface,NODE,NDP) result(vertices)
!
! NTYPE5 is 6-node PRISM
!
! So the Gambit convention is:
! 	Face 	Nodes
! 	1 	0,1,4,3
! 	2 	1,2,5,4
! 	3 	2,0,3,5
! 	4 	0,2,1
! 	5 	3,4,5
! NOTE: Gambit starts counting from 0, we here from 1!
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

end module utils
