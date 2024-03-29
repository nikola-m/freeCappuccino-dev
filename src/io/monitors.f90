module monitors

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    A module containing functions defining monitored variables during simulation.
!
!  Discussion:
!
!    It is often useful to monitor integral values during simulation to observe convergence process or to
!    follow temporal development of certain values.
!    To define a monitor you need a geometrical entity, such as a wall boundary defined by its name,
!    and integral quantity such as wall flux, force, force coefficient etc. recalculated and logged
!    every iteration and/or timestep.
!    Logged values are written in a dedicated file and can be plot e.g. using gnuplot during simulation run.
!
!
!  The following is the recap od the use of reference values (source: ANSYS Fluent User's guide):
!
!   * Force coefficients use the reference area, density, and velocity. In addition, the pressure force calculation 
!     uses the reference pressure.
!   * Moment coefficients use the reference length, area, density and velocity. In addition, the pressure force calculation 
!     uses the reference pressure.
!   * Reynolds number uses the reference length, density, and viscosity.
!   * Pressure and total pressure coefficients use the reference pressure, density, and velocity.
!   * Entropy uses the reference density, pressure, and temperature.
!   * Skin friction coefficient uses the reference density and velocity.
!   * Heat transfer coefficient uses the reference temperature.
!   * Turbomachinery efficiency calculations use the ratio of specific heats.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

  integer :: numMonVals = 0
  character( len = 30 ), dimension(:), allocatable :: monitored_value
  character( len = 30 ), dimension(:), allocatable :: monitored_location

  ! Directions to project lift and drag force
  ! real(dp), dimension(3) :: liftDir = (/ 0, 1, 0 /) 
  ! real(dp), dimension(3) :: dragDir = (/ 1, 0, 0 /)

  private

  public :: define_monitors, log_monitored_values

contains

subroutine define_monitors
!
!  Purpose: 
!
!    This is a driver subroutine for monitoring.
!
!  Discussion:
!
!    Read 'monitors' dictionary and define 
!    numMonVals, 
!    monitored_value(1,numMonVals)
!    monitored_location(1,numMonVals)
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    Feb 2021
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!
  use utils, only: get_unit

  implicit none

  integer :: i
  integer :: dfile
  logical :: file_exists

  ! Check if there is monitors file 
  inquire(FILE='monitors', EXIST=file_exists)

  if (file_exists) then

    call get_unit( dfile )

    open( unit = dfile, file='monitors' )
    rewind dfile

    read(dfile,*) numMonVals

    allocate( monitored_value(numMonVals), monitored_location(numMonVals) )

    do i=1,numMonVals
      read(dfile,*) monitored_value(i), monitored_location(i)
    enddo

    close( dfile )

  endif

end subroutine


subroutine log_monitored_values
!
!  Purpose: 
!
!    This is a driver subroutine for monitoring.
!
!  Discussion:
!
!    In the main routine we call this driver. 
!    The module has its own dictionary where we define number of monitored values,
!    we define which exactly are they (chosen from predefined given group, based on the 
!    monitoring subroutines we have here), and finally for each monitored value we have
!    a specific location, whether it's a boundary region, cell region, whole inner domain
!    (if we track volume weigthed averages), or certain points.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    Feb 2021
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!

  implicit none

  integer :: i

  do i=1,numMonVals

    select case( monitored_value(i) )

      case( 'surf_avg_shear' )
        call surf_avg_shear( monitored_location(i) )

      case( 'mean_ke_diss')
        call mean_ke_dissipation

      ! case( 'field_value' )
      !   call field_value( monitored_location(i) )

    end select

  enddo



end subroutine


subroutine surf_avg_shear( boundary_name )
!
!  Purpose: 
!
!    This subroutine computes and writes to monitor surface average value of Shear (Force).
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    27 December 2020
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Parameters:
!
!    Input, character, boundary_name, name of the boundary region for which we calculate wss.
!


  implicit none
!
! Parameters
!
  character (len = *), intent(in) :: boundary_name
!
! Local variables
!
  integer :: i, ib, iface, ijp, ijb, iWall
  real(dp) :: viss,are,nxf,nyf,nzf,vsol,upb,vpb,wpb
  real(dp) :: fshx, fshy, fshz, sumss, sumsf

  sumss = 0.
  sumsf = 0.

  iWall = 0

  ! Loop over boundaries
  do ib=1,numBoundaries
    
    if ( bctype(ib) == 'wall' .and. bcname(ib) /= boundary_name ) then
 
      iWall = iWall + 1

    elseif ( bctype(ib) == 'wall' .and. bcname(ib) == boundary_name ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        ! Diffusion coef. 
        viss=max(viscos,visw(iWall))
        vsol = viss*srdw(iWall)

        ! Face area
        are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

        ! Face normals
        nxf = arx(iface)/are
        nyf = ary(iface)/are
        nzf = arz(iface)/are


        ! Velocity difference vector components
        upb = u(ijp)-u(ijb)
        vpb = v(ijp)-v(ijb)
        wpb = w(ijp)-w(ijb)

        ! Shear force at wall projected in x, y and z direction.
        fshx = vsol*( (u(ijb)-u(ijp))*(1.-nxf**2) + vpb*nyf*nxf                  + wpb*nzf*nxf )
        fshy = vsol*( upb*nxf*nyf                 + (v(ijb)-v(ijp))*(1.-nyf**2)  + wpb*nzf*nyf )
        fshz = vsol*( upb*nxf*nzf                 + vpb*nyf*nzf                  + (w(ijb)-w(ijp))*(1.-nzf**2) )

        sumss = sumss + are * sqrt( fshx**2 + fshy**2 + fshz**2)
        sumsf = sumsf + are

      end do

    endif

  enddo

  write(*,'(3a,es15.7)') '  Shear at ', trim( boundary_name ) ,' =', sumss / sumsf

end subroutine


subroutine mean_ke_dissipation
!
!  Purpose: 
!
!    This subroutine computes and writes to monitor the dissipation of mean resolved kinetic energy.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    16. August 2022
!
!  Author:
!
!    Nikola Mirkov/Email: largeddysimulation@gmail.com
!
!  Parameters:
!
!    -
!
 
  use fieldManipulation, only: volumeWeightedAverage

  implicit none

!
! Local variables
!
  real(dp) :: K,Ko,diss


  Ko = volumeWeightedAverage( uo*uo+vo*vo+wo*wo )  

  K = volumeWeightedAverage( u*u+v*v+w*w )

  diss = -0.5*(K-Ko) / timestep


  write(*,'(1a,es15.7)') '  Mean resolved kinetic energy =', K
  write(*,'(1a,es15.7)') '  Dissipation of the mean resolved kinetic energy =', diss

end subroutine


end module



