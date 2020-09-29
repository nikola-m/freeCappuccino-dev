subroutine bcin

  use types
  use parameters
  use geometry
  use variables
  use mpi

  implicit none

  integer :: i,ib,ini,ino,iface
  ! integer :: input_unit
  ! integer :: totalNumInlets
  real(dp) :: uav,outare
  ! real(dp) :: flowen, flowte, flowed
  real(dp) :: are


  flomas = 0.0_dp 
  ! flomom = 0.0_dp
  ! flowen = 0.0_dp
  ! flowte = 0.0_dp
  ! flowed = 0.0_dp


  ! In the case we want to read recirculated profiles
  ! call get_unit ( input_unit )
  ! open ( unit = input_unit, file = 'inlet')
  ! rewind input_unit


  if(ninl.gt.0) then

    ! Loop over inlet boundaries
    do ib=1,numBoundaries
      
      if ( bctype(ib) == 'inlet' ) then

        do i=1,nfaces(ib)

          iface = startFace(ib) + i
          ini = iBndValueStart(ib) + i

          ! Example of generating inflow profile for one case (Ishihara hill)
          ! The numbers are then copied to 0/U file, under inlet \ nonuniform \ <list of values>
          ! write(input_unit,'(f9.6,1x,f3.1,1x,f3.1)') 5.5d0*(zf(iface)/0.2d0)**0.135,zero,zero  

          ! Example reading recirculated input
          ! read(input_unit,*) u(ini), v(ini), w(ini)!, te(ini), ed(ini) 

          vis(ini) = viscos

          flmass(iface) = den(ini)*(arx(iface)*u(ini)+ary(iface)*v(ini)+arz(iface)*w(ini))

          ! Face normal vector is faced outwards, while velocity vector at inlet
          ! is faced inwards. That means their scalar product will be negative,
          ! so minus signs here is to turn net mass influx - flomas, into positive value.
          flomas = flomas - flmass(iface) 
          ! flomom = flomom + abs(flmass(iface))*sqrt(u(ini)**2+v(ini)**2+w(ini)**2)
          ! if(lcal(ien)) flowen = flowen + abs(flmass(iface)*t(ini))
          ! flowte = flowte + abs(flmass(iface)*te(ini))
          ! flowed = flowed + abs(flmass(iface)*ed(ini))
          

        end do

      endif 

    enddo

    ! Correct turbulence at inlet for appropriate turbulence model
    if(lturb) call correct_turbulence_inlet()

  endif  

  ! close(input_unit)

  ! Outlet area
  outare = 0.0_dp

  if (nout > 0 ) then
    ! Loop over outlet boundaries
    do ib=1,numBoundaries    
      if ( bctype(ib) == 'outlet' ) then
        do i=1,nfaces(ib)
          iface = startFace(ib) + i
          outare = outare + sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
        end do
      endif 
    enddo
  endif  


  ! Global sum - MPI communication
  call global_sum( flomas )
  ! call global_sum( flomom )
  ! if(lcal(ien)) call global_sum( flowen )
  ! call global_sum( flowte )
  ! call global_sum( flowed )
  call global_sum( outare )

  ! Average velocity at outlet boundary
  uav = flomas/(densit*outare+small)


  if( myid .eq. 0 ) then
    write ( *, '(a)' ) '  Inlet boundary condition information:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,e12.6)' ) '  Mass inflow: ', flomas
    ! write ( *, '(a,e12.6)' ) '  Momentum inflow: ', flomom
    ! write ( *, '(a,e12.6)' ) '  Temperature inflow: ', flowen  
    ! write ( *, '(a,e12.6)' ) '  TKE inflow: ', flowte      
    ! write ( *, '(a,e12.6)' ) '  Dissipation inflow: ', flowed    
  endif


  ! Mass flow trough outlet faces using Uav velocity

  if(nout.gt.0) then

    ! Loop over outlet boundaries
    do ib=1,numBoundaries
      
      if ( bctype(ib) == 'outlet' ) then

        do i=1,nfaces(ib)

          iface = startFace(ib) + i
          ino = iBndValueStart(ib) + i

          are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

          u(ino) = uav * arx(iface)/are
          v(ino) = uav * ary(iface)/are
          w(ino) = uav * arz(iface)/are

          flmass(iface)=den(ino)*(arx(iface)*u(ino)+ary(iface)*v(ino)+arz(iface)*w(ino))

        end do

      endif 

    enddo

  endif  

end subroutine