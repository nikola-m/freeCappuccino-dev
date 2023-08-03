module nablap

use types
use interpolation, only: face_value_central, face_value_cds

implicit none

  character(len=10) :: pscheme = 'linear' ! Interpolation scheme for pressure: "linear", "central" or "weighted".

private

public :: gradp_and_sources
public :: pscheme ! parameter for input.nml

contains
    
!***********************************************************************
!
subroutine gradp_and_sources(p)
!
!***********************************************************************

  use geometry
  use parameters, only: small
  use variables, only: dPdxi
  use sparse_matrix, only: su,sv,sw,apu

  implicit none

  real(dp), dimension(numTotal), intent(in) :: p

  ! Local
  integer :: i,ijp,ijn,ijb,iface,istage
  real(dp) :: dfxe,dfye,dfze
  real(dp) :: pf
  real(dp) :: volr

!
!***********************************************************************
!
  ! Clear source arrays
  su = 0.0_dp
  sv = 0.0_dp
  sw = 0.0_dp  


  ! Accumulate  sum(-pf*Sf.x), sum(-pf*Sf.y), and sum(-pf*Sf.x) in source vectors su, sv, and sw.
  ! for each cell, by looping over inner faces at first.


  if ( pscheme == 'linear' .or. pscheme == 'central' ) then

    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      ! Value of the variable at cell-face center
      pf = face_value_cds( ijp, ijn, facint(i), p )


      ! Contribution =(interpolated mid-face value)x(area)
      dfxe = pf*arx(i)
      dfye = pf*ary(i)
      dfze = pf*arz(i)

      ! Accumulate contribution at cell center and neighbour.
      ! ***NOTE, we calculate negative Divergence, therefore opposite sign (minus in front of e.g. dfxe, etc.) 
      su(ijp) = su(ijp)-dfxe
      sv(ijp) = sv(ijp)-dfye
      sw(ijp) = sw(ijp)-dfze
       
      su(ijn) = su(ijn)+dfxe
      sv(ijn) = sv(ijn)+dfye
      sw(ijn) = sw(ijn)+dfze

    enddo

  elseif ( pscheme == 'weighted' ) then

    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      ! Pressure on face based on "Standard" interpolation in Fluent.
      ! This is weighted interpolation where weights are mass flows estimated at respective cell center
      pf = ( p(ijp)*Apu(ijp)+p(ijn)*Apu(ijn) ) / ( Apu(ijp) + Apu(ijn) + small )


      ! Contribution =(interpolated mid-face value)x(area)
      dfxe = pf*arx(i)
      dfye = pf*ary(i)
      dfze = pf*arz(i)

      ! Accumulate contribution at cell center and neighbour.
      ! ***NOTE, we calculate negative Divergence, therefore opposite sign (minus in front of e.g. dfxe, etc.) 
      su(ijp) = su(ijp)-dfxe
      sv(ijp) = sv(ijp)-dfye
      sw(ijp) = sw(ijp)-dfze
       
      su(ijn) = su(ijn)+dfxe
      sv(ijn) = sv(ijn)+dfye
      sw(ijn) = sw(ijn)+dfze

    enddo

  else

    write(*,'(a)') ' '
    write(*,'(a)') 'Fatal error: non-existing interpolation scheme for pressure!'
    stop

  endif


  ! From accumulation in su, sv, and sw we'll calculate the pressure gradient dPdxi.

 
  ! Update pressure at boundary faces and calculate gradient in two stages.
  ! Stage 1: by simple zero gradient extrapolation to boundaries;
  ! Stage 2: by linear extrapolation using newly created gradients to update p at wall boundries.
  do istage=1,2

    ! Update boundary pressure
    call bpres(p,istage)


    ! Recalculate sum(-pf*Sf) over inner faces in second stage only if we use 'central' scheme
    ! because in the first stage we didn't have necessary pressure gradients dPdxi.
    if ( istage == 2 .and. pscheme == 'central' ) then

      ! Clear source arrays
      su = 0.0_dp
      sv = 0.0_dp
      sw = 0.0_dp 

      do i=1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        pf = face_value_central( ijp, ijn, xf(i), yf(i), zf(i),  p, dPdxi ) 


        ! Contribution =(interpolated mid-face value)x(area)
        dfxe = pf*arx(i)
        dfye = pf*ary(i)
        dfze = pf*arz(i)

        ! Accumulate contribution at cell center and neighbour.
        ! ***NOTE, we calculate negative Divergence, therefore opposite sign (minus in front of e.g. dfxe, etc.) 
        su(ijp) = su(ijp)-dfxe
        sv(ijp) = sv(ijp)-dfye
        sw(ijp) = sw(ijp)-dfze
         
        su(ijn) = su(ijn)+dfxe
        sv(ijn) = sv(ijn)+dfye
        sw(ijn) = sw(ijn)+dfze

      enddo

    endif


    ! Take stored accumulation - note the minus sign because in su,sv and sw, we have negative divergence.
    dpdxi(1,1:numCells) = -su(1:numCells) 
    dpdxi(2,1:numCells) = -sv(1:numCells) 
    dpdxi(3,1:numCells) = -sw(1:numCells) 


    ! Contribution from boundaries.
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
    
      dpdxi(1,ijp) = dpdxi(1,ijp) + p(ijb)*arx(iface)
      dpdxi(2,ijp) = dpdxi(2,ijp) + p(ijb)*ary(iface)
      dpdxi(3,ijp) = dpdxi(3,ijp) + p(ijb)*arz(iface)

    enddo

    ! Finally calculate gradient components at cv-centers by dividing the accumulated sum by vol(inp)
    do ijp=1,numCells
      volr = 1.0_dp/vol(ijp)
      dpdxi(:,ijp) = dpdxi(:,ijp)*volr
    enddo

  end do


  ! Now when we have the correct p values at boundaries we'll add 
  ! boundary faces contribution to su,sv, and sw vectors.


  ! Contribution from boundaries.
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijb = numCells + i
  
    su(ijp) = su(ijp) - p(ijb)*arx(iface)
    sv(ijp) = sv(ijp) - p(ijb)*ary(iface)
    sw(ijp) = sw(ijp) - p(ijb)*arz(iface)

  enddo



end subroutine


end module