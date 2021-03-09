module nablap

use types
use gradients
use interpolation, only: face_value_central

implicit none

  logical :: ScndOrderPressIntrp  = .True. ! Second order (default), or massflow weighted interpolation of pressure to faces.

private

public :: surfaceIntegratePressure,surfaceIntegratePressureCrankNicolson,surfaceIntegratePressureCorr
public :: ScndOrderPressIntrp  ! parameter for input.nml

contains
    
!***********************************************************************
!
subroutine surfaceIntegratePressure 
!
!***********************************************************************
! 
!  -(nabla p) is a vector field: -( (nabla p)x . i + (nabla p)y . j + (nabla p)z . k )
!  Instead of computing the Grad(p) vector field in center and integrate volumetricaly
!  multiplying cell centroid value by Vol(ijp), we write it in divergence form.
!  Interpolate to face {-p e_i}_f * S_f
!  Source(u_i) = sum_over_cell_faces {-p}_f * S_fi
!  Interpolation to cell face centers done by central scheme.
!
!***********************************************************************
!
  use geometry
  use parameters, only: small
  use variables, only: p,dPdxi
  use sparse_matrix, only: su,sv,sw,apu

  implicit none

  ! Local
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: dfxe,dfye,dfze
  real(dp) :: pf

!
!***********************************************************************
!

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      if ( ScndOrderPressIntrp  ) then

        ! Value of the variable at cell-face center
        pf = face_value_central( ijp, ijn, xf(i), yf(i), zf(i),  p, dPdxi )

      else

        ! Pressure on face based on "Standard" interpolation in Fluent.
        ! This is weighted interpolation where weights are mass flows estimated at respective cell center
        pf = ( p(ijp)*Apu(ijp)+p(ijn)*Apu(ijn) ) / ( Apu(ijp) + Apu(ijn) + small )

        
      endif


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

    ! Contribution from boundaries

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
    
      su(ijp) = su(ijp) - p(ijb)*arx(iface)
      sv(ijp) = sv(ijp) - p(ijb)*ary(iface)
      sw(ijp) = sw(ijp) - p(ijb)*arz(iface)

    enddo

end subroutine


!***********************************************************************
!
subroutine surfaceIntegratePressureCorr
!
!***********************************************************************
! 
!  -(nabla p) is a vector field: -( (nabla p)x . i + (nabla p)y . j + (nabla p)z . k )
!  Instead of computing the Grad(p) vector field in center and integrate volumetricaly
!  multiplying cell centroid value by Vol(ijp), we write it in divergence form.
!  Interpolate to face {-p e_i}_f * S_f
!  Source(u_i) = sum_over_cell_faces {-p}_f * S_fi
!  Interpolation to cell face centers done by central scheme.
!
!***********************************************************************
!
  use geometry
  use parameters, only: small
  use variables, only: pp,dPdxi
  use sparse_matrix, only: su,sv,sw,apu

  implicit none

  ! Local
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: dfxe,dfye,dfze
  real(dp) :: pf

!
!***********************************************************************
!

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      if ( ScndOrderPressIntrp  ) then

        ! Value of the variable at cell-face center
        pf = face_value_central( ijp, ijn, xf(i), yf(i), zf(i),  pp, dPdxi )

      else

        ! Pressure on face based on "Standard" interpolation in Fluent.
        ! This is weighted interpolation where weights are mass flows estimated at respective cell center
        pf = ( pp(ijp)*Apu(ijp)+pp(ijn)*Apu(ijn) ) / ( Apu(ijp) + Apu(ijn) + small )

        
      endif


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

    ! Contribution from boundaries

    ! iPer = 0

    ! ! Loop over boundaries
    ! do ib=1,numBoundaries
      
    !   if ( bctype(ib) /= 'periodic' ) then

    !     do i=1,nfaces(ib)

    !       iface = startFace(ib) + i
    !       ijp = owner(iface)
    !       ijb = iBndValueStart(ib) + i

    !       su(ijp) = su(ijp) - pp(ijb)*arx(iface)
    !       sv(ijp) = sv(ijp) - pp(ijb)*ary(iface)
    !       sw(ijp) = sw(ijp) - pp(ijb)*arz(iface)

    !     end do


    !   elseif (  bctype(ib) == 'periodic' ) then

    !     iPer = iPer + 1

    !     ! Faces trough periodic boundaries
    !     do i=1,nfaces(ib)

    !       if = startFace(ib) + i
    !       ijp = owner(if)


    !       iftwin = startFaceTwin(iPer) + i
    !       ijn = owner(iftwin)

    !       ijb = iBndValueStart(ib) + i

    !       su(ijp) = su(ijp) - pp(ijb)*arx(iface)
    !       sv(ijp) = sv(ijp) - pp(ijb)*ary(iface)
    !       sw(ijp) = sw(ijp) - pp(ijb)*arz(iface)

    !       ijb = numCells + ( startFaceTwin(iPer) - numInnerFaces ) + i

    !       su(ijn) = su(ijn) - pp(ijb)*arx(iftwin)
    !       sv(ijn) = sv(ijn) - pp(ijb)*ary(iftwin)
    !       sw(ijn) = sw(ijn) - pp(ijb)*arz(iftwin)


    !     end do 

    !   endif 

    ! enddo 

    ! Contribution from boundaries
    
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
    
      su(ijp) = su(ijp) - pp(ijb)*arx(iface)
      sv(ijp) = sv(ijp) - pp(ijb)*ary(iface)
      sw(ijp) = sw(ijp) - pp(ijb)*arz(iface)

    enddo

end subroutine


!***********************************************************************
!
subroutine surfaceIntegratePressureCrankNicolson 
!
!***********************************************************************
! 
!  -(nabla p) is a vector field: -( (nabla p)x . i + (nabla p)y . j + (nabla p)z . k )
!  Instead of computing the Grad(p) vector field in center and integrate volumetricaly
!  multiplying cell centroid value by Vol(ijp), we write it in divergence form.
!  Interpolate to face {-p e_i}_f * S_f
!  Source(u_i) = sum_over_cell_faces {-p}_f * S_fi
!  Interpolation to cell face centers done by central scheme.
!
!***********************************************************************
!
  use geometry
  use parameters, only: small
  use variables, only: p,po,dPdxi
  use sparse_matrix, only: su,sv,sw,apu

  implicit none

  ! Local
  integer, parameter :: nipgrad = 2 
  integer :: i,ijp,ijn,ijb,iface,istage
  real(dp) :: dfxe,dfye,dfze
  real(dp) :: pf

!
!***********************************************************************
!

    ! Pressure gradient
    do istage=1,nipgrad
      ! Pressure at boundaries (for correct calculation of press. gradient)
      call bpres(po,istage)
      ! Calculate pressure gradient.
      call grad(po,dPdxi)
    end do

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      if ( ScndOrderPressIntrp  ) then

        ! Value of the variable at cell-face center
        pf = face_value_central( ijp, ijn, xf(i), yf(i), zf(i),  po, dPdxi )

      else

        ! Pressure on face based on "Standard" interpolation in Fluent.
        ! This is weighted interpolation where weights are mass flows estimated at respective cell center
        pf = ( po(ijp)*Apu(ijp)+po(ijn)*Apu(ijn) ) / ( Apu(ijp) + Apu(ijn) + small )

        
      endif

      ! For Crank-Nicolson - pressure source is split into two contributions from two consecutive timesteps
      ! each weighted by half.
      pf = 0.5_dp * pf 

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

    ! Contribution from boundaries

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i

      ! For Crank-Nicolson - pressure source is split into two contributions from two consecutive timesteps
      ! each weighted by half.
      pf = 0.5_dp * po(ijb)

      su(ijp) = su(ijp) - pf*arx(iface)
      sv(ijp) = sv(ijp) - pf*ary(iface)
      sw(ijp) = sw(ijp) - pf*arz(iface)

    enddo

  !
  ! > Pressure source from present timestep
  !
  
    ! Pressure gradient
    do istage=1,nipgrad
      ! Pressure at boundaries (for correct calculation of press. gradient)
      call bpres(p,istage)
      ! Calculate pressure gradient.
      call grad(p,dPdxi)
    end do

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      if ( ScndOrderPressIntrp  ) then

        ! Value of the variable at cell-face center
        pf = face_value_central( ijp, ijn, xf(i), yf(i), zf(i),  p, dPdxi )

      else

        ! Pressure on face based on "Standard" interpolation in Fluent.
        ! This is weighted interpolation where weights are mass flows estimated at respective cell center
        pf = ( p(ijp)*Apu(ijp)+p(ijn)*Apu(ijn) ) / ( Apu(ijp) + Apu(ijn) + small )

        
      endif

      ! For Crank-Nicolson - pressure source is split into two contributions from two consecutive timesteps
      ! each weighted by half.
      pf = 0.5_dp * pf 

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

    ! Contribution from boundaries

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i

      ! For Crank-Nicolson - pressure source is split into two contributions from two consecutive timesteps
      ! each weighted by half.
      pf = 0.5_dp * p(ijb)

      su(ijp) = su(ijp) - pf*arx(iface)
      sv(ijp) = sv(ijp) - pf*ary(iface)
      sw(ijp) = sw(ijp) - pf*arz(iface)

    enddo

end subroutine

end module