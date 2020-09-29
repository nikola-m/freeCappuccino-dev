    module fieldManipulation
!
!***********************************************************************
!
  use types
  use gradients
  use interpolation, only: face_value_w_option

  implicit none

  !****************************************************
  ! The inerpolations that we do for e.g. divergence is hardcoded here.
  ! Recommended are: 'central', 'cds', 'cdscorr'
  ! NOTE: 
  ! . scheme = 'cdscorr' or scheme = 'central' are good if we have specified ( interpolation_coeff_variant == 1 )
  !   in the "geometry" module. 
  ! . scheme = 'cds' or scheme = 'central' are good if we have specified ( interpolation_coeff_variant == 2 )
  !   in the "geometry" module.
  ! The "interpolation_coeff_variant" is set as parameter in the header of the "geometry" module and determines 
  ! how we calculate interpolation factor. Option "2" is typical for CFD codes, while "1" is something specific
  ! for freeCappuccino.
  !****************************************************
  character(len=10), parameter :: scheme = 'cds' 

  interface explDiv
    module procedure explDiv
    module procedure explDivMdot
  end interface

  public

  contains

!***********************************************************************
!
  pure function volumeWeightedAverage(U) result(wAvgU)
!
!***********************************************************************
    use geometry, only:numTotal,numCells,vol

    implicit none

!...Output
    real(dp) :: wAvgU

!...Input
    real(dp), dimension(numTotal), intent(in) :: U

!...Locals
    integer :: inp
    real(dp) :: sumvol
!
!***********************************************************************
!  
    sumvol = 0.0_dp
    wAvgU = 0.0_dp 

      do inp=1,numCells
          wAvgU = wAvgU + (Vol(inp)*U(inp))
          sumvol = sumvol + vol(inp)
      enddo
    
    wAvgU = wAvgU / sumvol

  end function


!***********************************************************************
!
  subroutine calcPressDiv 
!
!***********************************************************************
! 
!  -(nabla p) is a vector field: -( (nabla p)x . i + (nabla p)y . j + (nabla p)z . k )
!  Instead of computing the Grad(p) vector field in center and integrate volumetricaly
!  simply multiplying it by Vol(ijp), we write it in divergence form.
!  Interpolate to face {-(nabla p),i . e_i}_f * S_f
!  where ',i' is a partial derivative with respect to i, i={x,y,z},
!  '.' is a scalar product of vectors,
!  and {}_f is interpolation to face center.
!  And get {-(nabla p),i}_f * (S_f. e_i) = {-(nabla p),i}_f * S_fi
!  Source(u_i) = sum_over_cell_faces {-(nabla p),i}_f * S_fi
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
      call bpres(p,istage)
      ! Calculate pressure gradient.
      call grad(p,dPdxi)
    end do

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      ! ! Linear interpolation of rpessure to face
      ! call presFaceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
      !                   p, dPdxi,dfxe,dfye,dfze)


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

    ! Contribution from boundaries

    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
      call presFaceDivBoundary(arx(iface), ary(iface), arz(iface), p(ijb), su(ijp), sv(ijp), sw(ijp))

      ! su(ijp) = su(ijp) - p(ijb)*arx(iface)
      ! sv(ijp) = sv(ijp) - p(ijb)*ary(iface)
      ! sw(ijp) = sw(ijp) - p(ijb)*arz(iface)

    enddo

  return
  end


!***********************************************************************
!
  function explDiv(u,v,w) result(div)
!
!***********************************************************************
!  
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: div

!...Input
    real(dp), dimension(numTotal) :: u,v,w

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: dfxe
    real(dp), dimension(3,numTotal) :: dUdxi,dvdxi,dWdxi

    ! Calculate cell-centered gradient
    !call updateBoundary(u,'boundary_region_name', zeroGrad/value/nGrad)
    call grad(u,dUdxi)
    call grad(v,dVdxi)
    call grad(w,dWdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u,v,w, dUdxi,dvdxi,dWdxi, dfxe)

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp)+dfxe
      div(ijn) = div(ijn)-dfxe

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb),v(ijb),w(ijb), dfxe)   
      div(ijp) = div(ijp) + dfxe
    enddo

  return
  end



!***********************************************************************
!
  function explDivMdot(flmass,u,v,w) result(div)
!
!***********************************************************************
!  
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: div

!...Input
    real(dp), dimension(numTotal) :: flmass
    real(dp), dimension(numTotal) :: u,v,w

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: dfxe
    real(dp), dimension(3,numTotal) ::  dUdxi,dvdxi,dWdxi

    ! Calculate cell-centered gradient
    !call updateBoundary(u,'boundary_region_name', zeroGrad/value/nGrad)
    call grad(u,dUdxi)
    call grad(v,dVdxi)
    call grad(w,dWdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u,v,w,  dUdxi,dvdxi,dWdxi ,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp) + flmass(i) * dfxe

      div(ijn) = div(ijn) - flmass(i) * dfxe

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb),v(ijb),w(ijb), dfxe)
      div(ijp) = div(ijp) + flmass(iface) * dfxe
    enddo

  return
  end



!
!***********************************************************************
!
  subroutine presFaceDivInner(ijp,ijn, &
                              xfc,yfc,zfc,sx,sy,sz,fif, &
                              fi,df,dfxe,dfye,dfze)
!
!***********************************************************************
!
    use geometry, only: numTotal

    implicit none

    integer,  intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numTotal), intent(in) :: df
    real(dp), intent(out) :: dfxe,dfye,dfze

    real(dp) :: fie
!
!***********************************************************************
!

  ! Value of the variable at cell-face center
  fie = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, fi, df, scheme )

  ! (interpolated mid-face value)x(area)
  dfxe = fie*sx
  dfye = fie*sy
  dfze = fie*sz

  end subroutine


!
!***********************************************************************
!
  subroutine faceDivInner(ijp,ijn, &
                          xfc,yfc,zfc,sx,sy,sz,fif, &
                          u,v,w,du,dv,dw,dfxe)
!
!***********************************************************************
! 
    use geometry, only: numTotal

    implicit none

    integer,  intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: u,v,w
    real(dp), dimension(3,numTotal), intent(in) :: du,dv,dw
    real(dp), intent(out) :: dfxe

    real(dp) :: uf,vf,wf
!
!***********************************************************************
!

  ! Value of the variable at cell-face center
  uf = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, u, du, scheme )
  vf = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, v, dv, scheme )
  wf = face_value_w_option( ijp, ijn, xfc, yfc, zfc, fif, w, dw, scheme )

  ! (interpolated mid-face value)x(area)
  dfxe = uf*sx + vf*sy + wf*sz

  end subroutine


!***********************************************************************
!
  subroutine faceDivBoundary(sx,sy,sz,uf,vf,wf,dfx)
!
!***********************************************************************
!
    implicit none

    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: uf,vf,wf
    real(dp), intent(out)  :: dfx
!
!***********************************************************************
!
    dfx = uf*sx + vf*sy + wf*sz

  end subroutine

!***********************************************************************
!
  subroutine presFaceDivBoundary(sx,sy,sz,fi,dfx,dfy,dfz)
!
!***********************************************************************
!
    implicit none

    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fi
    real(dp), intent(inout)  :: dfx,dfy,dfz
!
!***********************************************************************
!
    dfx = dfx - fi*sx
    dfy = dfy - fi*sy
    dfz = dfz - fi*sz

  end subroutine



!***********************************************************************
!
  function average(u) result(aver)
!
!***********************************************************************
!                                       
!  average(u) = ( sum_{i=1}^{i=nf} (u)_f*sf) / (sum_{i=1}^{i=nf} sf)
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: aver

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: ui
    real(dp) :: are
    real(dp), dimension(3,numTotal) :: dUdxi
    real(dp), dimension(numCells) :: sumsf

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)


    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 ) 

      ui = face_value_w_option( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )

      ! Accumulate contribution at cell center and neighbour.
      aver(ijp) = aver(ijp)+ui*are
      aver(ijn) = aver(ijn)-ui*are

      sumsf(ijp) = sumsf(ijp) + are
      sumsf(ijn) = sumsf(ijn) + are

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i

      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      aver(ijp) = aver(ijp)+u(ijb)*are
      sumsf(ijp) = sumsf(ijp) + are

    enddo

    ! Divide by bounding faces total area
    do ijp = 1,numCells
      aver(ijp) = aver(ijp)/sumsf(ijp)
    enddo

  end function



!***********************************************************************
!
  function field_interpolate(u) result(ui)
!
!***********************************************************************
!  
!  volScalarField -> surfaceScalarField by interpolation.
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numFaces) :: ui

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp), dimension(3,numTotal) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      ui(i) = face_value_w_option( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )
    enddo

    ! Update boundaries?

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijb = numCells + i
      ui(iface) = u(ijb)
    enddo

  end function



!***********************************************************************
!
  function surfaceSum(u) result(ssum)
!
!***********************************************************************
!  
!  volScalarField -> volScalarField by interpolation + face summation.
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: ssum

!...Input
    real(dp), dimension(numTotal), intent(in) :: u

!...Local
    integer :: i,ijp,ijn,ijb,if
    real(dp) :: ui
    real(dp), dimension(3,numTotal) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      ui = face_value_w_option( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )

      ! Accumulate contribution at cell center and neighbour.
      ssum(ijp) = ssum(ijp)+ui

      ssum(ijn) = ssum(ijn)-ui

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      if = numInnerFaces + i
      ijp = owner(if)
      ijb = numCells + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

  end function



!***********************************************************************
!
  function surfaceIntegrate(u) result(ssum)
!
!***********************************************************************
!  
!  volScalarField -> volScalarField
!  Performs a sum of face values bounding each cell and dividing by the cell volume.
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: ssum

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,if
    real(dp) :: ui
    real(dp), dimension(3,numTotal) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      ui = face_value_w_option( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )

      ! Accumulate contribution at cell center and neighbour.
      ssum(ijp) = ssum(ijp)+ui

      ssum(ijn) = ssum(ijn)-ui

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      if = numInnerFaces + i
      ijp = owner(if)
      ijb = numCells + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    ! Divide by cell volume
    do ijp = 1,numCells
      ssum(ijp) = ssum(ijp)/vol(ijp)
    enddo

  end function



!***********************************************************************
!
  subroutine add_random_noise_to_field(Phi,percent)
!
!***********************************************************************
!
!     Level is Max perturbation aplitude in percent (%) from given field value.
!   
!     Example usage:
!       call add_random_noise_to_field(U,10)
!
  use geometry, only: numCells,numTotal
  use utils, only: init_random_seed

  implicit none
!
!***********************************************************************
!
  real(dp), dimension(numTotal), intent(inout) :: phi
  integer, intent(in) :: percent
  
  integer :: inp
  real(dp) :: level
  real(dp) :: perturb

  level = dble(percent)

  do inp=1,numCells

        ! Random number based fluctuation of mean profile            
        CALL init_random_seed()
        CALL RANDOM_NUMBER(perturb)
        
        ! perturb is now between 0. and 1., we want it to be from 0 to 2*amplitude
        ! e.g. perturb = 0.9+perturb/5. when Max perturbation is +/- 10% of mean profile
        perturb = ( 1.0_dp - level/100.0_dp ) + perturb * (2*level/100.0_dp)

        Phi(INP) = perturb*Phi(inp)

  enddo

  end subroutine

    
end module fieldManipulation
