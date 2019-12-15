    module fieldManipulation
!
!***********************************************************************
!
  use types
  use parameters
  use gradients

  implicit none


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
!  Calculates -Div(p)                                      
!  ExplDiv(u) = sum_{i=1}^{i=nf} (u)_f*sf or
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry
  use variables, only: p,dPdxi
  use sparse_matrix, only: su,sv,sw

  implicit none
!
!***********************************************************************
!


!...Local
    integer :: i,ijp,ijn,ijb,iface,istage
    real(dp) :: dfxe,dfye,dfze

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
      call presFaceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        p, dPdxi,dfxe,dfye,dfze)

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
    enddo

  return
  end


!***********************************************************************
!
  function explDiv(u) result(div)
!
!***********************************************************************
!  
!  Calculates explicit divergence of scalar field phi - Div(phi)                                       
!  explDiv(u) = sum_{i=1}^{i=nf} (u)_f*sf or
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
    real(dp), dimension(numCells) :: div

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: dfxe
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp)+dfxe
      div(ijn) = div(ijn)-dfxe

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), dfxe)
      div(ijp) = div(ijp) + dfxe
    enddo

  return
  end



!***********************************************************************
!
  function explDivMdot(flmass,u) result(div)
!
!***********************************************************************
!  
!  Calculates explicit divergence of scalar field mdot*phi - Div(mdot,phi)                                       
!  explDiv(u) = sum_{i=1}^{i=nf} mdot*(u)_f*sf or
!  Interpolation to cell face centers done by cds corrected scheme.
!  flmass - is surface field, no need for intepolation.
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
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: dfxe
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp) + flmass(i) * dfxe

      div(ijn) = div(ijn) - flmass(i) * dfxe

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), dfxe)
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
!     This routine calculates contribution to the explicit Divergence
!     of a scalar FI (pressure) arising from an inner cell face.
!
!***********************************************************************
!
    use types
    use parameters
    use geometry

    implicit none

    integer,    intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numCells), intent(in) :: df


    real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
    real(dp) :: fie,dfxe,dfye,dfze
    real(dp) :: fxn,fxp
!
!***********************************************************************
!
    fxn = fif 
    fxp = 1.0d0-fxn

    xi = xc(ijp)*fxp+xc(ijn)*fxn
    yi = yc(ijp)*fxp+yc(ijn)*fxn
    zi = zc(ijp)*fxp+zc(ijn)*fxn

    dfxi = df(ijp,1)*fxp+df(ijn,1)*fxn
    dfyi = df(ijp,2)*fxp+df(ijn,2)*fxn
    dfzi = df(ijp,3)*fxp+df(ijn,3)*fxn

    ! Value of the variable at cell-face center
    fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

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
                          fi,df,dfxe)
!
!***********************************************************************
! 
!     This routine calculates contribution to the explicit Divergence
!     of a scalar FI (pressure) arising from an inner cell face.
!
!***********************************************************************
!
    use types
    use parameters
    use geometry

    implicit none

    integer,  intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numCells), intent(in) :: df


    real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
    real(dp) :: fie,dfxe
    real(dp) :: fxn,fxp
    real(dp) :: are
!
!***********************************************************************
!
    fxn = fif 
    fxp = 1.0d0-fxn

    xi = xc(ijp)*fxp+xc(ijn)*fxn
    yi = yc(ijp)*fxp+yc(ijn)*fxn
    zi = zc(ijp)*fxp+zc(ijn)*fxn

    dfxi = df(ijp,1)*fxp+df(ijn,1)*fxn
    dfyi = df(ijp,2)*fxp+df(ijn,2)*fxn
    dfzi = df(ijp,3)*fxp+df(ijn,3)*fxn

    ! Value of the variable at cell-face center
    fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

    ! Face area
    are = sqrt(sx**2+sy**2+sz**2)

    ! (interpolated mid-face value)x(area)
    dfxe = fie*are

  end subroutine


!***********************************************************************
!
  subroutine faceDivBoundary(sx,sy,sz,fi,dfx)
!
!***********************************************************************
!
!  This routine calculates the contribution of a boundary cell face to 
!  explicit Divergence of scalar FI.
!
!***********************************************************************
!

    implicit none

    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fi
    real(dp), intent(out)  :: dfx
!
!***********************************************************************
!
    dfx = fi*sqrt(sx**2+sy**2+sz**2)

  end subroutine

!***********************************************************************
!
  subroutine presFaceDivBoundary(sx,sy,sz,fi,dfx,dfy,dfz)
!
!***********************************************************************
!
!  This routine calculates the contribution of a boundary cell face to 
!  explicit Divergence of scalar FI.
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
    integer :: i,ijp,ijn,iface
    real(dp) :: dfxe
    real(dp) :: are
    real(dp), dimension(3,numCells) :: dUdxi
    real(dp), dimension(numCells) :: sumsf

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      aver(ijp) = aver(ijp)+dfxe

      aver(ijn) = aver(ijn)-dfxe

      are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
      sumsf(ijn) = sumsf(ijn) + are

    enddo

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)

      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
    enddo

    ! Divide by bounding faces total area
    do ijp = 1,numCells
      aver(ijp) = aver(ijp)/sumsf(ijp)
    enddo

  end function



!***********************************************************************
!
  function Interpolate(u) result(ui)
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
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, ui(i) )
    enddo


    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijb = numCells + i
      ui(iface) = u(ijb)
    enddo

  return
  end



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
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,if
    real(dp) :: ui
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, ui )

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
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, ui )

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


!
!***********************************************************************
!
  subroutine faceInterpolateCdsCorr(ijp,ijn, &
                                    xfc,yfc,zfc,fif, &
                                    fi,df,fie)
!
!***********************************************************************
! 
! This routine calculates face interpolant using CDS corrected.
!
!***********************************************************************
!
    use types
    use parameters
    use geometry

    implicit none

    integer,    intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numCells), intent(in) :: df


    real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
    real(dp) :: fie
    real(dp) :: fxn,fxp
!
!***********************************************************************
!
    fxn = fif 
    fxp = 1.0d0-fxn

    xi = xc(ijp)*fxp+xc(ijn)*fxn
    yi = yc(ijp)*fxp+yc(ijn)*fxn
    zi = zc(ijp)*fxp+zc(ijn)*fxn

    dfxi = df(ijp,1)*fxp+df(ijn,1)*fxn
    dfyi = df(ijp,2)*fxp+df(ijn,2)*fxn
    dfzi = df(ijp,3)*fxp+df(ijn,3)*fxn

    ! Value of the variable at cell-face center
    fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

  end subroutine



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
