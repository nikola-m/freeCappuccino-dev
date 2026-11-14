    module fieldManipulation
!
!***********************************************************************
!
  use types
  use gradients
  use interpolation, only: face_value

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
  character(len=30), parameter :: scheme = 'central' 

  interface explDiv
    module procedure explDiv
    module procedure explDivFlmass
  end interface

  ! interface explDdt
  !   module procedure explDdt
  !   module procedure explDdtRho
  ! end interface

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
function explDivFlmass(flmass,u,v,w) result(div)
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
subroutine faceDivInner(ijp,ijn,xfc,yfc,zfc,sx,sy,sz,fif, &
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
  uf = face_value( ijp, ijn, xfc, yfc, zfc, fif, u, du, scheme )
  vf = face_value( ijp, ijn, xfc, yfc, zfc, fif, v, dv, scheme )
  wf = face_value( ijp, ijn, xfc, yfc, zfc, fif, w, dw, scheme )

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

      ui = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )

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
function fieldInterpolate(u) result(ui)
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

  ! Output
  real(dp), dimension(numFaces) :: ui

  ! Input
  real(dp), dimension(numTotal), intent(in) :: u

  ! Local
  integer :: i,ijp,ijn,ijb,iface
  real(dp), dimension(:,:), allocatable :: dUdxi

  allocate( dUdxi(3,numTotal) )

  ! Calculate cell-centered gradient
  call grad(u,dUdxi)

  ! Calculate terms integrated over surfaces

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    ui(i) = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )
  enddo

  ! Update boundaries?

  ! Contribution from boundaries
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijb = numCells + i
    ui(iface) = u(ijb)
  enddo

  deallocate(dUdxi)

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

      ui = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )

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

      ui = face_value( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, scheme )

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



! !***********************************************************************
! !
! function explDdtRho(den,u,uo,uoo) result(ddt)
! !
! !***********************************************************************
! !  
! !  volScalarField -> volScalarField
! !  Explicit time derivative in finite volume form.
! !
! !***********************************************************************
! !
!   use parameters
!   use geometry, only: vol, numCells, numTotal

!   implicit none
! !
! !***********************************************************************
! !

!   ! Output
!   real(dp), dimension(numCells) :: ddt

!   ! Input
!   real(dp), dimension(numTotal), intent(in) :: den
!   real(dp), dimension(numTotal), intent(in) :: u,uo
!   real(dp), dimension(numTotal), intent(in), optional :: uoo
!   ! Local
!   integer :: inp

!   ddt = 0.0_dp

!   !
!   ! > Loop over cells
!   !
!   do inp=1,numCells

!     if( bdf .or. cn ) then
 
!       ddt(inp) = den(inp)*vol(inp) * ( u(inp) - uo(inp) )

!     elseif( bdf2 ) then

!       ddt(inp) = den(inp)*vol(inp) * ( 1.5_dp*u(inp) - 2*uo(inp) + 0.5_dp*uoo(inp) )

!     endif

!   enddo

! end function


! !***********************************************************************
! !
! function explDdt(u,uo,uoo) result(ddt)
! !
! !***********************************************************************
! !  
! !  volScalarField -> volScalarField
! !  Explicit time derivative in finite volume form.
! !
! !***********************************************************************
! !
!   use parameters
!   use geometry, only: vol, numCells, numTotal

!   implicit none
! !
! !***********************************************************************
! !

!   ! Output
!   real(dp), dimension(numCells) :: ddt

!   ! Input
!   real(dp), dimension(numTotal), intent(in) :: u,uo
!   real(dp), dimension(numTotal), intent(in), optional :: uoo

!   ! Local
!   integer :: inp

!   ddt = 0.0_dp

!   !
!   ! > Loop over cells
!   !
!   do inp=1,numCells

!     if( bdf .or. cn ) then
 
!       ddt(inp) = vol(inp) * ( u(inp) - uo(inp) )

!     elseif( bdf2 ) then

!       ddt(inp) = vol(inp) * ( 1.5_dp*u(inp) - 2*uo(inp) + 0.5_dp*uoo(inp) )

!     endif

!   enddo

! end function

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
     
  CALL init_random_seed()

  do inp=1,numCells
    
    CALL RANDOM_NUMBER(perturb)
    ! perturb is now between 0. and 1., we want it to be from 0 to 2*amplitude
    ! e.g. perturb = 0.9+perturb/5. when Max perturbation is +/- 10% of mean profile
    perturb = ( 1.0_dp - level/100.0_dp ) + perturb * (2*level/100.0_dp)

    Phi(INP) = perturb*Phi(inp)

  enddo

end subroutine


subroutine update_inlet_mapping
! This was used to recirculate to inlet a slice from the inside of the domain
! so called mapping - we map back to inlet.
! Hard coded here so it is an ad hoc thing used for Jovic-Driver BFS case.
  use types
  ! use parameters, only: densit
  use geometry, only: numBoundaries,arx,ary,arz,startFace,iBndValueStart,nFaces,bctype
  use variables, only: u,v,w,flmass,den
  implicit none 

  integer :: i, ib, iface, ijn
  real(dp) :: flomas !, are, outare, flomas, uav

  flomas = 0.

  ! Initialize mass flow over inlet faces
  do ib=1,numBoundaries  
    if ( bctype(ib) == 'inlet' ) then

    ! Veca mreza
    u(startFace(ib)+1)=u(109)
    v(startFace(ib)+1)=v(109)
    u(startFace(ib)+2)=u(301)
    v(startFace(ib)+2)=v(301)
    u(startFace(ib)+3)=u(493)
    v(startFace(ib)+3)=v(493)
    u(startFace(ib)+4)=u(685)
    v(startFace(ib)+4)=v(685)
    u(startFace(ib)+5)=u(877)
    v(startFace(ib)+5)=v(877)
    u(startFace(ib)+6)=u(1069)
    v(startFace(ib)+6)=v(1069)
    u(startFace(ib)+7)=u(1261)
    v(startFace(ib)+7)=v(1261)
    u(startFace(ib)+8)=u(1453)
    v(startFace(ib)+8)=v(1453)
    u(startFace(ib)+9)=u(1645)
    v(startFace(ib)+9)=v(1645)
    u(startFace(ib)+10)=u(1837)
    v(startFace(ib)+10)=v(1837)
    u(startFace(ib)+11)=u(2029)
    v(startFace(ib)+11)=v(2029)
    u(startFace(ib)+12)=u(2221)
    v(startFace(ib)+12)=v(2221)
    u(startFace(ib)+13)=u(2413)
    v(startFace(ib)+13)=v(2413)
    u(startFace(ib)+14)=u(2605)
    v(startFace(ib)+14)=v(2605)
    u(startFace(ib)+15)=u(2797)
    v(startFace(ib)+15)=v(2797)
    u(startFace(ib)+16)=u(2989)
    v(startFace(ib)+16)=v(2989)
    u(startFace(ib)+17)=u(3181)
    v(startFace(ib)+17)=v(3181)
    u(startFace(ib)+18)=u(3373)
    v(startFace(ib)+18)=v(3373)
    u(startFace(ib)+19)=u(3565)
    v(startFace(ib)+19)=v(3565)
    u(startFace(ib)+20)=u(3757)
    v(startFace(ib)+20)=v(3757)
    u(startFace(ib)+21)=u(3949)
    v(startFace(ib)+21)=v(3949)
    u(startFace(ib)+22)=u(4141)
    v(startFace(ib)+22)=v(4141)
    u(startFace(ib)+23)=u(4333)
    v(startFace(ib)+23)=v(4333)
    u(startFace(ib)+24)=u(4525)
    v(startFace(ib)+24)=v(4525)
    u(startFace(ib)+25)=u(4717)
    v(startFace(ib)+25)=v(4717)
    u(startFace(ib)+26)=u(4909)
    v(startFace(ib)+26)=v(4909)
    u(startFace(ib)+27)=u(5101)
    v(startFace(ib)+27)=v(5101)
    u(startFace(ib)+28)=u(5293)
    v(startFace(ib)+28)=v(5293)
    u(startFace(ib)+29)=u(5485)
    v(startFace(ib)+29)=v(5485)
    u(startFace(ib)+30)=u(5677)
    v(startFace(ib)+30)=v(5677)
    u(startFace(ib)+31)=u(5869)
    v(startFace(ib)+31)=v(5869)
    u(startFace(ib)+32)=u(6061)
    v(startFace(ib)+32)=v(6061)
    u(startFace(ib)+33)=u(6253)
    v(startFace(ib)+33)=v(6253)
    u(startFace(ib)+34)=u(6445)
    v(startFace(ib)+34)=v(6445)
    u(startFace(ib)+35)=u(6637)
    v(startFace(ib)+35)=v(6637)
    u(startFace(ib)+36)=u(6829)
    v(startFace(ib)+36)=v(6829)
    u(startFace(ib)+37)=u(7021)
    v(startFace(ib)+37)=v(7021)
    u(startFace(ib)+38)=u(7213)
    v(startFace(ib)+38)=v(7213)
    u(startFace(ib)+39)=u(7405)
    v(startFace(ib)+39)=v(7405)

      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijn = iBndValueStart(ib) + i

        flmass(iface) = den(ijn)*(arx(iface)*u(ijn)+ary(iface)*v(ijn)+arz(iface)*w(ijn))

        flomas = flomas - flmass(iface)         
      end do
    endif 
  enddo


end subroutine
    
end module fieldManipulation
