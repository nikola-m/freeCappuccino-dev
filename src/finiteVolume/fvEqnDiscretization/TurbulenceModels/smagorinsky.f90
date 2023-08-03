module smagorinsky
!
! NOT FINISHED! Implementation of Smagorinsky SGS model
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

  real(dp), parameter :: C_Les = 0.17_dp ! 0.1-0.2; for channel flow even 0.065(Ferziger)   
  real(dp), parameter :: C_delta = 0.158_dp
  real(dp), parameter :: r13 = 1./3.0_dp
  real(dp), parameter :: Aplus = 26.0_dp

  private 

  public :: modify_viscosity_smagorinsky


contains


subroutine modify_viscosity_smagorinsky
!
! Update sub-grid scale (SGS) and effective viscosity.
!
  implicit none

  integer :: i,ib,inp
  integer :: iface, iwall, ijp,ijb
  real(dp) :: visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Utau,viscw
  real(dp) :: Vis_Smag,yplus,Cs

!
! Loop trough cells 
!

  do inp=1,numCells

    ! Store old value
    visold = vis(inp)


    ! Nondimensional distance
    yplus = Den(inp)*wallDistance(inp)*Utau/Viscos

    ! C_les coef with van Driest damping
    ! Cs = C_Les*(1.-Exp(-yplus/26.0))

    ! Smagorinsky SGS viscosity
    ! Vis_Smag = Den(Inp) * (Cs * Vol(Inp)**r13)**2 * magStrain(inp)


    ! Version implemented in OpenFOAM - merged ideas of Schumann(1991) and Van Driest:
    damping = (1.0_dp-Exp(-yplus/Aplus))
    Vis_Smag = Den(Inp) * ( min(C_Les*Vol(Inp)**r13, cappa*wallDistance(inp)/C_delta * damping)  )**2 * magStrain(inp)
               

    ! Update effective viscosity:
    Vis(inp) = Vis_Smag + Viscos

    ! Underelaxation
    vis(inp)= visold+urf(ivis)*(vis(inp)-visold) 

  enddo

  ! !
  ! ! Boundary faces 
  ! !

  ! iWall = 0

  ! do ib=1,numBoundaries

  !   if ( bctype(ib) == 'inlet' .or. bctype(ib) == 'outlet' .or. bctype(ib) == 'symmetry' ) then

  !     do i=1,nfaces(ib)

  !       iface = startFace(ib) + i
  !       ijp = owner(iface)
  !       ijb = iBndValueStart(ib) + i

  !       Vis(ijb) = Vis(ijp)

  !     end do

  !   elseif ( bctype(ib) == 'wall') then
  !   !
  !   ! Here we calculate modified viscosity stored for wall face, which will
  !   ! be used in expression for wall shear in momentum equation - calcuvw.
  !   !

  !     do i=1,nfaces(ib)

  !       iface = startFace(ib) + i
  !       ijp = owner(iface)
  !       ijb = iBndValueStart(ib) + i
  !       iWall = iWall + 1

  !       !
  !       ! Wall boundaries - update Visw and Ypl
  !       !

  !       ! Face area 
  !       are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

  !       ! Face normals
  !       nxf = arx(iface)/are
  !       nyf = ary(iface)/are
  !       nzf = arz(iface)/are

  !       ! Magnitude of a cell center velocity projected on boundary face normal
  !       Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

  !       ! Tangential velocity components 
  !       xtp = U(ijp)-Vnp*nxf
  !       ytp = V(ijp)-Vnp*nyf
  !       ztp = W(ijp)-Vnp*nzf

  !       ! Its magnitude
  !       Vtp = xtp*xtp+ytp*ytp+ztp*ztp

  !       ! Tangent direction
  !       xtp = xtp/vtp
  !       ytp = ytp/vtp
  !       ztp = ztp/vtp

  !       ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
  !       Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

  !       Tau(iWall) = viscos*Ut2/dnw(iWall)
  !       Utau = sqrt( Tau(iWall) / den(ijb) )
  !       ypl(iWall) = den(ijb)*Utau*dnw(iWall)/viscos

  !       ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
  !       ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
  !       ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

  !       viscw = zero

  !       ! if(ypl(iWall) > ctrans) then
  !       !   viscw = ypl(iWall)*viscos*cappa/log(Elog*ypl(iWall))
  !       ! endif

  !       visw(iWall) = max(viscos,viscw)
  !       vis(ijb) = visw(iWall)

  !     enddo

  !   endif 

  ! enddo


end subroutine


end module