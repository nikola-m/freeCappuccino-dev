module rheology

!********************************************************************************************************************************132 
!
!  Purpose: 
!
!    Constitutive relations for Non-Newtonian models.
!
!  Discussion:
!
!    Here we implement following models:
!    1) Power law constitutive Model (Ostwald)
!    2) Herschel-Bulkley constitutive regularized (Papanastasiou) Model
!    3) Bingham-Papanastasiou model
!    4) Carreau constitutive model
!    5) Casson model
!    6) Casson constitutive regularized (Papanastasiou) Model with consistency temperature dependance
!    7) Cross constitutive model
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
  use types
  use tensorFields
  use fvxGradient
  use variables, only: u,v,w,vis,visw
  use geometry, only: numCells,startFace,owner,iBndValueStart,numBoundaries,bctype,nfaces


  implicit none

  logical :: calcVis = .false. ! Logical switch - do we use this module

  ! Rheological model.
  ! Available options are listed below under 'case' switching statement in the
  ! 'modifyViscosityNonNewtonianFluid' routine.
  character( len=20 ) :: non_newtonian_model = 'PowerLaw' 

  ! n  : Power law index
  ! n < 1. ===> Shear thinning
  ! n > 1. ===> Shear thickening
  ! n = 1. and Tau_0  = 0 ===> Newtonian fluid.
  real(dp) :: npow = 0.4_dp    

  ! Consst : Consistency index
  real(dp) :: Consst = 10.0_dp    

  ! Lower limit for the shear rate magnitude 
  real(dp) :: shearmin = 0.01_dp 

  ! Minimum viscosity limit [kg/m-s]     
  real(dp) :: mumin = 0.01_dp 

  ! Maximum viscosity limit [kg/m-s]     
  real(dp) :: mumax = 1000.0_dp  

  ! Tau_0  : Yield stress    
  real(dp) :: Tau_0 = 0.0   

  ! Critical Shear Rate [1/s] 
  ! shearcr = 0.005 is a good value if yield stress is 18-60 Pa, 
  ! according to J. Gao, A. Fourie, Minerals Engineering 71 (2015) 120–132.
  real(dp) :: shearcr = 0.005

  ! megp : exponential growth parameter in Papanastasious regularization        
  real(dp) :: megp = 1000.0   

  ! Plastic viscosity [kg/m-s]  
  real(dp) :: muplastic = 0.0_dp 

  ! Zero shear rate viscosity [kg/m-s]
  real(dp) :: muzero =  300.0_dp  

  ! Infinity shear rate viscosity [kg/m-s]     
  real(dp) :: muinfty = 1000.0_dp     

  ! Time constant [s], also natural time 
  ! (i.e., inverse of the shear rate at which the fluid changes
  ! from Newtonian to power-law behavior)  
  real(dp) :: lamtime = 10.0_dp     

  ! Parameter governing the transition area in Carreau-Yasuda model (usually denoted 'a' in the literature)
  real(dp) :: acy = 2.0_dp

  ! Under-relaxation for recalculated viscosity  
  real(dp) :: urfVis = 1.0        

  real(dp), parameter :: sqrt2 = sqrt(2.0_dp)

  private 

  ! Params - change defaults through input.nml file
  public :: npow, Consst, shearmin, Tau_0, megp,  &
            muplastic, muzero, muinfty, lamtime, acy, &   
            mumin, mumax, shearcr, & 
            calcVis, urfVis, non_newtonian_model                           

  public :: modifyViscosityNonNewtonianFluid
  public :: YieldRegion

contains


subroutine modifyViscosityNonNewtonianFluid
!
!  Purpose: 
!
!    Driver subroutine for modification of viscosity using Non-Newtonian models.
!
!  Discussion:
!
!    Options are,
!    'PowerLaw'             : Power law constitutive Model (Ostwald)
!    'HerschelBulkley'      : Herschel-Bulkley constitutive regularized (Papanastasiou) Model
!    'Bingham'              : Bingham model
!    'BinghamPapanastasiou' : Bingham-Papanastasiou model
!    'Carreau'              : Carreau constitutive model
!    'CarreauYasuda'        : Carreau-Yasuda constitutive model
!    'CassonPapanastasiou'  : Casson constitutive regularized (Papanastasiou) Model with consistency temperature dependance
!    'CrossModel'           : Cross constitutive model
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!  Parameters:
!
!    Input, character, non_newtonian_model( 30 ), name of the non-newtonian model.
!

  implicit none


  select case( non_newtonian_model )

    case( 'PowerLaw' )

      call PowerLawModel

    case( 'HerschelBulkley' )

      call HerschelBulkleyModel

    case( 'Bingham' )

      call BinghamModel

    case( 'BinghamPapanastasiou' )

      call BinghamPapanastasiouModel

    case( 'Carreau' )

      call CarreauModel

    case( 'CarreauYasuda' )

      call CarreauYasudaModel

    case( 'Casson' )

      call CassonModel

    case( 'CassonPapanastasiou' )

      call CassonPapanastasiouModel

    case( 'CrossModel' )

      call CrossModel

  end select

end subroutine


subroutine PowerLawModel
!
!  Purpose: 
!
!    Power law constitutive Model (Ostwald).
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine 
!
  implicit none
!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
                    ( Consst * max(shear%mag(1:numCells), shearmin) ** (npow - 1.0_dp) ) &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - PowerLaw model."

end subroutine



subroutine HerschelBulkleyModel
!
!  Purpose: 
!
!    Herschel Bulkley Papanastasiou constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine, COMSOL theory, Fluent theory, etc.
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

 
  ! ! Now return to array based storage of variables
  ! vis(1:numCells) = urfVis * &
  !    min( mumin, ( Tau_0 + Consst * shear%mag(1:numCells) ** npow ) / ( shear%mag(1:numCells) + 1e-20 ) ) &
  !                  + (1.0_dp-urfVis)*vis(1:numCells)

  ! With Papansasstasiou regularization
  vis(1:numCells) = urfVis * &
     min( mumin, ( Tau_0*( 1.0_dp - exp( -megp*shear%mag(1:numCells) ) )  + &
                   Consst * shear%mag(1:numCells) ** npow ) / ( shear%mag(1:numCells) + 1e-30 ) ) &
                   + (1.0_dp-urfVis)*vis(1:numCells)
                 

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    else

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - Herschel-Bulkley model."

end subroutine



subroutine BinghamModel
!
!  Purpose: 
!
!    Bingham constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine, COMSOL theory, Fluent theory, etc.
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * ( muplastic + Tau_0/( shear%mag(1:numCells) + 1e-20 )  ) &
         + (1.0_dp-urfVis) * vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - Bingham model."

end subroutine



subroutine BinghamPapanastasiouModel
!
!  Purpose: 
!
!    Bingham Papanastasiou constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine, COMSOL theory, Fluent theory, etc.
!
!
!  Parameters for the model:
!    muplastic: Plastic viscosity [kg/m-s]
!    Tau_0  : Yield stress  
!    megp : exponential growth parameter in Papanastasious regularization    
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
   ( muplastic + Tau_0/( shear%mag(1:numCells) + 1e-20 ) * ( 1.0_dp - exp( -megp*shear%mag(1:numCells) ) ) ) &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - Bingham-Papanastasiou model."

end subroutine


subroutine CarreauModel
!
!  Purpose: 
!
!    Carreau constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Reference:
!
!    Non-Newtonian constitutive models implementation in COMSOL theory, Fluent theory, etc.
!
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
    ( muinfty + ( muzero - muinfty ) * ( 1.0_dp + ( lamtime * shear%mag(1:numCells) )**2 ) ** ( (npow-1.0_dp)/2.0_dp ) ) &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - Carreau model."

end subroutine

subroutine CarreauYasudaModel
!
!  Purpose: 
!
!    Carreau-Yasuda constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!
!  Reference:
!
!    Non-Newtonian constitutive models implementation in COMSOL theory, Fluent theory, etc.
!
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
    ( muinfty + ( muzero - muinfty ) * ( 1.0_dp + ( lamtime * shear%mag(1:numCells) )**acy ) ** ( (npow-1.0_dp)/acy ) ) &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) "  Updated effective viscosity: Carreau-Yasuda model."

end subroutine

subroutine CassonModel
!
!  Purpose: 
!
!    Casson constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine, COMSOL theory, Fluent theory, etc.
!
!  Parameters that need to be defined:
!    Tau_0  - Yield stress
!    muplastic - Plastic viscosity
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
    ( sqrt(muplastic) + sqrt( Tau_0/shear%mag(1:numCells) ) )**2 &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - Casson model."

end subroutine


subroutine CassonPapanastasiouModel
!
!  Purpose: 
!
!    Casson Papanastasiou constitutive model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine, COMSOL theory, Fluent theory, etc.
!
!  Parameters that need to be defined:
!    Tau_0  - Yield stress
!    megp - exponential growth parameter
!    muplastic - Plastic viscosity
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
    ( sqrt(muplastic) + sqrt(Tau_0/shear%mag(1:numCells)) * ( 1.0_dp - exp( -sqrt(megp*shear%mag(1:numCells)) ) ) )**2 &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo

  write(*,*) " Recalculated effective viscosity - Casson-Papanastasiou model."

end subroutine


subroutine CrossModel
!
!  Purpose: 
!
!    Cross constitutive Model.
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    09/12/2020
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Non-Newtonian constitutive models implementation by Seif-Eddine 
!

  implicit none

!
! Local variables
!
  integer :: i,iface,ijp,ijb,ib,iWall
  
  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: shear

  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  shear = sqrt2*(.mag.(.symm.( gradU ) ))

  ! Now return to array based storage of variables
  vis(1:numCells) = urfVis * &
                    muinfty + ( muzero - muinfty ) / (1.0_dp + ( lamtime * shear%mag(1:numCells) ) ** npow ) &
                   + (1.0_dp-urfVis)*vis(1:numCells)

  !
  ! Boundary faces 
  !
  iWall = 0

  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' .or. &
         bctype(ib) == 'outlet' .or. &
         bctype(ib) == 'symmetry' .or. &
         bctype(ib) == 'empty' .or. &
         bctype(ib) == 'pressure' ) then

      do i=1,nfaces(ib)

        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i

        Vis(ijb) = Vis(ijp)

      end do

    elseif ( bctype(ib) == 'wall') then
    !
    ! Here we calculate modified viscosity stored for wall face, which will
    ! be used in expression for wall shear in momentum equation - calcuvw.
    !

      do i=1,nfaces(ib)
        
        iface = startFace(ib) + i
        ijp = owner(iface)
        ijb = iBndValueStart(ib) + i
        iWall = iWall + 1

        visw(iWall) = max(vis(ijp),1e-20)
        vis(ijb) = visw(iWall)

      enddo

    endif 

  enddo
  
  write(*,*) " Recalculated effective viscosity - Cross model."

end subroutine


subroutine YieldRegion
!
!  Purpose: 
!
!    Calulates where did fluid yield, i.e. where did |\tau| > Tau_0
!
!  Discussion:
!
!    We use tensor field manipulation, not as efficient but very elegant.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license. 
!
!  Modified:
!
!    28/01/2021
!
!  Author:
!
!    Nikola Mirkov/Email: nikolamirkov@yahoo.com
!  
!
!  Reference:
!    
!    Papers by Mitsoulis.
!
  use utils, only: get_unit
  use output, only: vtu_write_scalar_field

  implicit none
!
! Local variables
!
  integer :: output_file
  integer :: i

  type(volVectorField) :: Uvec
  type(volTensorField) :: gradU
  type(volScalarField) :: tau


  Uvec = volVectorField( "VelocityVec",   &
                          U, &
                          V, &
                          W  &
                        )
  gradU = Grad( Uvec )
  tau = vis * sqrt2*(.mag.(.symm.( gradU ) ))

  ! Yielded region
  do i=1,numCells
    if ( tau%mag(i) >= Tau_0 ) then
      tau%mag(i) = 1.0
    else
      tau%mag(i) = 0.0
    endif
  enddo

  ! Write field to a Paraview .vtu file.
  call get_unit( output_file )
  open( unit = output_file, file='vtk/yield-region.vtu' )
  rewind output_file

  call vtu_write_scalar_field ( output_file, 'yield-region', tau%mag )

  close( output_file )

end subroutine


end module