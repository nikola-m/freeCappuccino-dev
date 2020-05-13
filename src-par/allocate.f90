!***********************************************************************
!
subroutine allocate_arrays
!
!***********************************************************************
!
  use parameters
  use geometry, only: numTotal,numCells,numFaces,nwal,npro
  use sparse_matrix, only: nnz
  use variables
  use hcoef
  use title_mod
  use statistics

  implicit none 
!
!***********************************************************************
!

  ! Variables

  ! Velocities 
  allocate( u( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: u" 

  allocate( v( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: v" 

  allocate( w( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: w" 

  allocate( uo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: uo" 

  allocate( vo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vo" 

  allocate( wo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: wo" 

  if( bdf2 .or. bdf3 ) then
    allocate( uoo( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uoo" 

    allocate( voo( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: voo" 

    allocate( woo( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: woo" 
  endif

  if( bdf3 ) then
    allocate( uooo( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uoo" 

    allocate( vooo( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: voo" 

    allocate( wooo( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: woo" 
  endif

  allocate( dUdxi( 3,numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dUdxi"

  allocate( dVdxi( 3,numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dVdxi"

  allocate( dWdxi( 3,numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dWdxi"

  if (ltransient) then
    allocate( u_aver( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: u_aver" 

    allocate( v_aver( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: v_aver" 

    allocate( w_aver( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: w_aver" 
  endif



  ! Pressure and pressure correction
  allocate( p( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: p"

  allocate(pp(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: pp"

  allocate( dPdxi( 3,numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dPdxi"


  if ( piso ) then
    
    allocate(po(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: po" 

    allocate(poo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: poo" 

    if (bdf3 ) then
      allocate(pooo(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: poo" 
    endif

  endif

  ! Turbulent K.E. and dissipation or any other turbulent scalar taking its place
  allocate( te( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: te" 

  allocate( ed( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ed" 

  allocate( teo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: teo" 

  allocate( edo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: edo" 

  if( bdf2 .or. bdf3 ) then

  allocate( teoo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: teoo" 

  allocate( edoo( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: edoo" 
  
  endif

  allocate( dTEdxi( 3,numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dTEdxi"

  allocate( dEDdxi( 3,numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dEDdxi"

  if (ltransient) then
    allocate( te_aver( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: te_aver" 
  endif


  ! Effective viscosity  
  allocate( vis( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vis" 



  ! Temperature
  if(lcal(ien)) then

    allocate( t( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: t"

    allocate( to( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: to"

    if( bdf2 .or. bdf3 ) then
            allocate( too( numTotal ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: too" 
    endif

    allocate( dTdxi( 3,numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dTdxi"

    if (ltransient) then
      allocate( t_aver( numTotal ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: t_aver" 
    endif

  endif

  ! Concentration
  if(lcal(icon)) then

    allocate( con( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: con"

    allocate( cono( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: cono"

    if( bdf .and. btime.gt.0.99 ) then
      allocate( conoo( numTotal ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: conoo" 
    endif

    allocate( dCondxi( 3,numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dCondxi"

    if (ltransient) then
      allocate( con_aver( numTotal ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: tt_aver" 
    endif

  endif

  ! Temperature variance and dissipation of temperature variance
  if(lcal(ivart)) then

    allocate( vart( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vart" 

    allocate( varto( numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: varto" 

    if( bdf2 .or. bdf3 ) then
      allocate( vartoo( numTotal ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: vartoo" 
    endif

    allocate( dVartdxi( 3,numTotal ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dVartdxi"

    if (ltransient) then
      allocate( tt_aver( numTotal ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: tt_aver" 
    endif

    ! allocate( edd( numTotal ), stat=ierr) 
    !   if(ierr /= 0)write(*,*)"allocation error: edd"

  endif


  ! Turbulent heat fluxes

  if(lcal(ien).and.lbuoy) then

    allocate( utt( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: utt" 

    allocate( vtt( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vtt" 

    allocate( wtt( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: wtt"

      if (ltransient) then

        allocate( ut_aver( numCells ), stat=ierr) 
          if(ierr /= 0)write(*,*)"allocation error: ut_aver" 

        allocate( vt_aver( numCells ), stat=ierr) 
          if(ierr /= 0)write(*,*)"allocation error: vt_aver" 

        allocate( wt_aver( numCells ), stat=ierr) 
          if(ierr /= 0)write(*,*)"allocation error: wt_aver"

      endif

  endif

  ! Density
  allocate( den( numTotal ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: den" 


  ! Mass flows trough cell faces
  allocate(flmass(numFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: flmass" 

  if (piso) then

    allocate(flmasso(numFaces),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: flmasso" 
      
    allocate(flmassoo(numFaces),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: flmassoo" 

    if (bdf3) then
      allocate(flmassooo(numFaces),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: flmassooo"
    endif

  endif


  ! Wall things
  allocate(  visw(nwal), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: visw"

  allocate(  ypl(nwal), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ypl"

  allocate(  tau(nwal), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: tau"

  if (ltransient) then
    allocate(  wss_aver(nwal), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: wss_aver"
  endif

  ! Turbulence production
  allocate( gen( numCells ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: gen"  

  allocate( magStrain( numCells ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: magStrain" 

  allocate( Vorticity( numCells ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: Vorticity" 



  ! Re_t used in low-re turbulent models
  ! allocate( ret( numCells ), stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ret" 


  ! Reynolds stresses
  if(lturb) then

    allocate( uu( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uu" 

    allocate( uv( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uv" 

    allocate( uw( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uw"

    allocate( vv( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vv"

    allocate( vw( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vw" 

    allocate( ww( numCells ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: ww" 

    if(lasm) then
      allocate( bij(5,numCells), stat=ierr)
        if(ierr/=0)write(*,*)"allocate error: bij"
     endif

    if(ltransient) then

      allocate( uu_aver( numCells ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: uu_aver"

      allocate( vv_aver( numCells ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: vv_aver"

      allocate( ww_aver( numCells ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: ww_aver" 

      allocate( uv_aver( numCells ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: uv_aver"

      allocate( uw_aver( numCells ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: uw_aver"

      allocate( vw_aver( numCells ), stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: vw_aver"


    endif

  endif 


  ! Coefficient arrays for PISO and PIMPLE
  if( piso ) then

    allocate( h( nnz ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: h" 

    allocate( hpr( npro ), stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: h" 

    allocate(rU(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: rU" 

    allocate(rV(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: rV" 

    allocate(rW(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: rW" 

  endif

  allocate( phimax( numCells ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: phimax" 

  allocate( phimin( numCells ), stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: phimin" 

end subroutine
