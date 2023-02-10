subroutine time_shift
  use types
  use parameters
  use geometry
  use variables
  use concentration, only: calcCon
  use temperature, only: calcT
  use energy, only: calcEn

  implicit none

  ! Transfer data from level n-2 to level n-3 for BDF3 time-stepping:
  if( bdf3 ) then
    uooo = uoo 
    vooo = voo 
    wooo = woo 
    flmassooo = flmassoo
    pooo = poo
  endif

  ! Transfer data from level n-1 to level n-2 for BDF3 or BDF2 time-stepping:
  if( bdf3 .or. bdf2 ) then
    uoo = uo 
    voo = vo 
    woo = wo 
    teoo = teo 
    edoo = edo
    if (calcT .or. calcEn ) Too = To 
    if ( calcEn ) Enoo = Eno   
    if ( calcCon ) Conoo = Cono 
    ! For consistent mass flux 
    ! if (piso) then 
    !   flmassoo = flmasso
    !   poo = po
    ! endif
  endif


  ! In all other cases (Crank-Nicolson, BDF) exclusively, 
  ! but also in the case of BDF2 and BDF3 as final shift, 
  ! transfer data from level n to level n-1:
  uo = u 
  vo = v 
  wo = w 
  teo = te 
  edo = ed 

  if ( CN ) po = p ! Only for Crank-Nicolson
  
  ! For buoyancy when not using Boussines and for pressure-based compressible solver, 
  if( (lbuoy .and. .not.boussinesq) .or. compressible ) deno = den  

  if ( calcT .or. calcEn ) To = T    
  if ( calcEn ) Eno = En      
  if ( calcCon ) Cono = Con 

  ! For consistent mass flux 
  ! if (piso) then
  !   flmasso = flmass
  !   po = p
  ! endif


end subroutine