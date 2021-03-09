subroutine time_shift
  use types
  use parameters
  use geometry
  use variables
  use concentration, only: calcCon
  use temperature, only: calcT
  use energy, only: calcEn

  implicit none

  if( bdf .or. cn ) then

    uo = u 
    vo = v 
    wo = w 

    teo = te 
    edo = ed 

    if ( calcT .or. calcEn ) To = T         
    if ( calcCon ) Cono = Con 

    if (CN) po = p

    ! For consistent mass flux 
    ! if (piso) then
    !   flmasso = flmass
    !   po = p
    ! endif

    deno = den
  
  elseif( bdf2 ) then

    uoo = uo 
    voo = vo 
    woo = wo 

    teoo = teo 
    edoo = edo

    if (calcT .or. calcEn ) Too = To 
    if ( calcCon ) Conoo = Cono 

    ! For consistent mass flux 
    ! if (piso) then 
    !   flmassoo = flmasso
    !   poo = po
    ! endif

    uo = u 
    vo = v 
    wo = w 

    teo = te 
    edo = ed 

    if (calcT .or. calcEn ) to = t         
    if ( calcCon ) Cono = con 
    
    ! For consistent mass flux     
    ! if (piso) then 
    !   flmasso = flmass
    !   po = p
    ! endif

  elseif( bdf3 ) then

    uooo = uoo 
    vooo = voo 
    wooo = woo 

    flmassooo = flmassoo
    pooo = poo

    uoo = uo 
    voo = vo 
    woo = wo 

    teoo = teo 
    edoo = edo
    if (calcT .or. calcEn ) Too = To 
    if ( calcCon ) Conoo = Cono 

    ! For consistent mass flux 
    ! if (piso) then 
    !   flmassoo = flmasso
    !   poo = po
    ! endif

    uo = u 
    vo = v 
    wo = w 
    teo = te 
    edo = ed 
    if (calcT .or. calcEn ) To = T          
    if ( calcCon ) Cono = con 

    ! For consistent mass flux 
    ! if (piso) then 
    !   flmasso = flmass
    !   po = p
    ! endif

  endif


end subroutine