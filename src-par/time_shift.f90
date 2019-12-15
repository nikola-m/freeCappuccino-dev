subroutine time_shift
  use types
  use parameters
  use geometry
  use variables

  implicit none

  if( (bdf .and. itime == 1) .or. (bdf2 .and. itime == 1) .or. (bdf3 .and. itime == 1) ) then

    uo = u 
    vo = v 
    wo = w 
    teo = te 
    edo = ed 
    if (lcal(ien)) to = t         
    if (lcal(ivart)) varto = vart 
    if (lcal(icon)) cono = con 
    if (piso) then 
      flmasso = flmass
      po = p
    endif
  
  elseif( (bdf .and. itime > 1) .or. (bdf2 .and. itime > 1) .or. (bdf3 .and. itime == 2) ) then

    uoo = uo 
    voo = vo 
    woo = wo 
    teoo = teo 
    edoo = edo
    if (lcal(ien)) too = to 
    if (lcal(ivart)) vartoo = varto 
    if (lcal(icon)) conoo = cono 
    if (piso) then 
      flmassoo = flmasso
      poo = po
    endif

    uo = u 
    vo = v 
    wo = w 
    teo = te 
    edo = ed 
    if (lcal(ien)) to = t         
    if (lcal(ivart)) varto = vart 
    if (lcal(icon)) cono = con 
    if (piso) then 
      flmasso = flmass
      po = p
    endif

  elseif( bdf3 .and. itime > 2) then

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
    if (lcal(ien)) too = to 
    if (lcal(ivart)) vartoo = varto 
    if (lcal(icon)) conoo = cono 
    if (piso) then 
      flmassoo = flmasso
      poo = po
    endif

    uo = u 
    vo = v 
    wo = w 
    teo = te 
    edo = ed 
    if (lcal(ien)) to = t         
    if (lcal(ivart)) varto = vart 
    if (lcal(icon)) cono = con 
    if (piso) then 
      flmasso = flmass
      po = p
    endif

  endif


end subroutine