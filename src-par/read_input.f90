!***********************************************************************
!
subroutine read_input_file
!
!***********************************************************************
!
! Open & Read and Process Input File
!
!***********************************************************************
  use types
  use parameters
  use gradients, only: lstsq, lstsq_qr, lstsq_dm, gauss, limiter
  use title_mod
  use utils
  use mpi

  implicit none

  integer :: i,imon
  character(len=2) :: trpn
  character(len=25) :: convective_scheme
  character(len=5) :: nproc_char
!
!***********************************************************************
!

! Root processor opens files
  if (myid .eq. 0) then

  open(unit=5,file=input_file)
  rewind 5

  read(5,'(a70)') title 
  read(5,*) lread,lwrite,ltest,lreadstat
  read(5,*) (lcal(i),i=1,nphi)
  read(5,*) monCell,pRefCell,MPoints
  read(5,*) slarge,sormax
  read(5,*) densit,viscos
  read(5,*) pranl,tref,beta
  read(5,*) lbuoy,gravx,gravy,gravz,boussinesq
  read(5,*) roughWall,EROUGH,ZZERO
  read(5,*) facnap,facflx
  read(5,*) ltransient,bdf,bdf2,bdf3,cn
  read(5,*) levm,lasm,lles,ldes
  read(5,*) lsgdh,lggdh,lafm
  read(5,*) TurbModel
  read(5,*) uin,vin,win,tein,edin,tin,vartin,conin
  read(5,*) convective_scheme
  read(5,*) limiter
  read(5,*) (gds(i),i=1,nphi)
  read(5,*) (urf(i),i=1,nphi)
  read(5,*) (sor(i),i=1,nphi)
  read(5,*) (nsw(i),i=1,nphi)
  read(5,*) numstep,timestep,nzapis,maxit
  read(5,*) lstsq, lstsq_qr, lstsq_dm, gauss
  read(5,*) npcor, nigrad
  read(5,*) simple,piso,ncorr
  read(5,*) const_mflux
  read(5,*) CoNumFix, CoNumFixValue

  close (5)

  ! Create an input file reading log:
  write(6,'(a)') '  Input file log: '
  write(6,'(a)') '---cut here-----------------------------------------------------------------------------'
  write(6,'(a70)') title
  write(6,'(4(l1,1x),5x,a)') lread,lwrite,ltest,lreadstat,'lread,lwrit,ltest,lreadstat'
  write(6,'(10(l1,1x),5x,a)') (lcal(i),i=1,nphi),'(lcal(i),i=1,nphi),ip=4,ite=5,ied=6,ien=7,ivis=8,ivart=9,icon=10'
  write(6,'(3(i3,1x),5x,a)') monCell,pRefCell,MPoints,'monCell,pRefCell,MPoints'
  write(6,'(2(es11.4,1x),5x,a)') slarge,sormax,'slarge,sormax'
  write(6,'(2(es11.4,1x),a)') densit,viscos,'densit,viscos'
  write(6,'(3(es11.4,1x),a)') pranl,tref,beta,'pranl,tref,beta'
  write(6,'(l1,1x,3f6.2,1x,l1,1x,a)') lbuoy,gravx,gravy,gravz,boussinesq,'lbuoy,gravx,gravy,gravz,boussinesq'
  write(6,'(L1,1x,f5.2,1x,es11.4,1x,a)') roughWall,erough,zzero,'roughWall,erough,zzero'
  write(6,'(2(f4.2,1x),a)') facnap,facflx,'facnap,facflx'
  write(6,'(5(l1,1x),a)') ltransient,bdf,bdf2,bdf3,cn,'ltransient,bdf,bdf2,bdf3,cn'
  write(6,'(4(l1,1x),a)') levm,lasm,lles,ldes,'levm,lasm,lles,ldes'
  write(6,'(3(l1,1x),a)') lsgdh,lggdh,lafm,'lsgdh,lggdh,lafm'
  write(6,'(i2,1x,a)') TurbModel, 'Turbulence Model'
  write(6,'(8(es11.4,1x),a)') uin,vin,win,tein,edin,tin,vartin,conin,'uin,vin,win,tein,edin,tin,vartin,conin'
  write(6,'(a,a)') convective_scheme, 'Convective scheme'
  write(6,'(a,1x,a)') limiter, 'Gradient limiter'
  write(6,'(10(f4.2,1x),a)') (gds(i),i=1,nphi),'(gds(i),i=1,nphi)'
  write(6,'(10(f4.2,1x),a)') (urf(i),i=1,nphi),'(urf(i),i=1,nphi)'
  write(6,'(10(es9.2,1x),a)') (sor(i),i=1,nphi),'(sor(i),i=1,nphi)'
  write(6,'(10(i3,1x),a)') (nsw(i),i=1,nphi),'(nsw(i),i=1,nphi)'
  write(6,'(i5,1x,es9.2,1x,i5,1x,i4,1x,a)') numstep,timestep,nzapis,maxit,'numstep,timestep,nzapis,maxit'
  write(6,'(4(L1,1x),a)') lstsq, lstsq_qr, lstsq_dm, gauss,'lstsq, lstsq_qr, lstsq_dm, gauss'
  write(6,'(i1,1x,i1,1x,a)') npcor, nigrad,'npcor, nigrad'
  write(6,'(2(l1,1x),i1,1x,a)') simple,piso,ncorr,'simple,piso,ncorr'
  write(6,'(1(L1,1x),5x,a)') const_mflux,'const_mflux'
  write(6,'(L1,es11.4,5x,a)') CoNumFix, CoNumFixValue,'CoNumFix, CoNumFixValue'
  write(6,'(a)') '---cut here-----------------------------------------------------------------------------'

  endif ! if(myid==0): end

! Broadcast input data to other processes
  call MPI_BCAST(title,70,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lread,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lwrite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ltest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lreadstat,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lcal,NPHI,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(monCell,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(pRefCell,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(MPOINTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  ! Treba naci kom procesu pripadaju monitoring tacke...

  call MPI_BCAST(slarge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(sormax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(densit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(viscos,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(pranl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lbuoy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(boussinesq,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(roughWall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(erough,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(zzero,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(facnap,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(facflx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(ltransient,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(bdf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(bdf2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)  
  call MPI_BCAST(bdf3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)   
  call MPI_BCAST(cn,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(levm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lasm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lles,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ldes,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lsgdh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lggdh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lafm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(turbmodel,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR) 

  call MPI_BCAST(uin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(vin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(win,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tein,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(edin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(vartin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(conin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(convective_scheme,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(limiter,20,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(GDS(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(URF(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(SOR(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(NSW(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(numstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(timestep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(nzapis,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(maxit,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lstsq,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lstsq_qr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lstsq_dm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gauss,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(npcor,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(nigrad,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(simple,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(piso,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ncorr,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(const_mflux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(conumfix,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(conumfixvalue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  ! Izgubi mu se pojam koji je process rank pa moram ovo da pozovem:
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr  ) 
  
  if (myid .eq. 0) then
    write(6,'(a)')' '
    write(6,'(a)')'  **Finished reading and broadcasting input data.**'
    write(6,'(a)')' '
  endif  


  !
  ! Turbulent flow computation condition:
  !
  lturb = levm.or.lasm.or.lles.or.ldes

  !
  ! Convective scheme:
  !
  select case(convective_scheme)

    case ('central')
      lcds = .true.
    case ('cds-corrected')
      lcdsc = .true.
    case ('linear')
      lluds = .true.
    case ('smart')
      lsmart = .true.
    case ('avl-smart')
      lavl = .true.
    case ('muscl')
      lmuscl = .true.
    case ('umist')
      lumist = .true.
    case ('koren')
      lkoren = .true.
    case ('charm')
      lcharm = .true.
    case ('ospre')
      lospre = .true.
    case ('central-f')
      lcds_flnt = .true.
    case ('linear-f')
      l2nd_flnt = .true.
    case('muscl-f')
      lmuscl_flnt = .true.
    case('cubic')
      lcubic = .true.      
    case default
      if (myid .eq. 0) write(*,'(a)') '  Using default convective scheme - 2nd order upwind.'
      l2nd_flnt = .true.

  end select

  ! Set value for flux_limiter logical
  if(lluds.or.lsmart.or.lavl.or.lmuscl.or.lumist.or.lkoren.or.lcharm.or.lospre) then
    flux_limiter = .true.
  else
    flux_limiter = .false.
  endif

  if (myid .eq. 0) then
    write(*,'(a)') ' '
    write(*,'(2a)') '  Convective scheme: ', adjustl(convective_scheme)
    write(*,'(a)') ' '
  
    !
    ! Gradient limiter:
    !
    if(adjustl(limiter) == 'Barth-Jespersen') then

      write(*,'(a)') '  Gradient limiter: Barth-Jespersen'

    elseif(adjustl(limiter) == 'Venkatakrishnan') then

      write(*,'(a)') '  Gradient limiter: Venkatakrishnan'

    elseif(adjustl(limiter) == 'mVenkatakrishnan') then

      write(*,'(a)') '  Gradient limiter: Wang modified Venkatakrishnan'

    elseif(adjustl(limiter) == 'MDL') then

      write(*,'(a)') '  Gradient limiter: Multidimensional'

    else !if(adjustl(limiter) == 'no-limit') then

      write(*,'(a)') '  Gradient limiter: no-limit'
      
    endif

    write(*,'(a)') ' '

    !
    ! Time stepping algorithm:
    !
    if( bdf ) then
      write(*,'(a)') '  Time stepping method: Euler Implicit'

    elseif(bdf2) then
        write(*,'(a)') '  Backward-Differentiation of 2nd order - BDF2'

    elseif(bdf2) then
        write(*,'(a)') '  Backward-Differentiation of 3nd order - BDF3'

    elseif( cn ) then
         write(*,'(a)') '  Time stepping method: Crank-Nicolson'   

    endif

  endif


  !
  ! Open files for data at monitoring points 
  !

  if( ltransient .and. mpoints>0 ) then

    ! nproc_char <- myid zapisan levo u vidu stringa.
    call i4_to_s_left ( myid, nproc_char )
    open(unit=89,file='processor'//trim(nproc_char)//'/monitoring_points')
    rewind 89

    read(89,*) mpoints

    do i=1,mpoints

      read(89,*) imon

      write(trpn,'(i0)') imon
      open(91+imon,file="transient_monitor_point_"//adjustl(trim(trpn)), access='append')
      if(.not.lreadstat) rewind(91+imon)

    end do

  end if


end subroutine