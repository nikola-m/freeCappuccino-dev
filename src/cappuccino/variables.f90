module variables
!%%%%%%%%%%%%%%
  use types


    real(dp), dimension(:), allocatable :: u,v,w              ! Velocity components     
    real(dp), dimension(:), allocatable :: uo, vo, wo         ! values from previous (n-1) timestep
    real(dp), dimension(:), allocatable :: uoo,voo,woo        ! values from n-2 time step 
    real(dp), dimension(:), allocatable :: uooo,vooo,wooo     ! values from n-3 time step
    real(dp),dimension(:,:), allocatable :: dUdxi,dVdxi,dWdxi ! Gradients


    real(dp), dimension(:), allocatable :: flmass             ! Mass fluxes trough inner faces
    real(dp), dimension(:), allocatable :: flmasso
    real(dp), dimension(:), allocatable :: flmassoo
    real(dp), dimension(:), allocatable :: flmassooo

    real(dp), dimension(:), allocatable :: p,pp               ! Pressure, Pressure correction,  
    real(dp), dimension(:), allocatable :: po
    real(dp), dimension(:), allocatable :: poo
    real(dp), dimension(:), allocatable :: pooo
    real(dp),dimension(:,:), allocatable :: dPdxi


    real(dp), dimension(:), allocatable :: te                 ! Turbulence kinetic energy, TKE Dissipation rate,
    real(dp), dimension(:), allocatable :: teo
    real(dp), dimension(:), allocatable :: teoo
    real(dp),dimension(:,:), allocatable :: dTEdxi

    real(dp), dimension(:), allocatable :: ed                 ! TKE Dissipation rate,
    real(dp), dimension(:), allocatable :: edo
    real(dp), dimension(:), allocatable :: edoo
    real(dp),dimension(:,:), allocatable :: dEDdxi

    real(dp), dimension(:), allocatable :: vis                ! Effective viscosity

    real(dp), dimension(:), allocatable :: den                ! Density
    real(dp), dimension(:), allocatable :: deno

    real(dp), dimension(:), allocatable :: gen                ! Turb. kin. energy generation (production)

    real(dp), dimension(:), allocatable :: T                  ! Temperature in K
    real(dp), dimension(:), allocatable :: to
    real(dp), dimension(:), allocatable :: too
    real(dp),dimension(:,:), allocatable :: dTdxi

    real(dp), dimension(:), allocatable :: con                ! Concentration
    real(dp), dimension(:), allocatable :: cono
    real(dp), dimension(:), allocatable :: conoo 
    real(dp),dimension(:,:), allocatable :: dCondxi

    real(dp), dimension(:), allocatable :: vart               ! Temperature variance
    real(dp), dimension(:), allocatable :: varto
    real(dp), dimension(:), allocatable :: vartoo
    real(dp),dimension(:,:), allocatable :: dVartdxi

    ! real(dp), dimension(:), allocatable :: edd              ! Dissipation of temperature variance
    ! real(dp), dimension(:), allocatable :: Ret              ! Turbulent Re number

    real(dp), dimension(:), allocatable :: uu,vv,ww,uv,uw,vw  ! Reynolds stress tensor components
    real(dp), dimension(:), allocatable :: utt,vtt,wtt        ! Turbulent heat fluxes
    real(dp), dimension(:,:), allocatable :: bij              ! Reynolds stress anisotropy tensor
    real(dp), dimension(:), allocatable :: visw               ! Effective visc. for boundary face,
    real(dp), dimension(:), allocatable :: ypl                ! y+ non-dimensional distance from wall
    real(dp), dimension(:), allocatable :: tau                ! tau - wall shear stress

    real(dp), dimension(:), allocatable :: magStrain          ! Strain magnitude
    real(dp), dimension(:), allocatable :: Vorticity          ! Vorticity magnitude

end module