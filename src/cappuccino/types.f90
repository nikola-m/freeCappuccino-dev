module types
!
! define precision for floating-point numbers
!
  ! double precision    
  integer, parameter :: dp = kind(1.0d0) 


  ! Field equation type
  type FieldEquation 
    logical  :: solve = .false.                      ! To activate the solution of this field in the main function. 
    real(dp) :: urf = 0.7_dp                         ! Under-relaxation factors.
    real(dp) :: gdsU = 1.0_dp                        ! Deferred correction factor.
    character ( len=30 ) :: cScheme = 'linearUpwind' ! Convection scheme - default is second order upwind.
    character ( len=12 ) :: dScheme = 'skewness'     ! Difussion scheme, i.e. the method for normal gradient at face skewness/offset.
    integer :: nrelax = 0                            ! Relaxation parameter non-orthogonal correction for face gradient minimal/orthogonal/over-relaxed, 1/0/-1.
    character ( len=12 ) :: lSolver = 'bicgstab'     ! Linear algebraic solver.
    integer  :: maxiter = 10                         ! Max number of iterations in linear solver.
    real(dp) :: tolAbs = 1e-10                       ! Absolute residual level.
    real(dp) :: tolRel = 0.01_dp                     ! Relative drop in residual to exit linear solver.
    real(dp) :: resor                                ! Residual norm - e.g. L1 norm of initial residual
  end type

end module types
