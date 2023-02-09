module vortexIdentification
!
! Purpose:
!   Module contains functions for extracting coherent vortex structures in the flow.
! Description:
!   We include several well known methods such as identification using Q-critera, 
!   $\lambda_2$ (Jeong&Hussain JFM 1995), and other.
!
!  Author: Nikola Mirkov
!  Email: nikolamirkov@yahoo.com
!
!  Modified:
!    Apr 10, 2020.
!
!  This is a part of freeCappuccino. 
!  The code is licenced under GPL licence.
!
use types
use geometry 
use variables, only: dUdxi,dVdxi,dWdxi

implicit none

real(dp), dimension(:), allocatable :: Qvortex
real(dp), dimension(:), allocatable :: lambda2
real(dp), dimension(:), allocatable :: sigma1,sigma2,sigma3

public 

contains

subroutine setQvortex
!
! Calculates so called Q criteria field, defined by Q = 1/2 * (S^2 - Omega^2).
! If Q > 0, vortical motion exists.
! Iso-surfaces of this field define coherent structures that approximate vortices in the flow field.
!
  implicit none

  integer :: inp
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz 
  real(dp) :: s11,s12,s13,s22,s23,s33,w12,w13,w23
  real(dp) :: magVorticitySq,magStrainSq

  if( .not.allocated( Qvortex ) ) then
    allocate( Qvortex(numCells) )
  endif

  do inp=1,numCells

    dudx = dudxi(1,inp)
    dudy = dudxi(2,inp)
    dudz = dudxi(3,inp)

    dvdx = dvdxi(1,inp)
    dvdy = dvdxi(2,inp)
    dvdz = dvdxi(3,inp)

    dwdx = dwdxi(1,inp)
    dwdy = dwdxi(2,inp)
    dwdz = dwdxi(3,inp)

    ! Find strain rate tensor
    ! [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
    s11=dudx
    s12=0.5*(dudy+dvdx)
    s13=0.5*(dudz+dwdx)
    s22=dvdy
    s23=0.5*(dvdz+dwdy) 
    s33=dwdz

    ! Find antisymmetric part of velocity gradient tensor
    ! [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
    w12=(dudy - dvdx)
    w13=(dudz - dwdx)
    w23=(dvdz - dwdy)


    ! Find strain rate squared s^2 = 2*sij*sij
    magStrainSq = 2*(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2))

    ! Find Vorticity mag squared Om^2 = 2*wij*wij
    magVorticitySq = (w12**2 + w23**2 + w13**2)

    Qvortex(inp) = 0.5_dp*( magVorticitySq - magStrainSq )

  enddo

end subroutine


subroutine setLambda2
!
! Compute second largest eigenvalue of (Sik*Skj + Wik*Wkj) tensor.
! Return the array of -$\lambda_2$ values.
! Iso-surfaces of this field define coherent structures that approximate vortices in the flow field.
!
! See Jeong and Hussain "On the identification of a vortex", JFM, 1995.
!
! This subroutine is based on Fluent UDF found on CFD-online forum:
! https://www.cfd-online.com/Forums/main/99674-lambda-2-criterion.html
!

  implicit none

  integer :: inp 
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz 
  real(dp) :: s11,s12,s13,s22,s23,s33,w12,w13,w23
  real(dp) :: a,b,c,d,i,j,k,m,n,p
  real(dp) :: x,y,z,tmp
  real(dp) :: P11,P12,P13,P22,P23,P33 
  real(dp) :: lambda(3)

  if( .not.allocated( lambda2 ) ) then
    allocate( lambda2(numCells) )
  endif


  do inp=1,numCells

    dudx = dudxi(1,inp)
    dudy = dudxi(2,inp)
    dudz = dudxi(3,inp)

    dvdx = dvdxi(1,inp)
    dvdy = dvdxi(2,inp)
    dvdz = dvdxi(3,inp)

    dwdx = dwdxi(1,inp)
    dwdy = dwdxi(2,inp)
    dwdz = dwdxi(3,inp)

    ! Find strain rate tensor
    ! [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
    s11=dudx
    s12=0.5*(dudy+dvdx)
    s13=0.5*(dudz+dwdx)
    s22=dvdy
    s23=0.5*(dvdz+dwdy) 
    s33=dwdz

    ! Find antisymmetric part of velocity gradient tensor
    ! [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
    w12=0.5*(dudy - dvdx)
    w13=0.5*(dudz - dwdx)
    w23=0.5*(dvdz - dwdy)


    P11=S11*S11+S12*S12+S13*S13-W12*W12-W13*W13
    P12=S12*(S11+S22)+S13*S23-W13*W23
    P13=S13*(S11+S33)+S12*S23+W12*W23
    P22=S12*S12+S22*S22+S23*S23-W12*W12-W23*W23
    P23=S23*(S22+S33)+S12*S13-W12*W13
    P33=S13*S13+S23*S23+S33*S33-W13*W13-W23*W23

    ! Coefficients of the characteristic polynomial

    ! a*lambda^3 + b*lambda^2 + c*lambda + d = 0

    a=-1.0
    b=P11+P22+P33
    c=P12*P12+P13*P13+P23*P23-P11*P22-P11*P33-P22*P33
    d=P11*P22*P33+2.0*P12*P13*P23-P12*P12*P33-P13*P13*P22-P23*P23*P11

    ! Resolution of the cubic equation, eigenvalues assumed to be real

    x=((3.0*c/a)-b*b/(a*a))/3.0
    y=(2.0*b*b*b/(a*a*a)-9.0*b*c/(a*a)+27.0*d/a)/27.0
    z=y*y/4.0+x*x*x/27.0

    i=sqrt(y*y/4.0-z)
    j=-i**(1.0/3.0)
    k=acos(-(y/(2.0*i)))
    m=cos(k/3.0)
    n=sqrt(3.0)*sin(k/3.0)
    p=b/(3.0*a)

    lambda(1)=2.0*j*m+p
    lambda(2)=-j*(m+n)+p
    lambda(3)=-j*(m-n)+p

    ! Ordering of the eigenvalues, lamda(1) >= lambda(2) >= lambda(3) >= 0,

    if(lambda(2)>lambda(1)) then
      tmp=lambda(2)
      lambda(2)=lambda(1)
      lambda(1)=tmp
    endif

    if(lambda(3)>lambda(2)) then

      tmp=lambda(3)
      lambda(3)=lambda(2)
      lambda(2)=tmp

      if(lambda(2)>lambda(1)) then 
        tmp=lambda(2)
        lambda(2)=lambda(1)
        lambda(1)=tmp
      endif

    endif

    ! Retrieval of the second eigenvalue

    lambda2(inp) = lambda(2)

  enddo

end subroutine



subroutine setSigma
!
! Compute singular values of (Sik*Skj + Wik*Wkj) tensor.
! We actually compute singular values of DT*D, where D i velocity gradient tensor.
! This subroutine is  based o the setLambda2 routine above.
!

  implicit none

  integer :: inp 
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz 
  real(dp) :: s11,s12,s13,s22,s23,s33,w12,w13,w23
  real(dp) :: a,b,c,d,i,j,k,m,n,p
  real(dp) :: x,y,z,tmp
  real(dp) :: P11,P12,P13,P22,P23,P33 
  real(dp) :: lambda(3)

  if( .not.allocated( sigma1 ) ) then
    allocate( sigma1(numCells),sigma2(numCells),sigma3(numCells) )
  endif


  do inp=1,numCells

    dudx = dudxi(1,inp)
    dudy = dudxi(2,inp)
    dudz = dudxi(3,inp)

    dvdx = dvdxi(1,inp)
    dvdy = dvdxi(2,inp)
    dvdz = dvdxi(3,inp)

    dwdx = dwdxi(1,inp)
    dwdy = dwdxi(2,inp)
    dwdz = dwdxi(3,inp)

    ! Find strain rate tensor
    ! [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
    s11=dudx
    s12=0.5*(dudy+dvdx)
    s13=0.5*(dudz+dwdx)
    s22=dvdy
    s23=0.5*(dvdz+dwdy) 
    s33=dwdz

    ! Find antisymmetric part of velocity gradient tensor
    ! [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
    w12=0.5*(dudy - dvdx)
    w13=0.5*(dudz - dwdx)
    w23=0.5*(dvdz - dwdy)


    P11=S11*S11+S12*S12+S13*S13-W12*W12-W13*W13
    P12=S12*(S11+S22)+S13*S23-W13*W23
    P13=S13*(S11+S33)+S12*S23+W12*W23
    P22=S12*S12+S22*S22+S23*S23-W12*W12-W23*W23
    P23=S23*(S22+S33)+S12*S13-W12*W13
    P33=S13*S13+S23*S23+S33*S33-W13*W13-W23*W23

    ! Coefficients of the characteristic polynomial

    ! a*lambda^3 + b*lambda^2 + c*lambda + d = 0

    a=-1.0
    b=P11+P22+P33
    c=P12*P12+P13*P13+P23*P23-P11*P22-P11*P33-P22*P33
    d=P11*P22*P33+2.0*P12*P13*P23-P12*P12*P33-P13*P13*P22-P23*P23*P11

    ! Resolution of the cubic equation, eigenvalues assumed to be real

    x=((3.0*c/a)-b*b/(a*a))/3.0
    y=(2.0*b*b*b/(a*a*a)-9.0*b*c/(a*a)+27.0*d/a)/27.0
    z=y*y/4.0+x*x*x/27.0

    i=sqrt(y*y/4.0-z)
    j=-i**(1.0/3.0)
    k=acos(-(y/(2.0*i)))
    m=cos(k/3.0)
    n=sqrt(3.0)*sin(k/3.0)
    p=b/(3.0*a)

    lambda(1)=2.0*j*m+p
    lambda(2)=-j*(m+n)+p
    lambda(3)=-j*(m-n)+p

    ! Ordering of the eigenvalues, lamda(1) >= lambda(2) >= lambda(3) >= 0,

    if(lambda(2)>lambda(1)) then
      tmp=lambda(2)
      lambda(2)=lambda(1)
      lambda(1)=tmp
    endif

    if(lambda(3)>lambda(2)) then

      tmp=lambda(3)
      lambda(3)=lambda(2)
      lambda(2)=tmp

      if(lambda(2)>lambda(1)) then 
        tmp=lambda(2)
        lambda(2)=lambda(1)
        lambda(1)=tmp
      endif

    endif

    ! Setting singular values

    sigma1(inp) = sqrt( lambda(1) )
    sigma2(inp) = sqrt( lambda(2) )
    sigma3(inp) = sqrt( lambda(3) )

    ! write(*,'(2x,3es11.4)') lambda(1),lambda(2),lambda(3)

  enddo

end subroutine


end module