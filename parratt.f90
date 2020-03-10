subroutine parratt_born(q,lambda,d,rho,beta,Rgen,M,N)
!***************************************************************************
!Subroutine to calculate Specular Reflectivity using Parrat Algorithm with 
!Born Approximation for roughness
!M = No. of data points
!N = No. of slabs
!lambda = wavelength
!d = thickness of each slab
!Rgen = generated reflectivtiy data
!q = change in wave vector
!***************************************************************************
double precision :: q(0:M), Rgen(0:M)
double precision :: d(0:N+1), rho(0:N+1), beta(0:N+1), sigma(0:N+1), qc(0:N+1)
double precision :: lambda
double complex :: X, fact1, fact2, r(0:N+1), k1, k2, fact
double precision, parameter :: re=2.814e-5, pi=3.14157

do j=0,N+1
   qc(j)=dsqrt(16.0d0*pi*re*rho(j))
enddo

do i = 1,M
   r(N+1)=dcmplx(0.0d0,0.0d0)
   do j=N,0,-1
      k1=cdsqrt(dcmplx(q(i)**2-qc(j)**2,-32.0d0*beta(j)*pi**2/lambda**2))
      k2=cdsqrt(dcmplx(q(i)**2-qc(j+1)**2,-32.0d0*beta(j+1)*pi**2/lambda**2))
      X=(k1-k2)*cdexp(-k1*k2*sigma(j+1)**2/2)/(k1+k2)
      fact1=dcmplx(dcos(dble(k2)*d(j+1)),dsin(dble(k2)*d(j+1)))
      fact2=dexp(-aimag(k2)*d(j+1))
      fact=fact1*fact2
      r(j)=(X+r(j+1)*fact)/(1.0+X*r(j+1)*fact)
   enddo
   Rgen(i)=cdabs(r(0))**2
enddo   
end subroutine parratt_born


subroutine parratt(q,lambda,d,rho,beta,Rgen,M,N)
!***************************************************************************
!
!M = No. of data points
!N = No. of slabs
!lambda = wavelength
!d = thickness of each slab
!Rgen = generated reflectivtiy data
!q = change in wave vector
!***************************************************************************
double precision :: q(0:M), Rgen(0:M)
double precision :: d(0:N+1), rho(0:N+1), beta(0:N+1), qc(0:N+1)
double precision :: lambda
double complex :: X, fact1, fact2, r(0:N+1), k1, k2, fact
double precision, parameter :: re=2.814e-5, pi=3.14157

do j=0,N+1
   qc(j)=dsqrt(16.0d0*pi*re*rho(j))
enddo

do i = 1,M
   r(N+1)=dcmplx(0.0d0,0.0d0)
   do j=N,0,-1
      k1=cdsqrt(dcmplx(q(i)**2-qc(j)**2,-32.0d0*beta(j)*pi**2/lambda**2))
      k2=cdsqrt(dcmplx(q(i)**2-qc(j+1)**2,-32.0d0*beta(j+1)*pi**2/lambda**2))
      X=(k1-k2)/(k1+k2)
      fact1=dcmplx(dcos(dble(k2)*d(j+1)),dsin(dble(k2)*d(j+1)))
      fact2=dexp(-aimag(k2)*d(j+1))
      fact=fact1*fact2
      r(j)=(X+r(j+1)*fact)/(1.0+X*r(j+1)*fact)
   enddo
   Rgen(i)=cdabs(r(0))**2
enddo   
end subroutine parratt
