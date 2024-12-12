!===============================================================================
!
! Cahn-Hilliard model for dynamic simulation of vapour-liquid interfaces
! in a one-dimensional ternary system of water, n-dodecane and 1-butanol.
!
! Author: Lukáš Šatura
! 2021
!
!===============================================================================
program main
!===============================================================================
use io
implicit none
!===============================================================================
! Variables
double precision, allocatable :: rho1_new(:), rho1_old(:), rho2_new(:), rho2_old(:),rho3_new(:), rho3_old(:), rhot_old(:) ! Density fields
double precision, allocatable :: y1(:),y2(:),y3(:)      ! Arrays of molar fractions
double precision, allocatable :: x(:)                   ! Array of spatial position
double precision :: factor                              ! Factor of numerical stability
double precision :: sig,igr                             ! Integrated surface tension (sig) with respective integral value (igr)

integer :: j                                            ! Iteration counter
integer :: i                                            ! Spatial index
integer :: p,No_plots                                   ! Plotting indeces
character(len=100) :: filename                          ! String representing output filename
character(len=10) :: iter                               ! String representing iteration number
double precision, allocatable :: mu1(:), dJrho1(:), d2rho1(:), mu_gen1(:), d2mu1(:), Jrho1(:)  !chemical potential, etc.
double precision, allocatable :: mu2(:), dJrho2(:), d2rho2(:), mu_gen2(:), d2mu2(:), Jrho2(:)
double precision, allocatable :: mu3(:), dJrho3(:), d2rho3(:), mu_gen3(:), d2mu3(:), Jrho3(:)
double precision, allocatable :: D12(:), D13(:), D23(:), phi1(:), phi2(:),phi3(:), phi12(:), phi13(:), phi23(:)
double precision, allocatable :: der1(:), der12(:), der2(:), der13(:), der3(:), der23(:), der(:), ampl(:)  ! Auxiliaries for surface tension integration der1 = (drho1/dx)^2
real :: time_begin, time_end                            ! Variables for measuring total CPU time
DOUBLE PRECISION, DIMENSION(2,2) :: MAT, MATINV
DOUBLE PRECISION, DIMENSION(2,1) :: mu_vecs, u_vecs
logical :: OK_FLAG
double precision, parameter :: Pi = 3.141592653589793238462643d0
!===============================================================================
call readpars
!===============================================================================
! Check numerical stability of the proposed scheme - CFL condition
!factor = Diff*Kap*rho_L*k/h**4
!write(*,*) 'factor of stability CFL = ',factor
write(*,*) 'spatial step h = ',h

if (factor.gt.0.5d0) then
  write(*,*) 'Numerical stability condition not fulfilled.'
  write(*,*) 'Adjust spatial and/or temporal step.'
  !stop
end if


! Allocate arrays for density fields
allocate(rho1_new(N), rho1_old(N), rho2_new(N), rho2_old(N),rho3_new(N), rho3_old(N), rhot_old(N),x(N))
allocate(y1(N),y2(N),y3(N))
allocate(mu1(N),dJrho1(N),d2rho1(N),mu_gen1(N),d2mu1(N),Jrho1(N))
allocate(mu2(N),dJrho2(N),d2rho2(N),mu_gen2(N),d2mu2(N),Jrho2(N))
allocate(mu3(N),dJrho3(N),d2rho3(N),mu_gen3(N),d2mu3(N),Jrho3(N))
allocate(D12(N),D13(N),D23(N),phi1(N), phi2(N),phi3(N), phi12(N), phi13(N), phi23(N))
allocate(der1(N-2),der12(N-2),der2(N-2),der13(N-2),der3(N-2),der23(N-2),der(N-2),ampl(N))
!===============================================================================
! Generate array of spatial position values (only for output purposes)
x(1) = 0.0d0
!x(1) = - dble(N/2)*h
do i=2,N
  x(i) = x(i-1) + h
end do
!======================= INITIAL CONDITIONS ====================================
if (init) then
  !if (per_BC) then
    ! Create random-noise density profile
    call random_number(ampl)
    do i=1,N
      rho1_old(i) = (rho_V(1)+rho_L(1))/(2.0d0) + 10.0d0*ampl(i)*cos(5.0d0*Pi*x(i))
      rho2_old(i) = (rho_V(2)+rho_L(2))/(2.0d0) + 10.0d0*ampl(i)*cos(5.0d0*Pi*x(i)-0.1d-8)
      rho3_old(i) = (rho_V(3)+rho_L(3))/(2.0d0) + 10.0d0*ampl(i)*cos(5.0d0*Pi*x(i)+0.1d-8)
    end do
!  else
    ! Create a linear density profile
    ! call linspace(rho1_old,rho_V(1)/rhoc(1),rho_L(1)/rhoc(1),N)
    ! call linspace(rho2_old,rho_V(2)/rhoc(2),rho_L(2)/rhoc(2),N)
    ! call linspace(rho3_old,rho_V(3)/rhoc(3),rho_L(3)/rhoc(3),N)
!  end if

  ! Save the initial density profile into a file
  open(100,file=''//trim(path)//'rho_ini.dat',status='unknown')
  do i=1,N
    write(100,'(E10.2,3(F10.2))') x(i), rho1_old(i), rho2_old(i), rho3_old(i)
    write(*,'(E10.2,3(F10.2))') x(i),  rho1_old(i), rho2_old(i), rho3_old(i)
  end do
  close(100)
  init_iter = 1
else
  write(iter,'(I10.10)') init_iter
  filename = 'rho_'//trim(iter)//'.dat'
  open(10,file=''//trim(path)//trim(filename)//'',status='old')
  do i=1,N
    read(10,'(E10.2,3(F10.2))') x(i), rho1_old(i), rho2_old(i), rho3_old(i)
    write(*,'(E10.2,3(F10.2))') x(i), rho1_old(i), rho2_old(i), rho3_old(i)
  end do
  close(10)
end if

!===============================================================================
! Start measuring time for computation
call CPU_TIME ( time_begin )

write(*,*) 'Iteration / Total'
!================================= MODEL =======================================
do j=init_iter,M

  if (per_BC) then
    ! Impose periodic BCs
    rho1_old(1:2) = rho1_old(N-3:N-2)
    rho1_old(N-1:N) = rho1_old(3:4)
    rho2_old(1:2) = rho2_old(N-3:N-2)
    rho2_old(N-1:N) = rho2_old(3:4)
    rho3_old(1:2) = rho3_old(N-3:N-2)
    rho3_old(N-1:N) = rho3_old(3:4)
  end if

  call potential(mu1,mu2,mu3,rho1_old,rho2_old,rho3_old,N)
  call laplace(d2rho1,rho1_old)
  call laplace(d2rho2,rho2_old)
  call laplace(d2rho3,rho3_old)

  do i=1,N
    mu_gen1(i) = mu1(i)/R/T - Kap(1)*d2rho1(i) - Kap(4)*d2rho2(i) - Kap(5)*d2rho3(i)
    mu_gen2(i) = mu2(i)/R/T - Kap(2)*d2rho2(i) - Kap(4)*d2rho1(i) - Kap(6)*d2rho3(i)
    mu_gen3(i) = mu3(i)/R/T - Kap(3)*d2rho3(i) - Kap(5)*d2rho1(i) - Kap(6)*d2rho2(i)
  end do

  call gradient(d2mu1,mu_gen1)
  call gradient(d2mu2,mu_gen2)
  call gradient(d2mu3,mu_gen3)

  do i=1,N
    phi1(i) = (R*T)/(Diff(1)*(rho1_old(i) + rho2_old(i)*Diff(1)/Diff(2) + rho3_old(i)*Diff(1)/Diff(3)))
    phi2(i) = (R*T)/(Diff(2)*(rho2_old(i) + rho1_old(i)*Diff(2)/Diff(1) + rho3_old(i)*Diff(2)/Diff(3)))
    phi3(i) = (R*T)/(Diff(3)*(rho3_old(i) + rho1_old(i)*Diff(3)/Diff(1) + rho2_old(i)*Diff(3)/Diff(2)))
    phi12(i)= dsqrt(phi1(i)*phi2(i))
    phi13(i)= dsqrt(phi1(i)*phi3(i))
    phi23(i)= dsqrt(phi3(i)*phi2(i))

    rhot_old(i) = rho1_old(i) + rho2_old(i) + rho3_old(i)
    y1(i) = rho1_old(i)/rhot_old(i)
    y2(i) = rho2_old(i)/rhot_old(i)
    y3(i) = rho3_old(i)/rhot_old(i)

    D12(i) = (R*T)/(rhot_old(i)*phi12(i))
    D13(i) = (R*T)/(rhot_old(i)*phi13(i))
    D23(i) = (R*T)/(rhot_old(i)*phi23(i))

    MAT(1,1) = y2(i)/D12(i) + y3(i)/D13(i) + y3(i)*Z(1)*rho1_old(i)/D13(i)/Z(3)/rho3_old(i)
    MAT(2,1) = -y1(i)/D12(i) + y3(i)/D23(i)*Z(1)*rho1_old(i)/Z(3)/rho3_old(i)
    MAT(1,2) = -y2(i)/D12(i) + y3(i)*Z(2)*rho2_old(i)/D13(i)/rho3_old(i)
    MAT(2,2) = y1(i)/D12(i) + y3(i)/D23(i) + y3(i)*Z(2)/Z(3)/D23(i)*rho2_old(i)/rho3_old(i)
    CALL M22INV (MAT, MATINV, OK_FLAG)

    mu_vecs(1,1) = d2mu1(i)
    mu_vecs(2,1) = d2mu2(i)

    u_vecs = matmul(MATINV,mu_vecs)

    Jrho1(i) = rho1_old(i)*u_vecs(1,1)
    Jrho2(i) = rho2_old(i)*u_vecs(2,1)
    Jrho3(i) = - (Z(1)*Jrho1(i) + Z(2)*Jrho2(i))/Z(3)
  end do

  call gradient(dJrho1,Jrho1)
  call gradient(dJrho2,Jrho2)
  call gradient(dJrho3,Jrho3)

  ! Integrate using explicit Euler scheme
  do i=1,N
    rho1_new(i) = rho1_old(i) + k*dJrho1(i)
    rho2_new(i) = rho2_old(i) + k*dJrho2(i)
    rho3_new(i) = rho3_old(i) + k*dJrho3(i)
  end do

  if (.not.per_BC) then
    ! Impose Dirichlet BCs
    rho1_new(N) = rho1_old(N) ! rho(N) = rho_L
    rho1_new(1) = rho1_old(1) ! rho(1) = rho_V
    rho2_new(N) = rho2_old(N) ! rho(N) = rho_V
    rho2_new(1) = rho2_old(1) ! rho(1) = rho_L
    rho3_new(N) = rho3_old(N) ! rho(N) = rho_V
    rho3_new(1) = rho3_old(1) ! rho(1) = rho_L
  end if

  ! Update the profile
  do i=1,N
    rho1_old(i) = rho1_new(i)
    rho2_old(i) = rho2_new(i)
    rho3_old(i) = rho3_new(i)
  end do
!=============================== PLOTTING ======================================
  ! Save output for selected iterations in variable 'plot_it'
 No_plots = size(plot_it)
  do p=1,No_plots
    if (j.eq.plot_it(p)) then
    ! Prepare output filename
      write(iter,'(I10.10)') j
      filename = 'rho_'//trim(iter)//'.dat'
      ! Save current profile into the file
      open(10,file=''//trim(path)//trim(filename)//'',status='unknown')
      write(*,*) 'j = ',j
      write(*,*) 'x','rho_old','mu'
        do i=1,N
          write(10,'(E10.2,3(F10.2))') x(i), rho1_old(i), rho2_old(i), rho3_old(i)

          write(*,'(E10.2,6(F10.2))') x(i), rho1_old(i), mu1(i), rho2_old(i), mu2(i), rho3_old(i), mu3(i)
        end do
      close(10)
    end if
  end do
  ! Print current progress on a screen
  if (mod(j,1000).eq.0) write(*,*) j, ' / ', M
end do
!===============================================================================
! End time measurement
CALL CPU_TIME ( time_end )

write(*,*) 'j = ',j
write(*,*) 'x','rho_old','mu'
do i=1,N
  write(*,'(E10.2,6(F10.2))') x(i), rho1_old(i), mu1(i), rho2_old(i), mu2(i), rho3_old(i), mu3(i)
end do

write(*,*) '========================================================'
write(*,'(A,F9.2,A)') 'Time of operation was ', (time_end - time_begin)/60, ' minutes.'
write(*,*) '========================================================'
!========================= SURFACE TENSION =====================================

do i=1,N-2
  der1(i)  = ((rho1_old(i+2) - rho1_old(i))/(2*h))**2
  der12(i) = ((rho1_old(i+2) - rho1_old(i))/(2*h))*((rho2_old(i+2) - rho2_old(i))/(2*h))
  der2(i)  = ((rho2_old(i+2) - rho2_old(i))/(2*h))**2
  der13(i)  = ((rho1_old(i+2) - rho1_old(i))/(2*h))*((rho3_old(i+2) - rho3_old(i))/(2*h))
  der3(i)  = ((rho3_old(i+2) - rho3_old(i))/(2*h))**2
  der23(i)  = ((rho3_old(i+2) - rho3_old(i))/(2*h))*((rho2_old(i+2) - rho2_old(i))/(2*h))
end do

der = Kap(1)*rhoc(1)**2*der1 + 2*Kap(4)*rhoc(1)*rhoc(2)*der12 + Kap(2)*rhoc(2)**2*der2 &
      + 2*Kap(5)*rhoc(1)*rhoc(3)*der13 +  Kap(3)*rhoc(3)**2*der3 + 2*Kap(6)*rhoc(2)*rhoc(3)*der23
do i=1,N-3
  igr=igr+h*dabs(der(i+1)-der(i))/2.0d0
end do

sig = R*T/L*igr

write(*,'(A,F9.2,A)') 'Integrated surface tension gamma:',real(sig*1.0d3),'mN•m-1 (DIVIDED BY THE NUMBER OF CREATED INTERFACES!!!).'
write(*,*) '========================================================'
write(*,'(A,/,(A,E,X,A,/))') 'EOS parameters water:','a = ',a(1),'J•m3/mol2','b = ',b(1),'m3/mol','Kappa = ',Kap(1),'m5/mol'
write(*,'(A,/,(A,E,X,A,/))') 'EOS parameters dodecane:','a = ',a(2),'J•m3/mol2','b = ',b(2),'m3/mol','Kappa = ',Kap(2),'m5/mol'
write(*,'(A,/,(A,E,X,A,/))') 'EOS parameters butanol:','a = ',a(3),'J•m3/mol2','b = ',b(3),'m3/mol','Kappa = ',Kap(3),'m5/mol'

open(20,file=''//trim(path)//'output.dat',status='unknown')
  write(20,'(A,I2.2)') 'Simulation number ',simNo
  write(20,'(A,/,A,I,/,A,I,/,(A,E,/))') 'Simulation data:','N = ',N,'M = ',M,'L = ',L,'k = ',k,'h = ',h,'T =',T
  write(20,'(A,F9.2,A)') 'Integrated surface tension gamma:',real(sig*1.0d3),'mN•m-1 (DIVIDED BY THE NUMBER OF CREATED INTERFACES!!!).'
  write(20,*) '========================================================'
  write(20,'(A,/,(A,E,X,A,/))') 'EOS parameters water:','a = ',a(1),'J•m3/mol2','b = ',b(1),'m3/mol','Kappa = ',Kap(1),'m5/mol'
  write(20,'(A,/,(A,E,X,A,/))') 'EOS parameters dodecane:','a = ',a(2),'J•m3/mol2','b = ',b(2),'m3/mol','Kappa = ',Kap(2),'m5/mol'
  write(20,'(A,/,(A,E,X,A,/))') 'EOS parameters butanol:','a = ',a(3),'J•m3/mol2','b = ',b(3),'m3/mol','Kappa = ',Kap(3),'m5/mol'
close(20)

! do i=1,N-2
!   write(*,*) der1,der12,der2
! end do

deallocate(rho1_new,rho1_old,x,mu1,dJrho1,d2rho1,mu_gen1,d2mu1,Jrho1)
deallocate(rho2_new,rho2_old,mu2,dJrho2,d2rho2,mu_gen2,d2mu2,Jrho2)
deallocate(rho3_new,rho3_old,mu3,dJrho3,d2rho3,mu_gen3,d2mu3,Jrho3, rhot_old)
deallocate(y1,y2,y3)
deallocate(phi1, phi2,phi3, phi12, phi13, phi23)
deallocate(D12,der1,der12,der2,der,ampl)
deallocate(D13,D23,der13,der23,der3)
stop
end program
!===============================================================================


!==================================SUBROUTTINES=================================


!===============================================================================
subroutine gradient(dy_dx, y)
!===============================================================================
use io
implicit none

integer :: i
double precision,intent(in) :: y(N)
double precision,intent(out) :: dy_dx(N)

!==ALL formulae are of 4th order accuracy==!

!$OMP PARALLEL DO
do i=3,(N-2)
  dy_dx(i) = (1.0d0*y(i-2) - 8.0d0*y(i-1) + 0.0d0*y(i) + 8.0d0*y(i+1) - 1.0d0*y(i+2))
end do

if (per_BC) then
  ! Imposing periodic BCs
  dy_dx(1:2) = dy_dx(N-3:N-2)
  dy_dx(N-1:N) = dy_dx(3:4)
else
  !dy_dx(N)   = 0.0d0    ! Imposing zero second derivative as a boundary condition
  dy_dx(N)   = ( 25.0d0*y(N) - 48.0d0*y(N-1) + 36.0d0*y(N-2) - 16.0d0*y(N-3) + 3.0d0*y(N-4))
  dy_dx(N-1) = ( 3.0d0*y(N)  + 10.0d0*y(N-1) - 18.0d0*y(N-2) + 6.0d0*y(N-3)  - 1.0d0*y(N-4))
  dy_dx(2)   = (-3.0d0*y(1)  - 10.0d0*y(2)   + 18.0d0*y(3)   - 6.0d0*y(4)    + 1.0d0*y(5))
  dy_dx(1)   = (-25.0d0*y(1) + 48.0d0*y(2)   - 36.0d0*y(3)   + 16.0d0*y(4)   - 3.0d0*y(5))
  !dy_dx(1)   = 0.0d0    ! Imposing zero second derivative as a boundary condition
end if

dy_dx = dy_dx / (12.0d0*h)

return
end subroutine
!===============================================================================
subroutine laplace(dy_dx, y)
!===============================================================================
use io
implicit none

integer :: i
double precision,intent(in) :: y(N)
double precision,intent(out) :: dy_dx(N)

!$OMP PARALLEL DO
do i=3,(N-2)
  dy_dx(i) = (-1.0d0*y(i-2) + 16.0d0*y(i-1) - 30.0d0*y(i) + 16.0d0*y(i+1) - 1.0d0*y(i+2))                      ! 4th order accuracy
end do

if (per_BC) then
  ! Imposing periodic BCs
  dy_dx(1:2) = dy_dx(N-3:N-2)
  dy_dx(N-1:N) = dy_dx(3:4)
else
  ! Explicit formulae for boundary points
  dy_dx(N) = 0.0d0      ! Imposing zero second derivative as a boundary condition
  ! dy_dx(N) = (35.0d0*y(N) - 104.0d0*y(N-1) + 114.0d0*y(N-2) - 56.0d0*y(N-3) + 11.0d0*y(N-4))                  ! 3rd order accuracy
  ! dy_dx(N) = (45.0d0*y(N) - 154.0d0*y(N-1) + 214.0d0*y(N-2) - 156.0d0*y(N-3) + 61.0d0*y(N-4) - 10.0d0*y(N-5)) ! 4th order accuracy

  dy_dx(N-1) = (10.0d0*y(N) - 15.0d0*y(N-1) - 4.0d0*y(N-2) + 14.0d0*y(N-3) - 6.0d0*y(N-4) + 1.0d0*y(N-5))      ! 4th order accuracy
  !dy_dx(N-1) = (11.0d0*y(N) - 20.0d0*y(N-1) + 6.0d0*y(N-2) + 4.0d0*y(N-3) - 1.0d0*y(N-4))                      ! 3rd order accuracy

  dy_dx(2)   = (10.0d0*y(1) - 15.0d0*y(2)   - 4.0d0*y(3)   + 14.0d0*y(4)   - 6.0d0*y(5)   + 1.0d0*y(6))        ! 4th order accuracy
  !dy_dx(2) = (11.0d0*y(1) - 20.0d0*y(2) + 6.0d0*y(3) + 4.0d0*y(4) - 1.0d0*y(5))                                ! 3rd order accuracy
  dy_dx(1) = 0.0d0      ! Imposing zero second derivative as a boundary condition
  !dy_dx(1) = (35.0d0*y(1) - 104.0d0*y(2) + 114.0d0*y(3) - 56.0d0*y(4) + 11.0d0*y(5))                          ! 3rd order accuracy
  !dy_dx(1) = (45.0d0*y(1) - 154.0d0*y(2) + 214.0d0*y(3) - 156.0d0*y(4) + 61.0d0*y(5) - 10.0d0*y(6))           ! 4th order accuracy
end if

dy_dx = dy_dx / (12.0d0*h**2)

return
end subroutine
!===============================================================================
subroutine linspace(x, x_start, x_end, x_len)
!===============================================================================
implicit none

double precision :: x_start, x_end, dx
integer :: x_len, k
double precision,dimension(x_len),intent(out) :: x

dx = (x_end - x_start) / dble(x_len - 1)

do k = 1, x_len
  x(k) = x_start + dble(k-1)*dx
end do

return
end subroutine
!===============================================================================
  SUBROUTINE M22INV (A, AINV, OK_FLAG)

  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN)  :: A
  DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: AINV
  LOGICAL, INTENT(OUT) :: OK_FLAG

  DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
  DOUBLE PRECISION :: DET
  DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR


  DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

  IF (ABS(DET) .LE. EPS) THEN
     AINV = 0.0D0
     OK_FLAG = .FALSE.
     RETURN
  END IF

  COFACTOR(1,1) = +A(2,2)
  COFACTOR(1,2) = -A(2,1)
  COFACTOR(2,1) = -A(1,2)
  COFACTOR(2,2) = +A(1,1)

  AINV = TRANSPOSE(COFACTOR) / DET

  OK_FLAG = .TRUE.

  RETURN

  END SUBROUTINE M22INV

!===============================================================================
SUBROUTINE M33INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV
!===============================================================================
subroutine potential(mu1,mu2,mu3,rho1_old,rho2_old,rho3_old)
!===============================================================================
use io
implicit none
! Subroutine for calculation of chemical potential acc. to Peng-Robinson's EoS
integer :: i
double precision,intent(in) :: rho1_old(N), rho2_old(N), rho3_old(N)
double precision,intent(out) :: mu1(N),mu2(N), mu3(N)
double precision :: rho1(N), rho2(N),rho3(N), rho(N), aa(N), bb(N)
call readpars

! Calculate potential profile for pertinent density profile

rho1 = rho1_old
rho2 = rho2_old
rho3 = rho3_old

!!$OMP PARALLEL PRIVATE(mu1,mu2), SHARED(rho,aa,bb) num_threads(2)
!!$OMP PARALLEL DO
do i=1,N
  rho(i) = rho1(i) + rho2(i) + rho3(i)
  aa(i) = a(1)*(rho1(i)**2)/(rho(i)**2) + a(2)*(rho2(i)**2)/(rho(i)**2) + a(3)*(rho3(i)**2)/(rho(i)**2) &
          + 2*dsqrt(a(1)*a(2))*(rho1(i)*rho2(i))/(rho(i)**2) + 2*dsqrt(a(1)*a(3))*(rho1(i)*rho3(i))/(rho(i)**2) &
          + 2*dsqrt(a(3)*a(2))*(rho3(i)*rho2(i))/(rho(i)**2)
  bb(i) = b(1)*rho1(i)/rho(i) + b(2)*rho2(i)/rho(i) + b(3)*rho3(i)/rho(i)
  mu1(i) = R*T*dlog(dabs(rho1(i)/(1.0d0 - bb(i)*rho(i)))) + R*T*b(1)*rho(i)/(1.0d0 - bb(i)*rho(i)) &
      + (1.0d0/2.0d0/dsqrt(2.0d0))*((2.0d0*a(1)*rho1(i) + 2.0d0*dsqrt(a(1)*a(2))*rho2(i) + 2.0d0*dsqrt(a(1)*a(3))*rho3(i) )/bb(i)/rho(i) - aa(i)*b(1)/bb(i)**2) &
      * dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*bb(i)*rho(i))/(1.0d0 + (1.0d0 + dsqrt(2.0d0))*bb(i)*rho(i)))) &
      - aa(i)*rho(i)*b(1)/bb(i)/(1.0d0 + 2.0d0*bb(i)*rho(i) - bb(i)**2*rho(i)**2)
  mu2(i) = R*T*dlog(dabs(rho2(i)/(1.0d0 - bb(i)*rho(i)))) + R*T*b(2)*rho(i)/(1.0d0 - bb(i)*rho(i)) &
      + (1.0d0/2.0d0/dsqrt(2.0d0))*((2.0d0*a(2)*rho2(i) + 2.0d0*dsqrt(a(1)*a(2))*rho1(i) + 2.0d0*dsqrt(a(3)*a(2))*rho3(i))/bb(i)/rho(i) - aa(i)*b(2)/bb(i)**2) &
      * dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*bb(i)*rho(i))/(1.0d0 + (1.0d0 + dsqrt(2.0d0))*bb(i)*rho(i)))) &
      - aa(i)*rho(i)*b(2)/bb(i)/(1.0d0 + 2.0d0*bb(i)*rho(i) - bb(i)**2*rho(i)**2)
  mu3(i) = R*T*dlog(dabs(rho3(i)/(1.0d0 - bb(i)*rho(i)))) + R*T*b(3)*rho(i)/(1.0d0 - bb(i)*rho(i)) &
      + (1.0d0/2.0d0/dsqrt(2.0d0))*((2.0d0*a(3)*rho3(i) + 2.0d0*dsqrt(a(1)*a(3))*rho1(i) + 2.0d0*dsqrt(a(3)*a(2))*rho2(i))/bb(i)/rho(i) - aa(i)*b(3)/bb(i)**2) &
      * dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*bb(i)*rho(i))/(1.0d0 + (1.0d0 + dsqrt(2.0d0))*bb(i)*rho(i)))) &
      - aa(i)*rho(i)*b(3)/bb(i)/(1.0d0 + 2.0d0*bb(i)*rho(i) - bb(i)**2*rho(i)**2)
end do
!!$OMP END PARALLEL
return
end subroutine
!===============================================================================
