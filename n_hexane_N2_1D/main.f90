!===============================================================================
!
! Cahn-Hilliard model for dynamic simulation of vapour-liquid interface
! in a one-dimensional binary system of n-hexane and nitrogen.
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
double precision, allocatable :: rho1_new(:), rho1_old(:), rho2_new(:), rho2_old(:) ! Density fields
double precision, allocatable :: x(:)                                               ! Array of spatial position
double precision :: factor                                                          ! Factor of numerical stability
double precision :: sig,igr                                                         ! Integrated surface tension (sig) with respective integral value (igr)

integer :: j                                                                        ! Iteration counter
integer :: i                                                                        ! Spatial index
integer :: p,No_plots                                                               ! Plotting indeces
character(len=100) :: filename                                                      ! String representing output filename
character(len=10) :: iter                                                           ! String representing iteration number
double precision, allocatable :: mu1(:), dJrho1(:), d2rho1(:), mu_gen1(:), d2mu1(:), Jrho1(:)  !chemical potential, etc.
double precision, allocatable :: mu2(:), dJrho2(:), d2rho2(:), mu_gen2(:), d2mu2(:), Jrho2(:)
double precision, allocatable :: D12(:), der1(:), der12(:), der2(:), der(:), ampl(:)  ! Auxiliaries for surface tension integration der1 = (drho1/dx)^2
real :: time_begin, time_end                                                        ! Variables for measuring total CPU time
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
allocate(rho1_new(N), rho1_old(N), rho2_new(N), rho2_old(N),x(N))
allocate(mu1(N),dJrho1(N),d2rho1(N),mu_gen1(N),d2mu1(N),Jrho1(N))
allocate(mu2(N),dJrho2(N),d2rho2(N),mu_gen2(N),d2mu2(N),Jrho2(N))
allocate(D12(N),der1(N-2),der12(N-2),der2(N-2),der(N-2),ampl(N))
!===============================================================================
! Generate array of spatial position values (only for output purposes)
x(1) = 0.0d0
!x(1) = - dble(N/2)*h
do i=2,N
  x(i) = x(i-1) + h
end do
!======================= INITIAL CONDITIONS ====================================
if (init) then
  if (per_BC) then
    ! Create random-noise density profile
    call random_number(ampl)
    do i=1,N
      rho1_old(i) = (rho_V(1)+rho_L(1))/(2.0d0) + ampl(i)/10.0d0*cos(5.0d0*Pi*x(i))
      rho2_old(i) = (rho_V(2)+rho_L(2))/(2.0d0) + ampl(i)/100.0d0*sin(5.0d0*Pi*x(i))
    end do
  else
    ! Create a linear density profile
    call linspace(rho1_old,rho_V(1),rho_L(1),N)
    call linspace(rho2_old,rho_V(2),rho_L(2),N)
  end if

  if (dmlss) then
    do i=1,N
      rho1_old(i) = rho1_old(i)/rhoc(1)
      rho2_old(i) = rho2_old(i)/rhoc(2)
    end do
  end if

  ! Save the initial density profile into a file
  open(10,file=''//trim(path)//'rho_ini.dat',status='unknown')
  do i=1,N
    write(10,'(E10.2,2(F10.2))') x(i), rho1_old(i), rho2_old(i)
    write(*,'(E10.2,2(F10.2))') x(i),  rho1_old(i), rho2_old(i)
  end do
  close(10)
  init_iter = 1
else
  write(iter,'(I10.10)') init_iter
  filename = 'rho_'//trim(iter)//'.dat'
  open(10,file=''//trim(path)//trim(filename)//'',status='old')
  do i=1,N
    read(10,'(E10.2,2(F10.2))') x(i), rho1_old(i), rho2_old(i)
    write(*,'(E10.2,2(F10.2))') x(i), rho1_old(i), rho2_old(i)
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
  end if

  call potential(mu1,mu2,rho1_old,rho2_old,N)
  call laplace(d2rho1,rho1_old)
  call laplace(d2rho2,rho2_old)

  if (dmlss) then
    do i=1,N
      mu_gen1(i) = mu1(i)/R/T - Kap(1)*rhoc(1)/L**2*d2rho1(i) - Kap(3)*rhoc(2)/L**2*d2rho2(i)
      mu_gen2(i) = mu2(i)/R/T - Kap(2)*rhoc(2)/L**2*d2rho2(i) - Kap(3)*rhoc(1)/L**2*d2rho1(i)
    end do
  else
    do i=1,N
      mu_gen1(i) = mu1(i)/R/T - Kap(1)*d2rho1(i) - Kap(3)*d2rho2(i)
      mu_gen2(i) = mu2(i)/R/T - Kap(2)*d2rho2(i) - Kap(3)*d2rho1(i)
    end do
  end if

  call gradient(d2mu1,mu_gen1)
  call gradient(d2mu2,mu_gen2)

  if (dmlss) then
    do i=1,N
      D12(i) = Diff(2)*rhoc(1)*rho1_old(i)/(rhoc(1)*rho1_old(i) + rhoc(2)*rho2_old(i)) &
             + Diff(1)*rhoc(2)*rho2_old(i)/(rhoc(1)*rho1_old(i) + rhoc(2)*rho2_old(i))
      Jrho1(i) = rho1_old(i)*D12(i)*d2mu1(i)
      Jrho2(i) = rho2_old(i)*D12(i)*d2mu2(i)
    end do
  else
    do i=1,N
      D12(i) = Diff(2)*rho1_old(i)/(rho1_old(i) + rho2_old(i)) &
             + Diff(1)*rho2_old(i)/(rho1_old(i) + rho2_old(i))
      Jrho1(i) = rho1_old(i)*D12(i)*d2mu1(i)
      Jrho2(i) = rho2_old(i)*D12(i)*d2mu2(i)
    end do
  end if

  call gradient(dJrho1,Jrho1)
  call gradient(dJrho2,Jrho2)

  ! Integrate using explicit Euler scheme
  do i=1,N
    rho1_new(i) = rho1_old(i) + k*dJrho1(i)
    rho2_new(i) = rho2_old(i) + k*dJrho2(i)
  end do

  if (.not.per_BC) then
    ! Impose Dirichlet BCs
    rho1_new(N) = rho1_old(N) ! rho(N) = rho_L
    rho1_new(1) = rho1_old(1) ! rho(1) = rho_V
    rho2_new(N) = rho2_old(N) ! rho(N) = rho_V
    rho2_new(1) = rho2_old(1) ! rho(1) = rho_L
  end if

  ! Update the profile
  do i=1,N
    rho1_old(i) = rho1_new(i)
    rho2_old(i) = rho2_new(i)
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
          write(10,'(E10.2,2(F10.2))') x(i), rho1_old(i), rho2_old(i)
          write(*,'(E10.2,4(F10.2))') x(i), rho1_old(i), mu1(i), rho2_old(i), mu2(i)
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
  write(*,'(E10.2,4(F10.2))') x(i), rho1_old(i), mu1(i), rho2_old(i), mu2(i)
end do

write(*,*) '========================================================'
write(*,'(A,F9.2,A)') 'Time of operation was ', (time_end - time_begin)/60, ' minutes.'
write(*,*) '========================================================'
!========================= SURFACE TENSION =====================================

do i=1,N-2
  der1(i)  = ((rho1_old(i+2) - rho1_old(i))/(2*h))**2
  der12(i) = ((rho1_old(i+2) - rho1_old(i))/(2*h))*((rho2_old(i+2) - rho2_old(i))/(2*h))
  der2(i)  = ((rho2_old(i+2) - rho2_old(i))/(2*h))**2
end do

if (dmlss) then
  der = Kap(1)*rhoc(1)**2*der1 + 2*Kap(3)*rhoc(1)*rhoc(2)*der12 + Kap(2)*rhoc(2)**2*der2
else
  der = Kap(1)*der1 + 2*Kap(3)*der12 + Kap(2)*der2
end if

do i=1,N-3
  igr=igr+h*(der(i+1) + der(i))/2.0d0
end do

if (dmlss) then
  sig = R*T/L*igr
else
  sig = R*T*igr
end if

write(*,'(A,F9.2,A)') 'Integrated surface tension gamma:',real(sig*1.0d3),'mN•m-1 (DIVIDED BY THE NUMBER OF CREATED INTERFACES!!!).'
write(*,*) '========================================================'
write(*,'(A,/,(A,E,X,A,/))') 'EOS parameters n-hexane:','a = ',a(1),'J•m3/mol2','b = ',b(1),'m3/mol','Kappa = ',Kap(1),'m5/mol'
write(*,'(A,/,(A,E,X,A,/))') 'EOS parameters N2:','a = ',a(2),'J•m3/mol2','b = ',b(2),'m3/mol','Kappa = ',Kap(2),'m5/mol'

open(20,file=''//trim(path)//'output.dat',status='unknown')
  write(20,'(A,I2.2)') 'Simulation number ',simNo
  write(20,'(A,/,A,I,/,A,I,/,(A,E,/))') 'Simulation data:','N = ',N,'M = ',M,'L = ',L,'k = ',k,'h = ',h,'T =',T,'p =',pressure
  write(20,'(A,F9.2,A)') 'Integrated surface tension gamma:',real(sig*1.0d3),'mN•m-1 (DIVIDED BY THE NUMBER OF CREATED INTERFACES!!!).'
  write(20,*) '========================================================'
  write(20,'(A,/,(A,E,X,A,/))') 'EOS parameters n-hexane:','a = ',a(1),'J•m3/mol2','b = ',b(1),'m3/mol','Kappa = ',Kap(1),'m5/mol'
  write(20,'(A,/,(A,E,X,A,/))') 'EOS parameters nitrogen:','a = ',a(2),'J•m3/mol2','b = ',b(2),'m3/mol','Kappa = ',Kap(2),'m5/mol'
close(20)


deallocate(rho1_new,rho1_old,x,mu1,dJrho1,d2rho1,mu_gen1,d2mu1,Jrho1)
deallocate(rho2_new,rho2_old,mu2,dJrho2,d2rho2,mu_gen2,d2mu2,Jrho2)
deallocate(D12,der1,der12,der2,der,ampl)
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
subroutine potential(mu1,mu2,rho1_old,rho2_old)
!===============================================================================
use io
implicit none
! Subroutine for calculation of chemical potential acc. to Peng-Robinson's EoS
integer :: i
double precision,intent(in) :: rho1_old(N), rho2_old(N)
double precision,intent(out) :: mu1(N),mu2(N)
double precision :: rho1(N), rho2(N), rho(N), aa(N), bb(N)
call readpars

! Calculate potential profile for pertinent density profile

if (dmlss) then
  rho1 = rho1_old * rhoc(1)
  rho2 = rho2_old * rhoc(2)
else
  rho1 = rho1_old
  rho2 = rho2_old
end if

!!$OMP PARALLEL PRIVATE(mu1,mu2), SHARED(rho,aa,bb) num_threads(2)
!!$OMP PARALLEL DO
do i=1,N
  rho(i) = rho1(i) + rho2(i)
  aa(i) = a(1)*(rho1(i)**2)/(rho(i)**2) + a(2)*(rho2(i)**2)/(rho(i)**2) + 2*dsqrt(a(1)*a(2))*(rho1(i)*rho2(i))/(rho(i)**2)
  bb(i) = b(1)*rho1(i)/rho(i) + b(2)*rho2(i)/rho(i)
  mu1(i) = R*T*dlog(dabs(rho1(i)/(1.0d0 - bb(i)*rho(i)))) + R*T*b(1)*rho(i)/(1.0d0 - bb(i)*rho(i)) &
      + (1.0d0/2.0d0/dsqrt(2.0d0))*((2.0d0*a(1)*rho1(i)+2.0d0*dsqrt(a(1)*a(2))*rho2(i))/bb(i)/rho(i) - aa(i)*b(1)/bb(i)**2) &
      * dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*bb(i)*rho(i))/(1.0d0 + (1.0d0 + dsqrt(2.0d0))*bb(i)*rho(i)))) &
      - aa(i)*rho(i)*b(1)/bb(i)/(1.0d0 + 2.0d0*bb(i)*rho(i) - bb(i)**2*rho(i)**2)
  mu2(i) = R*T*dlog(dabs(rho2(i)/(1.0d0 - bb(i)*rho(i)))) + R*T*b(2)*rho(i)/(1.0d0 - bb(i)*rho(i)) &
      + (1.0d0/2.0d0/dsqrt(2.0d0))*((2.0d0*a(2)*rho2(i)+2.0d0*dsqrt(a(1)*a(2))*rho1(i))/bb(i)/rho(i) - aa(i)*b(2)/bb(i)**2) &
      * dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*bb(i)*rho(i))/(1.0d0 + (1.0d0 + dsqrt(2.0d0))*bb(i)*rho(i)))) &
      - aa(i)*rho(i)*b(2)/bb(i)/(1.0d0 + 2.0d0*bb(i)*rho(i) - bb(i)**2*rho(i)**2)
end do
!!$OMP END PARALLEL
return
end subroutine
!===============================================================================
