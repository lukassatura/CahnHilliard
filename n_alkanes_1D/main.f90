!===============================================================================
!
! Cahn-Hilliard model for dynamic simulation of vapour-liquid interfaces
! in one-dimensional systems of non-polar compounds.
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
double precision, allocatable :: rho_new(:),rho_old(:),dJrho(:),d2rho(:),Jrho(:)  ! Density fields with respective derivatives
double precision, allocatable :: x(:)                                             ! Array of spatial position
double precision, allocatable :: mu(:),mu_gen(:),d2mu(:)                          ! Chemical potential with respective derivatives
double precision, allocatable :: ampl(:)                                          ! Random amplitude of noise-like initial condition for periodic BC
double precision, allocatable :: der(:)                                           ! Auxiliary for surface tension integration der = (drho/dx)^2
double precision :: factor                                                        ! Factor of numerical stability
double precision :: sig,igr                                                       ! Integrated surface tension (sig) with respective integral value (igr)
! Indeces
integer :: j                                                                      ! Iteration counter
integer :: i                                                                      ! Spatial index
integer :: p,No_plots                                                             ! Plotting indeces
! Miscellaneous
character(len=100) :: filename                                                    ! String representing output filename
character(len=10) :: iter                                                         ! String representing iteration number
real :: time_begin, time_end                                                      ! Variables for measuring total CPU time
double precision, parameter :: Pi = 3.141592653589793238462643d0
!===============================================================================
call readpars
!===============================================================================
!!! ZISTIŤ SKUTOČNÚ PODMIENKU STABILTY A DOPLNIŤ !!!
! Check numerical stability of the proposed scheme
factor = Kap*rhoc/L**2
write(*,*) 'factor of stability CFL = ',factor
write(*,*) 'spatial step h = ',h

if (factor.gt.0.5d0) then
  write(*,*) 'Numerical stability condition not fulfilled.'
  write(*,*) 'Adjust spatial and/or temporal step.'
  !stop
end if
!!! ========================================== !!!

! Allocate arrays of the variables
allocate(rho_new(N),rho_old(N),x(N),mu(N),dJrho(N),d2rho(N),mu_gen(N),d2mu(N),Jrho(N),der(N-2),ampl(N))

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
    if (dmlss) then
      do i=1,N
        rho_old(i) = (rho_V+rho_L)/(2*rhoc) + ampl(i)/10.0d0*cos(40.0d0*Pi*x(i))
      end do
    else
      do i=1,N
        rho_old(i) = (rho_V+rho_L)/(2.0d0) + 50.0d0*ampl(i)/1.0d0*cos(400.0d0*Pi*x(i))
      end do
    end if
  else
    if (dmlss) then
      call linspace(rho_old,rho_V/rhoc,rho_L/rhoc,N)
    else
      call linspace(rho_old,rho_V,rho_L,N)
    end if
    ! Create a linear density profile
  end if



  ! Save the initial density profile into a file
  open(10,file=''//trim(path)//'rho_ini.dat',status='unknown')
  do i=1,N
    write(10,*) x(i), rho_old(i)
    write(*,*) x(i), rho_old(i)
  end do
  close(10)
  init_iter = 1
else ! Initialisation from iteration 'init_iter'
  write(iter,'(I10.10)') init_iter
  filename = 'rho_'//trim(iter)//'.dat'
  open(10,file=''//trim(path)//trim(filename)//'',status='old')
  do i=1,N
    read(10,*) x(i), rho_old(i)
    write(*,*) x(i), rho_old(i)
  end do
  close(10)
end if
!===============================================================================
! Start measuring time of computation
call CPU_TIME ( time_begin )

write(*,*) 'Iteration / Total'

!================================= MODEL =======================================
do j=init_iter,M

  if (per_BC) then
    ! Impose periodic BCs
    rho_old(1:2) = rho_old(N-3:N-2)
    rho_old(N-1:N) = rho_old(3:4)
  end if

  call potential(mu,rho_old,N)
  call laplace(d2rho,rho_old)
  ! if (per_BC) then
  !   mu(1:2) = mu(N-3:N-2)
  !   mu(N-1:N) = mu(3:4)
  !   d2rho(1:2) = d2rho(N-3:N-2)
  !   d2rho(N-1:N) =d2rho(3:4)
  ! end if

  !$OMP PARALLEL DO
  if (dmlss) then
    do i=1,N
      mu_gen(i) = mu(i)/R/T - Kap*rhoc/L**2*d2rho(i)
    end do
  else
    do i=1,N
      mu_gen(i) = mu(i)/R/T - Kap*d2rho(i)
    end do
  end if

  call gradient(d2mu,mu_gen)

  !$OMP PARALLEL DO
  if (dmlss) then
    do i=1,N
      Jrho(i) = rho_old(i)*d2mu(i)
    end do
  else
    do i=1,N
      Jrho(i) = rho_old(i)*Diff*d2mu(i)
    end do
  end if

  call gradient(dJrho,Jrho)

  ! Integrate using explicit Euler scheme
  !$OMP PARALLEL DO
  do i=1,N
    rho_new(i) = rho_old(i) + k*dJrho(i)
  end do

  if (.not.per_BC) then
    ! Impose Dirichlet BCs
    rho_new(N) = rho_old(N)   ! rho(N) = rho_L
    rho_new(1) = rho_old(1)   ! rho(1) = rho_V
  end if

  ! Update the profile
  !$OMP PARALLEL DO
  do i=1,N
    rho_old(i) = rho_new(i)
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
          write(10,*) x(i), rho_old(i)
          write(*,*) x(i), rho_old(i)
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
write(*,*) '      x              rho_old                mu'
do i=1,N
  write(*,*) x(i), rho_old(i), mu(i)
end do

write(*,*) '========================================================'
write(*,'(A,F9.2,A)') 'Time of operation was ', (time_end - time_begin)/60, ' minutes.'
write(*,*) '========================================================'
!========================= SURFACE TENSION =====================================
!rho_old = rho_old*rhoc

do i=1,N-2
  der(i) = ((rho_old(i+2)-rho_old(i))/(2*h))**2
end do

do i=1,N-3
  igr=igr+h*(der(i+1)+der(i))/2.0d0                                           ! Trapezoidal numerical integration
end do

if (dmlss) then
  sig = Kap*R*T*rhoc**2/L*igr
else
  sig = Kap*R*T*igr
end if

write(*,'(A,F9.2,A)') 'Integrated surface tension gamma:',real(sig*1.0d3),'mN•m-1 (DIVIDED BY THE NUMBER OF CREATED INTERFACES!!!).'
write(*,*) '========================================================'
write(*,'(A,/,(A,E,X,A,/))') 'EOS parameters:','a = ',a,'J•m3/mol2','b = ',b,'m3/mol','Kappa = ',Kap,'m5/mol'

open(20,file=''//trim(path)//'output.dat',status='unknown')

  write(20,'(A,I2.2)') 'Simulation number ',simNo
  write(20,'(A,/,A,I,/,A,I,/,(A,E,/))') 'Simulation data:','N = ',N,'M = ',M,'L = ',L,'k = ',k,'h = ',h,'T =',T
  write(20,'(A,F9.2,A)') 'Integrated surface tension gamma:',real(sig*1.0d3),'mN•m-1 (DIVIDED BY THE NUMBER OF CREATED INTERFACES!!!).'
  write(20,*) '========================================================'
  write(20,'(A,/,(A,E,X,A,/))') 'EOS parameters:','a = ',a,'J•m3/mol2','b = ',b,'m3/mol','Kappa = ',Kap,'m5/mol'
close(20)

deallocate(rho_new,rho_old,x,mu,dJrho,d2rho,mu_gen,d2mu,Jrho,der,ampl)
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
  !dy_dx(N) = (35.0d0*y(N) - 104.0d0*y(N-1) + 114.0d0*y(N-2) - 56.0d0*y(N-3) + 11.0d0*y(N-4))                  ! 3rd order accuracy
  !dy_dx(N) = (45.0d0*y(N) - 154.0d0*y(N-1) + 214.0d0*y(N-2) - 156.0d0*y(N-3) + 61.0d0*y(N-4) - 10.0d0*y(N-5)) ! 4th order accuracy

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
subroutine potential(mu,rho_old)
!===============================================================================
use io
implicit none
! Subroutine for calculation of chemical potential acc. to Peng-Robinson's EoS
integer :: i
double precision,intent(in) :: rho_old(N)
double precision,intent(out) :: mu(N)
double precision :: rho(N)

if (dmlss) then
  rho = rho_old*rhoc
else
  rho = rho_old
end if

! Calculate potential profile for pertinent density profile
!$OMP PARALLEL DO
do i=1,N
  mu(i) = R*T*dlog(dabs(rho(i)/(1.0d0 - b*rho(i)))) &
          + R*T*b*rho(i)/(1.0d0 - b*rho(i)) &
          + (a/2.0d0/dsqrt(2.0d0)/b)*dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*b*rho(i)) &
          /(1.0d0 + (1.0d0 + dsqrt(2.0d0))*b*rho(i)))) &
          - a*rho(i)/(1.0d0 + 2.0d0*b*rho(i) - b**2*rho(i)**2)
end do

return
end subroutine
!===============================================================================
