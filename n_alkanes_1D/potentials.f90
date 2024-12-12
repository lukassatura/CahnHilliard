program main

use io
implicit none

double precision,allocatable :: rho(:,:),x(:),mu(:,:),d2rho(:),d1mu(:),Jrho(:,:)
integer                      :: i,j
integer                      :: file_num               ! Number of (pro)files
character(len=100)           :: filename, path_in, path_out, iter

call readpars
!N  = 100 ! Number of spatial nodes
file_num  = size(pack(plot_it,plot_it.le.M)) + 1 ! Number of evaluated profiles
path_in = ''                                     ! (Absolute) path to input files
path_out = './potentials/'                       ! (Absolute) path to output files
allocate(rho(N,file_num),x(N),mu(N,file_num),d2rho(N),d1mu(N),Jrho(N,file_num))

! Reading input files with density profiles
do j=1,file_num

  if (j.eq.1) then
    open(10,file=''//trim(path_in)//'rho_ini.dat',status='old')
  else
    write(iter,'(I10.10)') plot_it(j-1)
    filename = 'rho_'//trim(iter)//'.dat'
    open(10,file=''//trim(path_in)//trim(filename)//'',status='old')
  end if

  do i=1,N
    read(10,*) x(i), rho(i,j)
  end do
  close(10)
end do

rho = rho*rhoc

do j=1,file_num
  call laplace(d2rho,rho(:,j))
  do i=1,N
    mu(i,j) = R*T*dlog(dabs(rho(i,j)/(1.0d0 - b*rho(i,j)))) &
               + R*T*b*rho(i,j)/(1.0d0 - b*rho(i,j)) &
               + (a/2.0d0/dsqrt(2.0d0)/b)*dlog(dabs((1.0d0 + (1.0d0 - dsqrt(2.0d0))*b*rho(i,j)) &
               /(1.0d0 + (1.0d0 + dsqrt(2.0d0))*b*rho(i,j)))) &
               - a*rho(i,j)/(1.0d0 + 2.0d0*b*rho(i,j) - b**2*rho(i,j)**2)

    mu(i,j) = mu(i,j)/R/T - Kap/L**2*d2rho(i)
    !mu(i,j) = d2rho(i)
  end do

  call gradient(d1mu,mu(:,j),N,h)
  Jrho(:,j)  =  - rho(:,j)*d1mu
end do

mu = reshape(mu,(/N,file_num/))
Jrho = reshape(Jrho,(/N,file_num/))

! writing files with potential profiles
do j=1,file_num
  if (j.eq.1) then
    open(10,file=''//trim(path_out)//'mu_ini.dat',status='unknown')
  else
    write(iter,'(I10.10)') plot_it(j-1)
    filename = 'mu_'//trim(iter)//'.dat'
    open(10,file=''//trim(path_out)//trim(filename)//'',status='unknown')
  end if

  do i=1,N
    write(10,*) x(i), mu(i,j)
  end do
  close(10)
end do

! writing files with flux profiles
! do j=1,file_num
!
!   if (j.eq.1) then
!     open(10,file=''//trim(path_out)//'J_ini.dat',status='unknown')
!   else
!     write(iter,'(I10.10)') plot_it(j-1)
!     filename = 'J_'//trim(iter)//'.dat'
!     open(10,file=''//trim(path_out)//trim(filename)//'',status='unknown')
!   end if
!
!   do i=1,N
!     write(10,*) x(i), Jrho(i,j)
!   end do
! close(10)
! end do

deallocate(rho,x,mu,d2rho,d1mu,Jrho)
end program

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
!===============================================================
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
