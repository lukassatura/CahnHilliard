!===============================================================================
!
! Input module for Cahn-Hilliard model in main.f90 containing water,
! n-dodecane and 1-butanol.
!
! Author: Lukáš Šatura
! 2021
!
!===============================================================================
module io

implicit none

! Physical parameters
double precision :: rho_L(3),rho_V(3),rhoc(3)          ! MOLAR concentrations of bulk (L) and (V), and critical phases in [mol/m3]
double precision :: NA,kB,R                            ! Avogadro, Boltzmann, and the gas constants in base units
double precision :: T,Tr(3),Tc(3),pc(3),vc(3),w(3)     ! System temperature and critical parameters of the system
double precision :: a(3),b(3),Kap(6)                   ! Parameters of the Peng-Robinson EoS
double precision :: kappa(3), kappa_0(3), kappa_1(3), kappa_2(3), kappa_3(3), alpha(3)
double precision :: Z(3)                               ! Chain lengths in monomers (from Flory-Huggins)
double precision :: Diff(3)                            ! Self-diffusivity of the system [m2/s]
! Numerical parameters
double precision :: L                 ! Length of the space domain in [m]
double precision :: h                 ! Dimensionless spatial discretization step
double precision :: k                 ! Dimensionless temporal discretization step
integer          :: M,N               ! Number of iterations and number of spatial nodes
! Simulation parameters
logical          :: init              ! Logical switch to start from default (linear) init. conditions
integer          :: init_iter,simNo   ! The profile to star from if 'init' is .false., simNo is only an auxiliary info saved to output file
logical          :: per_BC            ! Logical switch for PERIODIC init. cond. (DIRICHLET BC applied if per_BC = .false.)
integer,dimension(60):: plot_it       ! Iterations to be saved
character(len=100) :: path            ! Path for file storage
!-------
contains
subroutine readpars
!============================ SET THESE VARIABLES =============================
path = ''
N = 100                                 ! Number of spatial nodes
M = 500000000                           ! Number of iterations
k = 1.0d-18         			              ! Temporal discretization step (s)
L = 1.0d-8                              ! Length of the space domain in [m]
h = L/dble(N)	                          ! Dimensionless spatial discretization step
init = .true.                           ! .false. starts from 'init_iter' profile
!init = .false.                          ! Logical switch to start from linear init. conditions
init_iter = 110000000                   ! Initial iteration in case 'init' is false
per_BC = .true.                         ! Logical switch for periodic boundary conditions
!per_BC = .false.                      ! Dirichlet BCs imposed if .false.
!==============================================================================
! Physical parameters
NA = 6.02214076d23 ! Precisely, acc. to 2019 redefinition of SI base units
kB = 1.380649d-23  ! Precisely, acc. to 2019 redefinition of SI base units
R = NA*kB
T = 298.15d0		      	! Temperature (K)
!p = 20.0d0*101325.0d0 ! Pressure (Pa)
simNo = 01

! Data from 'PRSV: An Improved Peng-Robinson Equation of State for Pure Compounds and Mixtures R. STRYJEK and J. H. VERA', 1986
! and 'PRSV2: A Cubic Equation of State for Accurate Vapor - Liquid Equilibria Calculations R. STRYJEK and J. H. VERA', 1986
!~~~~~~~~~ Water ~~~~~~~~~~~~
Tc(1) = 647.286d0;    ! Trange 274–623
pc(1) = 22089750d0;
w(1)  = 0.3438d0;
kappa_1(1) = -0.06635d0;
kappa_2(1) = 0.0199d0;
kappa_3(1) = 0.443d0;
!~~~~~~~~~ Dodecane ~~~~~~~~~~~~~~~~~~
Tc(2) = 658.2d0;      !Trange 312-520
pc(2) = 1823830d0;
w(2)  = 0.57508d0;
kappa_1(2) = 0.05426d0;
kappa_2(2) = 0.8744d0;
kappa_3(2) = 0.505d0;
!~~~~~~~~~~ 1-Butanol ~~~~~~~~~~~~~~~~~
Tc(3) = 562.98d0;     !Trange 352–399
pc(3) = 4412660d0;
w(3)  = 0.59022d0;
kappa_1(3) = 0.33431d0;
kappa_2(3) = -1.1748d0;
kappa_3(3) = 0.642d0;

! a = (/3.82713727d0, 0.087297706d0/)   		! Parameter a for hexane
! b = (/1.09064444d-4, 2.40809d-5/)  		 ! Parameter b for hexane
! Kap = (/1.50728d-22, 1.25598d-24, 1.3759d-23/) ! K1, K2, K12
Tr = T/Tc
kappa_0 = 0.378893d0 + 1.4897153d0*w - 0.17131848d0*w**2 + 0.0196544d0*w**3
kappa   = kappa_0 + (kappa_1 + kappa_2*(kappa_3 - Tr)*(1.0d0 - dsqrt(Tr)))*(1.0d0 + dsqrt(Tr))*(0.7d0 - Tr)
!kappa   = (0.37464d0 + 1.54226d0*w - 0.26992d0*w**2)
alpha = (1.0d0+kappa*(1-dsqrt(Tr)))**2
a = 0.457235d0*R**2*Tc**2/pc*alpha
b = 0.077796d0*R*Tc/pc

Kap(1:3) = 0.305d0*a*b**(2.0d0/3.0d0)*NA**(-2.0d0/3.0d0)/R/T
Kap(4) = dsqrt(Kap(1)*Kap(2))
Kap(5) = dsqrt(Kap(1)*Kap(3))
Kap(6) = dsqrt(Kap(3)*Kap(2))

!rho_L = (/55345.5454d0, 4379.73346d0, 10846.982d0/)  ! (/water, dodecane, butanol/)
!rho_V = (/1.3044684d0, 43.878d0, 102.1195d0/) ! (/water, dodecane, butanol/)
rho_L = (/4379.73346d0, 4379.73346d0, 4379.73346d0/)  ! (/water, dodecane, butanol/)
rho_V = (/4379.73346d0, 4379.73346d0, 4379.73346d0/) ! (/water, dodecane, butanol/)
!rhoc = (/2.702702702702703d3, 1.113585746102450d4/)   ! (/hexane, N2/)

! a = (/3.816361207131405d0, 0.087297647623448d0/)
! b = (/1.087247226613335d-4, 2.408093219033991d-5/)
! Kap = (/1.499910654123076d-22, 1.255977278004353d-24, 1.372535500675750d-23/)

Diff = (/2.299d-9, 0.9d-9,0.775d-10/)   			 ! Hexane, N2 Diffusivities
Z = (/1.0d0, 1.0d0, 1.0d0/)

plot_it = (/1,         50,        100,       200,       500,       1000,      &
            10000,     100000,    200000,    500000,    700000,    900000,    &
            1000000,   2000000,   5000000,   10000000,  20000000,  30000000, &
            40000000,  50000000,  60000000,  70000000,  80000000,  90000000,  &
            100000000, 110000000, 120000000, 130000000, 140000000, 150000000, &
            160000000, 170000000, 180000000, 190000000, 200000000, 210000000, &
            220000000, 230000000, 240000000, 250000000, 260000000, 270000000, &
            280000000, 290000000, 300000000, 350000000, &
            380000000, 400000000, 420000000, 450000000, 470000000, 490000000, &
            500000000, 510000000, 520000000, 550000000, 560000000, 580000000, &
            590000000, 600000000/)

return
end subroutine
!----------------------------------------------

end module
