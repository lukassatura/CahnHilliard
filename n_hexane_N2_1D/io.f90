!===============================================================================
!
! Input module for Cahn-Hilliard model in main.f90 containing hexane and nitrogen.
!
!
! Author: Lukáš Šatura
! 2021
!
!===============================================================================
module io

implicit none
!===============================================================================
! Physical parameters
!===============================================================================
double precision :: rho_L(2),rho_V(2),rhoc(2)                                     ! MOLAR concentrations of bulk (L) and (V), and critical phases in [mol/m3]
double precision :: NA,kB,R                                                       ! Avogadro, Boltzmann, and the gas constants in base units
double precision :: T,pressure,Tr(2),Tc(2),pc(2),vc(2),w(2)                       ! System temperature and critical parameters of the system
double precision :: a(2),b(2),Kap(3)                                              ! Parameters of the Peng-Robinson EoS
double precision :: kappa(2), kappa_0(2), kappa_1(2), kappa_2(2), kappa_3(2), alpha(2)
double precision :: Diff(2)                                                       ! Self-diffusivity of the system [m2/s]
! Numerical parameters
double precision :: L                 ! Length of the space domain in [m]
double precision :: h                 ! Dimensionless spatial discretization step
double precision :: k                 ! Dimensionless temporal discretization step
integer          :: M,N               ! Number of iterations and number of spatial nodes
! Simulation parameters
logical          :: init              ! Logical switch to start from default (linear) init. conditions
integer          :: init_iter, simNo  ! The profile to star from if 'init' is .false.
logical          :: per_BC,dmlss      ! Logical switch for PERIODIC init. cond. (DIRICHLET BC applied if per_BC = .false.), and for dimensionless formulation of CH eq.
integer,dimension(60):: plot_it       ! Iterations to be saved
character(len=100) :: path            ! Path for file storage
!-------
contains
subroutine readpars
!============================ SET THESE VARIABLES =============================
simNo = 01

path = ''
dmlss = .false.
!dmlss = .true.
N = 100                               ! Number of spatial nodes
M = 1000000                           ! Number of iterations
init = .true.                         ! .false. starts from 'init_iter' profile
!init = .false.                        ! Logical switch to start from linear init. conditions
init_iter = 500000000                   ! Initial iteration in case 'init' is false
per_BC = .true.                       ! Logical switch for periodic boundary conditions
per_BC = .false.                      ! Dirichlet BCs imposed if .false.
L = 1.0d-8                            ! Length of the space domain in [m]
if (dmlss) then
  k = 1.0d-8         			              ! Temporal discretization step (s)
  h = 1.d0/dble(N)	                    ! Dimensionless spatial discretization step
else
  k = 1.0d-21
  h = L/dble(N)
end if

!==============================================================================
! Physical parameters
NA = 6.02214076d23 ! Precisely, acc. to 2019 redefinition of SI base units
kB = 1.380649d-23  ! Precisely, acc. to 2019 redefinition of SI base units
R = NA*kB
T = 298.15d0		      	! Temperature (K)

!~~~~~~~~~~~~~HEXANE, Nitrogen
! Data from 'PRSV: An Improved Peng-Robinson Equation of State for Pure Compounds and Mixtures R. STRYJEK and J. H. VERA', 1986
! and 'PRSV2: A Cubic Equation of State for Accurate Vapor - Liquid Equilibria Calculations R. STRYJEK and J. H. VERA', 1986
Tc = (/507.3d0, 126.2d0/)
pc = (/3012360d0, 3400000d0/)
vc = (/370.0d-6, 89.8d-6/)
w = (/0.30075d0, 0.03726d0/)

kappa_1 = (/0.05104d0, 0.01996d0/)
kappa_2 = (/0.8634d0, 0.3162d0/)
kappa_3 = (/0.460d0, 0.535d0/)

! a = (/3.82713727d0, 0.087297706d0/)   	! Parameter a for hexane
! b = (/1.09064444d-4, 2.40809d-5/)  		 ! Parameter b for hexane
! Kap = (/1.50728d-22, 1.25598d-24, 1.3759d-23/) ! K1, K2, K12

Tr = T/Tc
kappa_0 = 0.378893d0 + 1.4897153d0*w - 0.17131848d0*w**2 + 0.0196544d0*w**3
kappa   = kappa_0 + (kappa_1 + kappa_2*(kappa_3 - Tr)*(1.0d0 - dsqrt(Tr)))*(1.0d0 + dsqrt(Tr))*(0.7d0 - Tr)
!kappa   = (0.37464d0 + 1.54226d0*w - 0.26992d0*w**2)
kappa(2) = kappa_0(2)
alpha = (1.0d0+kappa*(1-dsqrt(Tr)))**2
a = 0.457235d0*R**2*Tc**2/pc*alpha
b = 0.077796d0*R*Tc/pc

Kap(1:2) = 0.305d0*a*b**(2.0d0/3.0d0)*NA**(-2.0d0/3.0d0)/R/T
Kap(3) = dsqrt(Kap(1)*Kap(2))

rhoc = 1.0d0/vc

! rho_L = (/7.529714995075142d3, 0.312406981655924d3/)  ! (/hexane, N2/) !p = 20.0d0*101325.0d0 Pressure (Pa)
! rho_V = (/0.011640317416777d3, 0.815112143774716d3/) ! (/hexane, N2/)  !p = 20.0d0*101325.0d0 Pressure (Pa)
!rhoc = (/2.702702702702703d3, 1.113585746102450d4/)   ! (/hexane, N2/)
!rho_V = (/0.008455809316442d3, 0.032623527012931d3/)  !p = 101325.0d0 ! Pressure (Pa)
!rho_L = (/7.630985046361859d3, 0.012610015377510d3/)  !p = 101325.0d0 ! Pressure (Pa)

! rho_V = (/8.36559191585797d0, 32.7108383485089d0/) !p = 101325.0d0
! rho_L = (/7643.21881381869d0,	12.6247212317437d0/)
! pressure = 101325d0
! rho_V = (/8.51045732529446d0, 73.6065332194497d0/)	 !p = 2*101325.0d0
! rho_L = (/7637.88956614232d0, 28.3938488814009d0/)
! pressure = 2.0d0*101325d0
! rho_V = (/8.65748607225265d0, 114.538802401486d0/)		!p = 3*101325.0d0
! rho_L = (/7632.56231024064d0, 44.1612609197024d0/)
! pressure = 3*101325.0d0
! rho_V = (/8.95813142266327d0, 196.509838244133d0/)		!p = 5*01325.0d0
! rho_L = (/7621.91374381170d0, 75.6910308035023d0/)
! pressure = 5*101325.0d0
rho_V = (/9.42593920715926d0, 319.720495145512d0/)		!p = 8*101325.0d0
rho_L = (/7605.95564009330d0, 122.973408306950d0/)
pressure = 8*101325.0d0

Diff = (/4.21d-9, 2.12d-5/)   			 ! Hexane, N2 Diffusivities

plot_it = (/10,        50,        100,       200,       500,       1000,      &
            10000,     100000,    200000,    500000,    700000,    900000,    &
            1000000,   2000000,   5000000,   10000000,  20000000,  50000000,  &
            70000000,  100000000, 130000000, 150000000, 200000000, 210000000, 230000000, &
            250000000, 270000000, 280000000, 290000000, 300000000, 350000000, &
            380000000, 400000000, 420000000, 450000000, 470000000, 490000000, &
            500000000, 510000000, 520000000, 550000000, 560000000, 580000000, &
            590000000, 600000000, 610000000, 630000000, 650000000,            &
            680000000, 700000000, 720000000, 750000000, 770000000, 790000000, &
            800000000, 820000000, 850000000, 880000000, 900000000, 950000000/)

return
end subroutine
!----------------------------------------------

end module
