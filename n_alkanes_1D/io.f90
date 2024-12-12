!===============================================================================
!
! Input module for Cahn-Hilliard model in main.f90 containing selected
! non-polar compounds.
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
double precision :: rho_L,rho_V,rhoc  ! MOLAR concentrations of bulk (L) and (V), and critical phases in [mol/m3]
double precision :: NA,kB,R           ! Avogadro, Boltzmann, and the gas constants in base units
double precision :: T,Tc,Tr,pc,vc,w   ! System temperature and critical parameters of the system (w stands for omega in Peng-Rob. EoS)
double precision :: a,b,Kap,kappa, kappa_0, kappa_1, kappa_2, kappa_3, alpha     ! Parameters of the Peng-Robinson EoS
double precision :: Ai,Bi,di,Diff     ! Self-diffusivity of the system [m2/s]
! Numerical parameters
double precision :: L                 ! Length of the space domain in [m]
double precision :: h                 ! (Dimensionless) spatial discretization step
double precision :: k                 ! (Dimensionless) temporal discretization step
integer          :: M,N               ! Number of iterations and number of spatial nodes
! Simulation parameters
logical          :: init              ! Logical switch to start from default (linear) init. conditions
integer          :: init_iter,simNo   ! Profile to start from if 'init' is .false., simNo is only an auxiliary info saved to output file
logical          :: per_BC, dmlss     ! Logical switch for PERIODIC init. cond. (DIRICHLET BC applied if per_BC = .false.), and for dimensionless formulation of CH eq.
integer,dimension(60):: plot_it       ! Iterations to be saved
character(len=100) :: path            ! Path for file storage
!===============================================================================
contains
subroutine readpars
!============================ SET THESE VARIABLES =============================
path = ''                             ! Path to save output profiles
simNo = 01
dmlss = .false.
dmlss = .true.

N = 100                               ! Number of spatial nodes
M = 1000000                           ! Number of iterations
!k = 1.0d-8         			              ! Temporal discretization step (s)
L = 1.0d-8                            ! Length of the space domain in [m]
init = .true.                         ! .false. starts from 'init_iter' profile
!init = .false.                        ! Logical switch to start from linear init. conditions
init_iter = 20000000                   ! Initial iteration in case 'init' is false
per_BC = .false.                       ! Logical switch for periodic boundary conditions
!per_BC = .true.                      ! Dirichlet BCs imposed if .false.

if (dmlss) then
  k = 1.0d-8         			              ! Temporal discretization step (s)
  h = 1.d0/dble(N)	                    ! Dimensionless spatial discretization step
else
  k = 1.0d-17
  h = L/dble(N)
end if
!==============================================================================
! Physical parameters
!==============================================================================
NA = 6.02214076d23                    ! Precisely, acc. to 2019 redefinition of SI base units
kB = 1.380649d-23                     ! Precisely, acc. to 2019 redefinition of SI base units
R = NA*kB
T = 286.15d0    	                    ! Temperature (K)
!~~~~~~~~~~~~~HEPTANE~~~~~~~~~~~~~
! Tc = 540.3d0
! pc = 2740000d0
! vc = 432.0d-6
! w = 0.349d0
! rho_L = 6.707498357577608d3
! rho_V = 0.002547430256345d3
!~~~~~~~~~~~~~HEXANE~~~~~~~~~~~~~
! Data from 'PRSV: An Improved Peng-Robinson Equation of State for Pure Compounds and Mixtures R. STRYJEK and J. H. VERA', 1986
! and 'PRSV2: A Cubic Equation of State for Accurate Vapor - Liquid Equilibria Calculations R. STRYJEK and J. H. VERA', 1986
Tc = 507.3d0;
pc = 3012360d0;
w = 0.30075d0;
rho_V = 4.95556532139008d0
rho_L = 7758.96556288623d0
rhoc = 2.324990271027093d3

kappa_1 = 0.05104d0
kappa_2 = 0.8634d0
kappa_3 = 0.460d0
Diff = 4.21d-9   			 ! Hexane Diffusivity
!~~~~~~~~~~~~~CYCLOHEXANE~~~~~~~~~~~~~
! Tc = 553.8d0
! pc = 4070000d0
! vc = 308.0d-6
! w = 0.212d0
! rho_V = 0.005515050448729d3
! rho_L = 9.628572130732756d3
!~~~~~~~~~~~~~BENZENE~~~~~~~~~~~~~
! Tc = 562.1d0
! pc = 4890000.0d0
! vc = 259.0d-6
! w = 0.212d0
! rho_V = 0.000544434735709d4
! rho_L = 1.144840555288749d4
!==============================================================================
Tr = T/Tc
kappa_0 = 0.378893d0 + 1.4897153d0*w - 0.17131848d0*w**2 + 0.0196544d0*w**3                                  ! Peng-Robinson-Stryjek-Vera modification
kappa   = kappa_0 + (kappa_1 + kappa_2*(kappa_3 - Tr)*(1.0d0 - dsqrt(Tr)))*(1.0d0 + dsqrt(Tr))*(0.7d0 - Tr)  ! Peng-Robinson-Stryjek-Vera modification
!kappa   = (0.37464d0 + 1.54226d0*w - 0.26992d0*w**2)                                                        ! Peng-Robinson
alpha = (1.0d0+kappa*(1-dsqrt(Tr)))**2
a = 0.457235d0*R**2*Tc**2/pc*alpha
b = 0.077796d0*R*Tc/pc
di = 0.305d0*NA**(-2.0d0/3.0d0)

!Ai = - 10.0d-16/(1.2326d0 + 1.3757d0*w)
!Bi = + 10.0d-16/(0.9051d0 + 1.5410d0*w)
!di = Ai*(1 - Tr) + Bi           ! An alternative formula from literature

Kap = di*a*b**(2.0d0/3.0d0)/R/T  ! Kap is equivalent to the interaction parameter denoted also as 'c' in literature
!rhoc = 1.0d0/vc
plot_it = ( 100,      200,        500,      1000,       5000, &
            10000,     20000,    50000,      100000,   200000,     500000, &
            1000000,   2000000,  5000000,    10000000,   20000000, 50000000, &
            70000000, 100000000, 130000000, 150000000, 200000000, 230000000, &
            250000000, 270000000, 280000000, 290000000, 300000000, 350000000, &
            380000000, 400000000,420000000,  450000000, 470000000, 490000000,&
            500000000, 510000000,520000000, 550000000, 560000000,  580000000, &
            590000000, 600000000, 610000000, 630000000, 650000000, 670000000,   &
            680000000, 700000000,720000000, 750000000, 770000000, 790000000,&
            800000000, 820000000,850000000, 880000000, 900000000,950000000, 1000000000/)

return
end subroutine
!==============================================================================
end module
