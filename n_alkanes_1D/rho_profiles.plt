reset

set ytics nomirror
set xlabel 'Position x[m]'
set ylabel 'Hexane Density {/Symbol r}_1'


list = system('ls /rho_profiles/rho\_*.dat')

plot for [i in list] i using 1:2 title 'Hexane Iteration '.i with lines

pause -1
