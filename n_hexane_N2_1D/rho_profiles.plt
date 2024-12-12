reset

set ytics nomirror
set y2tics
set xlabel 'Position x[m]'

set ylabel 'N_2 Density {/Symbol r}_2'
set y2label 'Hexane Density {/Symbol r}_1'


list = system('ls rho\_*.dat')

plot for [i in list] i using 1:3 title 'N2 Iteration '.i axes x1y1 with lines
replot for [i in list] i using 1:2 title 'Hexane Iteration '.i axes x1y2 with lines



pause -1
