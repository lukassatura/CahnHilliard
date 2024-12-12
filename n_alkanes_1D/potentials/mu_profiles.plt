reset

set ytics nomirror
set xlabel 'Position x[m]'
set ylabel 'Hexane Potential {/Symbol m}_1'

list = system('ls mu*.dat')

plot for [i in list] i using 1:2 title 'Hexane Iteration '.i with lines

pause -1
