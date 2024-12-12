reset

set ytics nomirror
set y2tics
set xlabel 'Position x[m]'

set ylabel ' Density {/Symbol r}_2'
set y2label 'dodecane and butanol Density {/Symbol r}_1'


list = system('ls rho\_*.dat')

set multiplot layout 3,1 rowsfirst
# --- GRAPH a
set label 1 'water' at graph 0.92,0.9 font ',8'
plot for [i in list] i using 1:2 title '' with lines
# --- GRAPH b
set label 1 'dodecane' at graph 0.92,0.9 font ',8'
plot for [i in list] i using 1:3 title ''  with lines
# --- GRAPH c
set label 1 'butanol' at graph 0.92,0.9 font ',8'
plot for [i in list] i using 1:4 title ''  with lines
unset multiplot


pause -1
