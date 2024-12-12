reset

if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .14
if (!exists("MP_TOP"))    MP_TOP = 0.9
if (!exists("MP_xGAP"))   MP_xGAP = 0.1
if (!exists("MP_yGAP"))   MP_yGAP = 0.05

list = system('ls mu*.dat')

set multiplot layout 2,1 columnsfirst title "{/:Bold=15 Evolution of chemical potentials in a binary mixture}" \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP

unset title

set format y "%.1f"

set key box opaque
set ylabel '{/:Bold=12 N_2 Potential {/Symbol m}_2}'
unset xtics
set ytics font ",12"
plot for [i in list] i using 1:3 title "" with lines lt 1
set xtics font ",12" nomirror

set xlabel '{/:Bold=12 Position x[m]}'
set ylabel '{/:Bold=12 Hexane Potential {/Symbol m}_1}'
plot for [i in list] i using 1:2 title "" with lines lt 2

unset ylabel
unset ytics
unset xlabel

unset multiplot
pause -1
