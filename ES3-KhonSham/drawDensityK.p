
set datafile separator ','

set terminal pngcairo size 700, 550 font 'Verdana,15' enhanced
set output 'plots/K-density-interactiong.png'

set title "Electron density - K cluster"
set xlabel "r [a.u]"
set ylabel "{/Symbol r}(r)"
set key autotitle columnheader
set grid
set xrange[0:25]
set yrange[-0.0005:0.0035]

filename(n) = sprintf("dati/density%02i", n)
plot for [n=1:3] filename(2*n) using 2:3 with lines lw 3, \
     "dati/KData.csv" using 1:($2*3.15*10**-3) with lines dashtype 2  lw 1.5 lc "black" title "Van Giai"
#plot for [n=1:6] filename(n) using 2:4 with lines, \
#     for [n=1:6] filename(n) using 2:5 with lines
