
set datafile separator ','

set terminal pngcairo size 700, 550 font 'Verdana,15' enhanced
set output 'plots/Na-density-interactiong.png'

set title "Electron density - Na cluster"
set xlabel "r [a.u]"
set ylabel "{/Symbol r}(r)"
set key autotitle columnheader
set grid

filename(n) = sprintf("dati/density%02i", n)
plot for [n=1:3] filename(2*n) using 2:3 with lines lw 2, \
     "dati/NaData.csv" using 1:($2*3.76*10**-3) with lines title "Van Giai" lw 2
#plot for [n=1:6] filename(n) using 2:4 with lines, \
#     for [n=1:6] filename(n) using 2:5 with lines
