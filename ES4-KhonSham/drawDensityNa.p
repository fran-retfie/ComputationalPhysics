
set datafile separator ','

set terminal pngcairo size 700, 550 font 'Verdana,15' enhanced
set output 'plots/Na-density-interactiong.png'

set title "Electron density - Na cluster, N_e = 40"
set xlabel "r [a.u]"
set ylabel "{/Symbol r}(r)"
#set key autotitle columnheader
set grid

filename(n) = sprintf("dati/density%02i", n)
plot for [n=6:6] filename(n) using 2:3 with lines title "Our result" lw 2, \
     "dati/NaData.csv" using 1:($2*3.85*10**-3) with lines title "Van Giai" lw 2
#plot for [n=1:6] filename(n) using 2:4 with lines, \
#     for [n=1:6] filename(n) using 2:5 with lines
