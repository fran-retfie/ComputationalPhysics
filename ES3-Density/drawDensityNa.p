
set datafile separator ','

set terminal pngcairo size 700, 550 font 'Verdana,15' enhanced
set output 'plots/Na-density-noninteractiong.png'

set title "Independent electron model - Electron density - Na cluster"
set xlabel "r [a.u]"
set ylabel "{/Symbol r}(r)"
set key autotitle columnheader
set grid

filename(n) = sprintf("dati/density%02i", n)
plot for [n=1:4] filename(n) using 2:3 with lines lw 2 #title (sprintf("n = %i", n))

#pause -1
