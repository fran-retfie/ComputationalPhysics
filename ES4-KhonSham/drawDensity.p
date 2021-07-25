
set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

set title "Electron density"
set xlabel "r"
set ylabel "{/Symbol r}(r)"
set key autotitle columnheader

filename(n) = sprintf("dati/density%02i", n)
#plot for [n=1:4] filename(n) using 2:3 with lines #title (sprintf("n = %i", n))
plot for [n=1:4] filename(n) using 2:4 with lines, \
     for [n=1:4] filename(n) using 2:5 with lines #title (sprintf("n = %i", n))

pause -1
