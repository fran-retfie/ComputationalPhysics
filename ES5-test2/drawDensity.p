
set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

#set title "Electron density"
#set xlabel "r"
#set ylabel "{/Symbol r}(r)"
#set key autotitle columnheader

set grid
set key off
splot 'initialPos.csv' using 1:2:3:(sprintf("%d", $4)) with labels point pt 7 ps 3 lc rgb "blue" offset char 0.4,0.4,0.4

pause -1
