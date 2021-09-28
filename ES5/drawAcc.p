
#set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

#set title "Electron density"
#set xlabel "r"
#set ylabel "{/Symbol r}(r)"
#set key autotitle columnheader

set grid
set key off
plot 'ThermAcc.csv' using 1:2 with lines

pause -1
