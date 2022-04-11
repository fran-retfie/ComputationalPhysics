
#set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

set title "Total energy vs b"
set xlabel "b"
#set ylabel "Energy"
#set key autotitle columnheader

set grid

plot "dati/result.csv" using 1:3:($3-$5):($3+$5) with errorbars, "dati/result.csv" using 1:2:($2-$4):($2+$4) with errorbars

#plot "dati/result.csv" using 1:3 w points ps 2, "dati/result.csv" using 1:2 w points ps 2

pause -1
