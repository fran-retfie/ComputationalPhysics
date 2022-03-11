
set datafile separator ','

set terminal pngcairo size 700, 550 font 'Verdana,15' enhanced
set output 'plots/K-Spillout.png'

set title "Static polarizability of K cluster"
set xlabel "N^{-1/3}"
set ylabel "{/Symbol a}/N [a.u]"
set key autotitle columnheader
set grid

plot "dati/SpilloutK" using 4:5 lw 2, \
     "dati/JM.csv"    using 1:2 title "JM"   lw 2, \
     "dati/PHJM.csv"  using 1:2 title "PHJM" lw 2
