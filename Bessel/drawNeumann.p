
#set datafile separator ','
set terminal pngcairo size 1000, 800 enhanced dashed
set output 'plots/neumann.png'

set title "Bessel function"
set xlabel "x"
set ylabel "j(r)"
set yrange [0.5:-0.5]
set xrange [0:15]

plot for [l=2:6] "neumann.csv" using 1:l with lines linewidth 3 lc "orange" notitle ,\
     for [l=2:6] "neumannteo.csv" using 1:l with lines  dashtype 2 lc "black" notitle


pause -1
