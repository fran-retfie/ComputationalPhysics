
#set datafile separator ','
set terminal pngcairo size 1000, 800 enhanced dashed
set output 'plots/bessel.png'

set title "Bessel function"
set xlabel "x"
set ylabel "j(r)"
set xrange [0:15]
set yrange [-0.4:1.2]

plot for [l=2:6] "bessel.csv" using 1:l with lines linewidth 3 lc "orange" notitle ,\
     for [l=2:6] "besselteo.csv" using 1:l with lines  dashtype 2 lc "black" notitle

pause -1
