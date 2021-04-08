
set datafile separator ','

set xrange [0.5:5.5]
set yrange [0:5]
set xlabel "n"
set ylabel "E"
set grid xtics ytics mytics
set mytics 2

plot "energies.csv" using 1:3 with points pointtype 1

pause -1
